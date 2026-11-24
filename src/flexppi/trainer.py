"""Small transparent trainer for state pretraining and paired fine-tuning."""

from __future__ import annotations

import json
import random
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import torch
from torch.utils.data import DataLoader, Dataset, Subset, WeightedRandomSampler

from flexppi.config import ExperimentConfig
from flexppi.data.dataset import PairDataset, StateDataset, collate_pairs, collate_states
from flexppi.data.records import load_state_records
from flexppi.losses import paired_loss, state_loss
from flexppi.model.information_ablation import build_information_ablation
from flexppi.model.network import FlexPPI


def seed_everything(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)


def resolve_device(requested: str) -> torch.device:
    if requested == "auto":
        return torch.device("cuda" if torch.cuda.is_available() else "cpu")
    return torch.device(requested)


def move_to_device(value: Any, device: torch.device) -> Any:
    if torch.is_tensor(value):
        return value.to(device)
    if isinstance(value, dict):
        return {key: move_to_device(item, device) for key, item in value.items()}
    if isinstance(value, list):
        return [move_to_device(item, device) for item in value]
    return value


def _record_group(record: Any) -> str:
    metadata = getattr(record, "metadata", {})
    return str(
        metadata.get("cluster_id")
        or metadata.get("pdb_id")
        or getattr(record, "sample_id", getattr(record, "pair_id", "unknown"))
    )


def group_split(dataset: Dataset, val_fraction: float, seed: int) -> Tuple[Subset, Subset]:
    records = getattr(dataset, "records", None)
    if records is None or len(records) < 2 or val_fraction <= 0:
        indices = list(range(len(dataset)))
        return Subset(dataset, indices), Subset(dataset, [])
    group_to_indices: Dict[str, List[int]] = {}
    for index, record in enumerate(records):
        group_to_indices.setdefault(_record_group(record), []).append(index)
    groups = sorted(group_to_indices)
    random.Random(seed).shuffle(groups)
    val_group_count = max(1, int(round(len(groups) * val_fraction)))
    val_groups = set(groups[:val_group_count])
    train_indices, val_indices = [], []
    for group, indices in group_to_indices.items():
        (val_indices if group in val_groups else train_indices).extend(indices)
    if not train_indices:
        train_indices.append(val_indices.pop())
    return Subset(dataset, train_indices), Subset(dataset, val_indices)


def _metric_weights(config: ExperimentConfig) -> List[float]:
    return [float(config.train.metric_weights.get(name, 1.0)) for name in config.model.metrics]


def _make_dataset_and_collate(config: ExperimentConfig) -> Tuple[Dataset, Any]:
    if config.stage == "pretrain":
        if not config.data.state_manifest:
            raise ValueError("pretrain requires data.state_manifest")
        return StateDataset(
            config.data.state_manifest,
            config.data.max_residues,
            config.data.path_mappings,
        ), collate_states
    if config.stage == "finetune":
        if not config.data.pair_manifest:
            raise ValueError("finetune requires data.pair_manifest")
        return PairDataset(
            config.data.pair_manifest,
            config.data.max_residues,
            config.data.path_mappings,
        ), collate_pairs
    raise ValueError("stage must be pretrain or finetune")


def _batch_loss(
    config: ExperimentConfig,
    model: FlexPPI,
    batch: Dict[str, Any],
) -> Dict[str, torch.Tensor]:
    output = model(batch)
    weights = _metric_weights(config)
    if config.stage == "pretrain":
        loss = state_loss(output, batch, weights, config.model.metrics)
        return {"loss": loss, "state_loss": loss}
    return paired_loss(
        output,
        batch,
        weights,
        state_weight=config.train.state_weight,
        response_weight=config.train.response_weight,
        consistency_weight=config.train.consistency_weight,
        metric_names=config.model.metrics,
    )


def _record_groups(dataset: Dataset) -> set:
    records = getattr(dataset, "records", [])
    return {_record_group(record) for record in records}


def _fixed_state_validation_sets(
    config: ExperimentConfig,
    train_dataset: Dataset,
) -> Tuple[Dataset, Dataset]:
    if not config.data.validation_state_manifest:
        return group_split(train_dataset, config.train.val_fraction, config.seed)
    validation_dataset = StateDataset(
        config.data.validation_state_manifest,
        config.data.max_residues,
        config.data.path_mappings,
    )
    overlap = sorted(_record_groups(train_dataset) & _record_groups(validation_dataset))
    if overlap:
        raise ValueError(
            "Fixed train/validation manifests share biological groups: {}".format(
                ", ".join(overlap[:10])
            )
        )
    return train_dataset, validation_dataset


def _records_for_dataset(dataset: Dataset) -> Sequence[Any]:
    if isinstance(dataset, Subset):
        parent_records = _records_for_dataset(dataset.dataset)
        return [parent_records[index] for index in dataset.indices]
    return list(getattr(dataset, "records", []))


def _sampling_system(record: Any) -> str:
    metadata = getattr(record, "metadata", {})
    return str(
        metadata.get("biological_system")
        or metadata.get("pdb_id")
        or getattr(record, "sample_id", "unknown")
    )


def _source_system_balanced_weights(dataset: Dataset) -> List[float]:
    """Give each source equal mass, then each system equal mass within a source."""

    records = _records_for_dataset(dataset)
    if len(records) != len(dataset):
        raise ValueError("Balanced sampling requires a record-backed dataset")
    source_groups: Dict[str, Dict[str, List[int]]] = {}
    for index, record in enumerate(records):
        source = str(getattr(record, "dataset", "unknown"))
        system = _sampling_system(record)
        source_groups.setdefault(source, {}).setdefault(system, []).append(index)
    weights = [0.0] * len(records)
    source_count = max(len(source_groups), 1)
    for systems in source_groups.values():
        system_count = max(len(systems), 1)
        for indices in systems.values():
            record_weight = 1.0 / (source_count * system_count * len(indices))
            for index in indices:
                weights[index] = record_weight
    return weights


def _training_sampler(
    dataset: Dataset,
    strategy: str,
    samples_per_epoch: Optional[int],
    seed: int,
) -> Optional[WeightedRandomSampler]:
    if strategy == "uniform_records":
        return None
    if strategy != "source_system_balanced":
        raise ValueError("sampling_strategy must be uniform_records or source_system_balanced")
    weights = _source_system_balanced_weights(dataset)
    generator = torch.Generator()
    generator.manual_seed(seed)
    return WeightedRandomSampler(
        weights,
        num_samples=samples_per_epoch or len(dataset),
        replacement=True,
        generator=generator,
    )


def _rng_state() -> Dict[str, Any]:
    generator, keys, position, has_gaussian, cached_gaussian = np.random.get_state()
    state: Dict[str, Any] = {
        "random": random.getstate(),
        "numpy": (generator, keys.tolist(), position, has_gaussian, cached_gaussian),
        "torch": torch.get_rng_state(),
    }
    if torch.cuda.is_available():
        state["cuda"] = torch.cuda.get_rng_state_all()
    return state


def _restore_rng_state(state: Dict[str, Any]) -> None:
    random.setstate(state["random"])
    np.random.set_state(state["numpy"])
    torch.set_rng_state(state["torch"])
    if "cuda" in state and torch.cuda.is_available():
        torch.cuda.set_rng_state_all(state["cuda"])


@torch.no_grad()
def _validate(
    config: ExperimentConfig,
    model: FlexPPI,
    loader: DataLoader,
    device: torch.device,
) -> float:
    if len(loader) == 0:
        return float("nan")
    model.eval()
    losses = []
    for batch in loader:
        batch = move_to_device(batch, device)
        losses.append(float(_batch_loss(config, model, batch)["loss"].detach().cpu()))
    return float(np.mean(losses))


def _check_no_silent_drops(config: ExperimentConfig, dataset: Dataset) -> None:
    """Fail loudly instead of quietly discarding manifest records."""

    if config.stage != "pretrain":
        return
    expected = 0
    if isinstance(config.data.state_manifest, str):
        expected = len(load_state_records(config.data.state_manifest))
    else:
        for path in config.data.state_manifest or []:
            expected += len(load_state_records(path))
    if len(dataset) != expected:
        kept = {record.sample_id for record in dataset.records}
        wanted = set()
        if isinstance(config.data.state_manifest, str):
            wanted = {record.sample_id for record in load_state_records(config.data.state_manifest)}
        else:
            for path in config.data.state_manifest or []:
                wanted.update(record.sample_id for record in load_state_records(path))
        dropped = sorted(wanted - kept)
        raise RuntimeError(
            "StateDataset kept {} of {} manifest records; dropped: {}".format(
                len(dataset), expected, ", ".join(dropped[:20])
            )
        )
    missing_files = []
    for record in dataset.records:
        if not Path(record.structure_path).exists():
            missing_files.append("structure:{}".format(record.sample_id))
        if record.label_path and not Path(record.label_path).exists():
            missing_files.append("label:{}".format(record.sample_id))
    if missing_files:
        raise RuntimeError(
            "Manifest records reference missing files: {}".format(
                ", ".join(missing_files[:20])
            )
        )


def _latest_epoch_checkpoint(output_dir: Path) -> Optional[Tuple[int, Path]]:
    checkpoints = sorted(
        (int(name.name.split("_")[1].split(".")[0]), name)
        for name in output_dir.glob("epoch_*.pt")
        if name.name.split("_")[1].split(".")[0].isdigit()
    )
    return checkpoints[-1] if checkpoints else None


def train(config: ExperimentConfig) -> Dict[str, Any]:
    seed_everything(config.seed)
    device = resolve_device(config.device)
    dataset, collate = _make_dataset_and_collate(config)
    if len(dataset) == 0:
        raise RuntimeError(
            "No labelled samples survived filtering. For DynaRepo, compute labels "
            "and rerun index-dynarepo first."
        )
    _check_no_silent_drops(config, dataset)
    if config.stage == "pretrain":
        train_set, val_set = _fixed_state_validation_sets(config, dataset)
    else:
        train_set, val_set = group_split(dataset, config.train.val_fraction, config.seed)
    model = build_information_ablation(
        config.model, config.model.information_ablation
    ).to(device)
    if config.init_checkpoint:
        checkpoint = torch.load(config.init_checkpoint, map_location="cpu", weights_only=False)
        missing, unexpected = model.load_state_dict(checkpoint["model"], strict=False)
        if unexpected:
            raise RuntimeError("Unexpected checkpoint parameters: {}".format(unexpected))
        if missing:
            print("Initial checkpoint is missing newly introduced parameters: {}".format(missing))
    parameters = [parameter for parameter in model.parameters() if parameter.requires_grad]
    optimizer = torch.optim.AdamW(
        parameters,
        lr=config.train.learning_rate,
        weight_decay=config.train.weight_decay,
    )
    output_dir = Path(config.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    start_epoch = 1
    global_step = 0
    best_value = float("inf")
    history: List[Dict[str, Any]] = []
    if config.train.resume:
        latest = _latest_epoch_checkpoint(output_dir)
        if latest is not None:
            epoch, checkpoint_path = latest
            checkpoint = torch.load(
                str(checkpoint_path), map_location="cpu", weights_only=False
            )
            if checkpoint.get("epoch", 0) >= config.train.epochs:
                raise RuntimeError(
                    "Training already reached {} epochs in {}; refusing to resume a "
                    "finished run. Set train.resume: false or choose another output_dir."
                    .format(config.train.epochs, output_dir)
                )
            model.load_state_dict(checkpoint["model"])
            if "optimizer" in checkpoint:
                optimizer.load_state_dict(checkpoint["optimizer"])
            start_epoch = int(checkpoint.get("epoch", epoch)) + 1
            global_step = int(checkpoint.get("global_step", 0))
            best_value = float(checkpoint.get("best_value", float("inf")))
            history = list(checkpoint.get("history", []))
            if "rng" in checkpoint:
                _restore_rng_state(checkpoint["rng"])
            print(
                "Resumed {} from epoch={} step={} best_loss={} (history entries: {})".format(
                    output_dir, start_epoch - 1, global_step, best_value, len(history)
                ))
    with (output_dir / "resolved_config.json").open("w", encoding="utf-8") as handle:
        json.dump(config.to_dict(), handle, indent=2, ensure_ascii=False)

    for epoch in range(start_epoch, config.train.epochs + 1):
        train_sampler = _training_sampler(
            train_set,
            config.train.sampling_strategy,
            config.train.samples_per_epoch,
            config.seed + epoch,
        )
        train_loader = DataLoader(
            train_set,
            batch_size=config.data.batch_size,
            shuffle=train_sampler is None,
            sampler=train_sampler,
            num_workers=config.data.num_workers,
            collate_fn=collate,
        )
        model.train()
        epoch_losses = []
        for batch in train_loader:
            batch = move_to_device(batch, device)
            optimizer.zero_grad(set_to_none=True)
            loss_values = _batch_loss(config, model, batch)
            loss_values["loss"].backward()
            torch.nn.utils.clip_grad_norm_(parameters, config.train.grad_clip)
            optimizer.step()
            global_step += 1
            epoch_losses.append(float(loss_values["loss"].detach().cpu()))
            if global_step % config.train.log_every == 0:
                printable = {
                    key: round(float(value.detach().cpu()), 5)
                    for key, value in loss_values.items()
                }
                print("epoch={} step={} {}".format(epoch, global_step, printable))
        train_value = float(np.mean(epoch_losses))
        val_loader = DataLoader(
            val_set,
            batch_size=config.data.batch_size,
            shuffle=False,
            num_workers=config.data.num_workers,
            collate_fn=collate,
        )
        val_value = _validate(config, model, val_loader, device)
        selection_value = val_value if np.isfinite(val_value) else train_value
        row = {"epoch": epoch, "train_loss": train_value, "val_loss": val_value}
        history.append(row)
        print(row)
        checkpoint = {
            "epoch": epoch,
            "global_step": global_step,
            "model": model.state_dict(),
            "optimizer": optimizer.state_dict(),
            "config": config.to_dict(),
            "train_loss": train_value,
            "val_loss": val_value,
            "best_value": best_value,
            "history": history,
            "rng": _rng_state(),
        }
        if epoch % config.train.checkpoint_every == 0:
            torch.save(checkpoint, output_dir / "epoch_{:03d}.pt".format(epoch))
        if selection_value < best_value:
            best_value = selection_value
            checkpoint["best_value"] = best_value
            torch.save(checkpoint, output_dir / "best.pt")
        with (output_dir / "history.json").open("w", encoding="utf-8") as handle:
            json.dump(history, handle, indent=2)
    return {
        "device": str(device),
        "train_samples": len(train_set),
        "val_samples": len(val_set),
        "best_loss": best_value,
        "checkpoint": str((output_dir / "best.pt").resolve()),
    }
