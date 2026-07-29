"""Small transparent trainer for state pretraining and paired fine-tuning."""

from __future__ import annotations

import json
import random
from pathlib import Path
from typing import Any, Dict, List, Sequence, Tuple

import numpy as np
import torch
from torch.utils.data import DataLoader, Dataset, Subset

from flexppi.config import ExperimentConfig
from flexppi.data.dataset import PairDataset, StateDataset, collate_pairs, collate_states
from flexppi.losses import paired_loss, state_loss
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
        loss = state_loss(output, batch, weights)
        return {"loss": loss, "state_loss": loss}
    return paired_loss(
        output,
        batch,
        weights,
        state_weight=config.train.state_weight,
        response_weight=config.train.response_weight,
        consistency_weight=config.train.consistency_weight,
    )


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


def train(config: ExperimentConfig) -> Dict[str, Any]:
    seed_everything(config.seed)
    device = resolve_device(config.device)
    dataset, collate = _make_dataset_and_collate(config)
    if len(dataset) == 0:
        raise RuntimeError(
            "No labelled samples survived filtering. For DynaRepo, compute labels and rerun index-dynarepo first."
        )
    train_set, val_set = group_split(dataset, config.train.val_fraction, config.seed)
    train_loader = DataLoader(
        train_set,
        batch_size=config.data.batch_size,
        shuffle=True,
        num_workers=config.data.num_workers,
        collate_fn=collate,
    )
    val_loader = DataLoader(
        val_set,
        batch_size=config.data.batch_size,
        shuffle=False,
        num_workers=config.data.num_workers,
        collate_fn=collate,
    )
    model = FlexPPI(config.model).to(device)
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
    with (output_dir / "resolved_config.json").open("w", encoding="utf-8") as handle:
        json.dump(config.to_dict(), handle, indent=2, ensure_ascii=False)

    best_value = float("inf")
    history = []
    global_step = 0
    for epoch in range(1, config.train.epochs + 1):
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
        val_value = _validate(config, model, val_loader, device)
        selection_value = val_value if np.isfinite(val_value) else train_value
        row = {"epoch": epoch, "train_loss": train_value, "val_loss": val_value}
        history.append(row)
        print(row)
        checkpoint = {
            "epoch": epoch,
            "model": model.state_dict(),
            "optimizer": optimizer.state_dict(),
            "config": config.to_dict(),
            "train_loss": train_value,
            "val_loss": val_value,
        }
        if epoch % config.train.checkpoint_every == 0:
            torch.save(checkpoint, output_dir / "epoch_{:03d}.pt".format(epoch))
        if selection_value < best_value:
            best_value = selection_value
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
