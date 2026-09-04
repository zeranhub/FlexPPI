"""Leakage-aware, protein-level evaluation for flexibility and binding response."""

from __future__ import annotations

import math
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import torch
from torch.utils.data import ConcatDataset, DataLoader, Dataset, Subset

from flexppi.config import ExperimentConfig, ModelConfig
from flexppi.constants import METRIC_NAMES, RESTYPES
from flexppi.data.dataset import PairDataset, StateDataset, collate_pairs, collate_states
from flexppi.data.labels import load_label_file
from flexppi.data.records import PairRecord
from flexppi.data.structures import parse_pdb, target_mask
from flexppi.model.network import FlexPPI, gather_residues
from flexppi.trainer import group_split, move_to_device, resolve_device

DEFAULT_DELTA_THRESHOLDS = {"rmsf": 0.10, "chi_entropy": 0.10, "sasa": 2.0}


def _as_float(value: float) -> Optional[float]:
    return float(value) if np.isfinite(value) else None


def _rankdata(values: np.ndarray) -> np.ndarray:
    """Return average ranks for ties without requiring SciPy."""

    values = np.asarray(values)
    order = np.argsort(values, kind="mergesort")
    sorted_values = values[order]
    ranks = np.empty(len(values), dtype=np.float64)
    start = 0
    while start < len(values):
        stop = start + 1
        while stop < len(values) and sorted_values[stop] == sorted_values[start]:
            stop += 1
        ranks[order[start:stop]] = 0.5 * (start + stop - 1) + 1.0
        start = stop
    return ranks


def _correlation(first: np.ndarray, second: np.ndarray) -> Optional[float]:
    if len(first) < 2 or np.std(first) == 0 or np.std(second) == 0:
        return None
    return _as_float(np.corrcoef(first, second)[0, 1])


def _spearman(first: np.ndarray, second: np.ndarray) -> Optional[float]:
    return _correlation(_rankdata(first), _rankdata(second))


def _binary_auroc(scores: np.ndarray, labels: np.ndarray) -> Optional[float]:
    labels = np.asarray(labels, dtype=bool)
    positives = int(labels.sum())
    negatives = int((~labels).sum())
    if positives == 0 or negatives == 0:
        return None
    ranks = _rankdata(np.asarray(scores, dtype=np.float64))
    value = (ranks[labels].sum() - positives * (positives + 1) / 2.0) / (
        positives * negatives
    )
    return _as_float(value)


def _average_precision(scores: np.ndarray, labels: np.ndarray) -> Optional[float]:
    labels = np.asarray(labels, dtype=bool)
    positives = int(labels.sum())
    if positives == 0:
        return None
    order = np.argsort(-np.asarray(scores), kind="mergesort")
    ordered = labels[order].astype(np.float64)
    precision = np.cumsum(ordered) / np.arange(1, len(ordered) + 1)
    return _as_float(np.sum(precision * ordered) / positives)


def _bootstrap_mean_ci(
    values: Sequence[float],
    repetitions: int,
    generator: np.random.Generator,
) -> Optional[List[float]]:
    finite = np.asarray([value for value in values if np.isfinite(value)], dtype=np.float64)
    if repetitions <= 0 or len(finite) < 2:
        return None
    sampled = generator.choice(finite, size=(repetitions, len(finite)), replace=True).mean(axis=1)
    return [float(np.quantile(sampled, 0.025)), float(np.quantile(sampled, 0.975))]


def _regression_report(
    prediction: np.ndarray,
    target: np.ndarray,
    groups: np.ndarray,
    total_std: Optional[np.ndarray] = None,
    bootstrap_repetitions: int = 0,
    seed: int = 7,
) -> Dict[str, Any]:
    prediction = np.asarray(prediction, dtype=np.float64)
    target = np.asarray(target, dtype=np.float64)
    groups = np.asarray(groups, dtype=str)
    valid = np.isfinite(prediction) & np.isfinite(target)
    if total_std is not None:
        total_std = np.asarray(total_std, dtype=np.float64)
        valid &= np.isfinite(total_std) & (total_std > 0)
    prediction, target, groups = prediction[valid], target[valid], groups[valid]
    if total_std is not None:
        total_std = total_std[valid]
    if len(target) == 0:
        return {"count": 0}

    error = prediction - target
    pearson = _correlation(prediction, target)
    spearman = _spearman(prediction, target)
    target_variance = float(np.sum((target - target.mean()) ** 2))
    r2 = None if target_variance == 0 else 1.0 - float(np.sum(error ** 2)) / target_variance
    covariance = float(np.mean((prediction - prediction.mean()) * (target - target.mean())))
    ccc_denominator = float(
        np.var(prediction) + np.var(target) + (prediction.mean() - target.mean()) ** 2
    )
    ccc = None if ccc_denominator == 0 else 2.0 * covariance / ccc_denominator
    target_iqr = float(np.quantile(target, 0.75) - np.quantile(target, 0.25))
    report: Dict[str, Any] = {
        "count": int(len(target)),
        "system_count": int(len(np.unique(groups))),
        "mae": float(np.mean(np.abs(error))),
        "rmse": float(np.sqrt(np.mean(error ** 2))),
        "bias": float(np.mean(error)),
        "pearson": pearson,
        "spearman": spearman,
        "r2": _as_float(r2) if r2 is not None else None,
        "concordance_correlation": _as_float(ccc) if ccc is not None else None,
        "nmae_by_target_iqr": None
        if target_iqr == 0
        else float(np.mean(np.abs(error))) / target_iqr,
    }

    group_mae: List[float] = []
    group_rmse: List[float] = []
    group_pearson: List[float] = []
    group_spearman: List[float] = []
    for group in np.unique(groups):
        selected = groups == group
        local_error = error[selected]
        group_mae.append(float(np.mean(np.abs(local_error))))
        group_rmse.append(float(np.sqrt(np.mean(local_error ** 2))))
        local_pearson = _correlation(prediction[selected], target[selected])
        local_spearman = _spearman(prediction[selected], target[selected])
        if local_pearson is not None:
            group_pearson.append(local_pearson)
        if local_spearman is not None:
            group_spearman.append(local_spearman)
    report["macro_by_system"] = {
        "mae": float(np.mean(group_mae)),
        "rmse": float(np.mean(group_rmse)),
        "pearson": _as_float(np.mean(group_pearson)) if group_pearson else None,
        "spearman": _as_float(np.mean(group_spearman)) if group_spearman else None,
        "correlation_system_count": len(group_pearson),
    }
    if bootstrap_repetitions > 0:
        generator = np.random.default_rng(seed)
        report["macro_by_system"]["bootstrap_95ci"] = {
            "mae": _bootstrap_mean_ci(group_mae, bootstrap_repetitions, generator),
            "rmse": _bootstrap_mean_ci(group_rmse, bootstrap_repetitions, generator),
            "pearson": _bootstrap_mean_ci(group_pearson, bootstrap_repetitions, generator),
            "spearman": _bootstrap_mean_ci(group_spearman, bootstrap_repetitions, generator),
        }

    if total_std is not None:
        z_values = {"50": 0.67448975, "90": 1.64485363, "95": 1.95996398}
        variance = np.maximum(total_std ** 2, 1e-12)
        report["uncertainty"] = {
            "gaussian_nll": float(np.mean(0.5 * (error ** 2 / variance + np.log(variance)))),
            "mean_predicted_std": float(np.mean(total_std)),
            "abs_error_std_spearman": _spearman(np.abs(error), total_std),
            "interval_coverage": {
                level: float(np.mean(np.abs(error) <= multiplier * total_std))
                for level, multiplier in z_values.items()
            },
        }
    return report


def _retrieval_report(
    prediction: np.ndarray,
    target: np.ndarray,
    groups: np.ndarray,
    fraction: float = 0.20,
    use_magnitude: bool = False,
) -> Dict[str, Any]:
    prediction = np.abs(prediction) if use_magnitude else np.asarray(prediction)
    target = np.abs(target) if use_magnitude else np.asarray(target)
    groups = np.asarray(groups, dtype=str)
    labels: List[np.ndarray] = []
    scores: List[np.ndarray] = []
    recalls: List[float] = []
    group_aurocs: List[float] = []
    group_aps: List[float] = []
    for group in np.unique(groups):
        selected = groups == group
        local_target = target[selected]
        local_prediction = prediction[selected]
        if len(local_target) < 5:
            continue
        top_count = max(1, int(math.ceil(len(local_target) * fraction)))
        true_top = np.argsort(local_target)[-top_count:]
        predicted_top = np.argsort(local_prediction)[-top_count:]
        recalls.append(float(len(set(true_top) & set(predicted_top)) / top_count))
        local_labels = np.zeros(len(local_target), dtype=bool)
        local_labels[true_top] = True
        labels.append(local_labels)
        scores.append(local_prediction)
        local_auc = _binary_auroc(local_prediction, local_labels)
        local_ap = _average_precision(local_prediction, local_labels)
        if local_auc is not None:
            group_aurocs.append(local_auc)
        if local_ap is not None:
            group_aps.append(local_ap)
    if not labels:
        return {"fraction": fraction, "system_count": 0}
    pooled_labels = np.concatenate(labels)
    pooled_scores = np.concatenate(scores)
    return {
        "fraction": fraction,
        "random_expected_precision": fraction,
        "system_count": len(labels),
        "macro_top_fraction_recall": float(np.mean(recalls)),
        "macro_auroc": float(np.mean(group_aurocs)),
        "macro_average_precision": float(np.mean(group_aps)),
        "pooled_auroc": _binary_auroc(pooled_scores, pooled_labels),
        "pooled_average_precision": _average_precision(pooled_scores, pooled_labels),
    }


def _direction_report(
    prediction: np.ndarray,
    target: np.ndarray,
    threshold: float,
) -> Dict[str, Any]:
    if len(target) == 0:
        return {"dead_zone_threshold": threshold, "active_count": 0, "count": 0}
    prediction_class = np.where(prediction > threshold, 1, np.where(prediction < -threshold, -1, 0))
    target_class = np.where(target > threshold, 1, np.where(target < -threshold, -1, 0))
    active = target_class != 0
    class_f1: Dict[str, Optional[float]] = {}
    f1_values: List[float] = []
    for label, name in [(-1, "decrease"), (0, "neutral"), (1, "increase")]:
        true_positive = int(np.sum((prediction_class == label) & (target_class == label)))
        false_positive = int(np.sum((prediction_class == label) & (target_class != label)))
        false_negative = int(np.sum((prediction_class != label) & (target_class == label)))
        denominator = 2 * true_positive + false_positive + false_negative
        value = None if denominator == 0 else 2.0 * true_positive / denominator
        class_f1[name] = value
        if value is not None:
            f1_values.append(value)
    return {
        "dead_zone_threshold": threshold,
        "active_count": int(active.sum()),
        "active_sign_accuracy": None if not np.any(active) else float(
            np.mean(prediction_class[active] == target_class[active])
        ),
        "three_class_accuracy": float(np.mean(prediction_class == target_class)),
        "three_class_macro_f1": None if not f1_values else float(np.mean(f1_values)),
        "class_f1": class_f1,
        "target_class_counts": {
            "decrease": int(np.sum(target_class == -1)),
            "neutral": int(np.sum(target_class == 0)),
            "increase": int(np.sum(target_class == 1)),
        },
    }


def _concatenate(values: Sequence[np.ndarray], dtype: Any = np.float64) -> np.ndarray:
    if not values:
        return np.asarray([], dtype=dtype)
    return np.concatenate(values).astype(dtype, copy=False)


def _stratified_reports(
    prediction: np.ndarray,
    target: np.ndarray,
    groups: np.ndarray,
    total_std: np.ndarray,
    strata: np.ndarray,
    bootstrap_repetitions: int,
    seed: int,
) -> Dict[str, Any]:
    output: Dict[str, Any] = {}
    for value in sorted(np.unique(strata).tolist()):
        selected = strata == value
        output[str(value)] = _regression_report(
            prediction[selected], target[selected], groups[selected], total_std[selected],
            bootstrap_repetitions, seed,
        )
    return output


def _select_split(
    dataset: Dataset,
    config: ExperimentConfig,
    split: str,
    limit: Optional[int],
) -> Tuple[Dataset, Dict[str, Any], Subset]:
    train_set, val_set = group_split(dataset, config.train.val_fraction, config.seed)
    if split == "all":
        selected: Dataset = dataset
        role = "mixed_train_and_model_selection_validation"
    elif split == "train":
        selected = train_set
        role = "training"
    elif split == "val":
        selected = val_set
        role = "model_selection_validation"
    else:
        raise ValueError("split must be all, train, or val")
    total_before_limit = len(selected)
    if limit is not None:
        selected = Subset(selected, list(range(min(limit, len(selected)))))
    return selected, {
        "name": split,
        "role": role,
        "samples_before_limit": total_before_limit,
        "samples_evaluated": len(selected),
        "train_samples": len(train_set),
        "validation_samples": len(val_set),
        "warning": (
            "This trainer has no untouched test split. Validation selected the best checkpoint "
            "and is therefore optimistic; 'all' also includes training samples."
        ),
    }, train_set


def _training_label_means(
    dataset: StateDataset,
    train_set: Optional[Subset] = None,
) -> Dict[str, float]:
    sums = np.zeros(len(METRIC_NAMES), dtype=np.float64)
    counts = np.zeros(len(METRIC_NAMES), dtype=np.int64)
    indices = range(len(dataset)) if train_set is None else train_set.indices
    for index in indices:
        record = dataset.records[index]
        if not record.label_path:
            continue
        values, _, mask = load_label_file(record.label_path)
        for metric_index in range(len(METRIC_NAMES)):
            valid = mask[:, metric_index] & np.isfinite(values[:, metric_index])
            sums[metric_index] += float(values[valid, metric_index].sum())
            counts[metric_index] += int(valid.sum())
    return {
        name: float(sums[index] / counts[index])
        for index, name in enumerate(METRIC_NAMES)
        if counts[index] > 0
    }


def _training_amino_acid_means(
    dataset: StateDataset,
    metric_name: str,
) -> Dict[int, float]:
    """Fit amino-acid priors on training labels only."""

    metric_index = METRIC_NAMES.index(metric_name)
    sums: Dict[int, float] = {}
    counts: Dict[int, int] = {}
    for record in dataset.records:
        if not record.label_path:
            continue
        structure = parse_pdb(record.structure_path)
        values, _, mask = load_label_file(record.label_path)
        selected = (
            mask[:, metric_index]
            & np.isfinite(values[:, metric_index])
            & structure.frame_mask
            & target_mask(structure, record.target_chains)
        )
        for amino_acid in np.unique(structure.aatype[selected]):
            amino_acid = int(amino_acid)
            local = selected & (structure.aatype == amino_acid)
            sums[amino_acid] = sums.get(amino_acid, 0.0) + float(values[local, metric_index].sum())
            counts[amino_acid] = counts.get(amino_acid, 0) + int(local.sum())
    return {
        amino_acid: sums[amino_acid] / counts[amino_acid]
        for amino_acid in sums
        if counts[amino_acid] > 0
    }


def _amino_acid_name(index: int) -> str:
    return RESTYPES[index] if 0 <= index < len(RESTYPES) else "X"


def _amino_acid_corrected_report(
    prediction: np.ndarray,
    target: np.ndarray,
    groups: np.ndarray,
    total_std: np.ndarray,
    amino_acid_types: np.ndarray,
    training_priors: Mapping[int, float],
    fallback_prior: float,
    bootstrap_repetitions: int,
    seed: int,
) -> Dict[str, Any]:
    prior = np.asarray(
        [training_priors.get(int(index), fallback_prior) for index in amino_acid_types],
        dtype=np.float64,
    )
    residual_prediction = prediction - prior
    residual_target = target - prior
    type_names = np.asarray([_amino_acid_name(int(index)) for index in amino_acid_types])
    return {
        "definition": (
            "prediction and target minus the training-fold mean H-chi for that amino-acid type"
        ),
        "training_priors": {
            _amino_acid_name(index): float(value)
            for index, value in sorted(training_priors.items())
        },
        "residual_overall": _regression_report(
            residual_prediction,
            residual_target,
            groups,
            total_std,
            bootstrap_repetitions,
            seed,
        ),
        "within_amino_acid_type": _stratified_reports(
            prediction,
            target,
            groups,
            total_std,
            type_names,
            bootstrap_repetitions,
            seed,
        ),
    }


def _ca_interface_mask(record: PairRecord, cutoff: float = 8.0) -> np.ndarray:
    """Approximate interface residues by target/partner C-alpha distance."""

    structure = parse_pdb(record.holo.structure_path)
    target_chains = set(
        record.metadata.get("target_holo_chains", record.holo.target_chains)
    )
    target = np.asarray([chain in target_chains for chain in structure.chain_ids])
    partner = ~target
    ca_present = structure.backbone_mask[:, 1]
    target_indices = np.flatnonzero(target & ca_present)
    partner_indices = np.flatnonzero(partner & ca_present)
    interface = np.zeros(structure.num_residues, dtype=bool)
    if len(target_indices) and len(partner_indices):
        target_ca = structure.backbone[target_indices, 1]
        partner_ca = structure.backbone[partner_indices, 1]
        cutoff2 = cutoff ** 2
        for start in range(0, len(target_ca), 256):
            distance2 = np.sum(
                (target_ca[start:start + 256, None, :] - partner_ca[None, :, :]) ** 2,
                axis=-1,
            )
            interface[target_indices[start:start + 256]] = np.any(distance2 <= cutoff2, axis=1)
    return interface[np.asarray(record.holo_indices, dtype=np.int64)]


def _load_checkpoint_model(
    config: ExperimentConfig,
    checkpoint_path: str,
) -> Tuple[FlexPPI, ModelConfig, torch.device, Mapping[str, Any]]:
    device = resolve_device(config.device)
    checkpoint = torch.load(checkpoint_path, map_location="cpu", weights_only=False)
    checkpoint_model = checkpoint.get("config", {}).get("model")
    model_config = ModelConfig(**checkpoint_model) if checkpoint_model else config.model
    unknown_metrics = sorted(set(model_config.metrics) - set(METRIC_NAMES))
    if unknown_metrics:
        raise ValueError("Checkpoint uses unknown label metrics: {}".format(unknown_metrics))
    model = FlexPPI(model_config)
    model.load_state_dict(checkpoint["model"])
    model.to(device).eval()
    return model, model_config, device, checkpoint


def _state_evaluation(
    config: ExperimentConfig,
    model: FlexPPI,
    metric_names: Sequence[str],
    device: torch.device,
    split: str,
    limit: Optional[int],
    bootstrap_repetitions: int,
) -> Dict[str, Any]:
    if config.data.state_manifest is None:
        raise ValueError("state evaluation requires data.state_manifest")
    train_dataset = StateDataset(
        config.data.state_manifest,
        config.data.max_residues,
        config.data.path_mappings,
    )
    if config.data.validation_state_manifest:
        validation_dataset = StateDataset(
            config.data.validation_state_manifest,
            config.data.max_residues,
            config.data.path_mappings,
        )
        test_dataset = (
            StateDataset(
                config.data.test_state_manifest,
                config.data.max_residues,
                config.data.path_mappings,
            )
            if config.data.test_state_manifest
            else None
        )
        if split == "train":
            selected: Dataset = train_dataset
            role = "training"
        elif split == "val":
            selected = validation_dataset
            role = "model_selection_validation"
        elif split == "test":
            if test_dataset is None:
                raise ValueError("split=test requires data.test_state_manifest")
            selected = test_dataset
            role = "untouched_internal_test"
        elif split == "all":
            datasets: List[Dataset] = [train_dataset, validation_dataset]
            if test_dataset is not None:
                datasets.append(test_dataset)
            selected = ConcatDataset(datasets)
            role = "mixed_training_validation_and_test"
        else:
            raise ValueError("split must be all, train, val, or test")
        total_before_limit = len(selected)
        if limit is not None:
            selected = Subset(selected, list(range(min(limit, len(selected)))))
        split_info = {
            "name": split,
            "role": role,
            "samples_before_limit": total_before_limit,
            "samples_evaluated": len(selected),
            "train_samples": len(train_dataset),
            "validation_samples": len(validation_dataset),
            "test_samples": 0 if test_dataset is None else len(test_dataset),
            "warning": (
                "Validation selects checkpoints and is optimistic; the explicit test "
                "manifest is untouched during training. 'all' mixes statistical roles."
            ),
        }
        training_means = _training_label_means(train_dataset)
    else:
        selected, split_info, train_set = _select_split(
            train_dataset, config, split, limit
        )
        training_means = _training_label_means(train_dataset, train_set)
    loader = DataLoader(
        selected,
        batch_size=config.data.batch_size,
        num_workers=config.data.num_workers,
        collate_fn=collate_states,
    )
    accumulated: Dict[str, Dict[str, List[np.ndarray]]] = {
        name: {
            key: []
            for key in (
                "prediction",
                "target",
                "total_std",
                "group",
                "dataset",
                "source_group",
                "amino_acid_type",
            )
        }
        for name in metric_names
    }
    with torch.no_grad():
        for batch in loader:
            batch = move_to_device(batch, device)
            output = model.forward_state(batch)
            for output_index, name in enumerate(metric_names):
                label_index = METRIC_NAMES.index(name)
                for batch_index, sample_id in enumerate(batch["sample_id"]):
                    valid = (
                        batch["label_mask"][batch_index, :, label_index]
                        & batch["target_mask"][batch_index]
                        & output["mask"][batch_index]
                    )
                    if not torch.any(valid):
                        continue
                    count = int(valid.sum().item())
                    bucket = accumulated[name]
                    bucket["prediction"].append(
                        output["mean"][batch_index, :, output_index][valid].cpu().numpy()
                    )
                    bucket["target"].append(
                        batch["label_values"][batch_index, :, label_index][valid].cpu().numpy()
                    )
                    total_variance = (
                        torch.exp(output["log_variance"][batch_index, :, output_index])
                        + batch["label_variances"][batch_index, :, label_index]
                    )
                    bucket["total_std"].append(torch.sqrt(total_variance[valid]).cpu().numpy())
                    bucket["group"].append(np.repeat(str(sample_id), count))
                    bucket["dataset"].append(np.repeat(str(batch["dataset"][batch_index]), count))
                    bucket["source_group"].append(
                        np.repeat(str(batch["source_group"][batch_index]), count)
                    )
                    bucket["amino_acid_type"].append(
                        batch["aatype"][batch_index][valid].cpu().numpy()
                    )

    metrics: Dict[str, Any] = {}
    for name, bucket in accumulated.items():
        prediction = _concatenate(bucket["prediction"])
        target = _concatenate(bucket["target"])
        total_std = _concatenate(bucket["total_std"])
        groups = _concatenate(bucket["group"], str)
        datasets = _concatenate(bucket["dataset"], str)
        source_groups = _concatenate(bucket["source_group"], str)
        amino_acid_types = _concatenate(bucket["amino_acid_type"], np.int64)
        overall = _regression_report(
            prediction, target, groups, total_std, bootstrap_repetitions, config.seed
        )
        mean_baseline = training_means.get(name)
        baseline = None
        if mean_baseline is not None and len(target):
            baseline = _regression_report(
                np.full_like(target, mean_baseline), target, groups,
                bootstrap_repetitions=bootstrap_repetitions, seed=config.seed,
            )
            if baseline.get("mae", 0) > 0 and overall.get("mae") is not None:
                overall["mae_improvement_over_training_mean_pct"] = 100.0 * (
                    baseline["mae"] - overall["mae"]
                ) / baseline["mae"]
        metrics[name] = {
            "overall": overall,
            "training_mean_baseline": {"fitted_mean": mean_baseline, "metrics": baseline},
            "high_flexibility_retrieval": _retrieval_report(prediction, target, groups),
            "by_dataset": _stratified_reports(
                prediction,
                target,
                groups,
                total_std,
                datasets,
                bootstrap_repetitions,
                config.seed,
            )
            if len(target)
            else {},
            "by_source_group": _stratified_reports(
                prediction,
                target,
                groups,
                total_std,
                source_groups,
                bootstrap_repetitions,
                config.seed,
            )
            if len(target)
            else {},
        }
        if name == "chi_entropy" and len(target):
            amino_acid_priors = _training_amino_acid_means(train_dataset, "chi_entropy")
            metrics[name]["amino_acid_corrected"] = _amino_acid_corrected_report(
                prediction,
                target,
                groups,
                total_std,
                amino_acid_types,
                amino_acid_priors,
                float(training_means.get("chi_entropy", 0.0)),
                bootstrap_repetitions,
                config.seed,
            )
    return {"task": "absolute_state", "split": split_info, "metrics": metrics}


def _pair_evaluation(
    config: ExperimentConfig,
    model: FlexPPI,
    metric_names: Sequence[str],
    device: torch.device,
    split: str,
    limit: Optional[int],
    bootstrap_repetitions: int,
    delta_thresholds: Mapping[str, float],
    interface_cutoff: float,
) -> Dict[str, Any]:
    if config.data.pair_manifest is None:
        raise ValueError("pair evaluation requires data.pair_manifest")
    dataset = PairDataset(
        config.data.pair_manifest,
        config.data.max_residues,
        config.data.path_mappings,
    )
    selected, split_info, _ = _select_split(dataset, config, split, limit)
    interface_masks = {
        record.pair_id: _ca_interface_mask(record, interface_cutoff)
        for record in dataset.records
    }
    loader = DataLoader(
        selected,
        batch_size=config.data.batch_size,
        num_workers=config.data.num_workers,
        collate_fn=collate_pairs,
    )
    keys = (
        "direct", "state_subtraction", "target", "total_std", "group", "interface",
        "sequence_identity", "coverage",
    )
    accumulated: Dict[str, Dict[str, List[np.ndarray]]] = {
        name: {key: [] for key in keys} for name in metric_names
    }
    with torch.no_grad():
        for batch in loader:
            batch = move_to_device(batch, device)
            output = model.forward_pair(batch)
            apo_values = gather_residues(batch["apo"]["label_values"], batch["apo_indices"])
            holo_values = gather_residues(batch["holo"]["label_values"], batch["holo_indices"])
            apo_variances = gather_residues(
                batch["apo"]["label_variances"], batch["apo_indices"]
            )
            holo_variances = gather_residues(
                batch["holo"]["label_variances"], batch["holo_indices"]
            )
            apo_mask = gather_residues(batch["apo"]["label_mask"], batch["apo_indices"])
            holo_mask = gather_residues(batch["holo"]["label_mask"], batch["holo_indices"])
            for output_index, name in enumerate(metric_names):
                label_index = METRIC_NAMES.index(name)
                for batch_index, pair_id in enumerate(batch["pair_id"]):
                    valid = (
                        apo_mask[batch_index, :, label_index]
                        & holo_mask[batch_index, :, label_index]
                        & output["pair_mask"][batch_index]
                    )
                    if not torch.any(valid):
                        continue
                    count = int(valid.sum().item())
                    bucket = accumulated[name]
                    target = (
                        holo_values[batch_index, :, label_index]
                        - apo_values[batch_index, :, label_index]
                    )
                    total_variance = (
                        torch.exp(output["delta_log_variance"][batch_index, :, output_index])
                        + apo_variances[batch_index, :, label_index]
                        + holo_variances[batch_index, :, label_index]
                    )
                    bucket["direct"].append(
                        output["delta_mean"][batch_index, :, output_index][valid].cpu().numpy()
                    )
                    bucket["state_subtraction"].append(
                        output["state_difference"][batch_index, :, output_index][valid]
                        .cpu()
                        .numpy()
                    )
                    bucket["target"].append(target[valid].cpu().numpy())
                    bucket["total_std"].append(torch.sqrt(total_variance[valid]).cpu().numpy())
                    bucket["group"].append(np.repeat(str(pair_id), count))
                    pair_interface = torch.from_numpy(interface_masks[str(pair_id)]).to(
                        valid.device
                    )
                    bucket["interface"].append(pair_interface[valid].cpu().numpy())
                    bucket["sequence_identity"].append(
                        np.repeat(float(batch["sequence_identity"][batch_index].cpu()), count)
                    )
                    bucket["coverage"].append(
                        np.repeat(float(batch["coverage"][batch_index].cpu()), count)
                    )

    metrics: Dict[str, Any] = {}
    for name, bucket in accumulated.items():
        direct = _concatenate(bucket["direct"])
        subtraction = _concatenate(bucket["state_subtraction"])
        target = _concatenate(bucket["target"])
        total_std = _concatenate(bucket["total_std"])
        groups = _concatenate(bucket["group"], str)
        interface = _concatenate(bucket["interface"], bool)
        methods = {
            "response_head": direct,
            "state_subtraction": subtraction,
            "zero_change": np.zeros_like(target),
        }
        comparison: Dict[str, Any] = {}
        for method, prediction in methods.items():
            comparison[method] = {
                "overall": _regression_report(
                    prediction, target, groups,
                    total_std if method == "response_head" else None,
                    bootstrap_repetitions, config.seed,
                ),
                "direction": _direction_report(
                    prediction, target, float(delta_thresholds.get(name, 0.0))
                ),
                "large_effect_retrieval": _retrieval_report(
                    prediction, target, groups, use_magnitude=True
                ),
            }
        direct_mae = comparison["response_head"]["overall"].get("mae")
        for baseline_name in ("state_subtraction", "zero_change"):
            baseline_mae = comparison[baseline_name]["overall"].get("mae")
            if direct_mae is not None and baseline_mae:
                comparison["response_head"][
                    "mae_improvement_over_{}_pct".format(baseline_name)
                ] = 100.0 * (baseline_mae - direct_mae) / baseline_mae
        by_region: Dict[str, Any] = {}
        for region_name, selected_region in {
            "interface": interface,
            "non_interface": ~interface,
        }.items():
            by_region[region_name] = {
                method: _regression_report(
                    prediction[selected_region], target[selected_region], groups[selected_region],
                    total_std[selected_region] if method == "response_head" else None,
                    bootstrap_repetitions, config.seed,
                )
                for method, prediction in methods.items()
            }
        metrics[name] = {
            "method_comparison": comparison,
            "by_interface_region": by_region,
        }
    return {
        "task": "binding_response",
        "split": split_info,
        "interface_definition": "target/partner C-alpha distance <= {:.1f} Angstrom".format(
            interface_cutoff
        ),
        "delta_dead_zone_thresholds": dict(delta_thresholds),
        "metrics": metrics,
    }


def evaluate(
    config: ExperimentConfig,
    checkpoint_path: str,
    split: str = "all",
    limit: Optional[int] = None,
    bootstrap_repetitions: int = 200,
    delta_thresholds: Optional[Mapping[str, float]] = None,
    interface_cutoff: float = 8.0,
) -> Dict[str, Any]:
    """Evaluate a checkpoint with residue-pooled and equal-protein metrics.

    ``split='val'`` reconstructs the deterministic validation split used by the
    trainer.  It is useful for diagnosis but is not an untouched test set,
    because the same split selected ``best.pt``.
    """

    model, model_config, device, checkpoint = _load_checkpoint_model(config, checkpoint_path)
    if config.stage == "pretrain":
        task_report = _state_evaluation(
            config, model, model_config.metrics, device, split, limit, bootstrap_repetitions
        )
    elif config.stage == "finetune":
        thresholds = dict(DEFAULT_DELTA_THRESHOLDS)
        thresholds.update(delta_thresholds or {})
        task_report = _pair_evaluation(
            config, model, model_config.metrics, device, split, limit,
            bootstrap_repetitions, thresholds, interface_cutoff,
        )
    else:
        raise ValueError("stage must be pretrain or finetune")
    return {
        "checkpoint": checkpoint_path,
        "checkpoint_epoch": checkpoint.get("epoch"),
        "checkpoint_train_loss": checkpoint.get("train_loss"),
        "checkpoint_validation_loss": checkpoint.get("val_loss"),
        "device": str(device),
        "model_metrics": list(model_config.metrics),
        **task_report,
    }
