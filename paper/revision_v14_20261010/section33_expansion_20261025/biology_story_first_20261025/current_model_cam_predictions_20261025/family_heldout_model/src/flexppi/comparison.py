"""Paired, common-task comparison between FlexPPI and FastProtFlex.

FastProtFlex predicts a within-protein normalized flexibility score, whereas
FlexPPI predicts RMSF in Angstrom.  The primary comparison therefore uses
within-protein ranking, standardized error, and top-flexible-residue retrieval.
Absolute-scale error is reported only for FlexPPI and is never fabricated for
FastProtFlex.
"""

from __future__ import annotations

import csv
import json
import os
from pathlib import Path
import subprocess
import time
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import torch

from flexppi.config import load_config
from flexppi.baselines.fastprotflex import parse_heavy_protein_atoms
from flexppi.constants import METRIC_NAMES
from flexppi.data.dataset import StateDataset, collate_states
from flexppi.data.labels import load_label_file
from flexppi.data.structures import parse_pdb, target_mask
from flexppi.evaluation import _average_precision, _binary_auroc, _spearman
from flexppi.inference import _structure_example, load_model
from flexppi.trainer import group_split, move_to_device


def _finite_float(value: float) -> Optional[float]:
    return float(value) if np.isfinite(value) else None


def _zscore(values: np.ndarray) -> Optional[np.ndarray]:
    values = np.asarray(values, dtype=np.float64)
    if len(values) < 2:
        return None
    standard_deviation = float(np.std(values, ddof=0))
    if standard_deviation == 0 or not np.isfinite(standard_deviation):
        return None
    return (values - float(np.mean(values))) / standard_deviation


def _top_fraction_metrics(
    prediction: np.ndarray,
    target: np.ndarray,
    fraction: float = 0.20,
) -> Dict[str, Optional[float]]:
    count = len(target)
    if count < 2:
        return {"auroc": None, "average_precision": None, "precision_at_k": None}
    k = max(1, int(np.ceil(count * fraction)))
    target_top = np.zeros(count, dtype=bool)
    target_top[np.argsort(-target, kind="mergesort")[:k]] = True
    predicted_top = np.argsort(-prediction, kind="mergesort")[:k]
    return {
        "auroc": _binary_auroc(prediction, target_top),
        "average_precision": _average_precision(prediction, target_top),
        "precision_at_k": float(np.mean(target_top[predicted_top])),
    }


def _system_metrics(
    prediction: np.ndarray,
    target: np.ndarray,
    absolute_scale: bool,
) -> Dict[str, Optional[float]]:
    prediction = np.asarray(prediction, dtype=np.float64)
    target = np.asarray(target, dtype=np.float64)
    valid = np.isfinite(prediction) & np.isfinite(target)
    prediction, target = prediction[valid], target[valid]
    if len(target) < 2:
        return {}
    prediction_z, target_z = _zscore(prediction), _zscore(target)
    if prediction_z is None or target_z is None:
        return {}
    retrieval = _top_fraction_metrics(prediction, target)
    output: Dict[str, Optional[float]] = {
        "spearman": _spearman(prediction, target),
        "pearson": (
            None
            if np.std(prediction) == 0 or np.std(target) == 0
            else _finite_float(np.corrcoef(prediction, target)[0, 1])
        ),
        "z_mae": float(np.mean(np.abs(prediction_z - target_z))),
        "z_rmse": float(np.sqrt(np.mean((prediction_z - target_z) ** 2))),
        "top20_auroc": retrieval["auroc"],
        "top20_average_precision": retrieval["average_precision"],
        "top20_precision_at_k": retrieval["precision_at_k"],
    }
    if absolute_scale:
        output["absolute_mae_angstrom"] = float(np.mean(np.abs(prediction - target)))
        output["absolute_rmse_angstrom"] = float(
            np.sqrt(np.mean((prediction - target) ** 2))
        )
    return output


def _bootstrap_mean_ci(
    values: Sequence[float],
    repetitions: int,
    seed: int,
) -> Optional[List[float]]:
    array = np.asarray([value for value in values if np.isfinite(value)], dtype=np.float64)
    if repetitions <= 0 or len(array) < 2:
        return None
    generator = np.random.default_rng(seed)
    means = generator.choice(array, size=(repetitions, len(array)), replace=True).mean(axis=1)
    return [float(np.quantile(means, 0.025)), float(np.quantile(means, 0.975))]


def _method_summary(
    rows: Sequence[Mapping[str, Any]],
    method: str,
    bootstrap_repetitions: int,
    seed: int,
) -> Dict[str, Any]:
    metric_names = (
        "spearman",
        "pearson",
        "z_mae",
        "z_rmse",
        "top20_auroc",
        "top20_average_precision",
        "top20_precision_at_k",
        "absolute_mae_angstrom",
        "absolute_rmse_angstrom",
    )
    output: Dict[str, Any] = {}
    for metric in metric_names:
        values = [
            float(row[method][metric])
            for row in rows
            if metric in row.get(method, {}) and row[method][metric] is not None
        ]
        if values:
            output[metric] = {
                "macro_mean": float(np.mean(values)),
                "macro_median": float(np.median(values)),
                "protein_count": len(values),
                "protein_bootstrap_95_ci": _bootstrap_mean_ci(
                    values, bootstrap_repetitions, seed
                ),
            }
    return output


def _paired_comparison(
    rows: Sequence[Mapping[str, Any]],
    contender: str,
    reference: str,
    bootstrap_repetitions: int,
    seed: int,
) -> Dict[str, Any]:
    from scipy.stats import wilcoxon

    directions = {
        "spearman": 1.0,
        "pearson": 1.0,
        "z_mae": -1.0,
        "z_rmse": -1.0,
        "top20_auroc": 1.0,
        "top20_average_precision": 1.0,
        "top20_precision_at_k": 1.0,
    }
    output: Dict[str, Any] = {}
    for metric, direction in directions.items():
        differences = []
        for row in rows:
            contender_value = row.get(contender, {}).get(metric)
            reference_value = row.get(reference, {}).get(metric)
            if contender_value is None or reference_value is None:
                continue
            differences.append(direction * (float(contender_value) - float(reference_value)))
        if not differences:
            continue
        difference_array = np.asarray(differences, dtype=np.float64)
        nonzero = difference_array[np.abs(difference_array) > 1e-12]
        p_value = None
        if len(nonzero):
            try:
                p_value = float(wilcoxon(nonzero, alternative="two-sided").pvalue)
            except ValueError:
                pass
        output[metric] = {
            "advantage_definition": (
                "positive means {} is better than {}".format(contender, reference)
            ),
            "mean_advantage": float(np.mean(difference_array)),
            "median_advantage": float(np.median(difference_array)),
            "protein_bootstrap_95_ci": _bootstrap_mean_ci(
                difference_array.tolist(), bootstrap_repetitions, seed
            ),
            "win_rate": float(np.mean(difference_array > 1e-12)),
            "loss_rate": float(np.mean(difference_array < -1e-12)),
            "tie_rate": float(np.mean(np.abs(difference_array) <= 1e-12)),
            "paired_protein_count": len(difference_array),
            "wilcoxon_two_sided_p": p_value,
        }
    return output


def _timing_summary(values: Sequence[float]) -> Dict[str, Optional[float]]:
    array = np.asarray([value for value in values if np.isfinite(value)], dtype=np.float64)
    if len(array) == 0:
        return {"count": 0, "median_seconds": None}
    return {
        "count": int(len(array)),
        "mean_seconds": float(np.mean(array)),
        "median_seconds": float(np.median(array)),
        "p25_seconds": float(np.quantile(array, 0.25)),
        "p75_seconds": float(np.quantile(array, 0.75)),
        "p95_seconds": float(np.quantile(array, 0.95)),
    }


def _load_fast_outputs(
    prediction_path: Path,
    timing_path: Path,
) -> Tuple[
    Dict[str, np.ndarray],
    Dict[str, Dict[Tuple[str, str, str], float]],
    Dict[str, Dict[str, Any]],
]:
    predictions: Dict[str, List[Tuple[int, Tuple[str, str, str], float]]] = {}
    with prediction_path.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            predictions.setdefault(row["sample_id"], []).append(
                (
                    int(row["residue_order"]),
                    (
                        row["chain"] or "_",
                        row["residue_number"],
                        row["insertion_code"],
                    ),
                    float(row["fastprotflex_score"]),
                )
            )
    ordered = {
        sample_id: np.asarray(
            [value for _, _, value in sorted(values, key=lambda item: item[0])],
            dtype=np.float64,
        )
        for sample_id, values in predictions.items()
    }
    keyed = {
        sample_id: {
            key: value
            for _, key, value in sorted(values, key=lambda item: item[0])
        }
        for sample_id, values in predictions.items()
    }
    timings: Dict[str, Dict[str, Any]] = {}
    with timing_path.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            timings[row["sample_id"]] = dict(row)
    return ordered, keyed, timings


def _align_fast_scores(
    pdb_path: str,
    keyed_scores: Mapping[Tuple[str, str, str], float],
    expected_residue_count: int,
) -> np.ndarray:
    atoms = parse_heavy_protein_atoms(pdb_path)
    residue_keys = list(dict.fromkeys(atoms.residue_keys))
    if len(residue_keys) != expected_residue_count:
        raise ValueError(
            "standard residue-key count {} != FlexPPI residue count {}".format(
                len(residue_keys), expected_residue_count
            )
        )
    return np.asarray(
        [keyed_scores.get(key, np.nan) for key in residue_keys], dtype=np.float64
    )


def _run_fastprotflex(
    manifest_path: Path,
    prediction_path: Path,
    timing_path: Path,
    rscript_path: str,
    batch_script_path: str,
) -> None:
    environment = os.environ.copy()
    environment.setdefault(
        "R_LIBS_USER",
        str(Path(os.environ.get("LOCALAPPDATA", "")) / "R" / "win-library" / "4.6"),
    )
    subprocess.run(
        [
            rscript_path,
            batch_script_path,
            str(manifest_path),
            str(prediction_path),
            str(timing_path),
        ],
        check=True,
        env=environment,
    )


def _write_manifest(records: Sequence[Any], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample_id", "pdb_path"])
        for record in records:
            writer.writerow([record.sample_id, record.structure_path])


@torch.no_grad()
def compare_atlas_rmsf(
    config_path: str,
    pretrain_checkpoint: str,
    finetune_checkpoint: str,
    output_dir: str,
    rscript_path: str,
    fastprotflex_batch_script: str,
    limit: Optional[int] = None,
    bootstrap_repetitions: int = 2000,
    reuse_fastprotflex: bool = False,
    dataset_name: str = "atlas",
) -> Dict[str, Any]:
    """Run a paired state-validation comparison on the shared RMSF task."""

    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    dataset_name = dataset_name.lower()
    file_prefix = "".join(
        character if character.isalnum() or character in {"-", "_"} else "_"
        for character in dataset_name
    )
    manifest_path = output / (file_prefix + "_common_manifest.tsv")
    fast_prediction_path = output / (file_prefix + "_fastprotflex_predictions.csv")
    fast_timing_path = output / (file_prefix + "_fastprotflex_timings.csv")
    residue_path = output / (file_prefix + "_common_residue_predictions.csv")
    per_protein_path = output / (file_prefix + "_common_per_protein.csv")
    report_path = output / (file_prefix + "_comparison_report.json")

    config = load_config(config_path)
    dataset = StateDataset(
        config.data.state_manifest,
        config.data.max_residues,
        config.data.path_mappings,
    )
    train_set, validation_set = group_split(dataset, config.train.val_fraction, config.seed)
    validation_indices = [
        index
        for index in validation_set.indices
        if dataset.records[index].dataset.lower() == dataset_name
    ]
    if limit is not None:
        validation_indices = validation_indices[:limit]
    records = [dataset.records[index] for index in validation_indices]
    _write_manifest(records, manifest_path)

    if not (
        reuse_fastprotflex and fast_prediction_path.exists() and fast_timing_path.exists()
    ):
        _run_fastprotflex(
            manifest_path,
            fast_prediction_path,
            fast_timing_path,
            rscript_path,
            fastprotflex_batch_script,
        )
    fast_predictions, fast_keyed_predictions, fast_timings = _load_fast_outputs(
        fast_prediction_path, fast_timing_path
    )

    pretrain_model, pretrain_config, device = load_model(pretrain_checkpoint, config.device)
    finetune_model, finetune_config, finetune_device = load_model(
        finetune_checkpoint, config.device
    )
    if finetune_device != device:
        raise ValueError("Checkpoints resolved to different devices")
    pretrain_rmsf_index = pretrain_config.metrics.index("rmsf")
    finetune_rmsf_index = finetune_config.metrics.index("rmsf")
    dataset_code = 0.0 if dataset_name == "atlas" else 1.0

    # Warm up both networks once so one-time kernel setup is not assigned to a
    # biological sample. Checkpoint loading remains excluded for all methods.
    if records:
        warm_record = records[0]
        for model in (pretrain_model, finetune_model):
            warm_batch = collate_states([
                _structure_example(
                    warm_record.structure_path,
                    warm_record.target_chains,
                    warm_record.temperature,
                    warm_record.window_ns,
                    dataset_code,
                )
            ])
            model.forward_state(move_to_device(warm_batch, device))

    per_protein: List[Dict[str, Any]] = []
    residue_rows: List[Dict[str, Any]] = []
    skipped: List[Dict[str, str]] = []
    for sample_position, record in enumerate(records):
        fast_ordered = fast_predictions.get(record.sample_id)
        fast_keyed = fast_keyed_predictions.get(record.sample_id)
        timing = fast_timings.get(record.sample_id, {})
        if fast_ordered is None or fast_keyed is None or timing.get("status") != "ok":
            skipped.append({
                "sample_id": record.sample_id,
                "reason": timing.get("error", "missing FastProtFlex output"),
            })
            continue
        structure = parse_pdb(record.structure_path)
        target_values, _, label_mask = load_label_file(str(record.label_path))
        try:
            fast = _align_fast_scores(
                record.structure_path, fast_keyed, structure.num_residues
            )
        except ValueError as error:
            skipped.append({
                "sample_id": record.sample_id,
                "reason": str(error),
            })
            continue

        predictions: Dict[str, np.ndarray] = {}
        inference_seconds: Dict[str, float] = {}
        ordered_models = [
            ("flexppi_pretrain", pretrain_model, pretrain_rmsf_index),
            ("flexppi_finetune", finetune_model, finetune_rmsf_index),
        ]
        if sample_position % 2:
            ordered_models.reverse()
        for name, model, rmsf_index in ordered_models:
            started = time.perf_counter()
            batch = collate_states([
                _structure_example(
                    record.structure_path,
                    record.target_chains,
                    record.temperature,
                    record.window_ns,
                    dataset_code,
                )
            ])
            moved = move_to_device(batch, device)
            model_output = model.forward_state(moved)
            predictions[name] = model_output["mean"][0, :, rmsf_index].cpu().numpy()
            inference_seconds[name] = time.perf_counter() - started

        rmsf_label_index = METRIC_NAMES.index("rmsf")
        valid = (
            label_mask[:, rmsf_label_index]
            & structure.frame_mask
            & target_mask(structure, record.target_chains)
            & np.isfinite(target_values[:, rmsf_label_index])
            & np.isfinite(fast)
        )
        target = target_values[valid, rmsf_label_index].astype(np.float64)
        method_predictions = {
            "fastprotflex": fast[valid],
            "flexppi_pretrain": predictions["flexppi_pretrain"][valid],
            "flexppi_finetune": predictions["flexppi_finetune"][valid],
        }
        row: Dict[str, Any] = {
            "sample_id": record.sample_id,
            "residue_count": int(valid.sum()),
            "structure_residue_count": structure.num_residues,
            "fastprotflex_atom_count": int(timing["atom_count"]),
            "fastprotflex_end_to_end_seconds": float(timing["end_to_end_seconds"]),
            "flexppi_pretrain_end_to_end_seconds": inference_seconds["flexppi_pretrain"],
            "flexppi_finetune_end_to_end_seconds": inference_seconds["flexppi_finetune"],
        }
        for method, prediction in method_predictions.items():
            row[method] = _system_metrics(
                prediction, target, absolute_scale=method.startswith("flexppi")
            )
        per_protein.append(row)
        valid_indices = np.flatnonzero(valid)
        for compact_index, residue_index in enumerate(valid_indices):
            residue_rows.append({
                "sample_id": record.sample_id,
                "residue_index": int(residue_index),
                "chain": structure.chain_ids[residue_index],
                "author_residue_id": structure.author_residue_ids[residue_index],
                "target_rmsf_angstrom": float(target[compact_index]),
                "fastprotflex_score": float(method_predictions["fastprotflex"][compact_index]),
                "flexppi_pretrain_rmsf_angstrom": float(
                    method_predictions["flexppi_pretrain"][compact_index]
                ),
                "flexppi_finetune_rmsf_angstrom": float(
                    method_predictions["flexppi_finetune"][compact_index]
                ),
            })

    methods = ("fastprotflex", "flexppi_pretrain", "flexppi_finetune")
    method_summaries = {
        method: _method_summary(
            per_protein, method, bootstrap_repetitions, config.seed
        )
        for method in methods
    }
    comparisons = {
        method: _paired_comparison(
            per_protein,
            method,
            "fastprotflex",
            bootstrap_repetitions,
            config.seed,
        )
        for method in ("flexppi_pretrain", "flexppi_finetune")
    }
    length_strata: Dict[str, Any] = {}
    for stratum, lower, upper in (
        ("short_1_149", 1, 150),
        ("medium_150_299", 150, 300),
        ("long_300_plus", 300, None),
    ):
        selected = [
            row
            for row in per_protein
            if row["residue_count"] >= lower
            and (upper is None or row["residue_count"] < upper)
        ]
        length_strata[stratum] = {
            "protein_count": len(selected),
            "methods": {
                method: _method_summary(selected, method, 0, config.seed)
                for method in methods
            },
        }

    report: Dict[str, Any] = {
        "task": "shared_static_structure_to_rmsf",
        "primary_estimand": (
            "equal-protein macro performance on the intersection of evaluable C-alpha residues"
        ),
        "scale_note": (
            "FastProtFlex outputs a within-protein normalized score; direct shared metrics are "
            "rank correlation, within-protein z-score error, and top-20% retrieval. Absolute "
            "Angstrom error is valid only for FlexPPI."
        ),
        "dataset": {
            "name": "{} subset from FlexPPI deterministic validation split".format(
                dataset_name
            ),
            "combined_training_records": len(train_set),
            "combined_validation_records": len(validation_set),
            "dataset_validation_candidates": len(records),
            "paired_evaluable_proteins": len(per_protein),
            "paired_evaluable_residues": int(
                sum(row["residue_count"] for row in per_protein)
            ),
            "skipped": skipped,
            "warning": (
                (
                    "This is a same-domain model-selection diagnostic, not an untouched external "
                    "test: FlexPPI selected checkpoints using this validation split, and the "
                    "published FastProtFlex model was trained on an ATLAS subset."
                )
                if dataset_name == "atlas"
                else (
                    "This is not a symmetric generalization test: FlexPPI was pretrained on "
                    "DynaRepo records and selected using this validation partition, whereas "
                    "FastProtFlex was published as an ATLAS-trained model. Interpret differences "
                    "as the practical benefit of FlexPPI domain adaptation, not architecture-only "
                    "superiority."
                )
            ),
        },
        "methods": method_summaries,
        "paired_against_fastprotflex": comparisons,
        "by_protein_length": length_strata,
        "runtime": {
            "protocol": (
                "single-structure end-to-end prediction after package/checkpoint load; FlexPPI "
                "model order alternated by sample; CPU only"
            ),
            "fastprotflex": _timing_summary([
                row["fastprotflex_end_to_end_seconds"] for row in per_protein
            ]),
            "flexppi_pretrain": _timing_summary([
                row["flexppi_pretrain_end_to_end_seconds"] for row in per_protein
            ]),
            "flexppi_finetune": _timing_summary([
                row["flexppi_finetune_end_to_end_seconds"] for row in per_protein
            ]),
        },
        "artifacts": {
            "manifest": str(manifest_path.resolve()),
            "fastprotflex_predictions": str(fast_prediction_path.resolve()),
            "fastprotflex_timings": str(fast_timing_path.resolve()),
            "residue_predictions": str(residue_path.resolve()),
            "per_protein_metrics": str(per_protein_path.resolve()),
        },
    }

    with residue_path.open("w", encoding="utf-8", newline="") as handle:
        if residue_rows:
            writer = csv.DictWriter(handle, fieldnames=list(residue_rows[0]))
            writer.writeheader()
            writer.writerows(residue_rows)
    flat_rows: List[Dict[str, Any]] = []
    for row in per_protein:
        flat: Dict[str, Any] = {
            key: value for key, value in row.items() if not isinstance(value, dict)
        }
        for method in methods:
            for metric, value in row[method].items():
                flat[method + "__" + metric] = value
        flat_rows.append(flat)
    with per_protein_path.open("w", encoding="utf-8", newline="") as handle:
        if flat_rows:
            writer = csv.DictWriter(handle, fieldnames=list(flat_rows[0]))
            writer.writeheader()
            writer.writerows(flat_rows)
    with report_path.open("w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2, ensure_ascii=False)
    return report


def _classification_system_metrics(
    prediction: np.ndarray,
    labels: np.ndarray,
) -> Dict[str, Optional[float]]:
    prediction = np.asarray(prediction, dtype=np.float64)
    labels = np.asarray(labels, dtype=bool)
    valid = np.isfinite(prediction)
    prediction, labels = prediction[valid], labels[valid]
    if len(labels) == 0:
        return {}
    return {
        "auroc": _binary_auroc(prediction, labels),
        "average_precision": _average_precision(prediction, labels),
        "prevalence": float(np.mean(labels)),
    }


def _classification_summary(
    rows: Sequence[Mapping[str, Any]],
    method: str,
    bootstrap_repetitions: int,
    seed: int,
) -> Dict[str, Any]:
    output: Dict[str, Any] = {}
    for metric in ("auroc", "average_precision", "prevalence"):
        values = [
            float(row[method][metric])
            for row in rows
            if row.get(method, {}).get(metric) is not None
        ]
        if values:
            output[metric] = {
                "macro_mean": float(np.mean(values)),
                "macro_median": float(np.median(values)),
                "chain_count": len(values),
                "chain_bootstrap_95_ci": _bootstrap_mean_ci(
                    values, bootstrap_repetitions, seed
                ),
            }
    return output


def _paired_classification_comparison(
    rows: Sequence[Mapping[str, Any]],
    contender: str,
    reference: str,
    bootstrap_repetitions: int,
    seed: int,
) -> Dict[str, Any]:
    from scipy.stats import wilcoxon

    output: Dict[str, Any] = {}
    for metric in ("auroc", "average_precision"):
        differences = []
        for row in rows:
            left = row.get(contender, {}).get(metric)
            right = row.get(reference, {}).get(metric)
            if left is not None and right is not None:
                differences.append(float(left) - float(right))
        if not differences:
            continue
        values = np.asarray(differences, dtype=np.float64)
        nonzero = values[np.abs(values) > 1e-12]
        p_value = None
        if len(nonzero):
            try:
                p_value = float(wilcoxon(nonzero, alternative="two-sided").pvalue)
            except ValueError:
                pass
        output[metric] = {
            "advantage_definition": (
                "positive means {} is better than {}".format(contender, reference)
            ),
            "mean_advantage": float(np.mean(values)),
            "chain_bootstrap_95_ci": _bootstrap_mean_ci(
                values.tolist(), bootstrap_repetitions, seed
            ),
            "win_rate": float(np.mean(values > 1e-12)),
            "loss_rate": float(np.mean(values < -1e-12)),
            "tie_rate": float(np.mean(np.abs(values) <= 1e-12)),
            "paired_chain_count": len(values),
            "wilcoxon_two_sided_p": p_value,
        }
    return output


def compare_tcr3d_proxies(
    flexppi_residue_predictions: str,
    tcr3d_root: str,
    output_dir: str,
    rscript_path: str,
    fastprotflex_batch_script: str,
    limit: Optional[int] = None,
    bootstrap_repetitions: int = 2000,
    reuse_fastprotflex: bool = False,
) -> Dict[str, Any]:
    """Compare both models on TCR3D B-factor and CDR3 proxy tasks.

    These are deliberately described as proxies because TCR3D has neither MD
    trajectories nor experimental binding-induced RMSF changes.
    """

    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    manifest_path = output / "tcr3d_common_manifest.tsv"
    fast_prediction_path = output / "tcr3d_fastprotflex_predictions.csv"
    fast_timing_path = output / "tcr3d_fastprotflex_timings.csv"
    per_chain_path = output / "tcr3d_common_per_chain.csv"
    report_path = output / "tcr3d_comparison_report.json"

    flex_rows_by_sample: Dict[str, List[Dict[str, str]]] = {}
    with Path(flexppi_residue_predictions).open(
        "r", encoding="utf-8", newline=""
    ) as handle:
        for row in csv.DictReader(handle):
            flex_rows_by_sample.setdefault(row["sample_id"], []).append(dict(row))
    sample_ids = sorted(flex_rows_by_sample)
    if limit is not None:
        sample_ids = sample_ids[:limit]
        flex_rows_by_sample = {
            sample_id: flex_rows_by_sample[sample_id] for sample_id in sample_ids
        }

    path_by_stem = {
        path.stem: path
        for path in Path(tcr3d_root).joinpath("tcr_structures").rglob("*_aho.pdb")
    }
    missing_paths = [sample_id for sample_id in sample_ids if sample_id not in path_by_stem]
    available_ids = [sample_id for sample_id in sample_ids if sample_id in path_by_stem]
    with manifest_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample_id", "pdb_path"])
        for sample_id in available_ids:
            writer.writerow([sample_id, str(path_by_stem[sample_id])])

    if not (
        reuse_fastprotflex and fast_prediction_path.exists() and fast_timing_path.exists()
    ):
        _run_fastprotflex(
            manifest_path,
            fast_prediction_path,
            fast_timing_path,
            rscript_path,
            fastprotflex_batch_script,
        )
    fast_predictions, fast_keyed_predictions, fast_timings = _load_fast_outputs(
        fast_prediction_path, fast_timing_path
    )

    bfactor_rows: List[Dict[str, Any]] = []
    cdr3_rows: List[Dict[str, Any]] = []
    flat_rows: List[Dict[str, Any]] = []
    pooled_scores: Dict[str, List[float]] = {"fastprotflex": [], "flexppi_finetune": []}
    pooled_labels: List[bool] = []
    skipped: List[Dict[str, str]] = [
        {"sample_id": sample_id, "reason": "PDB path not found"}
        for sample_id in missing_paths
    ]
    for sample_id in available_ids:
        fast_ordered = fast_predictions.get(sample_id)
        fast_keyed = fast_keyed_predictions.get(sample_id)
        timing = fast_timings.get(sample_id, {})
        if fast_ordered is None or fast_keyed is None or timing.get("status") != "ok":
            skipped.append({
                "sample_id": sample_id,
                "reason": timing.get("error", "missing FastProtFlex output"),
            })
            continue
        rows = sorted(
            flex_rows_by_sample[sample_id], key=lambda row: int(row["residue_index"])
        )
        residue_indices = np.asarray([int(row["residue_index"]) for row in rows])
        try:
            structure_count = parse_pdb(str(path_by_stem[sample_id])).num_residues
            fast = _align_fast_scores(
                str(path_by_stem[sample_id]), fast_keyed, structure_count
            )
        except ValueError as error:
            skipped.append({
                "sample_id": sample_id,
                "reason": str(error),
            })
            continue
        if len(residue_indices) == 0 or int(residue_indices.max()) >= len(fast):
            skipped.append({
                "sample_id": sample_id,
                "reason": "FlexPPI residue indices exceed the aligned standard-residue output",
            })
            continue
        fast_selected = fast[residue_indices]
        flex_selected = np.asarray(
            [float(row["predicted_rmsf"]) for row in rows], dtype=np.float64
        )
        proxy = np.asarray(
            [
                float(row["bfactor_rmsf_proxy"])
                if row["bfactor_rmsf_proxy"] not in {"", "None"}
                else np.nan
                for row in rows
            ],
            dtype=np.float64,
        )
        proxy_valid = np.isfinite(proxy) & (proxy > 0)
        if int(proxy_valid.sum()) >= 2:
            bfactor_row: Dict[str, Any] = {
                "sample_id": sample_id,
                "residue_count": int(proxy_valid.sum()),
                "fastprotflex": _system_metrics(
                    fast_selected[proxy_valid], proxy[proxy_valid], absolute_scale=False
                ),
                "flexppi_finetune": _system_metrics(
                    flex_selected[proxy_valid], proxy[proxy_valid], absolute_scale=False
                ),
            }
            if bfactor_row["fastprotflex"] and bfactor_row["flexppi_finetune"]:
                bfactor_rows.append(bfactor_row)

        cdr_region = np.asarray([int(row["cdr_region"]) for row in rows], dtype=np.int8)
        cdr3_selected = (cdr_region == 0) | (cdr_region == 3)
        cdr3_labels = cdr_region[cdr3_selected] == 3
        if int(cdr3_labels.sum()) and int((~cdr3_labels).sum()):
            cdr3_row = {
                "sample_id": sample_id,
                "residue_count": int(cdr3_selected.sum()),
                "fastprotflex": _classification_system_metrics(
                    fast_selected[cdr3_selected], cdr3_labels
                ),
                "flexppi_finetune": _classification_system_metrics(
                    flex_selected[cdr3_selected], cdr3_labels
                ),
            }
            cdr3_rows.append(cdr3_row)
            pooled_scores["fastprotflex"].extend(fast_selected[cdr3_selected].tolist())
            pooled_scores["flexppi_finetune"].extend(flex_selected[cdr3_selected].tolist())
            pooled_labels.extend(cdr3_labels.tolist())
        flat_rows.append({
            "sample_id": sample_id,
            "residue_count": len(rows),
            "fastprotflex_end_to_end_seconds": float(timing["end_to_end_seconds"]),
            "fastprotflex_atom_count": int(timing["atom_count"]),
        })

    bfactor_methods = {
        method: _method_summary(
            bfactor_rows, method, bootstrap_repetitions, seed=7
        )
        for method in ("fastprotflex", "flexppi_finetune")
    }
    cdr3_methods: Dict[str, Any] = {}
    pooled_label_array = np.asarray(pooled_labels, dtype=bool)
    for method in ("fastprotflex", "flexppi_finetune"):
        method_summary = _classification_summary(
            cdr3_rows, method, bootstrap_repetitions, seed=7
        )
        pooled_score_array = np.asarray(pooled_scores[method], dtype=np.float64)
        method_summary["pooled"] = {
            "auroc": _binary_auroc(pooled_score_array, pooled_label_array),
            "average_precision": _average_precision(
                pooled_score_array, pooled_label_array
            ),
            "prevalence_and_random_average_precision": (
                None if len(pooled_label_array) == 0 else float(pooled_label_array.mean())
            ),
            "residue_count": int(len(pooled_label_array)),
        }
        cdr3_methods[method] = method_summary

    report: Dict[str, Any] = {
        "task": "tcr3d_external_structural_proxies",
        "scope_warning": (
            "TCR3D contains static crystallographic structures, not MD trajectories. "
            "B-factor ranking and CDR3 enrichment test transfer to orthogonal structural "
            "signals; they are not ground-truth RMSF or binding-induced flexibility change."
        ),
        "dataset": {
            "flexppi_evaluated_chain_candidates": len(sample_ids),
            "common_chain_files": len(available_ids),
            "bfactor_evaluable_chains": len(bfactor_rows),
            "cdr3_evaluable_chains": len(cdr3_rows),
            "skipped": skipped,
        },
        "bfactor_rank_proxy": {
            "primary_metric": "equal-chain macro Spearman",
            "methods": bfactor_methods,
            "paired_flexppi_against_fastprotflex": _paired_comparison(
                bfactor_rows,
                "flexppi_finetune",
                "fastprotflex",
                bootstrap_repetitions,
                seed=7,
            ),
        },
        "cdr3_vs_framework_retrieval": {
            "primary_metrics": "equal-chain macro AUROC and average precision",
            "methods": cdr3_methods,
            "paired_flexppi_against_fastprotflex": _paired_classification_comparison(
                cdr3_rows,
                "flexppi_finetune",
                "fastprotflex",
                bootstrap_repetitions,
                seed=7,
            ),
        },
        "runtime": {
            "fastprotflex": _timing_summary([
                row["fastprotflex_end_to_end_seconds"] for row in flat_rows
            ]),
            "note": (
                "The existing FlexPPI TCR3D artifact did not record per-chain timing; use the "
                "paired ATLAS timing experiment for direct runtime comparison."
            ),
        },
        "artifacts": {
            "manifest": str(manifest_path.resolve()),
            "fastprotflex_predictions": str(fast_prediction_path.resolve()),
            "fastprotflex_timings": str(fast_timing_path.resolve()),
            "per_chain_metrics": str(per_chain_path.resolve()),
        },
    }

    with per_chain_path.open("w", encoding="utf-8", newline="") as handle:
        rows_for_csv: List[Dict[str, Any]] = []
        bfactor_by_id = {row["sample_id"]: row for row in bfactor_rows}
        cdr3_by_id = {row["sample_id"]: row for row in cdr3_rows}
        for sample_id in sorted(set(bfactor_by_id) | set(cdr3_by_id)):
            row: Dict[str, Any] = {"sample_id": sample_id}
            if sample_id in bfactor_by_id:
                row["bfactor_residue_count"] = bfactor_by_id[sample_id]["residue_count"]
                for method in ("fastprotflex", "flexppi_finetune"):
                    for metric, value in bfactor_by_id[sample_id][method].items():
                        row["bfactor__{}__{}".format(method, metric)] = value
            if sample_id in cdr3_by_id:
                row["cdr3_residue_count"] = cdr3_by_id[sample_id]["residue_count"]
                for method in ("fastprotflex", "flexppi_finetune"):
                    for metric, value in cdr3_by_id[sample_id][method].items():
                        row["cdr3__{}__{}".format(method, metric)] = value
            rows_for_csv.append(row)
        if rows_for_csv:
            fieldnames = list(dict.fromkeys(
                key for row in rows_for_csv for key in row
            ))
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows_for_csv)
    with report_path.open("w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2, ensure_ascii=False)
    return report
