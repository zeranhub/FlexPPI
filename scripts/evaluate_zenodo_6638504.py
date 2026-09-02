"""Run external FlexPPI and biological analyses on Zenodo record 6638504.

This script deliberately separates exact-PDB-overlap systems from exact-PDB-
unseen systems.  The released FlexPPI state checkpoint saw DynaRepo simulations
for five of the eight complex PDB identifiers, although its paired-response
fine-tuning manifest contains none of the sixteen apo/holo chain-side pairs.

Outputs are residue-level source data, equal-complex macro metrics, biological
distance-bin summaries, remote high-effect residues, a JSON report, and a
compact audit figure.  RMSF and chi-entropy are evaluated; SASA is excluded
because it was not computed from these trajectories.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch
from scipy.spatial import cKDTree
from scipy.stats import pearsonr, spearmanr, wilcoxon

from flexppi.constants import METRIC_NAMES, THREE_TO_ONE
from flexppi.data.dataset import collate_pairs, collate_states
from flexppi.data.labels import load_label_file
from flexppi.data.pairing import match_apo_to_holo
from flexppi.data.records import PairRecord, load_pair_records
from flexppi.data.structures import ProteinStructure, parse_pdb
from flexppi.inference import _structure_example, load_model
from flexppi.trainer import move_to_device


EXACT_PDB_OVERLAP = {"1AK4", "1AY7", "1EMV", "1GCQ", "1PVH"}
EXACT_PDB_UNSEEN = {"1BRS", "2OOB", "3SGB"}
SOURCE_STABLE_INTERFACES = {"1GCQ", "2OOB"}
SOURCE_REORIENTATION = {"1AK4", "3SGB"}
METRICS = ("rmsf", "chi_entropy")
DEAD_ZONES = {"rmsf": 0.10, "chi_entropy": 0.10}
REGION_ORDER = ("interface_0_5A", "near_5_12A", "remote_gt_12A")


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--pair-manifest",
        default="artifacts/zenodo_6638504/data/zenodo_6638504_pairs.jsonl",
    )
    parser.add_argument("--checkpoint", default="artifacts/finetune/best.pt")
    parser.add_argument("--native-evaluation", default="")
    parser.add_argument("--output-dir", default="artifacts/zenodo_6638504/evaluation")
    parser.add_argument("--device", default="auto")
    parser.add_argument("--bootstrap", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=6638504)
    return parser.parse_args()


def safe_correlation(x: Sequence[float], y: Sequence[float], method: str) -> float:
    first = np.asarray(x, dtype=float)
    second = np.asarray(y, dtype=float)
    valid = np.isfinite(first) & np.isfinite(second)
    first, second = first[valid], second[valid]
    if len(first) < 3 or np.unique(first).size < 2 or np.unique(second).size < 2:
        return math.nan
    value = spearmanr(first, second).statistic if method == "spearman" else pearsonr(first, second).statistic
    return float(value) if np.isfinite(value) else math.nan


def zscore(values: Sequence[float]) -> np.ndarray:
    array = np.asarray(values, dtype=float)
    std = float(np.nanstd(array))
    if not np.isfinite(std) or std < 1e-12:
        return np.zeros_like(array)
    return (array - float(np.nanmean(array))) / std


def bootstrap_mean_ci(
    values: Sequence[float], repetitions: int, seed: int
) -> Tuple[float, float]:
    array = np.asarray(values, dtype=float)
    array = array[np.isfinite(array)]
    if not len(array):
        return math.nan, math.nan
    generator = np.random.default_rng(seed)
    samples = generator.choice(array, size=(repetitions, len(array)), replace=True).mean(axis=1)
    return float(np.quantile(samples, 0.025)), float(np.quantile(samples, 0.975))


def top_fraction_precision(prediction: Sequence[float], target: Sequence[float], fraction: float = 0.2) -> float:
    prediction = np.asarray(prediction, dtype=float)
    target = np.asarray(target, dtype=float)
    valid = np.isfinite(prediction) & np.isfinite(target)
    prediction, target = prediction[valid], target[valid]
    if len(target) < 5:
        return math.nan
    if np.unique(prediction).size < 2:
        return math.nan
    count = max(1, int(math.ceil(fraction * len(target))))
    predicted = set(np.argsort(np.abs(prediction))[-count:].tolist())
    observed = set(np.argsort(np.abs(target))[-count:].tolist())
    return len(predicted & observed) / count


def regression_metrics(
    prediction: Sequence[float],
    target: Sequence[float],
    dead_zone: float,
    prediction_std: Optional[Sequence[float]] = None,
    target_variance: Optional[Sequence[float]] = None,
) -> Dict[str, float]:
    prediction = np.asarray(prediction, dtype=float)
    target = np.asarray(target, dtype=float)
    valid = np.isfinite(prediction) & np.isfinite(target)
    prediction, target = prediction[valid], target[valid]
    result = {
        "residue_count": int(len(target)),
        "spearman": safe_correlation(prediction, target, "spearman"),
        "pearson": safe_correlation(prediction, target, "pearson"),
        "mae": float(np.mean(np.abs(prediction - target))) if len(target) else math.nan,
        "rmse": float(np.sqrt(np.mean((prediction - target) ** 2))) if len(target) else math.nan,
        "bias": float(np.mean(prediction - target)) if len(target) else math.nan,
        "top20_precision": top_fraction_precision(prediction, target),
    }
    predicted_class = np.where(prediction > dead_zone, 1, np.where(prediction < -dead_zone, -1, 0))
    target_class = np.where(target > dead_zone, 1, np.where(target < -dead_zone, -1, 0))
    active = target_class != 0
    result["active_count"] = int(active.sum())
    result["active_sign_accuracy"] = (
        float(np.mean(predicted_class[active] == target_class[active]))
        if active.any() else math.nan
    )
    result["three_class_accuracy"] = (
        float(np.mean(predicted_class == target_class)) if len(target) else math.nan
    )
    if prediction_std is not None and target_variance is not None:
        std = np.asarray(prediction_std, dtype=float)[valid]
        variance = np.asarray(target_variance, dtype=float)[valid]
        total_std = np.sqrt(np.maximum(std * std + variance, 1e-12))
        result["interval_coverage_95"] = float(
            np.mean(np.abs(prediction - target) <= 1.96 * total_std)
        )
        result["uncertainty_error_spearman"] = safe_correlation(
            total_std, np.abs(prediction - target), "spearman"
        )
    else:
        result["interval_coverage_95"] = math.nan
        result["uncertainty_error_spearman"] = math.nan
    return result


def _pair_output(
    model: torch.nn.Module,
    apo_pdb: str,
    holo_pdb: str,
    target_chain: str,
    device: torch.device,
) -> Dict[str, np.ndarray]:
    apo_structure = parse_pdb(apo_pdb)
    holo_structure = parse_pdb(holo_pdb)
    apo_indices, holo_indices, _, identity, coverage = match_apo_to_holo(
        apo_structure, holo_structure
    )
    keep = [
        position for position, index in enumerate(holo_indices)
        if holo_structure.chain_ids[index] == target_chain
    ]
    apo_indices = [apo_indices[position] for position in keep]
    holo_indices = [holo_indices[position] for position in keep]
    example = {
        "pair_id": "zenodo_6638504",
        "apo": _structure_example(
            apo_pdb, apo_structure.unique_chains(), 300.0, 100.0, 1.0
        ),
        "holo": _structure_example(holo_pdb, [target_chain], 300.0, 100.0, 1.0),
        "apo_indices": torch.tensor(apo_indices, dtype=torch.long),
        "holo_indices": torch.tensor(holo_indices, dtype=torch.long),
        "sequence_identity": identity,
        "coverage": coverage,
        "metadata": {"target_chain": target_chain},
    }
    batch = move_to_device(collate_pairs([example]), device)
    with torch.no_grad():
        output = model.forward_pair(batch)
    mask = output["pair_mask"][0].detach().cpu().numpy().astype(bool)
    ai = batch["apo_indices"][0].detach().cpu().numpy()[mask]
    hi = batch["holo_indices"][0].detach().cpu().numpy()[mask]
    return {
        "apo_indices": ai,
        "holo_indices": hi,
        "apo_mean": output["apo"]["mean"][0, batch["apo_indices"][0]][mask].detach().cpu().numpy(),
        "holo_mean": output["holo"]["mean"][0, batch["holo_indices"][0]][mask].detach().cpu().numpy(),
        "delta_mean": output["delta_mean"][0][mask].detach().cpu().numpy(),
        "delta_std": torch.exp(0.5 * output["delta_log_variance"][0][mask]).detach().cpu().numpy(),
        "state_difference": output["state_difference"][0][mask].detach().cpu().numpy(),
    }


def _state_output(
    model: torch.nn.Module,
    pdb_path: str,
    target_chain: str,
    device: torch.device,
) -> Dict[str, np.ndarray]:
    structure = parse_pdb(pdb_path)
    batch = move_to_device(
        collate_states([_structure_example(pdb_path, [target_chain], 300.0, 100.0, 1.0)]),
        device,
    )
    with torch.no_grad():
        output = model.forward_state(batch)
    indices = torch.nonzero(
        batch["target_mask"][0] & output["mask"][0], as_tuple=False
    ).squeeze(-1)
    return {
        "indices": indices.detach().cpu().numpy(),
        "mean": output["mean"][0, indices].detach().cpu().numpy(),
        "std": torch.exp(0.5 * output["log_variance"][0, indices]).detach().cpu().numpy(),
        "author_residue_id": np.asarray(
            [structure.author_residue_ids[int(index)] for index in indices.detach().cpu().tolist()]
        ),
    }


def elastic_network_scores(path: str, cutoff: float = 10.0) -> Dict[str, np.ndarray]:
    structure = parse_pdb(path)
    coordinates = np.asarray(structure.backbone[:, 1], dtype=float)
    distance = np.linalg.norm(coordinates[:, None, :] - coordinates[None, :, :], axis=-1)
    adjacency = (distance > 0.0) & (distance <= cutoff)
    degree = adjacency.sum(axis=1).astype(float)
    kirchhoff = np.diag(degree) - adjacency.astype(float)
    covariance = np.linalg.pinv(kirchhoff, hermitian=True, rcond=1e-8)
    gnm = np.maximum(np.diag(covariance), 0.0)
    inverse_contact = 1.0 / np.maximum(degree, 1.0)
    terminal = np.zeros(structure.num_residues, dtype=float)
    for chain in structure.unique_chains():
        indices = structure.indices_for_chain(chain)
        position = np.arange(len(indices), dtype=float)
        terminal[indices] = -np.minimum(position, position[::-1])
    return {"gnm": gnm, "inverse_contact": inverse_contact, "terminal": terminal}


def residue_partner_distances(path: str, target_chain: str) -> Dict[str, float]:
    target: Dict[str, List[np.ndarray]] = {}
    partner: List[np.ndarray] = []
    with Path(path).open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith(("ATOM  ", "HETATM")) or len(line) < 54:
                continue
            residue_name = line[17:20].strip().upper()
            if residue_name not in THREE_TO_ONE:
                continue
            atom_name = line[12:16].strip().upper().lstrip("0123456789")
            if atom_name.startswith("H"):
                continue
            coordinate = np.asarray(
                [float(line[30:38]), float(line[38:46]), float(line[46:54])], dtype=float
            )
            chain = line[21:22].strip() or "_"
            author_id = line[22:26].strip() + line[26:27].strip()
            if chain == target_chain:
                target.setdefault(author_id, []).append(coordinate)
            else:
                partner.append(coordinate)
    if not target or not partner:
        raise ValueError("Could not identify target and partner heavy atoms in {}".format(path))
    tree = cKDTree(np.asarray(partner))
    return {
        author_id: float(np.min(tree.query(np.asarray(coordinates), k=1)[0]))
        for author_id, coordinates in target.items()
    }


def overlap_class(pdb_id: str) -> str:
    return "exact_pdb_overlap" if pdb_id in EXACT_PDB_OVERLAP else "exact_pdb_unseen"


def region_from_distance(distance: float) -> str:
    if distance <= 5.0:
        return REGION_ORDER[0]
    if distance <= 12.0:
        return REGION_ORDER[1]
    return REGION_ORDER[2]


def build_residue_table(
    records: Sequence[PairRecord],
    model: torch.nn.Module,
    device: torch.device,
) -> pd.DataFrame:
    rows: List[Dict[str, Any]] = []
    network_cache: Dict[str, Dict[str, np.ndarray]] = {}
    for position, record in enumerate(records, start=1):
        pdb_id = str(record.metadata["pdb_id"]).upper()
        side = str(record.metadata["side"])
        target_only_path = str(record.metadata["holo_target_only_pdb"])
        full = _pair_output(
            model, record.apo.structure_path, record.holo.structure_path, side, device
        )
        removed = _pair_output(
            model, record.apo.structure_path, target_only_path, side, device
        )
        removed_state = _state_output(model, target_only_path, side, device)
        if not np.array_equal(full["apo_indices"], np.asarray(record.apo_indices)):
            raise ValueError("{} apo model/manifest mapping differs".format(record.pair_id))
        if not np.array_equal(full["holo_indices"], np.asarray(record.holo_indices)):
            raise ValueError("{} holo model/manifest mapping differs".format(record.pair_id))
        if len(removed_state["indices"]) != len(full["holo_indices"]):
            raise ValueError("{} partner-removed residue count differs".format(record.pair_id))

        apo_values, apo_variances, apo_mask = load_label_file(str(record.apo.label_path))
        holo_values, holo_variances, holo_mask = load_label_file(str(record.holo.label_path))
        ai, hi = full["apo_indices"], full["holo_indices"]
        apo_values, apo_variances, apo_mask = apo_values[ai], apo_variances[ai], apo_mask[ai]
        holo_values, holo_variances, holo_mask = holo_values[hi], holo_variances[hi], holo_mask[hi]
        structure = parse_pdb(record.holo.structure_path)
        distances = residue_partner_distances(record.holo.structure_path, side)

        if record.apo.structure_path not in network_cache:
            network_cache[record.apo.structure_path] = elastic_network_scores(record.apo.structure_path)
        if record.holo.structure_path not in network_cache:
            network_cache[record.holo.structure_path] = elastic_network_scores(record.holo.structure_path)
        apo_network = network_cache[record.apo.structure_path]
        holo_network = network_cache[record.holo.structure_path]
        apo_gnm = apo_network["gnm"][ai]
        holo_gnm = holo_network["gnm"][hi]
        apo_contact = apo_network["inverse_contact"][ai]
        holo_contact = holo_network["inverse_contact"][hi]
        apo_terminal = apo_network["terminal"][ai]
        holo_terminal = holo_network["terminal"][hi]

        for local_index, (apo_index, holo_index) in enumerate(zip(ai, hi)):
            author_id = structure.author_residue_ids[int(holo_index)]
            distance = distances[author_id]
            row: Dict[str, Any] = {
                "pdb_id": pdb_id,
                "side": side,
                "pair_id": record.pair_id,
                "overlap_class": overlap_class(pdb_id),
                "apo_residue_index": int(apo_index),
                "holo_residue_index": int(holo_index),
                "author_residue_id": author_id,
                "amino_acid": structure.sequence[int(holo_index)],
                "partner_distance_angstrom": distance,
                "distance_region": region_from_distance(distance),
                "is_interface_5A": distance <= 5.0,
                "apo_gnm": float(apo_gnm[local_index]),
                "holo_gnm": float(holo_gnm[local_index]),
                "delta_gnm_z": float(zscore(holo_gnm)[local_index] - zscore(apo_gnm)[local_index]),
                "apo_inverse_contact": float(apo_contact[local_index]),
                "holo_inverse_contact": float(holo_contact[local_index]),
                "delta_inverse_contact_z": float(
                    zscore(holo_contact)[local_index] - zscore(apo_contact)[local_index]
                ),
                "apo_terminal_proximity": float(apo_terminal[local_index]),
                "holo_terminal_proximity": float(holo_terminal[local_index]),
                "interface_rigidification": float(-math.exp(-distance / 5.0)),
            }
            for metric_index, metric in enumerate(METRIC_NAMES):
                valid = bool(apo_mask[local_index, metric_index] and holo_mask[local_index, metric_index])
                row["true_apo_" + metric] = (
                    float(apo_values[local_index, metric_index]) if apo_mask[local_index, metric_index] else math.nan
                )
                row["true_holo_" + metric] = (
                    float(holo_values[local_index, metric_index]) if holo_mask[local_index, metric_index] else math.nan
                )
                row["true_delta_" + metric] = (
                    float(holo_values[local_index, metric_index] - apo_values[local_index, metric_index])
                    if valid else math.nan
                )
                row["true_delta_variance_" + metric] = (
                    float(holo_variances[local_index, metric_index] + apo_variances[local_index, metric_index])
                    if valid else math.nan
                )
                row["flexppi_apo_" + metric] = float(full["apo_mean"][local_index, metric_index])
                row["flexppi_holo_full_" + metric] = float(full["holo_mean"][local_index, metric_index])
                row["flexppi_holo_partner_removed_" + metric] = float(
                    removed_state["mean"][local_index, metric_index]
                )
                row["flexppi_delta_response_" + metric] = float(full["delta_mean"][local_index, metric_index])
                row["flexppi_delta_response_std_" + metric] = float(full["delta_std"][local_index, metric_index])
                row["flexppi_delta_state_subtraction_" + metric] = float(
                    full["state_difference"][local_index, metric_index]
                )
                row["flexppi_delta_partner_removed_response_" + metric] = float(
                    removed["delta_mean"][local_index, metric_index]
                )
                row["flexppi_holo_context_effect_" + metric] = float(
                    full["holo_mean"][local_index, metric_index]
                    - removed_state["mean"][local_index, metric_index]
                )
            rows.append(row)
        print("[{}/{}] {} side {}: {} residues".format(
            position, len(records), pdb_id, side, len(ai)
        ), flush=True)
    return pd.DataFrame(rows)


def per_side_metrics(frame: pd.DataFrame) -> pd.DataFrame:
    rows: List[Dict[str, Any]] = []
    for (pdb_id, side), group in frame.groupby(["pdb_id", "side"], sort=True):
        category = str(group["overlap_class"].iloc[0])
        for metric in METRICS:
            dead_zone = DEAD_ZONES[metric]
            response_target = group["true_delta_" + metric].to_numpy(dtype=float)
            response_methods: Mapping[str, np.ndarray] = {
                "flexppi_response_head": group["flexppi_delta_response_" + metric].to_numpy(dtype=float),
                "flexppi_state_subtraction": group["flexppi_delta_state_subtraction_" + metric].to_numpy(dtype=float),
                "flexppi_partner_removed_response": group[
                    "flexppi_delta_partner_removed_response_" + metric
                ].to_numpy(dtype=float),
                "flexppi_holo_context_effect": group[
                    "flexppi_holo_context_effect_" + metric
                ].to_numpy(dtype=float),
                "interface_rigidification": group["interface_rigidification"].to_numpy(dtype=float),
                "zero_change": np.zeros(len(group), dtype=float),
            }
            if metric == "rmsf":
                response_methods = {
                    **response_methods,
                    "gnm_state_subtraction_z": group["delta_gnm_z"].to_numpy(dtype=float),
                    "inverse_contact_subtraction_z": group[
                        "delta_inverse_contact_z"
                    ].to_numpy(dtype=float),
                }
            for method, prediction in response_methods.items():
                result = regression_metrics(
                    prediction,
                    response_target,
                    dead_zone,
                    prediction_std=(
                        group["flexppi_delta_response_std_" + metric].to_numpy(dtype=float)
                        if method == "flexppi_response_head" else None
                    ),
                    target_variance=(
                        group["true_delta_variance_" + metric].to_numpy(dtype=float)
                        if method == "flexppi_response_head" else None
                    ),
                )
                rows.append({
                    "pdb_id": pdb_id,
                    "side": side,
                    "overlap_class": category,
                    "task": "binding_response",
                    "metric": metric,
                    "method": method,
                    **result,
                })

            state_specs: List[Tuple[str, str, Mapping[str, np.ndarray]]] = [
                (
                    "apo_state", "true_apo_" + metric,
                    {"flexppi": group["flexppi_apo_" + metric].to_numpy(dtype=float)},
                ),
                (
                    "holo_state", "true_holo_" + metric,
                    {
                        "flexppi_full_complex": group["flexppi_holo_full_" + metric].to_numpy(dtype=float),
                        "flexppi_partner_removed": group[
                            "flexppi_holo_partner_removed_" + metric
                        ].to_numpy(dtype=float),
                    },
                ),
            ]
            if metric == "rmsf":
                state_specs[0][2].update({
                    "gnm": group["apo_gnm"].to_numpy(dtype=float),
                    "inverse_contact": group["apo_inverse_contact"].to_numpy(dtype=float),
                    "terminal_proximity": group["apo_terminal_proximity"].to_numpy(dtype=float),
                })
                state_specs[1][2].update({
                    "gnm": group["holo_gnm"].to_numpy(dtype=float),
                    "inverse_contact": group["holo_inverse_contact"].to_numpy(dtype=float),
                    "terminal_proximity": group["holo_terminal_proximity"].to_numpy(dtype=float),
                })
            for task, target_column, methods in state_specs:
                target = group[target_column].to_numpy(dtype=float)
                for method, prediction in methods.items():
                    result = regression_metrics(prediction, target, dead_zone=0.0)
                    rows.append({
                        "pdb_id": pdb_id,
                        "side": side,
                        "overlap_class": category,
                        "task": task,
                        "metric": metric,
                        "method": method,
                        **result,
                    })
    return pd.DataFrame(rows)


def complex_macro(frame: pd.DataFrame) -> pd.DataFrame:
    numeric = [
        "residue_count", "spearman", "pearson", "mae", "rmse", "bias",
        "top20_precision", "active_count", "active_sign_accuracy",
        "three_class_accuracy", "interval_coverage_95", "uncertainty_error_spearman",
    ]
    result = (
        frame.groupby(
            ["pdb_id", "overlap_class", "task", "metric", "method"], sort=True
        )[numeric]
        .mean()
        .reset_index()
    )
    result["source_interface_behavior"] = np.where(
        result["pdb_id"].isin(SOURCE_STABLE_INTERFACES),
        "single_stable_interface_cluster",
        "multiple_interface_substates",
    )
    result["source_reorientation_behavior"] = np.where(
        result["pdb_id"].isin(SOURCE_REORIENTATION),
        "reported_partner_reorientation",
        "no_major_reorientation_reported",
    )
    return result


def summary_metrics(
    per_complex: pd.DataFrame, repetitions: int, seed: int
) -> pd.DataFrame:
    rows: List[Dict[str, Any]] = []
    subsets = {
        "all_8": set(EXACT_PDB_OVERLAP | EXACT_PDB_UNSEEN),
        "exact_pdb_unseen_3": set(EXACT_PDB_UNSEEN),
        "exact_pdb_overlap_5": set(EXACT_PDB_OVERLAP),
    }
    value_columns = [
        "spearman", "pearson", "mae", "rmse", "top20_precision",
        "active_sign_accuracy", "three_class_accuracy", "interval_coverage_95",
        "uncertainty_error_spearman",
    ]
    group_columns = ["task", "metric", "method"]
    for subset_index, (subset, pdb_ids) in enumerate(subsets.items()):
        selected = per_complex[per_complex["pdb_id"].isin(pdb_ids)]
        for group_key, group in selected.groupby(group_columns, sort=True):
            task, metric, method = group_key
            for value_index, column in enumerate(value_columns):
                values = group[column].to_numpy(dtype=float)
                values = values[np.isfinite(values)]
                low, high = bootstrap_mean_ci(
                    values, repetitions, seed + subset_index * 1000 + value_index
                )
                rows.append({
                    "subset": subset,
                    "complex_count": int(group["pdb_id"].nunique()),
                    "task": task,
                    "metric": metric,
                    "method": method,
                    "measure": column,
                    "macro_mean": float(np.mean(values)) if len(values) else math.nan,
                    "macro_median": float(np.median(values)) if len(values) else math.nan,
                    "bootstrap_95ci_low": low,
                    "bootstrap_95ci_high": high,
                })
    return pd.DataFrame(rows)


def binding_effects(frame: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    side_rows: List[Dict[str, Any]] = []
    for (pdb_id, side), group in frame.groupby(["pdb_id", "side"], sort=True):
        for metric in METRICS:
            threshold = DEAD_ZONES[metric]
            for region in REGION_ORDER:
                selected = group[group["distance_region"].eq(region)]
                values = selected["true_delta_" + metric].to_numpy(dtype=float)
                values = values[np.isfinite(values)]
                side_rows.append({
                    "pdb_id": pdb_id,
                    "side": side,
                    "overlap_class": overlap_class(pdb_id),
                    "metric": metric,
                    "distance_region": region,
                    "residue_count": int(len(values)),
                    "mean_delta": float(np.mean(values)) if len(values) else math.nan,
                    "mean_absolute_delta": float(np.mean(np.abs(values))) if len(values) else math.nan,
                    "fraction_decreased": float(np.mean(values < -threshold)) if len(values) else math.nan,
                    "fraction_increased": float(np.mean(values > threshold)) if len(values) else math.nan,
                })
        coupled = group[["true_delta_rmsf", "true_delta_chi_entropy"]].dropna()
        side_rows.append({
            "pdb_id": pdb_id,
            "side": side,
            "overlap_class": overlap_class(pdb_id),
            "metric": "rmsf_chi_delta_coupling",
            "distance_region": "all",
            "residue_count": int(len(coupled)),
            "mean_delta": safe_correlation(
                coupled["true_delta_rmsf"], coupled["true_delta_chi_entropy"], "spearman"
            ),
            "mean_absolute_delta": math.nan,
            "fraction_decreased": math.nan,
            "fraction_increased": math.nan,
        })
    side = pd.DataFrame(side_rows)
    numeric = [
        "residue_count", "mean_delta", "mean_absolute_delta",
        "fraction_decreased", "fraction_increased",
    ]
    complex_frame = (
        side.groupby(
            ["pdb_id", "overlap_class", "metric", "distance_region"], sort=True
        )[numeric]
        .mean()
        .reset_index()
    )
    return side, complex_frame


def biological_summary(
    per_complex: pd.DataFrame, repetitions: int, seed: int
) -> pd.DataFrame:
    rows: List[Dict[str, Any]] = []
    subsets = {
        "all_8": EXACT_PDB_OVERLAP | EXACT_PDB_UNSEEN,
        "exact_pdb_unseen_3": EXACT_PDB_UNSEEN,
        "exact_pdb_overlap_5": EXACT_PDB_OVERLAP,
    }
    for subset_index, (subset, pdb_ids) in enumerate(subsets.items()):
        selected = per_complex[per_complex["pdb_id"].isin(pdb_ids)]
        for key, group in selected.groupby(["metric", "distance_region"], sort=True):
            metric, region = key
            for value_index, column in enumerate(
                ["mean_delta", "mean_absolute_delta", "fraction_decreased", "fraction_increased"]
            ):
                values = group[column].to_numpy(dtype=float)
                values = values[np.isfinite(values)]
                low, high = bootstrap_mean_ci(
                    values, repetitions, seed + subset_index * 1000 + value_index
                )
                rows.append({
                    "subset": subset,
                    "complex_count": int(group["pdb_id"].nunique()),
                    "metric": metric,
                    "distance_region": region,
                    "measure": column,
                    "macro_mean": float(np.mean(values)) if len(values) else math.nan,
                    "bootstrap_95ci_low": low,
                    "bootstrap_95ci_high": high,
                })
    return pd.DataFrame(rows)


def top_remote_residues(frame: pd.DataFrame, count: int = 5) -> pd.DataFrame:
    rows: List[pd.DataFrame] = []
    remote = frame[frame["partner_distance_angstrom"] > 12.0]
    for (pdb_id, side), group in remote.groupby(["pdb_id", "side"], sort=True):
        for metric in METRICS:
            column = "true_delta_" + metric
            selected = group.dropna(subset=[column]).copy()
            selected["absolute_effect"] = selected[column].abs()
            selected = selected.nlargest(count, "absolute_effect").copy()
            selected["metric"] = metric
            selected["true_delta"] = selected[column]
            selected["flexppi_delta"] = selected["flexppi_delta_response_" + metric]
            selected["window_standard_deviation"] = np.sqrt(
                selected["true_delta_variance_" + metric]
            )
            selected["effect_to_window_sd"] = (
                selected["true_delta"] / selected["window_standard_deviation"].replace(0.0, np.nan)
            )
            selected["sign_correct"] = (
                np.sign(selected["true_delta"]) == np.sign(selected["flexppi_delta"])
            )
            rows.append(selected[[
                "pdb_id", "side", "overlap_class", "metric", "author_residue_id",
                "amino_acid", "partner_distance_angstrom", "true_delta",
                "flexppi_delta", "absolute_effect", "window_standard_deviation",
                "effect_to_window_sd", "sign_correct",
            ]])
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def paired_comparisons(per_complex: pd.DataFrame) -> pd.DataFrame:
    specifications = [
        ("binding_response", "rmsf", "mae", "flexppi_response_head", "zero_change"),
        ("binding_response", "rmsf", "spearman", "flexppi_response_head", "flexppi_state_subtraction"),
        ("binding_response", "rmsf", "spearman", "flexppi_response_head", "interface_rigidification"),
        ("binding_response", "chi_entropy", "mae", "flexppi_response_head", "zero_change"),
        ("binding_response", "chi_entropy", "spearman", "flexppi_response_head", "interface_rigidification"),
        ("holo_state", "rmsf", "spearman", "flexppi_full_complex", "flexppi_partner_removed"),
        ("holo_state", "chi_entropy", "spearman", "flexppi_full_complex", "flexppi_partner_removed"),
    ]
    rows: List[Dict[str, Any]] = []
    subsets = {
        "all_8": EXACT_PDB_OVERLAP | EXACT_PDB_UNSEEN,
        "exact_pdb_unseen_3": EXACT_PDB_UNSEEN,
        "exact_pdb_overlap_5": EXACT_PDB_OVERLAP,
    }
    for subset, pdb_ids in subsets.items():
        for task, metric, measure, first, second in specifications:
            selected = per_complex[
                per_complex["pdb_id"].isin(pdb_ids)
                & per_complex["task"].eq(task) & per_complex["metric"].eq(metric)
                & per_complex["method"].isin([first, second])
            ]
            pivot = selected.pivot(index="pdb_id", columns="method", values=measure).dropna()
            if first not in pivot or second not in pivot:
                continue
            difference = pivot[first] - pivot[second]
            if measure in {"mae", "rmse"}:
                difference = -difference
                orientation = "positive means first has lower error"
            else:
                orientation = "positive means first is better"
            try:
                p_value = float(wilcoxon(difference).pvalue) if len(difference) >= 2 else math.nan
            except ValueError:
                p_value = 1.0
            rows.append({
                "subset": subset,
                "task": task,
                "metric": metric,
                "measure": measure,
                "first_method": first,
                "second_method": second,
                "complex_count": int(len(difference)),
                "mean_first_advantage": float(difference.mean()),
                "median_first_advantage": float(difference.median()),
                "win_rate": float((difference > 0).mean()),
                "wilcoxon_two_sided_p": p_value,
                "orientation": orientation,
            })
    return pd.DataFrame(rows)


def spatial_prior_analysis(residue: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    rows: List[Dict[str, Any]] = []
    for (pdb_id, side), group in residue.groupby(["pdb_id", "side"], sort=True):
        for metric in METRICS:
            rows.append({
                "pdb_id": pdb_id,
                "side": side,
                "overlap_class": overlap_class(pdb_id),
                "metric": metric,
                "truth_vs_partner_distance_spearman": safe_correlation(
                    group["true_delta_" + metric], group["partner_distance_angstrom"], "spearman"
                ),
                "flexppi_vs_partner_distance_spearman": safe_correlation(
                    group["flexppi_delta_response_" + metric],
                    group["partner_distance_angstrom"], "spearman"
                ),
                "flexppi_vs_truth_spearman": safe_correlation(
                    group["flexppi_delta_response_" + metric], group["true_delta_" + metric], "spearman"
                ),
            })
    side = pd.DataFrame(rows)
    numeric = [
        "truth_vs_partner_distance_spearman",
        "flexppi_vs_partner_distance_spearman",
        "flexppi_vs_truth_spearman",
    ]
    complex_frame = (
        side.groupby(["pdb_id", "overlap_class", "metric"], sort=True)[numeric]
        .mean()
        .reset_index()
    )
    return side, complex_frame


def source_paper_strata(per_complex: pd.DataFrame) -> pd.DataFrame:
    selected = per_complex[
        per_complex["task"].eq("binding_response")
        & per_complex["metric"].eq("rmsf")
        & per_complex["method"].isin(["flexppi_response_head", "interface_rigidification"])
    ][[
        "pdb_id", "overlap_class", "source_interface_behavior",
        "source_reorientation_behavior", "method", "spearman", "mae",
    ]]
    return selected.sort_values(["source_interface_behavior", "pdb_id", "method"])


def _weighted_linear_fit(
    features: np.ndarray,
    target: np.ndarray,
    groups: Sequence[str],
    ridge: float = 1e-6,
) -> np.ndarray:
    features = np.asarray(features, dtype=float)
    target = np.asarray(target, dtype=float)
    groups = np.asarray(groups, dtype=str)
    valid = np.isfinite(target) & np.all(np.isfinite(features), axis=1)
    features, target, groups = features[valid], target[valid], groups[valid]
    design = np.column_stack([features, np.ones(len(features), dtype=float)])
    counts = pd.Series(groups).value_counts().to_dict()
    weights = np.asarray([1.0 / counts[group] for group in groups], dtype=float)
    weights /= weights.mean()
    penalty = np.eye(design.shape[1], dtype=float) * ridge
    penalty[-1, -1] = 0.0
    return np.linalg.solve(
        design.T @ (weights[:, None] * design) + penalty,
        design.T @ (weights * target),
    )


def leave_one_complex_out_calibration(
    residue: pd.DataFrame,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    feature_sets = {
        "affine_flexppi_response": ["flexppi_delta_response_{metric}"],
        "affine_state_subtraction": ["flexppi_delta_state_subtraction_{metric}"],
        "calibrated_interface_distance": ["interface_rigidification"],
        "hybrid_flexppi_plus_distance": [
            "flexppi_delta_response_{metric}", "interface_rigidification"
        ],
    }
    rows: List[Dict[str, Any]] = []
    for metric in METRICS:
        target_column = "true_delta_" + metric
        for held_out in sorted(residue["pdb_id"].unique()):
            train = residue[~residue["pdb_id"].eq(held_out)].dropna(subset=[target_column])
            test = residue[residue["pdb_id"].eq(held_out)].dropna(subset=[target_column])
            for model_name, templates in feature_sets.items():
                columns = [template.format(metric=metric) for template in templates]
                coefficients = _weighted_linear_fit(
                    train[columns].to_numpy(dtype=float),
                    train[target_column].to_numpy(dtype=float),
                    train["pdb_id"].astype(str).to_numpy(),
                )
                prediction = (
                    test[columns].to_numpy(dtype=float) @ coefficients[:-1] + coefficients[-1]
                )
                side_values: List[Dict[str, float]] = []
                for side, indices in test.groupby("side").groups.items():
                    local_positions = test.index.get_indexer(indices)
                    result = regression_metrics(
                        prediction[local_positions],
                        test.loc[indices, target_column].to_numpy(dtype=float),
                        DEAD_ZONES[metric],
                    )
                    side_values.append(result)
                numeric = side_values[0].keys()
                averaged: Dict[str, float] = {}
                for key in numeric:
                    values = np.asarray([value[key] for value in side_values], dtype=float)
                    finite = values[np.isfinite(values)]
                    averaged[key] = float(np.mean(finite)) if len(finite) else math.nan
                rows.append({
                    "held_out_pdb": held_out,
                    "overlap_class": overlap_class(held_out),
                    "metric": metric,
                    "calibration_model": model_name,
                    "training_pdb_count": 7,
                    "feature_columns": ";".join(columns),
                    "coefficients": ";".join("{:.9g}".format(value) for value in coefficients),
                    **averaged,
                })
    per_complex = pd.DataFrame(rows)
    summary_rows: List[Dict[str, Any]] = []
    subsets = {
        "all_8": EXACT_PDB_OVERLAP | EXACT_PDB_UNSEEN,
        "exact_pdb_unseen_3": EXACT_PDB_UNSEEN,
        "exact_pdb_overlap_5": EXACT_PDB_OVERLAP,
    }
    for subset, pdb_ids in subsets.items():
        selected = per_complex[per_complex["held_out_pdb"].isin(pdb_ids)]
        for key, group in selected.groupby(["metric", "calibration_model"], sort=True):
            metric, model_name = key
            for measure in [
                "spearman", "mae", "rmse", "bias", "top20_precision",
                "active_sign_accuracy", "three_class_accuracy",
            ]:
                values = group[measure].to_numpy(dtype=float)
                values = values[np.isfinite(values)]
                summary_rows.append({
                    "subset": subset,
                    "complex_count": int(group["held_out_pdb"].nunique()),
                    "metric": metric,
                    "calibration_model": model_name,
                    "measure": measure,
                    "macro_mean": float(np.mean(values)) if len(values) else math.nan,
                    "macro_median": float(np.median(values)) if len(values) else math.nan,
                })
    return per_complex, pd.DataFrame(summary_rows)


def make_figure(
    path_prefix: Path,
    residue: pd.DataFrame,
    per_complex: pd.DataFrame,
    effects: pd.DataFrame,
) -> None:
    plt.rcParams.update({"font.size": 8.0, "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(2, 2, figsize=(8.2, 6.3), constrained_layout=True)

    ax = axes[0, 0]
    methods = ["flexppi_response_head", "flexppi_state_subtraction", "interface_rigidification"]
    labels = ["FlexPPI response", "State subtraction", "Interface heuristic"]
    colors = ["#2457A6", "#C65D21", "#777777"]
    for index, (method, label, color) in enumerate(zip(methods, labels, colors)):
        values = per_complex[
            per_complex["task"].eq("binding_response")
            & per_complex["metric"].eq("rmsf")
            & per_complex["method"].eq(method)
        ]["spearman"].dropna().to_numpy()
        ax.scatter(np.full(len(values), index) + np.linspace(-0.10, 0.10, len(values)), values,
                   s=18, alpha=0.75, color=color)
        ax.plot([index - 0.22, index + 0.22], [np.mean(values)] * 2, color="black", lw=1.4)
    ax.axhline(0.0, color="0.5", linestyle="--", lw=0.8)
    ax.set_xticks(range(len(labels)), labels, rotation=15, ha="right")
    ax.set_ylabel("Per-complex Spearman vs MD ΔRMSF")
    ax.set_title("a  Binding-response ranking")

    ax = axes[0, 1]
    rmsf_effect = effects[
        effects["metric"].eq("rmsf") & effects["distance_region"].isin(REGION_ORDER)
    ]
    for index, region in enumerate(REGION_ORDER):
        values = rmsf_effect[rmsf_effect["distance_region"].eq(region)]["mean_delta"].dropna()
        ax.scatter(np.full(len(values), index) + np.linspace(-0.10, 0.10, len(values)), values,
                   s=18, alpha=0.75, color="#2D8C61")
        ax.plot([index - 0.22, index + 0.22], [values.mean()] * 2, color="black", lw=1.4)
    ax.axhline(0.0, color="0.5", linestyle="--", lw=0.8)
    ax.set_xticks(range(3), ["0–5 Å", "5–12 Å", ">12 Å"])
    ax.set_ylabel("MD ΔRMSF (Å)")
    ax.set_title("b  Binding effect by partner distance")

    ax = axes[1, 0]
    selected = per_complex[
        per_complex["task"].eq("holo_state")
        & per_complex["metric"].eq("rmsf")
        & per_complex["method"].isin(["flexppi_full_complex", "flexppi_partner_removed"])
    ]
    pivot = selected.pivot(index="pdb_id", columns="method", values="spearman").dropna()
    ax.scatter(pivot["flexppi_partner_removed"], pivot["flexppi_full_complex"],
               c=["#C65D21" if pdb in EXACT_PDB_UNSEEN else "#2457A6" for pdb in pivot.index], s=32)
    limits = [float(np.nanmin(pivot.to_numpy()) - 0.05), float(np.nanmax(pivot.to_numpy()) + 0.05)]
    ax.plot(limits, limits, color="0.4", linestyle="--", lw=0.9)
    for pdb_id, row in pivot.iterrows():
        ax.annotate(pdb_id, (row["flexppi_partner_removed"], row["flexppi_full_complex"]), fontsize=6)
    ax.set(xlim=limits, ylim=limits, xlabel="Partner removed Spearman", ylabel="Full complex Spearman")
    ax.set_title("c  Does partner context improve bound RMSF?")

    ax = axes[1, 1]
    valid = residue[["true_delta_rmsf", "flexppi_delta_response_rmsf", "is_interface_5A"]].dropna()
    ax.scatter(valid.loc[~valid["is_interface_5A"], "true_delta_rmsf"],
               valid.loc[~valid["is_interface_5A"], "flexppi_delta_response_rmsf"],
               s=5, alpha=0.18, color="#777777", label="non-interface")
    ax.scatter(valid.loc[valid["is_interface_5A"], "true_delta_rmsf"],
               valid.loc[valid["is_interface_5A"], "flexppi_delta_response_rmsf"],
               s=7, alpha=0.35, color="#C65D21", label="interface")
    ax.axhline(0.0, color="0.6", lw=0.7)
    ax.axvline(0.0, color="0.6", lw=0.7)
    ax.set(xlabel="MD ΔRMSF (Å)", ylabel="FlexPPI ΔRMSF (Å)")
    ax.legend(frameon=False, fontsize=7)
    ax.set_title("d  Residue-level calibration")

    path_prefix.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(path_prefix) + ".png", dpi=300, bbox_inches="tight")
    fig.savefig(str(path_prefix) + ".pdf", bbox_inches="tight")
    plt.close(fig)


def _json_ready(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): _json_ready(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_json_ready(item) for item in value]
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        return None if not np.isfinite(value) else float(value)
    return value


def selected_summary(summary: pd.DataFrame) -> List[Dict[str, Any]]:
    wanted = summary[
        summary["measure"].isin(["spearman", "mae", "top20_precision", "active_sign_accuracy"])
        & (
            summary["method"].isin([
                "flexppi_response_head", "flexppi_state_subtraction", "zero_change",
                "flexppi_full_complex", "flexppi_partner_removed", "gnm", "inverse_contact",
            ])
        )
    ]
    return wanted.to_dict(orient="records")


def main() -> None:
    args = parse_arguments()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    records = load_pair_records(args.pair_manifest)
    if {str(record.metadata["pdb_id"]).upper() for record in records} != EXACT_PDB_OVERLAP | EXACT_PDB_UNSEEN:
        raise ValueError("Pair manifest does not contain the expected eight PDB systems")
    model, model_config, device = load_model(args.checkpoint, args.device)
    residue = build_residue_table(records, model, device)
    per_side = per_side_metrics(residue)
    per_complex = complex_macro(per_side)
    summary = summary_metrics(per_complex, args.bootstrap, args.seed)
    effect_side, effect_complex = binding_effects(residue)
    effect_summary = biological_summary(effect_complex, args.bootstrap, args.seed + 10000)
    remote = top_remote_residues(residue)
    comparisons = paired_comparisons(per_complex)
    spatial_side, spatial_complex = spatial_prior_analysis(residue)
    source_strata = source_paper_strata(per_complex)
    loco_calibration, loco_summary = leave_one_complex_out_calibration(residue)

    outputs: Mapping[str, pd.DataFrame] = {
        "residue_predictions.csv": residue,
        "per_side_metrics.csv": per_side,
        "per_complex_metrics.csv": per_complex,
        "summary_metrics.csv": summary,
        "binding_effects_per_side.csv": effect_side,
        "binding_effects_per_complex.csv": effect_complex,
        "biological_summary.csv": effect_summary,
        "top_remote_residues.csv": remote,
        "paired_comparisons.csv": comparisons,
        "spatial_prior_per_side.csv": spatial_side,
        "spatial_prior_per_complex.csv": spatial_complex,
        "source_paper_strata.csv": source_strata,
        "loco_calibration_per_complex.csv": loco_calibration,
        "loco_calibration_summary.csv": loco_summary,
    }
    for name, value in outputs.items():
        value.to_csv(output_dir / name, index=False)
    make_figure(output_dir / "zenodo_6638504_audit", residue, per_complex, effect_complex)

    native = None
    if args.native_evaluation and Path(args.native_evaluation).is_file():
        with Path(args.native_evaluation).open("r", encoding="utf-8") as handle:
            native = json.load(handle)
    report = {
        "experiment": "FlexPPI external apo/holo evaluation on Zenodo 6638504",
        "dataset": {
            "pdb_systems": sorted(EXACT_PDB_OVERLAP | EXACT_PDB_UNSEEN),
            "complex_count": 8,
            "paired_chain_sides": 16,
            "residue_pairs": int(len(residue)),
            "production_trajectory_total_microseconds": 8.4,
            "label_windows_ns": 100,
            "label_window_count": 84,
            "rmsf_alignment": "each protein chain fitted independently",
            "interface_definition": "minimum inter-partner heavy-atom distance <= 5 A in starting holo structure",
        },
        "leakage_audit": {
            "exact_pdb_overlap_with_state_pretraining": sorted(EXACT_PDB_OVERLAP),
            "exact_pdb_unseen": sorted(EXACT_PDB_UNSEEN),
            "exact_pair_overlap_with_pair_finetuning": [],
            "primary_generalization_subset": "exact_pdb_unseen",
            "homology_cluster_exclusion_performed": False,
        },
        "model": {
            "checkpoint": str(Path(args.checkpoint).resolve()),
            "device": str(device),
            "metrics": list(model_config.metrics),
            "parameter_updates_during_evaluation": 0,
        },
        "statistics": {
            "independent_unit": "PDB complex after averaging its A and B sides",
            "bootstrap_repetitions": args.bootstrap,
            "paired_test": "two-sided Wilcoxon signed-rank across the eight complexes",
            "dead_zones": DEAD_ZONES,
        },
        "selected_summary": selected_summary(summary),
        "paired_comparisons": comparisons.to_dict(orient="records"),
        "spatial_prior_per_complex": spatial_complex.to_dict(orient="records"),
        "source_paper_strata": source_strata.to_dict(orient="records"),
        "leave_one_complex_out_calibration": loco_summary.to_dict(orient="records"),
        "biological_summary": effect_summary.to_dict(orient="records"),
        "native_flexppi_evaluation": native,
    }
    with (output_dir / "report.json").open("w", encoding="utf-8") as handle:
        json.dump(_json_ready(report), handle, indent=2, ensure_ascii=False)
    print(json.dumps({
        "device": str(device),
        "complexes": int(residue["pdb_id"].nunique()),
        "paired_sides": int(residue[["pdb_id", "side"]].drop_duplicates().shape[0]),
        "residue_pairs": int(len(residue)),
        "output_dir": str(output_dir.resolve()),
    }, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
