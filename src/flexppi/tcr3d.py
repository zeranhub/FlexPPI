"""Independent, proxy-based TCR3D evaluation from crystallographic structures.

TCR3D does not contain MD trajectories.  Consequently, crystallographic
B-factors and CDR annotations are reported as *orthogonal proxy tests* and are
never presented as ground-truth MD RMSF or binding-induced flexibility change.
"""

from __future__ import annotations

import csv
import json
import math
import re
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import torch

from flexppi.constants import THREE_TO_ONE
from flexppi.data.dataset import collate_states
from flexppi.data.structures import parse_pdb
from flexppi.evaluation import (
    _average_precision,
    _binary_auroc,
    _regression_report,
)
from flexppi.inference import _structure_example, load_model
from flexppi.trainer import move_to_device


CHAIN_FILE_PATTERN = re.compile(
    r"^(?P<pdb>[0-9A-Za-z]{4})_(?P<chain>[^_]+)_(?P<tcr>[AB])_aho\.pdb$",
    flags=re.IGNORECASE,
)


def _split_chains(value: str) -> List[str]:
    return [part.strip() for part in value.split(",") if part.strip()]


def _load_chain_metadata(path: str) -> Dict[Tuple[str, str, str], Dict[str, str]]:
    output: Dict[Tuple[str, str, str], Dict[str, str]] = {}
    with Path(path).open("r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            pdb_id = str(row.get("pdb_id", "")).lower()
            tcr_type = str(row.get("tcr_type", "")).lower()
            chains = set()
            for column in ("cdr1_chains", "cdr2_chains", "cdr3_chains"):
                chains.update(_split_chains(str(row.get(column, ""))))
            for chain in chains:
                output[(pdb_id, tcr_type, chain)] = dict(row)
    return output


def _training_pdb_ids(manifest_paths: Sequence[str]) -> set:
    identifiers = set()
    for manifest_path in manifest_paths:
        with Path(manifest_path).open("r", encoding="utf-8") as handle:
            for line in handle:
                if not line.strip():
                    continue
                row = json.loads(line)
                states = [row]
                if "apo" in row and "holo" in row:
                    states = [row["apo"], row["holo"]]
                    if row.get("metadata", {}).get("pdb_id"):
                        identifiers.add(str(row["metadata"]["pdb_id"]).lower())
                for state in states:
                    for pdb_id in state.get("metadata", {}).get("pdb_ids", []):
                        identifiers.add(str(pdb_id).lower())
                    if str(state.get("dataset", "")).lower() == "atlas":
                        identifiers.add(str(state.get("sample_id", ""))[:4].lower())
    return identifiers


def _ca_bfactors(path: str) -> np.ndarray:
    order: List[Tuple[str, str, str]] = []
    values: Dict[Tuple[str, str, str], float] = {}
    saw_model = False
    with Path(path).open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            record = line[:6].strip()
            if record == "MODEL":
                if saw_model:
                    break
                saw_model = True
                continue
            if record == "ENDMDL":
                break
            if record not in {"ATOM", "HETATM"} or len(line) < 66:
                continue
            if line[16:17] not in {" ", "A", "1"}:
                continue
            if line[17:20].strip().upper() not in THREE_TO_ONE:
                continue
            key = (line[21:22].strip() or "_", line[22:26].strip(), line[26:27].strip())
            if key not in values:
                order.append(key)
                values[key] = float("nan")
            if line[12:16].strip().upper() == "CA" and not np.isfinite(values[key]):
                try:
                    values[key] = float(line[60:66])
                except ValueError:
                    pass
    return np.asarray([values[key] for key in order], dtype=np.float32)


def _cdr_regions(sequence: str, metadata: Mapping[str, str]) -> np.ndarray:
    regions = np.zeros(len(sequence), dtype=np.int8)
    for region_index, column in enumerate(
        ("cdr1_sequences", "cdr2_sequences", "cdr3_sequences"), start=1
    ):
        motif = str(metadata.get(column, "")).strip().upper().replace(" ", "")
        if not motif:
            continue
        start = sequence.find(motif)
        if start >= 0:
            regions[start:start + len(motif)] = region_index
    return regions


def _classification_report(scores: np.ndarray, labels: np.ndarray, groups: np.ndarray) -> Dict[str, Any]:
    valid = np.isfinite(scores)
    scores, labels, groups = scores[valid], labels[valid].astype(bool), groups[valid]
    macro_auc: List[float] = []
    macro_ap: List[float] = []
    for group in np.unique(groups):
        selected = groups == group
        auc = _binary_auroc(scores[selected], labels[selected])
        ap = _average_precision(scores[selected], labels[selected])
        if auc is not None:
            macro_auc.append(auc)
        if ap is not None:
            macro_ap.append(ap)
    positive = scores[labels]
    negative = scores[~labels]
    return {
        "count": int(len(scores)),
        "positive_count": int(labels.sum()),
        "prevalence_and_random_ap": None if len(labels) == 0 else float(labels.mean()),
        "pooled_auroc": _binary_auroc(scores, labels),
        "pooled_average_precision": _average_precision(scores, labels),
        "macro_chain_auroc": None if not macro_auc else float(np.mean(macro_auc)),
        "macro_chain_average_precision": None if not macro_ap else float(np.mean(macro_ap)),
        "evaluable_chain_count": len(macro_auc),
        "mean_score_positive": None if len(positive) == 0 else float(positive.mean()),
        "mean_score_framework": None if len(negative) == 0 else float(negative.mean()),
        "positive_minus_framework": None if len(positive) == 0 or len(negative) == 0 else float(
            positive.mean() - negative.mean()
        ),
    }


def _within_group_zscore(values: np.ndarray, groups: np.ndarray) -> np.ndarray:
    output = np.full_like(values, np.nan, dtype=np.float64)
    for group in np.unique(groups):
        selected = (groups == group) & np.isfinite(values)
        if selected.sum() < 2:
            continue
        standard_deviation = float(np.std(values[selected]))
        if standard_deviation == 0:
            continue
        output[selected] = (values[selected] - float(np.mean(values[selected]))) / standard_deviation
    return output


def _candidate_chain_files(root: str) -> List[Path]:
    base = Path(root)
    return sorted(
        path for path in base.rglob("*_aho.pdb")
        if CHAIN_FILE_PATTERN.match(path.name) and path.parent.name in {"a_chain", "b_chain"}
    )


def evaluate_tcr3d(
    checkpoint_path: str,
    tcr3d_root: str,
    device_name: str = "auto",
    exclude_manifests: Optional[Sequence[str]] = None,
    limit: Optional[int] = None,
    dataset_code: float = 0.0,
    alternate_dataset_code: Optional[float] = 1.0,
    temperature: float = 300.0,
    window_ns: float = 100.0,
    bootstrap_repetitions: int = 200,
) -> Tuple[Dict[str, Any], List[Dict[str, Any]]]:
    root = Path(tcr3d_root).resolve()
    chain_table = root / "tcr_chain_data.tsv"
    structure_root = root / "tcr_structures"
    metadata = _load_chain_metadata(str(chain_table))
    excluded_ids = _training_pdb_ids(exclude_manifests or [])
    candidates = _candidate_chain_files(str(structure_root))
    exact_overlap_count = sum(
        bool(CHAIN_FILE_PATTERN.match(path.name).group("pdb").lower() in excluded_ids)
        for path in candidates
    )
    candidates = [
        path for path in candidates
        if CHAIN_FILE_PATTERN.match(path.name).group("pdb").lower() not in excluded_ids
    ]
    if limit is not None:
        candidates = candidates[:limit]

    model, model_config, device = load_model(checkpoint_path, device_name)
    rows: List[Dict[str, Any]] = []
    skipped_no_metadata = 0
    skipped_length_mismatch = 0
    with torch.no_grad():
        for path in candidates:
            match = CHAIN_FILE_PATTERN.match(path.name)
            if match is None:
                continue
            pdb_id = match.group("pdb").lower()
            chain_id = match.group("chain")
            chain_type = "alpha" if match.group("tcr").upper() == "A" else "beta"
            row_metadata = metadata.get((pdb_id, chain_type, chain_id))
            if row_metadata is None:
                skipped_no_metadata += 1
                continue
            structure = parse_pdb(str(path))
            bfactor = _ca_bfactors(str(path))
            if len(bfactor) != structure.num_residues:
                skipped_length_mismatch += 1
                continue
            cdr_region = _cdr_regions(structure.sequence, row_metadata)
            primary_batch = collate_states(
                [_structure_example(
                    str(path), structure.unique_chains(), temperature, window_ns, dataset_code
                )]
            )
            primary_output = model.forward_state(move_to_device(primary_batch, device))
            primary_mean = primary_output["mean"][0].cpu().numpy()
            primary_std = torch.exp(0.5 * primary_output["log_variance"][0]).cpu().numpy()
            valid_geometry = primary_output["mask"][0].cpu().numpy()
            alternate_mean = None
            if alternate_dataset_code is not None:
                alternate_batch = collate_states(
                    [_structure_example(
                        str(path), structure.unique_chains(), temperature, window_ns,
                        alternate_dataset_code,
                    )]
                )
                alternate_output = model.forward_state(move_to_device(alternate_batch, device))
                alternate_mean = alternate_output["mean"][0].cpu().numpy()
            sample_id = path.stem
            for residue_index in range(structure.num_residues):
                if not valid_geometry[residue_index]:
                    continue
                value: Dict[str, Any] = {
                    "sample_id": sample_id,
                    "pdb_id": pdb_id,
                    "chain": chain_id,
                    "tcr_type": chain_type,
                    "tcr_name": str(row_metadata.get("tcr_name", "")),
                    "complex_type": str(row_metadata.get("complex_type", "")),
                    "residue_index": residue_index,
                    "author_residue_id": structure.author_residue_ids[residue_index],
                    "amino_acid": structure.sequence[residue_index],
                    "cdr_region": int(cdr_region[residue_index]),
                    "ca_bfactor": float(bfactor[residue_index]),
                    "bfactor_rmsf_proxy": float(
                        math.sqrt(max(float(bfactor[residue_index]), 0.0) / (8.0 * math.pi ** 2))
                    ) if np.isfinite(bfactor[residue_index]) else None,
                }
                for metric_index, name in enumerate(model_config.metrics):
                    value["predicted_" + name] = float(primary_mean[residue_index, metric_index])
                    value["predicted_" + name + "_std"] = float(
                        primary_std[residue_index, metric_index]
                    )
                    if alternate_mean is not None:
                        value["alternate_condition_" + name] = float(
                            alternate_mean[residue_index, metric_index]
                        )
                rows.append(value)

    groups = np.asarray([row["sample_id"] for row in rows], dtype=str)
    cdr = np.asarray([row["cdr_region"] for row in rows], dtype=np.int8)
    report: Dict[str, Any] = {
        "task": "tcr3d_external_proxy",
        "checkpoint": checkpoint_path,
        "device": str(device),
        "scope_warning": (
            "TCR3D has no MD trajectories. B-factor ranking and CDR enrichment are independent "
            "proxy tests, not MD-RMSF or binding-response ground truth."
        ),
        "decontamination": {
            "exact_pdb_ids_loaded_from_manifests": len(excluded_ids),
            "candidate_chain_files_removed_for_exact_pdb_overlap": exact_overlap_count,
            "homology_filter": "not_performed; run MMseqs2/structure clustering before publication",
        },
        "dataset": {
            "candidate_chain_files_after_exclusion_and_limit": len(candidates),
            "evaluated_chains": int(len(np.unique(groups))) if len(groups) else 0,
            "evaluated_residues": len(rows),
            "skipped_no_metadata": skipped_no_metadata,
            "skipped_length_mismatch": skipped_length_mismatch,
        },
        "primary_condition": {
            "dataset_code": dataset_code,
            "temperature": temperature,
            "window_ns": window_ns,
        },
        "cdr_retrieval": {},
    }
    for metric_name in model_config.metrics:
        scores = np.asarray([row["predicted_" + metric_name] for row in rows], dtype=np.float64)
        metric_report: Dict[str, Any] = {
            "all_cdr_vs_framework": _classification_report(scores, cdr > 0, groups)
        }
        for region in (1, 2, 3):
            selected = (cdr == 0) | (cdr == region)
            metric_report["cdr{}_vs_framework".format(region)] = _classification_report(
                scores[selected], cdr[selected] == region, groups[selected]
            )
        report["cdr_retrieval"][metric_name] = metric_report

    if rows and "rmsf" in model_config.metrics:
        prediction = np.asarray([row["predicted_rmsf"] for row in rows], dtype=np.float64)
        proxy = np.asarray([
            np.nan if row["bfactor_rmsf_proxy"] is None else row["bfactor_rmsf_proxy"]
            for row in rows
        ])
        valid = np.isfinite(proxy) & (proxy > 0)
        prediction_z = _within_group_zscore(prediction[valid], groups[valid])
        proxy_z = _within_group_zscore(proxy[valid], groups[valid])
        z_valid = np.isfinite(prediction_z) & np.isfinite(proxy_z)
        report["bfactor_rank_proxy"] = {
            "raw_scale": _regression_report(
                prediction[valid], proxy[valid], groups[valid],
                bootstrap_repetitions=bootstrap_repetitions,
            ),
            "within_chain_zscore": _regression_report(
                prediction_z[z_valid], proxy_z[z_valid], groups[valid][z_valid],
                bootstrap_repetitions=bootstrap_repetitions,
            ),
            "interpretation": (
                "Use within-chain Spearman as the primary statistic; crystal B-factor scale is "
                "confounded by resolution, packing, refinement, and TLS treatment."
            ),
        }

    if rows and alternate_dataset_code is not None:
        sensitivity: Dict[str, Any] = {}
        for metric_name in model_config.metrics:
            primary = np.asarray([row["predicted_" + metric_name] for row in rows])
            alternate = np.asarray([row["alternate_condition_" + metric_name] for row in rows])
            sensitivity[metric_name] = {
                "mean_absolute_change": float(np.mean(np.abs(primary - alternate))),
                "mean_signed_change": float(np.mean(alternate - primary)),
                "pearson_between_conditions": _regression_report(
                    alternate, primary, groups
                ).get("pearson"),
            }
        report["dataset_condition_sensitivity"] = {
            "alternate_dataset_code": alternate_dataset_code,
            "metrics": sensitivity,
        }
    return report, rows


def write_tcr3d_outputs(
    report: Mapping[str, Any],
    rows: Iterable[Mapping[str, Any]],
    output_dir: str,
) -> Dict[str, str]:
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    report_path = output / "tcr3d_report.json"
    residue_path = output / "tcr3d_residue_predictions.csv"
    with report_path.open("w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2, ensure_ascii=False)
    rows = list(rows)
    with residue_path.open("w", encoding="utf-8", newline="") as handle:
        if rows:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    return {"report": str(report_path.resolve()), "residue_predictions": str(residue_path.resolve())}
