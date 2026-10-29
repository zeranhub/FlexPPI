"""Chain-context experiments for FlexPPI versus FastProtFlex.

The central comparison is deliberately chain-aware.  FastProtFlex scores from
separately processed chains are normalized independently, so metrics are first
computed within each protein chain and only then averaged within a complex.
This avoids treating a score of +1 in chain A as directly comparable to +1 in
chain B.
"""

from __future__ import annotations

from collections import Counter, defaultdict
import csv
import json
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence

import numpy as np
from scipy.spatial import cKDTree

from flexppi.baselines.fastprotflex import parse_heavy_protein_atoms
from flexppi.comparison import (
    _align_fast_scores,
    _load_fast_outputs,
    _method_summary,
    _paired_comparison,
    _run_fastprotflex,
    _system_metrics,
)
from flexppi.config import load_config
from flexppi.constants import METRIC_NAMES
from flexppi.data.dataset import StateDataset
from flexppi.data.labels import load_label_file
from flexppi.data.structures import parse_pdb, target_mask
from flexppi.trainer import group_split


COMPLEX_METHODS = (
    "fastprotflex_chain",
    "fastprotflex_complex",
    "flexppi_pretrain_complex",
    "flexppi_finetune_complex",
)
SINGLE_METHODS = ("fastprotflex", "flexppi_pretrain", "flexppi_finetune")


def _load_flex_residue_predictions(path: str) -> Dict[str, Dict[int, Dict[str, float]]]:
    output: Dict[str, Dict[int, Dict[str, float]]] = defaultdict(dict)
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            output[row["sample_id"]][int(row["residue_index"])] = {
                "target": float(row["target_rmsf_angstrom"]),
                "pretrain": float(row["flexppi_pretrain_rmsf_angstrom"]),
                "finetune": float(row["flexppi_finetune_rmsf_angstrom"]),
            }
    return output


def _interface_mask(structure: Any, cutoff: float) -> np.ndarray:
    mask = np.zeros(structure.num_residues, dtype=bool)
    ca_present = structure.backbone_mask[:, 1]
    indices = np.flatnonzero(ca_present)
    if len(indices) < 2:
        return mask
    coordinates = structure.backbone[indices, 1].astype(np.float64)
    for left_position, right_position in cKDTree(coordinates).query_pairs(cutoff):
        left = int(indices[left_position])
        right = int(indices[right_position])
        if structure.chain_ids[left] != structure.chain_ids[right]:
            mask[left] = True
            mask[right] = True
    return mask


def _aggregate_chains_to_complexes(
    rows: Sequence[Mapping[str, Any]],
    methods: Sequence[str],
) -> List[Dict[str, Any]]:
    grouped: Dict[str, List[Mapping[str, Any]]] = defaultdict(list)
    for row in rows:
        grouped[str(row["sample_id"])].append(row)
    output: List[Dict[str, Any]] = []
    for sample_id, chain_rows in grouped.items():
        complex_row: Dict[str, Any] = {
            "sample_id": sample_id,
            "chain_count": len(chain_rows),
            "residue_count": int(sum(int(row["residue_count"]) for row in chain_rows)),
        }
        for method in methods:
            metric_names = set().union(
                *(set(row.get(method, {})) for row in chain_rows)
            )
            complex_row[method] = {}
            for metric in metric_names:
                values = [
                    float(row[method][metric])
                    for row in chain_rows
                    if row.get(method, {}).get(metric) is not None
                ]
                if values:
                    complex_row[method][metric] = float(np.mean(values))
        output.append(complex_row)
    return sorted(output, key=lambda row: row["sample_id"])


def _summaries(
    rows: Sequence[Mapping[str, Any]],
    methods: Sequence[str],
    bootstrap_repetitions: int,
    seed: int,
) -> Dict[str, Any]:
    return {
        method: _method_summary(rows, method, bootstrap_repetitions, seed)
        for method in methods
    }


def _flat_chain_rows(
    rows: Sequence[Mapping[str, Any]],
    methods: Sequence[str],
) -> List[Dict[str, Any]]:
    output: List[Dict[str, Any]] = []
    for row in rows:
        flat = {
            key: value
            for key, value in row.items()
            if not isinstance(value, dict)
        }
        for method in methods:
            for metric, value in row.get(method, {}).items():
                flat[method + "__" + metric] = value
        output.append(flat)
    return output


def _write_rows(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        if not rows:
            return
        fieldnames = list(dict.fromkeys(key for row in rows for key in row))
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def compare_chain_context(
    config_path: str,
    atlas_report_path: str,
    dynarepo_flex_predictions: str,
    dynarepo_fast_complex_predictions: str,
    dynarepo_fast_complex_timings: str,
    output_dir: str,
    rscript_path: str,
    fastprotflex_batch_script: str,
    interface_cutoff: float = 8.0,
    bootstrap_repetitions: int = 2000,
    reuse_fastprotflex_chain: bool = False,
) -> Dict[str, Any]:
    """Evaluate single-chain and complex-context settings separately."""

    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    chain_manifest_path = output / "fastprotflex_independent_chains_manifest.tsv"
    chain_prediction_path = output / "fastprotflex_independent_chain_predictions.csv"
    chain_timing_path = output / "fastprotflex_independent_chain_timings.csv"
    per_chain_path = output / "complex_per_chain_metrics.csv"
    per_complex_path = output / "complex_per_system_metrics.csv"
    report_path = output / "chain_context_comparison_report.json"

    config = load_config(config_path)
    dataset = StateDataset(
        config.data.state_manifest,
        config.data.max_residues,
        config.data.path_mappings,
    )
    _, validation = group_split(dataset, config.train.val_fraction, config.seed)
    records = [dataset.records[index] for index in validation.indices]
    records_by_id = {record.sample_id: record for record in records}

    chain_ids_by_sample: Dict[str, List[str]] = {}
    chain_count_audit: Counter = Counter()
    source_group_audit: Counter = Counter()
    for record in records:
        structure = parse_pdb(record.structure_path)
        chains = structure.unique_chains()
        chain_ids_by_sample[record.sample_id] = chains
        structure_class = "single_chain" if len(chains) == 1 else "protein_complex"
        chain_count_audit[(record.dataset, structure_class, len(chains))] += 1
        source_group_audit[(record.dataset, structure_class, record.source_group)] += 1

    dynarepo_complex_records = [
        record for record in records
        if record.dataset.lower() == "dynarepo"
        and len(chain_ids_by_sample[record.sample_id]) >= 2
    ]
    with chain_manifest_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample_id", "pdb_path", "chain_filter"])
        for record in dynarepo_complex_records:
            for chain in chain_ids_by_sample[record.sample_id]:
                writer.writerow([
                    record.sample_id + "__chain__" + chain,
                    record.structure_path,
                    chain,
                ])
    if not (
        reuse_fastprotflex_chain
        and chain_prediction_path.exists()
        and chain_timing_path.exists()
    ):
        _run_fastprotflex(
            chain_manifest_path,
            chain_prediction_path,
            chain_timing_path,
            rscript_path,
            fastprotflex_batch_script,
        )

    _, fast_complex_keyed, fast_complex_timings = _load_fast_outputs(
        Path(dynarepo_fast_complex_predictions),
        Path(dynarepo_fast_complex_timings),
    )
    _, fast_chain_keyed, fast_chain_timings = _load_fast_outputs(
        chain_prediction_path, chain_timing_path
    )
    flex_predictions = _load_flex_residue_predictions(dynarepo_flex_predictions)

    single_rows: List[Dict[str, Any]] = []
    complex_chain_rows: List[Dict[str, Any]] = []
    interface_chain_rows: List[Dict[str, Any]] = []
    non_interface_chain_rows: List[Dict[str, Any]] = []
    skipped: List[Dict[str, str]] = []
    rmsf_index = METRIC_NAMES.index("rmsf")

    for record in [r for r in records if r.dataset.lower() == "dynarepo"]:
        structure = parse_pdb(record.structure_path)
        chains = chain_ids_by_sample[record.sample_id]
        target_values, _, label_mask = load_label_file(str(record.label_path))
        full_keyed = fast_complex_keyed.get(record.sample_id)
        sample_flex = flex_predictions.get(record.sample_id)
        if full_keyed is None or sample_flex is None:
            skipped.append({"sample_id": record.sample_id, "reason": "missing cached prediction"})
            continue
        try:
            fast_complex = _align_fast_scores(
                record.structure_path, full_keyed, structure.num_residues
            )
        except ValueError as error:
            skipped.append({"sample_id": record.sample_id, "reason": str(error)})
            continue
        flex_pretrain = np.full(structure.num_residues, np.nan, dtype=np.float64)
        flex_finetune = np.full(structure.num_residues, np.nan, dtype=np.float64)
        cached_target = np.full(structure.num_residues, np.nan, dtype=np.float64)
        for residue_index, values in sample_flex.items():
            if residue_index < structure.num_residues:
                flex_pretrain[residue_index] = values["pretrain"]
                flex_finetune[residue_index] = values["finetune"]
                cached_target[residue_index] = values["target"]
        target = target_values[:, rmsf_index].astype(np.float64)
        target_valid = (
            label_mask[:, rmsf_index]
            & structure.frame_mask
            & target_mask(structure, record.target_chains)
            & np.isfinite(target)
            & np.isfinite(cached_target)
        )
        interface = _interface_mask(structure, interface_cutoff)

        if len(chains) == 1:
            valid = (
                target_valid
                & np.isfinite(fast_complex)
                & np.isfinite(flex_pretrain)
                & np.isfinite(flex_finetune)
            )
            if int(valid.sum()) >= 2:
                single_rows.append({
                    "sample_id": record.sample_id,
                    "residue_count": int(valid.sum()),
                    "fastprotflex": _system_metrics(
                        fast_complex[valid], target[valid], absolute_scale=False
                    ),
                    "flexppi_pretrain": _system_metrics(
                        flex_pretrain[valid], target[valid], absolute_scale=True
                    ),
                    "flexppi_finetune": _system_metrics(
                        flex_finetune[valid], target[valid], absolute_scale=True
                    ),
                })
            continue

        heavy_atoms = parse_heavy_protein_atoms(record.structure_path)
        standard_keys = list(dict.fromkeys(heavy_atoms.residue_keys))
        if len(standard_keys) != structure.num_residues:
            skipped.append({
                "sample_id": record.sample_id,
                "reason": "standard residue key count differs from FlexPPI residues",
            })
            continue
        for chain in chains:
            job_id = record.sample_id + "__chain__" + chain
            chain_keyed = fast_chain_keyed.get(job_id)
            chain_timing = fast_chain_timings.get(job_id, {})
            if chain_keyed is None or chain_timing.get("status") != "ok":
                skipped.append({
                    "sample_id": record.sample_id,
                    "reason": "missing independent FastProtFlex prediction for chain " + chain,
                })
                continue
            fast_chain = np.asarray(
                [chain_keyed.get(key, np.nan) for key in standard_keys], dtype=np.float64
            )
            selected_chain = np.asarray(
                [chain_id == chain for chain_id in structure.chain_ids], dtype=bool
            )
            valid = (
                target_valid
                & selected_chain
                & np.isfinite(fast_chain)
                & np.isfinite(fast_complex)
                & np.isfinite(flex_pretrain)
                & np.isfinite(flex_finetune)
            )
            if int(valid.sum()) < 2:
                continue
            method_predictions = {
                "fastprotflex_chain": fast_chain,
                "fastprotflex_complex": fast_complex,
                "flexppi_pretrain_complex": flex_pretrain,
                "flexppi_finetune_complex": flex_finetune,
            }
            chain_row: Dict[str, Any] = {
                "sample_id": record.sample_id,
                "chain": chain,
                "residue_count": int(valid.sum()),
                "interface_residue_count": int(np.sum(valid & interface)),
                "non_interface_residue_count": int(np.sum(valid & ~interface)),
                "chain_count_in_complex": len(chains),
                "fastprotflex_chain_seconds": float(chain_timing["end_to_end_seconds"]),
                "fastprotflex_complex_seconds": float(
                    fast_complex_timings[record.sample_id]["end_to_end_seconds"]
                ),
            }
            for method, prediction in method_predictions.items():
                chain_row[method] = _system_metrics(
                    prediction[valid],
                    target[valid],
                    absolute_scale=method.startswith("flexppi"),
                )
            complex_chain_rows.append(chain_row)

            for region_name, region_mask, destination in (
                ("interface", interface, interface_chain_rows),
                ("non_interface", ~interface, non_interface_chain_rows),
            ):
                region_valid = valid & region_mask
                if int(region_valid.sum()) < 5:
                    continue
                region_row: Dict[str, Any] = {
                    "sample_id": record.sample_id,
                    "chain": chain,
                    "region": region_name,
                    "residue_count": int(region_valid.sum()),
                }
                for method, prediction in method_predictions.items():
                    region_row[method] = _system_metrics(
                        prediction[region_valid],
                        target[region_valid],
                        absolute_scale=method.startswith("flexppi"),
                    )
                destination.append(region_row)

    complex_rows = _aggregate_chains_to_complexes(complex_chain_rows, COMPLEX_METHODS)
    interface_rows = _aggregate_chains_to_complexes(interface_chain_rows, COMPLEX_METHODS)
    non_interface_rows = _aggregate_chains_to_complexes(
        non_interface_chain_rows, COMPLEX_METHODS
    )
    with Path(atlas_report_path).open("r", encoding="utf-8") as handle:
        atlas_report = json.load(handle)

    report: Dict[str, Any] = {
        "task": "single_chain_and_complex_context_comparison",
        "dataset_audit": {
            "validation_structure_counts": [
                {
                    "dataset": dataset_name,
                    "structure_class": structure_class,
                    "protein_chain_count": chain_count,
                    "system_count": count,
                }
                for (dataset_name, structure_class, chain_count), count
                in sorted(chain_count_audit.items())
            ],
            "source_group_counts": [
                {
                    "dataset": dataset_name,
                    "structure_class": structure_class,
                    "source_group": source_group,
                    "system_count": count,
                }
                for (dataset_name, structure_class, source_group), count
                in sorted(source_group_audit.items())
            ],
            "skipped": skipped,
        },
        "single_chain": {
            "atlas": {
                "system_count": atlas_report["dataset"]["paired_evaluable_proteins"],
                "residue_count": atlas_report["dataset"]["paired_evaluable_residues"],
                "methods": atlas_report["methods"],
                "paired_against_fastprotflex": atlas_report[
                    "paired_against_fastprotflex"
                ],
                "interpretation": (
                    "All ATLAS validation structures contain exactly one protein chain."
                ),
            },
            "dynarepo_single_protein_chain": {
                "system_count": len(single_rows),
                "residue_count": int(sum(row["residue_count"] for row in single_rows)),
                "methods": _summaries(
                    single_rows, SINGLE_METHODS, bootstrap_repetitions, config.seed
                ),
                "paired_against_fastprotflex": {
                    method: _paired_comparison(
                        single_rows,
                        method,
                        "fastprotflex",
                        bootstrap_repetitions,
                        config.seed,
                    )
                    for method in ("flexppi_pretrain", "flexppi_finetune")
                },
                "interpretation": (
                    "These systems contain one protein chain but may originate from protein-"
                    "ligand or monomer MD sources; both predictors receive protein atoms only."
                ),
            },
        },
        "protein_complex": {
            "protocol": (
                "FastProtFlex-chain sees one chain at a time. FastProtFlex-complex and FlexPPI "
                "see the full multichain PDB. Metrics are computed within chains, averaged "
                "across chains inside each complex, then macro-averaged across complexes."
            ),
            "system_count": len(complex_rows),
            "chain_count": len(complex_chain_rows),
            "residue_count": int(sum(row["residue_count"] for row in complex_rows)),
            "interface_definition": (
                "a residue has a C-alpha atom within {:.1f} Angstrom of another protein chain"
            ).format(interface_cutoff),
            "equal_complex_methods": _summaries(
                complex_rows, COMPLEX_METHODS, bootstrap_repetitions, config.seed
            ),
            "equal_chain_methods": _summaries(
                complex_chain_rows, COMPLEX_METHODS, bootstrap_repetitions, config.seed
            ),
            "primary_paired_against_independent_chain_fastprotflex": {
                method: _paired_comparison(
                    complex_rows,
                    method,
                    "fastprotflex_chain",
                    bootstrap_repetitions,
                    config.seed,
                )
                for method in (
                    "fastprotflex_complex",
                    "flexppi_pretrain_complex",
                    "flexppi_finetune_complex",
                )
            },
            "by_region": {
                "interface": {
                    "evaluable_complexes": len(interface_rows),
                    "methods": _summaries(
                        interface_rows,
                        COMPLEX_METHODS,
                        bootstrap_repetitions,
                        config.seed,
                    ),
                    "paired_against_independent_chain_fastprotflex": {
                        method: _paired_comparison(
                            interface_rows,
                            method,
                            "fastprotflex_chain",
                            bootstrap_repetitions,
                            config.seed,
                        )
                        for method in (
                            "fastprotflex_complex",
                            "flexppi_pretrain_complex",
                            "flexppi_finetune_complex",
                        )
                    },
                },
                "non_interface": {
                    "evaluable_complexes": len(non_interface_rows),
                    "methods": _summaries(
                        non_interface_rows,
                        COMPLEX_METHODS,
                        bootstrap_repetitions,
                        config.seed,
                    ),
                    "paired_against_independent_chain_fastprotflex": {
                        method: _paired_comparison(
                            non_interface_rows,
                            method,
                            "fastprotflex_chain",
                            bootstrap_repetitions,
                            config.seed,
                        )
                        for method in (
                            "fastprotflex_complex",
                            "flexppi_pretrain_complex",
                            "flexppi_finetune_complex",
                        )
                    },
                },
            },
        },
        "artifacts": {
            "fastprotflex_chain_manifest": str(chain_manifest_path.resolve()),
            "fastprotflex_chain_predictions": str(chain_prediction_path.resolve()),
            "fastprotflex_chain_timings": str(chain_timing_path.resolve()),
            "per_chain_metrics": str(per_chain_path.resolve()),
            "per_complex_metrics": str(per_complex_path.resolve()),
        },
    }

    _write_rows(
        per_chain_path,
        _flat_chain_rows(complex_chain_rows, COMPLEX_METHODS),
    )
    _write_rows(per_complex_path, _flat_chain_rows(complex_rows, COMPLEX_METHODS))
    with report_path.open("w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2, ensure_ascii=False)
    return report
