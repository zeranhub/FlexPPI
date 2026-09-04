"""Select one physiology-range ThermoPCD condition per biological system.

The model released from the performance-first run is structure-only, so the
same PDB cannot appear several times with contradictory temperature-specific
labels.  This inventory therefore prefers 310 K and falls back to 311 K only
when 310 K is unavailable.  Replicate trajectories at the selected temperature
are grouped into one record and must later be reduced to one mean/variance NPZ.
No trajectory is processed by this inventory command.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Tuple

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from flexppi.data.records import StateRecord, write_jsonl  # noqa: E402
from flexppi.data.structures import parse_pdb  # noqa: E402


def _split_field(value: str) -> List[str]:
    return [item.strip() for item in value.split("|") if item.strip()]


def _preferred_structure(directory: Path, pdb_id: str) -> Optional[Path]:
    candidates = sorted(directory.glob("*.pdb"))
    if not candidates:
        return None
    prefix = pdb_id.lower()
    model = [path for path in candidates if path.name.lower() == prefix + "-m.pdb"]
    if model:
        return model[0]
    zero_ns = [
        path
        for path in candidates
        if path.name.lower().startswith(prefix) and "0ns" in path.name.lower()
    ]
    if zero_ns:
        return zero_ns[0]
    prefixed = [path for path in candidates if path.name.lower().startswith(prefix)]
    return prefixed[0] if prefixed else candidates[0]


def _sequence_digest(structure_path: Path) -> Tuple[str, int, List[str]]:
    structure = parse_pdb(str(structure_path))
    chain_sequences = [structure.chain_sequence(chain) for chain in structure.unique_chains()]
    digest = hashlib.sha256("|".join(chain_sequences).encode("utf-8")).hexdigest()
    return digest, structure.num_residues, structure.unique_chains()


def _trajectory_paths(row: Mapping[str, str]) -> List[Path]:
    directory = Path(row["path"])
    return [directory / name for name in _split_field(row["production_xtc_names"])]


def prepare_candidates(
    local_entries_path: Path,
    output_manifest: Path,
    qc_table: Path,
    summary_path: Path,
    label_dir: Path,
    max_residues: int = 1024,
) -> Mapping[str, object]:
    with local_entries_path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = [dict(row) for row in csv.DictReader(handle)]
    production_complexes = [
        row
        for row in rows
        if row.get("entry_kind") == "complex" and row.get("has_production_xtc") == "True"
    ]
    by_pdb: Dict[str, List[Dict[str, str]]] = defaultdict(list)
    for row in production_complexes:
        by_pdb[row["pdb_id"].strip().upper()].append(row)

    selected_by_pdb: Dict[str, List[Dict[str, str]]] = {}
    chosen_temperature: Dict[str, int] = {}
    for pdb_id, pdb_rows in by_pdb.items():
        exact = defaultdict(list)
        for row in pdb_rows:
            temperature = row.get("temperature_k", "").strip()
            if temperature in {"310", "311"}:
                exact[int(temperature)].append(row)
        if exact.get(310):
            chosen_temperature[pdb_id] = 310
            selected_by_pdb[pdb_id] = exact[310]
        elif exact.get(311):
            chosen_temperature[pdb_id] = 311
            selected_by_pdb[pdb_id] = exact[311]

    qc_rows = []
    records: List[StateRecord] = []
    ready_systems = 0
    pending_systems = 0
    for pdb_id in sorted(by_pdb):
        selected = selected_by_pdb.get(pdb_id, [])
        selected_ids = {id(row) for row in selected}
        per_entry_details: Dict[int, Dict[str, object]] = {}
        for row in by_pdb[pdb_id]:
            detail: Dict[str, object] = {
                "structure_path": "",
                "sequence_digest": "",
                "residue_count": 0,
                "error": "",
            }
            if id(row) in selected_ids:
                structure_path = _preferred_structure(Path(row["path"]), pdb_id)
                try:
                    if structure_path is None:
                        raise FileNotFoundError("no_pdb_structure")
                    digest, residue_count, _ = _sequence_digest(structure_path)
                    detail.update(
                        {
                            "structure_path": str(structure_path.resolve()),
                            "sequence_digest": digest,
                            "residue_count": residue_count,
                        }
                    )
                except (OSError, ValueError) as exc:
                    detail["error"] = str(exc)
            per_entry_details[id(row)] = detail

        selected_details = [per_entry_details[id(row)] for row in selected]
        selected_digests = {
            str(detail["sequence_digest"])
            for detail in selected_details
            if detail["sequence_digest"]
        }
        errors = [str(detail["error"]) for detail in selected_details if detail["error"]]
        residue_counts = {
            int(detail["residue_count"])
            for detail in selected_details
            if int(detail["residue_count"])
        }
        trajectories = [path for row in selected for path in _trajectory_paths(row)]
        missing_trajectories = [str(path) for path in trajectories if not path.exists()]
        ready = bool(selected) and not errors and not missing_trajectories
        ready = ready and len(selected_digests) == 1 and len(residue_counts) == 1
        ready = ready and next(iter(residue_counts), max_residues + 1) <= max_residues

        if ready:
            first_structure = Path(str(selected_details[0]["structure_path"]))
            structure = parse_pdb(str(first_structure))
            temperature = chosen_temperature[pdb_id]
            sample_id = "thermopcd_{}_{}k".format(pdb_id.lower(), temperature)
            expected_label = label_dir / (sample_id + ".npz")
            label_path = str(expected_label.resolve()) if expected_label.exists() else None
            metadata = {
                "pdb_id": pdb_id,
                "pdb_ids": [pdb_id],
                "biological_system": "pdb:{}".format(pdb_id),
                "residue_count": structure.num_residues,
                "selected_temperature_k": temperature,
                "temperature_policy": "prefer_310_else_311",
                "source_directories": [row["path"] for row in selected],
                "seed_groups": sorted(
                    {seed for row in selected for seed in _split_field(row.get("seed", ""))}
                ),
                "trajectory_count": len(trajectories),
                "replicate_reduction": "single_system_mean_and_sample_variance",
                "quality_tier": "ready_for_uniform_label_extraction",
                "expected_label_path": str(expected_label.resolve()),
            }
            records.append(
                StateRecord(
                    sample_id=sample_id,
                    dataset="thermopcd",
                    structure_path=str(first_structure.resolve()),
                    label_path=label_path,
                    trajectory_paths=[str(path.resolve()) for path in trajectories],
                    temperature=float(temperature),
                    window_ns=max(
                        [
                            float(value)
                            for row in selected
                            for value in _split_field(row.get("duration_ns", ""))
                        ]
                        or [0.0]
                    ),
                    target_chains=structure.unique_chains(),
                    source_group="thermopcd_complex",
                    metadata=metadata,
                )
            )
            if label_path:
                ready_systems += 1
            else:
                pending_systems += 1

        for row in by_pdb[pdb_id]:
            detail = per_entry_details[id(row)]
            if id(row) in selected_ids:
                if ready:
                    status = "selected_main_pending_labels"
                    reason = "preferred_physiological_temperature"
                else:
                    status = "excluded_qc"
                    reasons = []
                    if errors:
                        reasons.append("structure_error")
                    if missing_trajectories:
                        reasons.append("missing_trajectory")
                    if len(selected_digests) > 1:
                        reasons.append("incompatible_sequences_across_replicates")
                    if len(residue_counts) > 1:
                        reasons.append("inconsistent_residue_counts")
                    if residue_counts and max(residue_counts) > max_residues:
                        reasons.append("length_gt_{}".format(max_residues))
                    reason = ";".join(reasons) or "unresolved_qc"
            else:
                raw_temperature = row.get("temperature_k", "")
                if "|" in raw_temperature:
                    reason = "ambiguous_merged_temperature_metadata"
                elif raw_temperature in {"312", "313"}:
                    reason = "reserved_temperature_ood"
                elif raw_temperature == "309":
                    reason = "outside_main_temperature_policy"
                elif pdb_id in chosen_temperature:
                    reason = "nonpreferred_physiological_temperature"
                else:
                    reason = "no_exact_310_or_311_condition"
                status = "not_selected"
            qc_rows.append(
                {
                    "pdb_id": pdb_id,
                    "temperature_k": row.get("temperature_k", ""),
                    "seed": row.get("seed", ""),
                    "duration_ns": row.get("duration_ns", ""),
                    "directory": row["path"],
                    "trajectory_count": row.get("production_xtc_count", ""),
                    "structure_path": detail["structure_path"],
                    "residue_count": detail["residue_count"],
                    "sequence_digest": detail["sequence_digest"],
                    "error": detail["error"],
                    "status": status,
                    "reason": reason,
                }
            )

    output_manifest.parent.mkdir(parents=True, exist_ok=True)
    write_jsonl(str(output_manifest), (record.to_dict() for record in records))
    qc_table.parent.mkdir(parents=True, exist_ok=True)
    with qc_table.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(qc_rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(qc_rows)
    summary: Mapping[str, object] = {
        "production_complex_directories": len(production_complexes),
        "production_complex_pdb_ids": len(by_pdb),
        "selected_systems": len(records),
        "selected_systems_with_labels": ready_systems,
        "selected_systems_pending_uniform_labels": pending_systems,
        "selected_temperature_counts": {
            str(temperature): sum(1 for record in records if int(record.temperature) == temperature)
            for temperature in (310, 311)
        },
        "pre_qc_temperature_choice_counts": {
            str(temperature): sum(
                1 for value in chosen_temperature.values() if value == temperature
            )
            for temperature in (310, 311)
        },
        "selection_policy": (
            "One record per PDB; prefer exact 310 K, otherwise exact 311 K; merge trajectories "
            "only when protein sequences and residue counts agree; reserve 312/313 K for OOD."
        ),
        "output_manifest": str(output_manifest.resolve()),
    }
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--local-entries",
        default=str(ROOT / "artifacts/thermopcd_audit/local_entries.csv"),
    )
    parser.add_argument(
        "--output-manifest",
        default=str(ROOT / "artifacts/data/thermopcd_performance_first_candidates.jsonl"),
    )
    parser.add_argument(
        "--qc-table",
        default=str(ROOT / "data/splits/state_rmsf_hchi_v2_performance_first/thermopcd_qc.tsv"),
    )
    parser.add_argument(
        "--summary",
        default=str(
            ROOT / "data/splits/state_rmsf_hchi_v2_performance_first/thermopcd_qc_summary.json"
        ),
    )
    parser.add_argument("--label-dir", default=str(ROOT / "artifacts/data/thermopcd_labels_v2"))
    parser.add_argument("--max-residues", type=int, default=1024)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    summary = prepare_candidates(
        Path(args.local_entries),
        Path(args.output_manifest),
        Path(args.qc_table),
        Path(args.summary),
        Path(args.label_dir),
        args.max_residues,
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
