"""ATLAS indexing using its precomputed three-replica RMSF tables."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

from flexppi.constants import METRIC_NAMES
from flexppi.data.pairing import global_sequence_alignment
from flexppi.data.records import StateRecord, write_jsonl
from flexppi.data.structures import parse_pdb


def _read_rmsf(path: Path) -> Tuple[str, np.ndarray]:
    sequence: List[str] = []
    values: List[List[float]] = []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"seq", "RMSF_R1", "RMSF_R2", "RMSF_R3"}
        if not required.issubset(set(reader.fieldnames or [])):
            raise ValueError("Unexpected ATLAS RMSF columns in {}".format(path))
        for row in reader:
            sequence.append(row["seq"].strip().upper() or "X")
            values.append([float(row["RMSF_R1"]), float(row["RMSF_R2"]), float(row["RMSF_R3"])])
    return "".join(sequence), np.asarray(values, dtype=np.float32)


def prepare_atlas(
    root: str,
    output_dir: str,
    limit: Optional[int] = None,
) -> Dict[str, int]:
    """Create an absolute-state manifest and compact NPZ label files."""

    dataset_root = Path(root).resolve()
    output_root = Path(output_dir).resolve()
    label_root = output_root / "atlas_labels"
    label_root.mkdir(parents=True, exist_ok=True)
    records: List[StateRecord] = []
    skipped = 0

    entries = sorted(path for path in dataset_root.iterdir() if path.is_dir())
    if limit is not None:
        entries = entries[:limit]
    for entry in entries:
        sample_id = entry.name
        pdb_path = entry / (sample_id + ".pdb")
        rmsf_path = entry / (sample_id + "_RMSF.tsv")
        if not pdb_path.exists() or not rmsf_path.exists():
            skipped += 1
            continue
        try:
            label_path = label_root / (sample_id + ".npz")
            trajectories = sorted(str(path.resolve()) for path in entry.glob(sample_id + "_R*.xtc"))
            if label_path.exists():
                with np.load(label_path, allow_pickle=False) as cached:
                    residue_count = int(cached["values"].shape[0])
                records.append(
                    StateRecord(
                        sample_id=sample_id,
                        dataset="atlas",
                        structure_path=str(pdb_path.resolve()),
                        label_path=str(label_path),
                        trajectory_paths=trajectories,
                        temperature=300.0,
                        window_ns=100.0,
                        target_chains=[],
                        source_group="atlas_absolute",
                        metadata={
                            "residue_count": residue_count,
                            "replicates": 3,
                            "resumed_from_cached_label": True,
                        },
                    )
                )
                continue
            structure = parse_pdb(str(pdb_path))
            label_sequence, replicas = _read_rmsf(rmsf_path)
            values = np.full((structure.num_residues, len(METRIC_NAMES)), np.nan, dtype=np.float32)
            variances = np.zeros_like(values)
            mask = np.zeros_like(values, dtype=bool)
            alignment = global_sequence_alignment(structure.sequence, label_sequence)
            for structure_idx, label_idx in zip(alignment.indices_a, alignment.indices_b):
                row = replicas[label_idx]
                values[structure_idx, 0] = float(row.mean())
                variances[structure_idx, 0] = float(row.var(ddof=1))
                mask[structure_idx, 0] = True
            np.savez_compressed(
                label_path,
                values=values,
                variances=variances,
                mask=mask,
                metric_names=np.asarray(METRIC_NAMES),
            )
            records.append(
                StateRecord(
                    sample_id=sample_id,
                    dataset="atlas",
                    structure_path=str(pdb_path.resolve()),
                    label_path=str(label_path),
                    trajectory_paths=trajectories,
                    temperature=300.0,
                    window_ns=100.0,
                    target_chains=structure.unique_chains(),
                    source_group="atlas_absolute",
                    metadata={
                        "residue_count": structure.num_residues,
                        "alignment_identity": alignment.identity,
                        "alignment_coverage": alignment.coverage,
                        "replicates": 3,
                    },
                )
            )
        except (OSError, ValueError, KeyError):
            skipped += 1

    manifest_path = output_root / "atlas_states.jsonl"
    write_jsonl(str(manifest_path), (record.to_dict() for record in records))
    return {"indexed": len(records), "skipped": skipped, "manifest": str(manifest_path)}
