"""DynaRepo inventory and test53 apo/holo pair construction."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List, Optional

from flexppi.data.pairing import match_apo_to_holo
from flexppi.data.records import PairRecord, StateRecord, write_jsonl
from flexppi.data.structures import parse_pdb


def _source_group(description: str) -> str:
    lower = description.lower()
    if "ppis from test53" in lower:
        return "test53_complex"
    if "single chains from test53" in lower:
        return "test53_single"
    if "dynabench" in lower:
        return "dynabench_complex"
    if "pdbbind" in lower:
        return "pdbbind_complex"
    if "antigen" in lower:
        return "antigen_single"
    return "other"


def _read_metadata(path: Path) -> Dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def index_dynarepo(
    root: str,
    output_dir: str,
    label_dir: Optional[str] = None,
    limit: Optional[int] = None,
    minimum_identity: float = 0.90,
    minimum_coverage: float = 0.80,
) -> Dict[str, Any]:
    dataset_root = Path(root).resolve()
    output_root = Path(output_dir).resolve()
    labels = Path(label_dir).resolve() if label_dir else output_root / "dynarepo_labels"
    entries = sorted(path for path in dataset_root.iterdir() if path.is_dir())
    if limit is not None:
        entries = entries[:limit]

    states: List[StateRecord] = []
    by_pdb_complex: Dict[str, List[StateRecord]] = {}
    by_pdb_single: Dict[str, List[StateRecord]] = {}
    skipped = 0
    group_counts: Dict[str, int] = {}

    for entry in entries:
        metadata_path = entry / "dynarepo_metadata.json"
        structure_path = entry / "structure.pdb"
        trajectory_path = entry / "trajectory.xtc"
        if not metadata_path.exists() or not structure_path.exists():
            skipped += 1
            continue
        try:
            metadata = _read_metadata(metadata_path)
            structure = parse_pdb(str(structure_path))
        except (OSError, ValueError, json.JSONDecodeError):
            skipped += 1
            continue
        accession = str(metadata.get("accession", entry.name))
        pdb_ids = [str(value).upper() for value in metadata.get("pdb_ids", [])]
        description = str(metadata.get("description", ""))
        group = _source_group(description)
        group_counts[group] = group_counts.get(group, 0) + 1
        expected_label = labels / (accession + ".npz")
        record = StateRecord(
            sample_id=accession,
            dataset="dynarepo",
            structure_path=str(structure_path.resolve()),
            label_path=str(expected_label) if expected_label.exists() else None,
            trajectory_paths=[str(trajectory_path.resolve())] if trajectory_path.exists() else [],
            temperature=float(metadata.get("temperature", 310.0)),
            window_ns=100.0,
            target_chains=structure.unique_chains(),
            source_group=group,
            metadata={
                "residue_count": structure.num_residues,
                "pdb_ids": pdb_ids,
                "frames": int(metadata.get("frames", 0)),
                "declared_replicates": len(metadata.get("mds", [])),
                "local_trajectory_count": int(trajectory_path.exists()),
            },
        )
        states.append(record)
        for pdb_id in pdb_ids:
            if group == "test53_complex":
                by_pdb_complex.setdefault(pdb_id, []).append(record)
            elif group == "test53_single":
                by_pdb_single.setdefault(pdb_id, []).append(record)

    pairs: List[PairRecord] = []
    rejected_pairs = 0
    rejected_examples: List[Dict[str, Any]] = []
    for pdb_id in sorted(set(by_pdb_complex) & set(by_pdb_single)):
        for holo_record in by_pdb_complex[pdb_id]:
            for apo_record in by_pdb_single[pdb_id]:
                try:
                    apo_structure = parse_pdb(apo_record.structure_path)
                    holo_structure = parse_pdb(holo_record.structure_path)
                    apo_idx, holo_idx, holo_chains, identity, coverage = match_apo_to_holo(
                        apo_structure, holo_structure
                    )
                except (OSError, ValueError) as exc:
                    rejected_pairs += 1
                    rejected_examples.append(
                        {
                            "pdb_id": pdb_id,
                            "apo": apo_record.sample_id,
                            "holo": holo_record.sample_id,
                            "reason": str(exc),
                        }
                    )
                    continue
                if identity < minimum_identity or coverage < minimum_coverage:
                    rejected_pairs += 1
                    rejected_examples.append(
                        {
                            "pdb_id": pdb_id,
                            "apo": apo_record.sample_id,
                            "holo": holo_record.sample_id,
                            "identity": identity,
                            "coverage": coverage,
                            "reason": "below_threshold",
                        }
                    )
                    continue
                apo_copy = StateRecord.from_dict(apo_record.to_dict())
                holo_copy = StateRecord.from_dict(holo_record.to_dict())
                apo_copy.target_chains = apo_structure.unique_chains()
                holo_copy.target_chains = holo_chains
                pairs.append(
                    PairRecord(
                        pair_id="{}__{}__{}".format(pdb_id, apo_copy.sample_id, holo_copy.sample_id),
                        apo=apo_copy,
                        holo=holo_copy,
                        apo_indices=apo_idx,
                        holo_indices=holo_idx,
                        sequence_identity=identity,
                        coverage=coverage,
                        metadata={"pdb_id": pdb_id, "target_holo_chains": holo_chains},
                    )
                )

    state_manifest = output_root / "dynarepo_states.jsonl"
    pair_manifest = output_root / "dynarepo_pairs.jsonl"
    write_jsonl(str(state_manifest), (record.to_dict() for record in states))
    write_jsonl(str(pair_manifest), (record.to_dict() for record in pairs))
    return {
        "indexed": len(states),
        "skipped": skipped,
        "pairs": len(pairs),
        "rejected_pairs": rejected_pairs,
        "rejected_examples": rejected_examples,
        "groups": group_counts,
        "state_manifest": str(state_manifest),
        "pair_manifest": str(pair_manifest),
    }
