"""Summarize local FlexPPI datasets without changing training artifacts.

The report deliberately distinguishes raw files, state records, independent PDB
identifiers, residue labels, and strict apo/holo pairs.  Those counts are not
interchangeable and have very different statistical meaning.
"""

from __future__ import annotations

import argparse
import csv
import json
import random
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Sequence, Set

import numpy as np


METRICS = ("rmsf", "chi_entropy", "sasa")


def read_jsonl(path: Path) -> List[Dict[str, Any]]:
    with path.open("r", encoding="utf-8") as handle:
        return [json.loads(line) for line in handle if line.strip()]


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def local_label_path(root: Path, record: Mapping[str, Any]) -> Path:
    dataset = str(record["dataset"]).lower()
    sample_id = str(record["sample_id"])
    if dataset == "atlas":
        return root / "artifacts" / "data" / "atlas_labels" / f"{sample_id}.npz"
    if dataset == "dynarepo":
        return root / "artifacts" / "data" / "dynarepo_labels" / f"{sample_id}.npz"
    if dataset == "zenodo_6638504":
        return root / "artifacts" / "zenodo_6638504" / "data" / "labels" / f"{sample_id}.npz"
    supplied = Path(str(record.get("label_path", "")))
    return supplied


def pdb_ids_for_state(record: Mapping[str, Any]) -> Set[str]:
    metadata = record.get("metadata", {})
    values = metadata.get("pdb_ids") or ([metadata["pdb_id"]] if metadata.get("pdb_id") else [])
    if values:
        return {str(value).upper() for value in values}
    sample_id = str(record.get("sample_id", ""))
    if len(sample_id) >= 4 and sample_id[:4].isalnum():
        return {sample_id[:4].upper()}
    return set()


def label_counts(root: Path, records: Sequence[Mapping[str, Any]]) -> Dict[str, Any]:
    residue_labels = Counter({metric: 0 for metric in METRICS})
    states_with_metric = Counter({metric: 0 for metric in METRICS})
    missing_files: List[str] = []
    for record in records:
        path = local_label_path(root, record)
        if not path.exists():
            missing_files.append(str(record.get("sample_id")))
            continue
        with np.load(path, allow_pickle=False) as label:
            names = [str(value) for value in label["metric_names"].tolist()]
            mask = np.asarray(label["mask"], dtype=bool)
        for metric in METRICS:
            if metric not in names:
                continue
            count = int(mask[:, names.index(metric)].sum())
            residue_labels[metric] += count
            if count:
                states_with_metric[metric] += 1
    return {
        "residue_labels": dict(residue_labels),
        "states_with_metric": dict(states_with_metric),
        "missing_label_files": missing_files,
    }


def state_summary(root: Path, records: Sequence[Mapping[str, Any]], limit: int) -> Dict[str, Any]:
    kept = [record for record in records if int(record.get("metadata", {}).get("residue_count", 0)) <= limit]
    all_pdb = set().union(*(pdb_ids_for_state(record) for record in records)) if records else set()
    kept_pdb = set().union(*(pdb_ids_for_state(record) for record in kept)) if kept else set()
    return {
        "manifest_states": len(records),
        "unique_sample_ids": len({str(record["sample_id"]) for record in records}),
        "unique_pdb_ids": len(all_pdb),
        "total_residues_from_metadata": sum(int(record.get("metadata", {}).get("residue_count", 0)) for record in records),
        "source_groups": dict(sorted(Counter(str(record.get("source_group", "unknown")) for record in records).items())),
        "temperatures_k": dict(sorted(Counter(str(record.get("temperature")) for record in records).items())),
        "length_limit": limit,
        "states_after_length_limit": len(kept),
        "pdb_ids_after_length_limit": len(kept_pdb),
        "residues_after_length_limit": sum(int(record.get("metadata", {}).get("residue_count", 0)) for record in kept),
        "excluded_by_length": len(records) - len(kept),
        "all_label_coverage": label_counts(root, records),
        "kept_label_coverage": label_counts(root, kept),
    }


def trainer_group(record: Mapping[str, Any], is_pair: bool = False) -> str:
    metadata = record.get("metadata", {})
    return str(metadata.get("cluster_id") or metadata.get("pdb_id") or record.get("pair_id" if is_pair else "sample_id", "unknown"))


def current_split(records: Sequence[Mapping[str, Any]], fraction: float, seed: int, is_pair: bool = False) -> Dict[str, Any]:
    grouped: Dict[str, List[str]] = defaultdict(list)
    for record in records:
        grouped[trainer_group(record, is_pair=is_pair)].append(str(record.get("pair_id" if is_pair else "sample_id")))
    groups = sorted(grouped)
    random.Random(seed).shuffle(groups)
    validation_count = max(1, int(round(len(groups) * fraction)))
    validation_groups = set(groups[:validation_count])
    train_items = [record for record in records if trainer_group(record, is_pair=is_pair) not in validation_groups]
    validation_items = [record for record in records if trainer_group(record, is_pair=is_pair) in validation_groups]
    if is_pair:
        record_pdb_ids = lambda record: {str(record.get("metadata", {}).get("pdb_id", "")).upper()}
    else:
        record_pdb_ids = pdb_ids_for_state
    train_pdb = set().union(*(record_pdb_ids(record) for record in train_items)) if train_items else set()
    validation_pdb = set().union(*(record_pdb_ids(record) for record in validation_items)) if validation_items else set()
    leaked_pdb = sorted(train_pdb & validation_pdb)
    return {
        "group_count": len(groups),
        "train_records": len(train_items),
        "validation_records": len(validation_items),
        "exact_pdb_ids_shared_between_train_and_validation": leaked_pdb,
        "exact_pdb_leakage_count": len(leaked_pdb),
        "validation_groups": sorted(validation_groups),
    }


def pair_summary(root: Path, records: Sequence[Mapping[str, Any]], limit: int) -> Dict[str, Any]:
    kept = [
        record for record in records
        if int(record["apo"].get("metadata", {}).get("residue_count", 0)) <= limit
        and int(record["holo"].get("metadata", {}).get("residue_count", 0)) <= limit
    ]
    mapped = lambda items: sum(len(record.get("apo_indices", [])) for record in items)
    pdbs = lambda items: {str(record.get("metadata", {}).get("pdb_id", "")).upper() for record in items}
    return {
        "manifest_pairs": len(records),
        "unique_pdb_ids": len(pdbs(records)),
        "mapped_target_residues": mapped(records),
        "length_limit": limit,
        "pairs_after_length_limit": len(kept),
        "pdb_ids_after_length_limit": len(pdbs(kept)),
        "mapped_target_residues_after_length_limit": mapped(kept),
        "excluded_pair_ids": [str(record["pair_id"]) for record in records if record not in kept],
        "current_split": current_split(kept, 0.2, 7, is_pair=True),
    }


def exact_overlap_rows(sets: Mapping[str, Set[str]]) -> List[Dict[str, Any]]:
    names = list(sets)
    rows: List[Dict[str, Any]] = []
    for index, left in enumerate(names):
        for right in names[index + 1 :]:
            overlap = sorted(sets[left] & sets[right])
            rows.append({"left": left, "right": right, "count": len(overlap), "pdb_ids": "|".join(overlap)})
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--output-dir", type=Path, default=None)
    args = parser.parse_args()
    root = args.root.resolve()
    output = (args.output_dir or root / "artifacts" / "dataset_inventory").resolve()
    output.mkdir(parents=True, exist_ok=True)

    atlas = read_jsonl(root / "artifacts" / "data" / "atlas_states.jsonl")
    dyna = read_jsonl(root / "artifacts" / "data" / "dynarepo_states.jsonl")
    dyna_pairs = read_jsonl(root / "artifacts" / "data" / "dynarepo_pairs.jsonl")
    zenodo_states = read_jsonl(root / "artifacts" / "zenodo_6638504" / "data" / "zenodo_6638504_states.jsonl")
    zenodo_pairs = read_jsonl(root / "artifacts" / "zenodo_6638504" / "data" / "zenodo_6638504_pairs.jsonl")

    atlas_kept = [record for record in atlas if int(record["metadata"].get("residue_count", 0)) <= 768]
    dyna_kept = [record for record in dyna if int(record["metadata"].get("residue_count", 0)) <= 768]
    pretrain_kept = atlas_kept + dyna_kept

    thermopcd_summary = json.loads((root / "artifacts" / "thermopcd_audit" / "summary.json").read_text(encoding="utf-8"))
    thermopcd_refit = json.loads((root / "artifacts" / "thermopcd_audit" / "refit_pair_summary.json").read_text(encoding="utf-8"))
    with (root / "artifacts" / "thermopcd_audit" / "pdb_inventory.csv").open("r", encoding="utf-8-sig", newline="") as handle:
        thermopcd_pdb = {row["pdb_id"].upper() for row in csv.DictReader(handle)}

    tcr_manifest = root / "artifacts" / "evaluation" / "fastprotflex_comparison" / "tcr3d" / "tcr3d_common_manifest.tsv"
    with tcr_manifest.open("r", encoding="utf-8-sig", newline="") as handle:
        tcr_rows = list(csv.DictReader(handle, delimiter="\t"))
    tcr_pdb = {row["sample_id"][:4].upper() for row in tcr_rows}

    itrace_preparation = json.loads((root / "artifacts" / "data" / "itrace_context_diagnostic" / "preparation_report.json").read_text(encoding="utf-8"))
    itrace_states = read_jsonl(root / "artifacts" / "data" / "itrace_context_diagnostic" / "itrace_states.jsonl")
    itrace_pdb = set().union(*(pdb_ids_for_state(record) for record in itrace_states))

    datasets = {
        "atlas": state_summary(root, atlas, 768),
        "dynarepo_states": state_summary(root, dyna, 768),
        "dynarepo_pairs": pair_summary(root, dyna_pairs, 1024),
        "current_pretrain_effective": {
            "states": len(pretrain_kept),
            "residues_from_metadata": sum(int(record["metadata"].get("residue_count", 0)) for record in pretrain_kept),
            "current_split": current_split(pretrain_kept, 0.1, 7),
        },
        "zenodo_6638504": {
            **state_summary(root, zenodo_states, 1024),
            "strict_pairs": len(zenodo_pairs),
            "pair_pdb_ids": len({str(record["metadata"]["pdb_id"]).upper() for record in zenodo_pairs}),
            "mapped_target_residues": sum(len(record.get("apo_indices", [])) for record in zenodo_pairs),
        },
        "thermopcd": {"inventory": thermopcd_summary, "refit_pilot": thermopcd_refit},
        "tcr3d": {
            "common_chain_files": len(tcr_rows),
            "unique_pdb_ids": len(tcr_pdb),
            "role": "external crystallographic proxy; no MD trajectories",
        },
        "itrace": {
            **itrace_preparation,
            "unique_pdb_ids": len(itrace_pdb),
            "role": "analysis-only known-failure stress diagnostic; not untouched test",
        },
    }

    sets = {
        "atlas": set().union(*(pdb_ids_for_state(record) for record in atlas)),
        "dynarepo_states": set().union(*(pdb_ids_for_state(record) for record in dyna)),
        "dynarepo_pairs": {str(record["metadata"]["pdb_id"]).upper() for record in dyna_pairs},
        "zenodo_6638504": {str(record["metadata"]["pdb_id"]).upper() for record in zenodo_pairs},
        "thermopcd": thermopcd_pdb,
        "tcr3d": tcr_pdb,
        "itrace": itrace_pdb,
    }
    overlaps = exact_overlap_rows(sets)
    with (output / "exact_pdb_overlaps.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=("left", "right", "count", "pdb_ids"))
        writer.writeheader()
        writer.writerows(overlaps)

    exclusions: List[Dict[str, Any]] = []
    for dataset_name, records, limit in (("atlas", atlas, 768), ("dynarepo_states", dyna, 768)):
        for record in records:
            count = int(record["metadata"].get("residue_count", 0))
            if count > limit:
                exclusions.append({"dataset": dataset_name, "record_id": record["sample_id"], "residue_count": count, "limit": limit})
    for record in dyna_pairs:
        apo_count = int(record["apo"]["metadata"].get("residue_count", 0))
        holo_count = int(record["holo"]["metadata"].get("residue_count", 0))
        if apo_count > 1024 or holo_count > 1024:
            exclusions.append({"dataset": "dynarepo_pairs", "record_id": record["pair_id"], "residue_count": max(apo_count, holo_count), "limit": 1024})
    with (output / "length_filter_exclusions.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=("dataset", "record_id", "residue_count", "limit"))
        writer.writeheader()
        writer.writerows(exclusions)

    write_json(output / "dataset_summary.json", datasets)
    print(json.dumps(datasets, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
