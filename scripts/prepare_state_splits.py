"""Build leakage-aware fixed splits for absolute RMSF and chi-entropy training.

The script treats a complete state record as the indivisible statistical unit.
All records sharing a PDB identifier are joined, and records are also joined
when any of their protein chains belong to the same MMseqs2 sequence cluster.
External manifests participate in clustering but never enter training; every
training record connected to an external sequence cluster is excluded.

Only absolute state manifests are accepted.  No apo/holo pair manifest or
binding-response label is read by this script.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import shlex
import subprocess
import sys
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Mapping, Sequence, Set, Tuple

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from flexppi.constants import METRIC_NAMES  # noqa: E402
from flexppi.data.labels import load_label_file  # noqa: E402
from flexppi.data.records import (  # noqa: E402
    StateRecord,
    load_state_records,
    remap_state_record,
    write_jsonl,
)
from flexppi.data.structures import parse_pdb  # noqa: E402


TRAINING_METRICS = ("rmsf", "chi_entropy")
SPLIT_NAMES = ("train", "validation", "test")
PDB_PATTERN = re.compile(r"^[0-9][A-Za-z0-9]{3}$")


class UnionFind:
    def __init__(self, size: int):
        self.parent = list(range(size))
        self.rank = [0] * size

    def find(self, value: int) -> int:
        while self.parent[value] != value:
            self.parent[value] = self.parent[self.parent[value]]
            value = self.parent[value]
        return value

    def union(self, first: int, second: int) -> None:
        first_root = self.find(first)
        second_root = self.find(second)
        if first_root == second_root:
            return
        if self.rank[first_root] < self.rank[second_root]:
            first_root, second_root = second_root, first_root
        self.parent[second_root] = first_root
        if self.rank[first_root] == self.rank[second_root]:
            self.rank[first_root] += 1


@dataclass
class Candidate:
    key: str
    role: str
    record: StateRecord
    original_record: Dict[str, object]
    manifest: str
    pdb_ids: List[str] = field(default_factory=list)
    residue_count: int = 0
    chain_sequences: List[Tuple[str, str]] = field(default_factory=list)
    metric_counts: Dict[str, int] = field(default_factory=dict)
    sequence_clusters: Set[str] = field(default_factory=set)
    exclusion_reason: str = ""
    biological_group: str = ""
    split: str = ""


def _parse_path_mappings(values: Sequence[str]) -> Dict[str, str]:
    mappings: Dict[str, str] = {}
    for value in values:
        if "=" not in value:
            raise ValueError("--path-map must use SOURCE=TARGET")
        source, target = value.split("=", 1)
        mappings[source.strip()] = target.strip()
    return mappings


def _default_path_mappings() -> Dict[str, str]:
    candidates = {
        "/ranze/repos/FlexPPI/artifacts/data": str(ROOT / "artifacts" / "data"),
        "/ranze/data/ATLAS_extracted": r"E:\data\protein_md_datasets\ATLAS_extracted",
        "/ranze/data/dynarepo/dynarepo_non_nucleic_full": (
            r"D:\data\dynarepo_non_nucleic_full"
        ),
    }
    return {source: target for source, target in candidates.items() if Path(target).exists()}


def _pdb_ids(record: StateRecord) -> List[str]:
    raw: List[object] = []
    metadata = record.metadata or {}
    for name in ("pdb_id", "pdb_ids"):
        value = metadata.get(name)
        if isinstance(value, (list, tuple)):
            raw.extend(value)
        elif value:
            raw.append(value)
    sample_prefix = record.sample_id.split("_", 1)[0]
    if PDB_PATTERN.fullmatch(sample_prefix):
        raw.append(sample_prefix)
    output = []
    for value in raw:
        normalised = str(value).strip().upper()
        if PDB_PATTERN.fullmatch(normalised) and normalised not in output:
            output.append(normalised)
    return output


def _length_bin(length: int) -> str:
    if length <= 256:
        return "short"
    if length <= 512:
        return "medium"
    if length <= 768:
        return "long"
    return "very_long"


def _stable_digest(*values: object) -> str:
    text = "|".join(str(value) for value in values)
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _load_candidates(
    manifests: Sequence[str],
    role: str,
    mappings: Mapping[str, str],
    max_residues: int,
) -> List[Candidate]:
    candidates: List[Candidate] = []
    for manifest_index, manifest in enumerate(manifests):
        for record_index, original in enumerate(load_state_records(manifest)):
            record = remap_state_record(original, mappings)
            key = "{}:{:02d}:{:06d}:{}".format(
                role, manifest_index, record_index, record.sample_id
            )
            candidate = Candidate(
                key=key,
                role=role,
                record=record,
                original_record=original.to_dict(),
                manifest=str(Path(manifest).resolve()),
                pdb_ids=_pdb_ids(record),
            )
            try:
                structure_path = Path(record.structure_path)
                if not structure_path.exists():
                    raise FileNotFoundError("missing_structure")
                structure = parse_pdb(str(structure_path))
                candidate.residue_count = structure.num_residues
                if role == "training" and structure.num_residues > max_residues:
                    candidate.exclusion_reason = "length_gt_{}".format(max_residues)
                    candidates.append(candidate)
                    continue
                for chain_id in structure.unique_chains():
                    sequence = structure.chain_sequence(chain_id)
                    if len(sequence) >= 20:
                        candidate.chain_sequences.append((chain_id, sequence))
                if not candidate.chain_sequences:
                    raise ValueError("no_chain_with_at_least_20_residues")
                if role == "training":
                    if not record.label_path or not Path(record.label_path).exists():
                        raise FileNotFoundError("missing_label")
                    values, _, mask = load_label_file(str(record.label_path))
                    if values.shape[0] != structure.num_residues:
                        raise ValueError("structure_label_length_mismatch")
                    for metric in TRAINING_METRICS:
                        metric_index = METRIC_NAMES.index(metric)
                        valid = mask[:, metric_index] & np.isfinite(values[:, metric_index])
                        candidate.metric_counts[metric] = int(valid.sum())
                    if candidate.metric_counts["rmsf"] == 0:
                        raise ValueError("no_rmsf_labels")
            except Exception as exc:
                reason = str(exc).strip() or type(exc).__name__
                candidate.exclusion_reason = re.sub(r"\s+", "_", reason)[:160]
            candidates.append(candidate)
    return candidates


def _write_fasta_and_index(
    candidates: Sequence[Candidate],
    fasta_path: Path,
    index_path: Path,
) -> Dict[str, int]:
    sequence_to_candidate: Dict[str, int] = {}
    with fasta_path.open("w", encoding="utf-8", newline="\n") as fasta, index_path.open(
        "w", encoding="utf-8", newline=""
    ) as index_handle:
        writer = csv.writer(index_handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sequence_id", "candidate_key", "role", "sample_id", "chain_id", "length"])
        for candidate_index, candidate in enumerate(candidates):
            if candidate.exclusion_reason:
                continue
            for chain_index, (chain_id, sequence) in enumerate(candidate.chain_sequences):
                sequence_id = "SEQ{:06d}_{:02d}".format(candidate_index, chain_index)
                fasta.write(">{}\n{}\n".format(sequence_id, sequence))
                writer.writerow(
                    [
                        sequence_id,
                        candidate.key,
                        candidate.role,
                        candidate.record.sample_id,
                        chain_id,
                        len(sequence),
                    ]
                )
                sequence_to_candidate[sequence_id] = candidate_index
    return sequence_to_candidate


def _windows_to_wsl(path: Path) -> str:
    resolved = path.resolve()
    drive, tail = os.path.splitdrive(str(resolved))
    if not drive:
        return str(resolved).replace("\\", "/")
    return "/mnt/host/{}/{}".format(
        drive[0].lower(), tail.lstrip("\\/").replace("\\", "/")
    )


def _run_mmseqs(
    fasta_path: Path,
    output_dir: Path,
    mmseqs_path: Path,
    minimum_identity: float,
    minimum_coverage: float,
    threads: int,
) -> Tuple[Path, str]:
    # WSL creates a ``latest`` symbolic link inside the MMseqs2 work directory.
    # Python 3.9 on Windows cannot reliably remove that link, so the cache lives
    # under the already Git-ignored artifacts tree instead of a TemporaryDirectory.
    cache_root = ROOT / "artifacts" / "mmseqs_work"
    cache_root.mkdir(parents=True, exist_ok=True)
    temp_dir = Path(tempfile.mkdtemp(prefix="state_split_", dir=str(cache_root)))
    result_prefix = temp_dir / "result"
    work_dir = temp_dir / "work"
    arguments = [
        "easy-cluster",
        str(fasta_path),
        str(result_prefix),
        str(work_dir),
        "--min-seq-id",
        str(minimum_identity),
        "-c",
        str(minimum_coverage),
        "--cov-mode",
        "0",
        "--cluster-mode",
        "0",
        "-s",
        "7.5",
        "--threads",
        str(threads),
        "-v",
        "1",
    ]
    if os.name == "nt":
        linux_arguments = [
            arguments[0],
            _windows_to_wsl(Path(arguments[1])),
            _windows_to_wsl(Path(arguments[2])),
            _windows_to_wsl(Path(arguments[3])),
            *arguments[4:],
        ]
        command = "{} {}".format(
            shlex.quote(_windows_to_wsl(mmseqs_path)),
            " ".join(shlex.quote(value) for value in linux_arguments),
        )
        subprocess.run(["wsl.exe", "sh", "-lc", command], check=True)
        version_command = "{} version".format(shlex.quote(_windows_to_wsl(mmseqs_path)))
        version = subprocess.check_output(
            ["wsl.exe", "sh", "-lc", version_command], text=True
        ).strip().splitlines()[-1]
    else:
        subprocess.run([str(mmseqs_path), *arguments], check=True)
        version = subprocess.check_output(
            [str(mmseqs_path), "version"], text=True
        ).strip().splitlines()[-1]
    generated = Path(str(result_prefix) + "_cluster.tsv")
    if not generated.exists():
        raise FileNotFoundError("MMseqs2 did not create {}".format(generated))
    destination = output_dir / "sequence_clusters.tsv"
    destination.write_bytes(generated.read_bytes())
    return destination, version


def _apply_sequence_clusters(
    cluster_path: Path,
    sequence_to_candidate: Mapping[str, int],
    candidates: Sequence[Candidate],
    union_find: UnionFind,
) -> None:
    representative_first: Dict[str, int] = {}
    with cluster_path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            columns = line.rstrip("\n").split("\t")
            if len(columns) < 2:
                raise ValueError("Malformed MMseqs2 TSV line {}".format(line_number))
            representative, member = columns[:2]
            if member not in sequence_to_candidate:
                continue
            member_index = sequence_to_candidate[member]
            candidates[member_index].sequence_clusters.add(representative)
            first_index = representative_first.setdefault(representative, member_index)
            union_find.union(first_index, member_index)


def _union_pdb_identifiers(candidates: Sequence[Candidate], union_find: UnionFind) -> None:
    first_by_pdb: Dict[str, int] = {}
    for index, candidate in enumerate(candidates):
        if candidate.exclusion_reason:
            continue
        for pdb_id in candidate.pdb_ids:
            if pdb_id in first_by_pdb:
                union_find.union(first_by_pdb[pdb_id], index)
            else:
                first_by_pdb[pdb_id] = index


def _component_stratum(indices: Sequence[int], candidates: Sequence[Candidate]) -> str:
    training = [candidates[index] for index in indices if candidates[index].role == "training"]
    datasets = sorted({candidate.record.dataset for candidate in training})
    hchi = any(candidate.metric_counts.get("chi_entropy", 0) > 0 for candidate in training)
    median_length = int(np.median([candidate.residue_count for candidate in training]))
    return "{}|hchi={}|{}".format(
        "+".join(datasets), int(hchi), _length_bin(median_length)
    )


def _assign_splits(
    components: Sequence[Tuple[str, List[int]]],
    candidates: Sequence[Candidate],
    fractions: Mapping[str, float],
    seed: int,
) -> Dict[str, str]:
    by_stratum: Dict[str, List[Tuple[str, List[int]]]] = {}
    for group_id, indices in components:
        by_stratum.setdefault(_component_stratum(indices, candidates), []).append(
            (group_id, indices)
        )
    assignments: Dict[str, str] = {}
    for stratum, groups in sorted(by_stratum.items()):
        groups = sorted(groups, key=lambda item: _stable_digest(seed, stratum, item[0]))
        total = sum(len(indices) for _, indices in groups)
        targets = {name: total * fractions[name] for name in SPLIT_NAMES}
        current = {name: 0 for name in SPLIT_NAMES}
        for group_id, indices in groups:
            size = len(indices)
            best_split = min(
                SPLIT_NAMES,
                key=lambda proposed: (
                    sum(
                        (
                            (
                                current[name]
                                + (size if name == proposed else 0)
                                - targets[name]
                            )
                            / max(targets[name], 1.0)
                        )
                        ** 2
                        for name in SPLIT_NAMES
                    ),
                    SPLIT_NAMES.index(proposed),
                ),
            )
            assignments[group_id] = best_split
            current[best_split] += size
    return assignments


def _record_with_split(candidate: Candidate, version: str) -> Dict[str, object]:
    output = dict(candidate.original_record)
    metadata = dict(output.get("metadata", {}))
    metadata.update(
        {
            "pdb_id": candidate.pdb_ids[0] if candidate.pdb_ids else None,
            "pdb_ids": candidate.pdb_ids,
            "cluster_id": candidate.biological_group,
            "sequence_clusters": sorted(candidate.sequence_clusters),
            "split": candidate.split,
            "split_version": version,
            "training_metrics": list(TRAINING_METRICS),
        }
    )
    if candidate.exclusion_reason:
        metadata["exclusion_reason"] = candidate.exclusion_reason
    output["metadata"] = metadata
    return output


def _write_master_table(path: Path, candidates: Sequence[Candidate]) -> None:
    columns = [
        "role",
        "dataset",
        "sample_id",
        "source_group",
        "pdb_ids",
        "residue_count",
        "rmsf_labels",
        "chi_entropy_labels",
        "chain_count",
        "sequence_clusters",
        "biological_group",
        "split",
        "exclusion_reason",
        "source_manifest",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for candidate in candidates:
            writer.writerow(
                {
                    "role": candidate.role,
                    "dataset": candidate.record.dataset,
                    "sample_id": candidate.record.sample_id,
                    "source_group": candidate.record.source_group,
                    "pdb_ids": ";".join(candidate.pdb_ids),
                    "residue_count": candidate.residue_count,
                    "rmsf_labels": candidate.metric_counts.get("rmsf", 0),
                    "chi_entropy_labels": candidate.metric_counts.get("chi_entropy", 0),
                    "chain_count": len(candidate.chain_sequences),
                    "sequence_clusters": ";".join(sorted(candidate.sequence_clusters)),
                    "biological_group": candidate.biological_group,
                    "split": candidate.split,
                    "exclusion_reason": candidate.exclusion_reason,
                    "source_manifest": candidate.manifest,
                }
            )


def _validate_no_leakage(candidates: Sequence[Candidate]) -> Dict[str, object]:
    group_splits: Dict[str, Set[str]] = {}
    pdb_splits: Dict[str, Set[str]] = {}
    sequence_splits: Dict[str, Set[str]] = {}
    for candidate in candidates:
        if candidate.role != "training" or candidate.split not in SPLIT_NAMES:
            continue
        group_splits.setdefault(candidate.biological_group, set()).add(candidate.split)
        for pdb_id in candidate.pdb_ids:
            pdb_splits.setdefault(pdb_id, set()).add(candidate.split)
        for sequence_cluster in candidate.sequence_clusters:
            sequence_splits.setdefault(sequence_cluster, set()).add(candidate.split)
    violations = {
        "biological_group": sorted(key for key, value in group_splits.items() if len(value) > 1),
        "pdb_id": sorted(key for key, value in pdb_splits.items() if len(value) > 1),
        "sequence_cluster": sorted(
            key for key, value in sequence_splits.items() if len(value) > 1
        ),
    }
    if any(violations.values()):
        raise RuntimeError("Leakage validation failed: {}".format(violations))
    return {"passed": True, "violations": violations}


def prepare_splits(args: argparse.Namespace) -> Dict[str, object]:
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    mappings = _default_path_mappings()
    mappings.update(_parse_path_mappings(args.path_map))
    training = _load_candidates(args.manifest, "training", mappings, args.max_residues)
    external = _load_candidates(args.external_manifest, "external", mappings, args.max_residues)
    candidates = training + external

    fasta_path = output_dir / "all_chains.fasta"
    sequence_index_path = output_dir / "sequence_index.tsv"
    sequence_to_candidate = _write_fasta_and_index(
        candidates, fasta_path, sequence_index_path
    )
    cluster_path, mmseqs_version = _run_mmseqs(
        fasta_path,
        output_dir,
        Path(args.mmseqs),
        args.minimum_identity,
        args.minimum_coverage,
        args.threads,
    )

    union_find = UnionFind(len(candidates))
    _union_pdb_identifiers(candidates, union_find)
    _apply_sequence_clusters(cluster_path, sequence_to_candidate, candidates, union_find)

    components: Dict[int, List[int]] = {}
    for index, candidate in enumerate(candidates):
        if not candidate.exclusion_reason:
            components.setdefault(union_find.find(index), []).append(index)

    eligible_components: List[Tuple[str, List[int]]] = []
    for indices in components.values():
        digest = _stable_digest(*(sorted(candidates[index].key for index in indices)))[:12]
        group_id = "biogroup_{}".format(digest)
        contains_external = any(candidates[index].role == "external" for index in indices)
        training_indices = [
            index for index in indices if candidates[index].role == "training"
        ]
        for index in indices:
            candidates[index].biological_group = group_id
        if contains_external:
            for index in training_indices:
                candidates[index].split = "excluded"
                candidates[index].exclusion_reason = "homologous_to_external_manifest"
        elif training_indices:
            eligible_components.append((group_id, training_indices))

    fractions = {
        "train": args.train_fraction,
        "validation": args.validation_fraction,
        "test": round(1.0 - args.train_fraction - args.validation_fraction, 12),
    }
    if min(fractions.values()) <= 0:
        raise ValueError("All train/validation/test fractions must be positive")
    assignments = _assign_splits(eligible_components, candidates, fractions, args.seed)
    for group_id, indices in eligible_components:
        for index in indices:
            candidates[index].split = assignments[group_id]

    for candidate in candidates:
        if candidate.role == "external" and not candidate.split:
            candidate.split = "external"
        if candidate.role == "training" and candidate.exclusion_reason and not candidate.split:
            candidate.split = "excluded"

    version = args.version
    for split in SPLIT_NAMES:
        write_jsonl(
            str(output_dir / (split + ".jsonl")),
            (
                _record_with_split(candidate, version)
                for candidate in candidates
                if candidate.role == "training" and candidate.split == split
            ),
        )
    write_jsonl(
        str(output_dir / "excluded.jsonl"),
        (
            _record_with_split(candidate, version)
            for candidate in candidates
            if candidate.role == "training" and candidate.split == "excluded"
        ),
    )
    _write_master_table(output_dir / "master_cluster_split.tsv", candidates)
    leakage = _validate_no_leakage(candidates)

    split_summary: Dict[str, object] = {}
    for split in (*SPLIT_NAMES, "excluded"):
        selected = [
            candidate
            for candidate in candidates
            if candidate.role == "training" and candidate.split == split
        ]
        split_summary[split] = {
            "records": len(selected),
            "biological_groups": len({candidate.biological_group for candidate in selected}),
            "pdb_ids": len({pdb for candidate in selected for pdb in candidate.pdb_ids}),
            "rmsf_labels": sum(candidate.metric_counts.get("rmsf", 0) for candidate in selected),
            "chi_entropy_labels": sum(
                candidate.metric_counts.get("chi_entropy", 0) for candidate in selected
            ),
            "datasets": {
                dataset: sum(candidate.record.dataset == dataset for candidate in selected)
                for dataset in sorted({candidate.record.dataset for candidate in selected})
            },
        }
    summary: Dict[str, object] = {
        "version": version,
        "purpose": "absolute RMSF and chi-entropy training only",
        "pair_manifests_used": False,
        "sasa_used_for_training": False,
        "seed": args.seed,
        "max_residues": args.max_residues,
        "mmseqs": {
            "version": mmseqs_version,
            "minimum_sequence_identity": args.minimum_identity,
            "minimum_bidirectional_coverage": args.minimum_coverage,
            "sensitivity": 7.5,
            "cluster_mode": 0,
        },
        "source_manifests": [
            {"path": str(Path(path).resolve()), "sha256": _file_sha256(Path(path))}
            for path in args.manifest
        ],
        "external_manifests": [
            {"path": str(Path(path).resolve()), "sha256": _file_sha256(Path(path))}
            for path in args.external_manifest
        ],
        "split_fractions": fractions,
        "splits": split_summary,
        "external_records": len(external),
        "leakage_validation": leakage,
    }
    with (output_dir / "summary.json").open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(summary, handle, indent=2, ensure_ascii=False)
        handle.write("\n")
    return summary


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest",
        action="append",
        default=[],
        help="absolute state manifest; repeat for multiple sources",
    )
    parser.add_argument(
        "--external-manifest",
        action="append",
        default=[],
        help="absolute-state external benchmark used only for homology exclusion",
    )
    parser.add_argument(
        "--output-dir",
        default="data/splits/state_rmsf_hchi_v1",
    )
    parser.add_argument(
        "--mmseqs",
        default="artifacts/tools/mmseqs/bin/mmseqs",
        help="MMseqs2 executable; a Linux binary is invoked through WSL on Windows",
    )
    parser.add_argument("--path-map", action="append", default=[])
    parser.add_argument("--max-residues", type=int, default=1024)
    parser.add_argument("--minimum-identity", type=float, default=0.30)
    parser.add_argument("--minimum-coverage", type=float, default=0.70)
    parser.add_argument("--train-fraction", type=float, default=0.80)
    parser.add_argument("--validation-fraction", type=float, default=0.10)
    parser.add_argument("--seed", type=int, default=7)
    parser.add_argument("--threads", type=int, default=max(1, min(os.cpu_count() or 1, 16)))
    parser.add_argument("--version", default="state_rmsf_hchi_v1")
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    if not args.manifest:
        parser.error("at least one --manifest is required")
    print(json.dumps(prepare_splits(args), indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
