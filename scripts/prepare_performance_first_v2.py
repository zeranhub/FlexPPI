"""Build the performance-first RMSF/H-chi pool and five system-level folds.

Ready and pending-label systems both receive a stable fold.  Only records with
validated RMSF labels are written to train/test manifests.  ITrace is evaluated
as unseen complete TCR-pMHC combinations, so its ubiquitous individual chain
families do not create a transitive mega-component.  Exact PDB identifiers,
complete-complex sequence signatures, and duplicate label files remain
indivisible for every data source.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Set, Tuple

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
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
from scripts.prepare_state_splits import UnionFind, _run_mmseqs  # noqa: E402

VERSION = "state_rmsf_hchi_v2_performance_first"
TRAINING_METRICS = ("rmsf", "chi_entropy")
PDB_PATTERN = re.compile(r"^[0-9][A-Za-z0-9]{3}$")


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _digest(*parts: object) -> str:
    return hashlib.sha256("|".join(str(part) for part in parts).encode("utf-8")).hexdigest()


def _default_path_mappings() -> Dict[str, str]:
    candidates = {
        "/ranze/repos/FlexPPI/artifacts/data": str(ROOT / "artifacts/data"),
        "/ranze/data/ATLAS_extracted": r"E:\data\protein_md_datasets\ATLAS_extracted",
        "/ranze/data/dynarepo/dynarepo_non_nucleic_full": (r"D:\data\dynarepo_non_nucleic_full"),
    }
    return {source: target for source, target in candidates.items() if Path(target).exists()}


def _parse_path_mappings(values: Sequence[str]) -> Dict[str, str]:
    output = {}
    for value in values:
        if "=" not in value:
            raise ValueError("--path-map must use SOURCE=TARGET")
        source, target = value.split("=", 1)
        output[source.strip()] = target.strip()
    return output


def _pdb_ids(record: StateRecord) -> List[str]:
    raw: List[object] = []
    for field_name in ("pdb_id", "pdb_ids"):
        value = record.metadata.get(field_name)
        if isinstance(value, (list, tuple)):
            raw.extend(value)
        elif value:
            raw.append(value)
    for token in (record.sample_id, record.sample_id.split("_", 1)[0]):
        if PDB_PATTERN.fullmatch(token):
            raw.append(token)
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


@dataclass
class Candidate:
    key: str
    record: StateRecord
    original_record: Dict[str, object]
    source_manifest: str
    pdb_ids: List[str] = field(default_factory=list)
    residue_count: int = 0
    chain_sequences: List[Tuple[str, str]] = field(default_factory=list)
    complex_signature: str = ""
    biological_system: str = ""
    metric_counts: Dict[str, int] = field(default_factory=dict)
    label_sha256: str = ""
    label_status: str = ""
    quality_tier: str = ""
    ppi_classification: str = "not_audited"
    sequence_clusters: Set[str] = field(default_factory=set)
    biological_group: str = ""
    fold: Optional[int] = None
    exclusion_reason: str = ""

    @property
    def active(self) -> bool:
        return not self.exclusion_reason and self.label_status in {
            "ready",
            "pending_uniform_labels",
        }

    @property
    def ready(self) -> bool:
        return self.active and self.label_status == "ready"


def _load_pdbbind_audit(path: Optional[Path]) -> Dict[str, str]:
    if path is None or not path.exists():
        return {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        return {
            row["sample_id"]: row["classification"]
            for row in csv.DictReader(handle, delimiter="\t")
        }


def _load_candidates(
    manifests: Sequence[Path],
    mappings: Mapping[str, str],
    max_residues: int,
    pdbbind_audit: Mapping[str, str],
) -> List[Candidate]:
    candidates = []
    for manifest_index, manifest in enumerate(manifests):
        for record_index, original in enumerate(load_state_records(str(manifest))):
            record = remap_state_record(original, mappings)
            candidate = Candidate(
                key="{:02d}:{:06d}:{}".format(manifest_index, record_index, record.sample_id),
                record=record,
                original_record=original.to_dict(),
                source_manifest=str(manifest.resolve()),
                pdb_ids=_pdb_ids(record),
                ppi_classification=pdbbind_audit.get(record.sample_id, "not_audited"),
            )
            try:
                structure_path = Path(record.structure_path)
                if not structure_path.exists():
                    raise FileNotFoundError("missing_structure")
                structure = parse_pdb(str(structure_path))
                candidate.residue_count = structure.num_residues
                if structure.num_residues > max_residues:
                    raise ValueError("length_gt_{}".format(max_residues))
                candidate.chain_sequences = [
                    (chain, structure.chain_sequence(chain))
                    for chain in structure.unique_chains()
                    if len(structure.chain_sequence(chain)) >= 20
                ]
                if not candidate.chain_sequences:
                    raise ValueError("no_chain_with_at_least_20_residues")
                complete_sequences = sorted(sequence for _, sequence in candidate.chain_sequences)
                candidate.complex_signature = "complex_" + _digest(*complete_sequences)[:16]
                explicit_system = str(record.metadata.get("biological_system") or "").strip()
                if explicit_system:
                    candidate.biological_system = explicit_system
                elif candidate.pdb_ids:
                    candidate.biological_system = "pdb:" + "+".join(candidate.pdb_ids)
                else:
                    candidate.biological_system = "sample:{}:{}".format(
                        record.dataset, record.sample_id
                    )

                label_path = Path(record.label_path) if record.label_path else None
                if label_path is None or not label_path.exists():
                    if record.dataset.lower() == "thermopcd":
                        candidate.label_status = "pending_uniform_labels"
                        candidate.quality_tier = "protocol_ready_label_pending"
                    else:
                        raise FileNotFoundError("missing_label")
                else:
                    values, _, mask = load_label_file(str(label_path))
                    if values.shape[0] != structure.num_residues:
                        raise ValueError("structure_label_length_mismatch")
                    for name in TRAINING_METRICS:
                        metric_index = METRIC_NAMES.index(name)
                        valid = mask[:, metric_index] & np.isfinite(values[:, metric_index])
                        candidate.metric_counts[name] = int(valid.sum())
                    if candidate.metric_counts.get("rmsf", 0) == 0:
                        raise ValueError("no_rmsf_labels")
                    candidate.label_sha256 = _sha256_file(label_path)
                    candidate.label_status = "ready"
                    candidate.quality_tier = "training_ready"
            except (OSError, ValueError, KeyError) as exc:
                candidate.exclusion_reason = re.sub(r"\s+", "_", str(exc).strip())[:160]
                candidate.label_status = candidate.label_status or "excluded"
                candidate.quality_tier = "excluded"
            candidates.append(candidate)
    return candidates


def _write_sequences(
    candidates: Sequence[Candidate], fasta_path: Path, index_path: Path
) -> Dict[str, int]:
    sequence_to_candidate = {}
    with (
        fasta_path.open("w", encoding="utf-8", newline="\n") as fasta,
        index_path.open("w", encoding="utf-8", newline="") as index_handle,
    ):
        writer = csv.writer(index_handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sequence_id", "candidate_key", "dataset", "sample_id", "chain", "length"])
        for candidate_index, candidate in enumerate(candidates):
            if not candidate.active:
                continue
            for chain_index, (chain, sequence) in enumerate(candidate.chain_sequences):
                sequence_id = "V2SEQ{:06d}_{:02d}".format(candidate_index, chain_index)
                fasta.write(">{}\n{}\n".format(sequence_id, sequence))
                writer.writerow(
                    [
                        sequence_id,
                        candidate.key,
                        candidate.record.dataset,
                        candidate.record.sample_id,
                        chain,
                        len(sequence),
                    ]
                )
                sequence_to_candidate[sequence_id] = candidate_index
    return sequence_to_candidate


def _apply_clusters(
    cluster_path: Path,
    sequence_to_candidate: Mapping[str, int],
    candidates: Sequence[Candidate],
) -> None:
    with cluster_path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            columns = line.rstrip("\n").split("\t")
            if len(columns) < 2:
                raise ValueError("Malformed MMseqs2 cluster row {}".format(line_number))
            representative, member = columns[:2]
            candidate_index = sequence_to_candidate.get(member)
            if candidate_index is None:
                continue
            candidate = candidates[candidate_index]
            candidate.sequence_clusters.add(representative)


def _union_by_key(
    candidates: Sequence[Candidate], union_find: UnionFind, values: Sequence[Sequence[str]]
) -> None:
    first_by_value: Dict[str, int] = {}
    for index, candidate_values in enumerate(values):
        if not candidates[index].active:
            continue
        for value in candidate_values:
            if not value:
                continue
            if value in first_by_value:
                union_find.union(first_by_value[value], index)
            else:
                first_by_value[value] = index


def _build_groups(candidates: Sequence[Candidate], union_find: UnionFind) -> List[List[int]]:
    _union_by_key(candidates, union_find, [candidate.pdb_ids for candidate in candidates])
    _union_by_key(
        candidates,
        union_find,
        [[candidate.biological_system] for candidate in candidates],
    )
    _union_by_key(
        candidates,
        union_find,
        [[candidate.complex_signature] for candidate in candidates],
    )
    _union_by_key(
        candidates,
        union_find,
        [[candidate.label_sha256] for candidate in candidates],
    )
    # A single-chain protein is indivisible from its directly homologous family.
    # For a complex, an individual common partner (antibody framework, MHC, etc.)
    # is not enough to merge the whole system.  Complexes are joined only when
    # their complete set of chain-family clusters matches.  ITrace uses an even
    # narrower complete-combination definition in the primary CV.
    single_chain_cluster_values = []
    complete_complex_cluster_values = []
    for candidate in candidates:
        is_itrace = candidate.record.dataset.lower() == "itrace"
        if not is_itrace and len(candidate.chain_sequences) == 1:
            single_chain_cluster_values.append(sorted(candidate.sequence_clusters))
        else:
            single_chain_cluster_values.append([])
        if not is_itrace and len(candidate.chain_sequences) > 1:
            signature = "complex_cluster_set:" + ";".join(sorted(candidate.sequence_clusters))
            complete_complex_cluster_values.append([signature])
        else:
            complete_complex_cluster_values.append([])
    _union_by_key(candidates, union_find, single_chain_cluster_values)
    _union_by_key(candidates, union_find, complete_complex_cluster_values)
    components: Dict[int, List[int]] = defaultdict(list)
    for index, candidate in enumerate(candidates):
        if candidate.active:
            components[union_find.find(index)].append(index)
    groups = []
    for indices in components.values():
        group_id = "biogroup_" + _digest(*(sorted(candidates[index].key for index in indices)))[:12]
        for index in indices:
            candidates[index].biological_group = group_id
        groups.append(indices)
    return groups


def _assign_folds(
    groups: Sequence[Sequence[int]], candidates: Sequence[Candidate], folds: int, seed: int
) -> None:
    datasets = sorted({candidate.record.dataset for candidate in candidates if candidate.active})
    target = {
        dataset: sum(
            candidate.record.dataset == dataset for candidate in candidates if candidate.active
        )
        / folds
        for dataset in datasets
    }
    current = [{dataset: 0 for dataset in datasets} for _ in range(folds)]
    current_total = [0] * folds
    active_count = sum(candidate.active for candidate in candidates)
    target_total = active_count / folds
    ordered = sorted(
        groups,
        key=lambda indices: (
            -len(indices),
            _digest(seed, *(sorted(candidates[index].key for index in indices))),
        ),
    )
    for indices in ordered:
        additions = Counter(candidates[index].record.dataset for index in indices)

        def allocation_cost(proposed: int) -> Tuple[float, float, str]:
            dataset_cost = 0.0
            total_cost = 0.0
            for fold in range(folds):
                for dataset in datasets:
                    value = current[fold][dataset]
                    if fold == proposed:
                        value += additions[dataset]
                    dataset_cost += ((value - target[dataset]) / max(target[dataset], 1.0)) ** 2
                fold_total = current_total[fold] + (len(indices) if fold == proposed else 0)
                total_cost += ((fold_total - target_total) / max(target_total, 1.0)) ** 2
            return (
                dataset_cost,
                total_cost,
                _digest(seed, proposed, candidates[indices[0]].biological_group),
            )

        chosen = min(range(folds), key=allocation_cost)
        for index in indices:
            candidates[index].fold = chosen
        for dataset, count in additions.items():
            current[chosen][dataset] += count
        current_total[chosen] += len(indices)


def _output_record(candidate: Candidate) -> Dict[str, object]:
    output = dict(candidate.original_record)
    metadata = dict(output.get("metadata", {}))
    metadata.update(
        {
            "pdb_id": candidate.pdb_ids[0] if candidate.pdb_ids else None,
            "pdb_ids": candidate.pdb_ids,
            "biological_system": candidate.biological_system,
            "cluster_id": candidate.biological_group,
            "complex_signature": candidate.complex_signature,
            "sequence_clusters": sorted(candidate.sequence_clusters),
            "fold": candidate.fold,
            "split_version": VERSION,
            "training_metrics": list(TRAINING_METRICS),
            "label_status": candidate.label_status,
            "quality_tier": candidate.quality_tier,
            "ppi_classification": candidate.ppi_classification,
        }
    )
    if candidate.exclusion_reason:
        metadata["exclusion_reason"] = candidate.exclusion_reason
    output["metadata"] = metadata
    return output


def _write_master(path: Path, candidates: Sequence[Candidate]) -> None:
    columns = [
        "dataset",
        "sample_id",
        "pdb_id",
        "biological_system",
        "source_group",
        "temperature",
        "seed",
        "trajectory_window_ns",
        "residue_count",
        "rmsf_labels",
        "hchi_labels",
        "quality_tier",
        "label_status",
        "sequence_clusters",
        "complex_signature",
        "biological_group",
        "fold",
        "ppi_classification",
        "exclusion_reason",
        "structure_path",
        "label_path",
        "source_manifest",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for candidate in candidates:
            seed = candidate.record.metadata.get(
                "seed_groups", candidate.record.metadata.get("seed", "")
            )
            if isinstance(seed, list):
                seed = ";".join(str(value) for value in seed)
            writer.writerow(
                {
                    "dataset": candidate.record.dataset,
                    "sample_id": candidate.record.sample_id,
                    "pdb_id": ";".join(candidate.pdb_ids),
                    "biological_system": candidate.biological_system,
                    "source_group": candidate.record.source_group,
                    "temperature": candidate.record.temperature,
                    "seed": seed,
                    "trajectory_window_ns": candidate.record.window_ns,
                    "residue_count": candidate.residue_count,
                    "rmsf_labels": candidate.metric_counts.get("rmsf", 0),
                    "hchi_labels": candidate.metric_counts.get("chi_entropy", 0),
                    "quality_tier": candidate.quality_tier,
                    "label_status": candidate.label_status,
                    "sequence_clusters": ";".join(sorted(candidate.sequence_clusters)),
                    "complex_signature": candidate.complex_signature,
                    "biological_group": candidate.biological_group,
                    "fold": "" if candidate.fold is None else candidate.fold,
                    "ppi_classification": candidate.ppi_classification,
                    "exclusion_reason": candidate.exclusion_reason,
                    "structure_path": candidate.original_record.get("structure_path", ""),
                    "label_path": candidate.original_record.get("label_path", ""),
                    "source_manifest": candidate.source_manifest,
                }
            )


def _write_family_inventory(path: Path, candidates: Sequence[Candidate]) -> None:
    cluster_members: Dict[str, List[Candidate]] = defaultdict(list)
    for candidate in candidates:
        if candidate.active:
            for cluster in candidate.sequence_clusters:
                cluster_members[cluster].append(candidate)
    columns = [
        "sequence_cluster",
        "records",
        "ready_records",
        "biological_groups",
        "datasets",
        "single_chain_records",
        "complex_records",
        "itrace_records",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for cluster, members in sorted(cluster_members.items()):
            writer.writerow(
                {
                    "sequence_cluster": cluster,
                    "records": len(members),
                    "ready_records": sum(member.ready for member in members),
                    "biological_groups": len({member.biological_group for member in members}),
                    "datasets": ";".join(sorted({member.record.dataset for member in members})),
                    "single_chain_records": sum(
                        len(member.chain_sequences) == 1 for member in members
                    ),
                    "complex_records": sum(len(member.chain_sequences) > 1 for member in members),
                    "itrace_records": sum(
                        member.record.dataset.lower() == "itrace" for member in members
                    ),
                }
            )


def _write_leave_one_source_out(
    output_dir: Path, candidates: Sequence[Candidate]
) -> Mapping[str, object]:
    ppi_groups = {"dynabench_complex", "test53_complex", "pdbbind_complex"}
    definitions = {
        "atlas": lambda candidate: candidate.record.dataset.lower() == "atlas",
        "dynarepo": lambda candidate: candidate.record.dataset.lower() == "dynarepo",
        "itrace": lambda candidate: candidate.record.dataset.lower() == "itrace",
        "thermopcd": lambda candidate: candidate.record.dataset.lower() == "thermopcd",
        "dynarepo_ppi": lambda candidate: (
            candidate.record.dataset.lower() == "dynarepo"
            and candidate.record.source_group in ppi_groups
            and (
                candidate.record.source_group != "pdbbind_complex"
                or candidate.ppi_classification == "high_confidence_ppi"
            )
        ),
    }
    summary = {}
    for name, predicate in definitions.items():
        test_candidates = [
            candidate for candidate in candidates if candidate.active and predicate(candidate)
        ]
        test_groups = {candidate.biological_group for candidate in test_candidates}
        train_candidates = [
            candidate
            for candidate in candidates
            if candidate.ready
            and not predicate(candidate)
            and candidate.biological_group not in test_groups
        ]
        ready_test = [candidate for candidate in test_candidates if candidate.ready]
        pending_test = [
            candidate
            for candidate in test_candidates
            if candidate.label_status == "pending_uniform_labels"
        ]
        write_jsonl(
            str(output_dir / "loso_{}_train.jsonl".format(name)),
            (_output_record(candidate) for candidate in train_candidates),
        )
        write_jsonl(
            str(output_dir / "loso_{}_test.jsonl".format(name)),
            (_output_record(candidate) for candidate in ready_test),
        )
        if pending_test:
            write_jsonl(
                str(output_dir / "loso_{}_pending_test.jsonl".format(name)),
                (_output_record(candidate) for candidate in pending_test),
            )
        summary[name] = {
            "train_ready_records": len(train_candidates),
            "test_ready_records": len(ready_test),
            "test_pending_records": len(pending_test),
            "overlapping_groups_removed_from_train": sum(
                candidate.ready
                and not predicate(candidate)
                and candidate.biological_group in test_groups
                for candidate in candidates
            ),
        }
    (output_dir / "leave_one_source_out_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    return summary


def _violations(candidates: Sequence[Candidate]) -> Mapping[str, List[str]]:
    mappings: Dict[str, Dict[str, Set[int]]] = {
        "biological_group": defaultdict(set),
        "pdb_id": defaultdict(set),
        "complex_signature": defaultdict(set),
        "label_sha256": defaultdict(set),
        "single_chain_sequence_cluster": defaultdict(set),
        "complete_complex_cluster_signature": defaultdict(set),
    }
    for candidate in candidates:
        if not candidate.active or candidate.fold is None:
            continue
        mappings["biological_group"][candidate.biological_group].add(candidate.fold)
        for pdb_id in candidate.pdb_ids:
            mappings["pdb_id"][pdb_id].add(candidate.fold)
        mappings["complex_signature"][candidate.complex_signature].add(candidate.fold)
        if candidate.label_sha256:
            mappings["label_sha256"][candidate.label_sha256].add(candidate.fold)
        if candidate.record.dataset.lower() != "itrace":
            if len(candidate.chain_sequences) == 1:
                for cluster in candidate.sequence_clusters:
                    mappings["single_chain_sequence_cluster"][cluster].add(candidate.fold)
            elif len(candidate.chain_sequences) > 1:
                cluster_signature = ";".join(sorted(candidate.sequence_clusters))
                mappings["complete_complex_cluster_signature"][cluster_signature].add(
                    candidate.fold
                )
    return {
        name: sorted(value for value, folds in mapping.items() if len(folds) > 1)
        for name, mapping in mappings.items()
    }


def _summarise(
    candidates: Sequence[Candidate], manifests: Sequence[Path], args: argparse.Namespace
) -> Mapping[str, object]:
    fold_rows = {}
    for fold in range(args.folds):
        fold_candidates = [candidate for candidate in candidates if candidate.fold == fold]
        fold_rows[str(fold)] = {
            "ready_records": sum(candidate.ready for candidate in fold_candidates),
            "pending_records": sum(
                candidate.label_status == "pending_uniform_labels" for candidate in fold_candidates
            ),
            "biological_groups": len({candidate.biological_group for candidate in fold_candidates}),
            "datasets": dict(
                sorted(Counter(candidate.record.dataset for candidate in fold_candidates).items())
            ),
            "rmsf_labels": sum(
                candidate.metric_counts.get("rmsf", 0) for candidate in fold_candidates
            ),
            "hchi_labels": sum(
                candidate.metric_counts.get("chi_entropy", 0) for candidate in fold_candidates
            ),
        }
    violations = _violations(candidates)
    return {
        "version": VERSION,
        "purpose": "structure-only absolute RMSF and chi-entropy performance-first training",
        "seed": args.seed,
        "folds": args.folds,
        "max_residues": args.max_residues,
        "primary_cv_protocol": (
            "each fold is test once and belongs to training in the other four runs"
        ),
        "validation_protocol": (
            "fixed epoch budget; no outer-test fold used for checkpoint selection"
        ),
        "itrace_protocol": (
            "complete-combination generalization; chain-family clusters are descriptive "
            "only and do "
            "not join ITrace records"
        ),
        "other_complex_protocol": (
            "single-chain systems are grouped by chain homology; multi-chain systems are "
            "grouped by "
            "their complete set of chain-family clusters, not by a shared common partner alone"
        ),
        "mmseqs": {
            "version": getattr(args, "mmseqs_version", "unknown"),
            "minimum_sequence_identity": args.minimum_identity,
            "minimum_bidirectional_coverage": args.minimum_coverage,
            "cluster_mode": 0,
            "sensitivity": 7.5,
        },
        "source_manifests": [
            {"path": str(path.resolve()), "sha256": _sha256_file(path)} for path in manifests
        ],
        "records": {
            "total": len(candidates),
            "ready": sum(candidate.ready for candidate in candidates),
            "pending_uniform_labels": sum(
                candidate.label_status == "pending_uniform_labels" for candidate in candidates
            ),
            "excluded": sum(bool(candidate.exclusion_reason) for candidate in candidates),
        },
        "ready_datasets": dict(
            sorted(
                Counter(
                    candidate.record.dataset for candidate in candidates if candidate.ready
                ).items()
            )
        ),
        "fold_summary": fold_rows,
        "leakage_validation": {
            "passed": not any(violations.values()),
            "violations": violations,
            "itrace_sequence_cluster_scope": "exempt_from_primary_combination_cv_only",
        },
    }


def prepare(args: argparse.Namespace) -> Mapping[str, object]:
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    manifests = [Path(path) for path in args.manifest]
    mappings = _default_path_mappings()
    mappings.update(_parse_path_mappings(args.path_map))
    audit = _load_pdbbind_audit(Path(args.pdbbind_audit) if args.pdbbind_audit else None)
    candidates = _load_candidates(manifests, mappings, args.max_residues, audit)

    fasta_path = output_dir / "all_chains.fasta"
    sequence_index = output_dir / "sequence_index.tsv"
    sequence_to_candidate = _write_sequences(candidates, fasta_path, sequence_index)
    cluster_path, mmseqs_version = _run_mmseqs(
        fasta_path,
        output_dir,
        Path(args.mmseqs),
        args.minimum_identity,
        args.minimum_coverage,
        args.threads,
    )
    args.mmseqs_version = mmseqs_version
    union_find = UnionFind(len(candidates))
    _apply_clusters(cluster_path, sequence_to_candidate, candidates)
    groups = _build_groups(candidates, union_find)
    _assign_folds(groups, candidates, args.folds, args.seed)

    _write_master(output_dir / "master_manifest.tsv", candidates)
    _write_family_inventory(output_dir / "family_cluster_inventory.tsv", candidates)
    write_jsonl(
        str(output_dir / "all_ready.jsonl"),
        (_output_record(candidate) for candidate in candidates if candidate.ready),
    )
    write_jsonl(
        str(output_dir / "pending_uniform_labels.jsonl"),
        (
            _output_record(candidate)
            for candidate in candidates
            if candidate.label_status == "pending_uniform_labels"
        ),
    )
    write_jsonl(
        str(output_dir / "excluded.jsonl"),
        (_output_record(candidate) for candidate in candidates if candidate.exclusion_reason),
    )
    for test_fold in range(args.folds):
        write_jsonl(
            str(output_dir / "fold_{}_train.jsonl".format(test_fold)),
            (
                _output_record(candidate)
                for candidate in candidates
                if candidate.ready and candidate.fold != test_fold
            ),
        )
        write_jsonl(
            str(output_dir / "fold_{}_test.jsonl".format(test_fold)),
            (
                _output_record(candidate)
                for candidate in candidates
                if candidate.ready and candidate.fold == test_fold
            ),
        )

    _write_leave_one_source_out(output_dir, candidates)

    summary = _summarise(candidates, manifests, args)
    if not summary["leakage_validation"]["passed"]:
        raise RuntimeError("Leakage validation failed: {}".format(summary["leakage_validation"]))
    (output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    return summary


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", action="append", required=True)
    parser.add_argument("--output-dir", default=str(ROOT / "data/splits" / VERSION))
    parser.add_argument("--pdbbind-audit", default=None)
    parser.add_argument("--mmseqs", required=True)
    parser.add_argument("--path-map", action="append", default=[])
    parser.add_argument("--max-residues", type=int, default=1024)
    parser.add_argument("--minimum-identity", type=float, default=0.30)
    parser.add_argument("--minimum-coverage", type=float, default=0.70)
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--seed", type=int, default=7)
    parser.add_argument("--threads", type=int, default=8)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    summary = prepare(args)
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
