"""Create training-ready ITrace label files with RMSF and H-chi targets.

The diagnostic ITrace NPZ files contain RMSF only.  This script copies those
labels into a new, versioned directory and fills the ``chi_entropy`` column
from the three-replica consensus table produced by
``scripts/extract_itrace_hchi.py`` and ``scripts/build_itrace_hchi_benchmark.py``.
The diagnostic labels are never modified in place.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path
from typing import Dict, Mapping, Tuple

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from flexppi.constants import METRIC_NAMES  # noqa: E402
from flexppi.data.labels import load_label_file, save_label_file  # noqa: E402
from flexppi.data.records import load_state_records, write_jsonl  # noqa: E402
from flexppi.data.structures import parse_pdb  # noqa: E402

ResidueKey = Tuple[str, str, str, str]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_hchi_consensus(path: Path) -> Dict[ResidueKey, Tuple[float, float, int]]:
    """Return PDB/chain/author-id/residue-name keyed H-chi statistics."""

    output: Dict[ResidueKey, Tuple[float, float, int]] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        required = {
            "pdb_id",
            "chain",
            "residue_number",
            "insertion_code",
            "amino_acid_3",
            "target_chi_entropy",
            "target_chi_entropy_std",
            "target_chi_replica_count",
        }
        missing = sorted(required - set(reader.fieldnames or []))
        if missing:
            raise ValueError("Consensus CSV is missing columns: {}".format(", ".join(missing)))
        for row in reader:
            value_text = row["target_chi_entropy"].strip()
            if not value_text:
                continue
            author_id = row["residue_number"].strip() + row["insertion_code"].strip()
            key = (
                row["pdb_id"].strip().lower(),
                row["chain"].strip() or "_",
                author_id,
                row["amino_acid_3"].strip().upper(),
            )
            if key in output:
                raise ValueError("Duplicate ITrace H-chi residue key: {}".format(key))
            standard_deviation = float(row["target_chi_entropy_std"] or 0.0)
            replica_count = int(float(row["target_chi_replica_count"] or 0.0))
            output[key] = (float(value_text), standard_deviation**2, replica_count)
    return output


def integrate_labels(
    manifest_path: Path,
    consensus_path: Path,
    output_label_dir: Path,
    output_manifest: Path,
    report_path: Path,
) -> Mapping[str, object]:
    consensus = load_hchi_consensus(consensus_path)
    output_label_dir.mkdir(parents=True, exist_ok=True)
    records = load_state_records(str(manifest_path))
    updated_records = []
    systems = []
    used_keys = set()
    total_rmsf = 0
    total_hchi = 0
    hchi_index = METRIC_NAMES.index("chi_entropy")

    for record in records:
        pdb_id = str(record.metadata.get("pdb_id") or record.sample_id).lower()
        structure = parse_pdb(record.structure_path)
        if not record.label_path:
            raise ValueError("ITrace record has no RMSF label path: {}".format(record.sample_id))
        values, variances, mask = load_label_file(record.label_path)
        if values.shape[0] != structure.num_residues:
            raise ValueError(
                "Structure/label length mismatch for {}: {} versus {}".format(
                    record.sample_id, structure.num_residues, values.shape[0]
                )
            )

        local_hchi = 0
        replica_count_histogram: Dict[str, int] = {}
        for index, (chain, author_id, residue_name) in enumerate(
            zip(structure.chain_ids, structure.author_residue_ids, structure.residue_names)
        ):
            key = (pdb_id, chain, author_id, residue_name)
            statistics = consensus.get(key)
            if statistics is None:
                continue
            value, variance, replica_count = statistics
            values[index, hchi_index] = value
            variances[index, hchi_index] = variance
            mask[index, hchi_index] = True
            local_hchi += 1
            used_keys.add(key)
            histogram_key = str(replica_count)
            replica_count_histogram[histogram_key] = (
                replica_count_histogram.get(histogram_key, 0) + 1
            )

        if local_hchi == 0:
            raise ValueError("No H-chi residues mapped for {}".format(record.sample_id))
        output_label = output_label_dir / (record.sample_id + ".npz")
        save_label_file(str(output_label), values, variances, mask)
        copied = record.to_dict()
        copied["label_path"] = str(output_label.resolve())
        copied["source_group"] = "itrace_complex"
        metadata = dict(copied.get("metadata", {}))
        metadata.update(
            {
                "analysis_only": False,
                "historical_diagnostic_use": True,
                "hchi_labels": local_hchi,
                "hchi_replica_count_histogram": replica_count_histogram,
                "label_protocol": "rmsf_existing_chainfit+hchi_120deg_joint_rotamer_3replica",
                "training_role": "performance_first_pool",
            }
        )
        copied["metadata"] = metadata
        updated_records.append(copied)
        rmsf_count = int(mask[:, METRIC_NAMES.index("rmsf")].sum())
        total_rmsf += rmsf_count
        total_hchi += local_hchi
        systems.append(
            {
                "sample_id": record.sample_id,
                "pdb_id": pdb_id,
                "residues": structure.num_residues,
                "rmsf_labels": rmsf_count,
                "hchi_labels": local_hchi,
                "output_label_sha256": _sha256(output_label),
            }
        )

    write_jsonl(str(output_manifest), updated_records)
    unused_by_pdb: Dict[str, int] = {}
    for key in set(consensus) - used_keys:
        unused_by_pdb[key[0]] = unused_by_pdb.get(key[0], 0) + 1
    report: Mapping[str, object] = {
        "version": "itrace_rmsf_hchi_training_v2",
        "source_manifest": str(manifest_path.resolve()),
        "source_manifest_sha256": _sha256(manifest_path),
        "source_consensus": str(consensus_path.resolve()),
        "source_consensus_sha256": _sha256(consensus_path),
        "output_manifest": str(output_manifest.resolve()),
        "systems": len(updated_records),
        "rmsf_labels": total_rmsf,
        "hchi_labels": total_hchi,
        "consensus_hchi_rows": len(consensus),
        "unused_consensus_rows_by_pdb": dict(sorted(unused_by_pdb.items())),
        "system_details": systems,
    }
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    return report


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest",
        default=str(ROOT / "artifacts/data/itrace_context_diagnostic/itrace_states.jsonl"),
    )
    parser.add_argument(
        "--consensus",
        default=str(ROOT / "output/itrace_hchi_benchmark/raw/hchi_consensus_per_residue.csv"),
    )
    parser.add_argument(
        "--output-label-dir",
        default=str(ROOT / "artifacts/data/itrace_performance_first_labels"),
    )
    parser.add_argument(
        "--output-manifest",
        default=str(ROOT / "artifacts/data/itrace_performance_first_states.jsonl"),
    )
    parser.add_argument(
        "--report",
        default=str(ROOT / "artifacts/data/itrace_performance_first_integration.json"),
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    report = integrate_labels(
        Path(args.manifest),
        Path(args.consensus),
        Path(args.output_label_dir),
        Path(args.output_manifest),
        Path(args.report),
    )
    print(
        json.dumps(
            {key: report[key] for key in ("systems", "rmsf_labels", "hchi_labels")}, indent=2
        )
    )


if __name__ == "__main__":
    main()
