"""Audit DynaRepo ``pdbbind_complex`` records for protein-protein content.

PDBbind provenance alone does not prove that a record is a PPI.  This audit
uses only protein residues in the supplied structure and reports a conservative
geometric classification.  A high-confidence PPI needs two chains of at least
20 residues and a cross-chain C-alpha interface with at least five contacts
and at least three residues on each side within 8 Angstrom.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Dict, List, Mapping, Tuple

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from flexppi.data.records import load_state_records, remap_path  # noqa: E402
from flexppi.data.structures import ProteinStructure, parse_pdb  # noqa: E402


def _default_mappings() -> Dict[str, str]:
    mappings = {
        "/ranze/data/dynarepo/dynarepo_non_nucleic_full": (r"D:\data\dynarepo_non_nucleic_full")
    }
    return {source: target for source, target in mappings.items() if Path(target).exists()}


def _chain_contact_summary(
    structure: ProteinStructure,
    minimum_chain_residues: int,
    cutoff: float,
) -> Tuple[List[str], List[Mapping[str, object]]]:
    chain_indices = {
        chain: structure.indices_for_chain(chain) for chain in structure.unique_chains()
    }
    qualifying = sorted(
        chain for chain, indices in chain_indices.items() if len(indices) >= minimum_chain_residues
    )
    pairs: List[Mapping[str, object]] = []
    cutoff2 = cutoff**2
    for first_position, first_chain in enumerate(qualifying):
        first_indices = chain_indices[first_chain]
        first_indices = first_indices[structure.backbone_mask[first_indices, 1]]
        first_ca = structure.backbone[first_indices, 1]
        for second_chain in qualifying[first_position + 1 :]:
            second_indices = chain_indices[second_chain]
            second_indices = second_indices[structure.backbone_mask[second_indices, 1]]
            second_ca = structure.backbone[second_indices, 1]
            if len(first_ca) == 0 or len(second_ca) == 0:
                continue
            distance2 = np.sum(
                (first_ca[:, None, :] - second_ca[None, :, :]) ** 2,
                axis=-1,
            )
            contact = distance2 <= cutoff2
            contact_pairs = int(contact.sum())
            first_interface = int(np.any(contact, axis=1).sum())
            second_interface = int(np.any(contact, axis=0).sum())
            if contact_pairs:
                pairs.append(
                    {
                        "chains": "{}:{}".format(first_chain, second_chain),
                        "ca_contact_pairs": contact_pairs,
                        "first_interface_residues": first_interface,
                        "second_interface_residues": second_interface,
                        "minimum_ca_distance": float(np.sqrt(distance2.min())),
                    }
                )
    return qualifying, pairs


def classify_structure(
    structure: ProteinStructure,
    minimum_chain_residues: int = 20,
    cutoff: float = 8.0,
    minimum_contacts: int = 5,
    minimum_interface_residues: int = 3,
) -> Mapping[str, object]:
    chain_lengths = {
        chain: int(len(structure.indices_for_chain(chain))) for chain in structure.unique_chains()
    }
    qualifying, pairs = _chain_contact_summary(structure, minimum_chain_residues, cutoff)
    accepted_pairs = [
        pair
        for pair in pairs
        if int(pair["ca_contact_pairs"]) >= minimum_contacts
        and int(pair["first_interface_residues"]) >= minimum_interface_residues
        and int(pair["second_interface_residues"]) >= minimum_interface_residues
    ]
    if len(qualifying) < 2:
        classification = "not_ppi_single_protein_chain"
    elif accepted_pairs:
        classification = "high_confidence_ppi"
    else:
        classification = "needs_manual_review"
    return {
        "classification": classification,
        "chain_lengths": chain_lengths,
        "qualifying_chains": qualifying,
        "contacting_chain_pairs": pairs,
        "accepted_chain_pairs": accepted_pairs,
    }


def run_audit(
    manifest: Path,
    output_table: Path,
    output_summary: Path,
    minimum_chain_residues: int = 20,
    cutoff: float = 8.0,
    minimum_contacts: int = 5,
    minimum_interface_residues: int = 3,
) -> Mapping[str, object]:
    mappings = _default_mappings()
    rows = []
    counts: Dict[str, int] = {}
    for record in load_state_records(str(manifest)):
        if record.source_group != "pdbbind_complex":
            continue
        structure_path = remap_path(record.structure_path, mappings)
        try:
            structure = parse_pdb(structure_path)
            result = classify_structure(
                structure,
                minimum_chain_residues,
                cutoff,
                minimum_contacts,
                minimum_interface_residues,
            )
            error = ""
        except (OSError, ValueError) as exc:
            result = {
                "classification": "audit_error",
                "chain_lengths": {},
                "qualifying_chains": [],
                "contacting_chain_pairs": [],
                "accepted_chain_pairs": [],
            }
            error = str(exc)
        classification = str(result["classification"])
        counts[classification] = counts.get(classification, 0) + 1
        pdb_ids = [str(value).upper() for value in record.metadata.get("pdb_ids", [])]
        accepted = result["accepted_chain_pairs"]
        all_pairs = result["contacting_chain_pairs"]
        rows.append(
            {
                "sample_id": record.sample_id,
                "pdb_ids": ";".join(pdb_ids),
                "residue_count": sum(result["chain_lengths"].values()),
                "protein_chain_count": len(result["chain_lengths"]),
                "qualifying_chain_count": len(result["qualifying_chains"]),
                "qualifying_chains": ";".join(result["qualifying_chains"]),
                "accepted_chain_pairs": ";".join(str(pair["chains"]) for pair in accepted),
                "maximum_ca_contact_pairs": max(
                    [int(pair["ca_contact_pairs"]) for pair in all_pairs], default=0
                ),
                "classification": classification,
                "manual_review_required": classification == "needs_manual_review",
                "error": error,
                "structure_path": structure_path,
            }
        )

    output_table.parent.mkdir(parents=True, exist_ok=True)
    columns = list(rows[0]) if rows else []
    with output_table.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    summary: Mapping[str, object] = {
        "records": len(rows),
        "classification_counts": dict(sorted(counts.items())),
        "criteria": {
            "minimum_chain_residues": minimum_chain_residues,
            "ca_contact_cutoff_angstrom": cutoff,
            "minimum_ca_contact_pairs": minimum_contacts,
            "minimum_interface_residues_per_chain": minimum_interface_residues,
        },
        "interpretation": (
            "Automated protein-chain geometry audit. High-confidence records may be called "
            "PPI for inventory purposes; needs_manual_review records require curation before "
            "supporting a PPI-count claim. All categories may still provide absolute protein "
            "dynamics labels."
        ),
    }
    output_summary.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", default=str(ROOT / "artifacts/data/dynarepo_states.jsonl"))
    parser.add_argument(
        "--output-table",
        default=str(
            ROOT / "data/splits/state_rmsf_hchi_v2_performance_first/pdbbind_ppi_audit.tsv"
        ),
    )
    parser.add_argument(
        "--output-summary",
        default=str(
            ROOT / "data/splits/state_rmsf_hchi_v2_performance_first/pdbbind_ppi_audit_summary.json"
        ),
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    summary = run_audit(Path(args.manifest), Path(args.output_table), Path(args.output_summary))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
