"""Prepare the complete, label-blind template x candidate cross product only.

No checkpoint is loaded and no model predictions are generated. The source CSV
is read only for candidate identities and sequences; its closest-template columns
are deliberately not consulted. All listed sensitivity templates receive every
candidate in their own panel, including the common panel parent.
"""
from __future__ import annotations

from collections import Counter
import csv
from datetime import datetime, timezone
import json
from pathlib import Path

from local_cpu_inference import import_runtime, prepare_fixed_backbone, sha256, write_json

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "structures"
EXPECTED_COUNTS = {"1G4_NYESO": 14, "A6_TAX": 16}
PANEL_TO_SYSTEM = {"1G4_NYESO": "nyeso", "A6_TAX": "tax"}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    runtime = import_runtime()
    np, parse_pdb = runtime["np"], runtime["parse_pdb"]
    label_path = ROOT / "experimental_labels" / "SPR_candidates.csv"
    templates_path = ROOT / "experimental_labels" / "templates.json"
    protocol_path = ROOT / "evaluation_review" / "fixed_analysis_protocol.json"
    with label_path.open(encoding="utf-8-sig", newline="") as handle:
        source_candidates = list(csv.DictReader(handle))
    assert Counter(row["panel_id"] for row in source_candidates) == EXPECTED_COUNTS
    assert len({row["candidate_id"] for row in source_candidates}) == 30
    templates = json.loads(templates_path.read_text(encoding="utf-8"))
    protocol = json.loads(protocol_path.read_text(encoding="utf-8"))
    template_by_id = {row["pdb_id"]: row for row in templates}
    systems = protocol["systems"]
    rows = []
    template_inventory = []
    for panel_id in ["1G4_NYESO", "A6_TAX"]:
        system = systems[PANEL_TO_SYSTEM[panel_id]]
        template_ids = [system["primary_template"]] + system["sensitivity_templates"]
        assert template_ids == (["2BNQ", "2BNR"] if panel_id == "1G4_NYESO" else ["1AO7", "1QSE", "1QSF"])
        candidates = [row for row in source_candidates if row["panel_id"] == panel_id]
        assert {row["panel_reference_peptide"] for row in candidates} == {system["parent_peptide"]}
        parent = [row for row in candidates if row["peptide_sequence"] == system["parent_peptide"]]
        assert len(parent) == 1
        for template_id in template_ids:
            supplement = template_id not in template_by_id
            if not supplement:
                t = template_by_id[template_id]
                raw_pdb = Path(t["pdb_path"])
                assert sha256(raw_pdb) == t["pdb_sha256"]
                peptide_chain = t["peptide_chain"]
            else:
                assert template_id == "1QSE", "Unexpected template missing from templates.json"
                raw_pdb = ROOT / "experimental_labels/raw/github/mutation_effects/tax/pdbs/1qse.pdb"
                peptide_chain = "C"
            st = parse_pdb(str(raw_pdb))
            # Molecular roles are verified using complete per-template sequences.
            # Their letter mapping is recorded, rather than assumed in later data.
            assert 270 <= len(st.chain_sequence("A")) <= 280
            assert 95 <= len(st.chain_sequence("B")) <= 105
            assert len(st.chain_sequence(peptide_chain)) == 9
            assert len(st.chain_sequence("D")) > 100 and len(st.chain_sequence("E")) > 100
            assert st.chain_sequence("A").startswith("GSHSMRY")
            assert st.chain_sequence("B").startswith("MIQRTPKI") or st.chain_sequence("B").startswith("IQRTPKI")
            if not supplement:
                assert st.chain_sequence(peptide_chain) == t["peptide_sequence"]
                for chain, meta in t["chain_sequences"].items():
                    assert st.chain_sequence(chain) == meta["sequence"], (template_id, chain)
            peptide_indices = st.indices_for_chain(peptide_chain)
            assert st.frame_mask[peptide_indices].all()
            roles = {"hla_heavy": "A", "beta2m": "B", "peptide": peptide_chain,
                     "tcr_alpha": "D", "tcr_beta": "E"}
            inventory = {"template_pdb": template_id, "panel_id": panel_id,
                "primary": template_id == system["primary_template"],
                "sensitivity": template_id != system["primary_template"],
                "raw_pdb": str(raw_pdb), "raw_pdb_sha256": sha256(raw_pdb),
                "template_peptide_sequence": st.chain_sequence(peptide_chain),
                "peptide_chain": peptide_chain, "roles": roles,
                "chain_sequences": {ch: st.chain_sequence(ch) for ch in st.unique_chains()},
                "chain_residue_counts": {ch: len(st.chain_sequence(ch)) for ch in st.unique_chains()},
                "residues": st.num_residues, "geometry_valid_count": int(st.frame_mask.sum()),
                "supplemented_from_existing_raw_pdb": supplement,
                "supplement_reason": "1QSE is predeclared in protocol and raw source inventory but omitted from templates.json." if supplement else None,
                "backbone_float32_sha256": __import__("hashlib").sha256(st.backbone.tobytes()).hexdigest()}
            template_inventory.append(inventory)
            for candidate in candidates:
                unique_id = template_id + "__" + candidate["candidate_id"]
                folder = OUT / "fixed_backbone" / unique_id
                folder.mkdir(parents=True, exist_ok=True)
                prepared_pdb = folder / "model_input_fixed_backbone.pdb"
                prep = prepare_fixed_backbone(raw_pdb, {peptide_chain: candidate["peptide_sequence"]}, prepared_pdb, runtime)
                after = parse_pdb(str(prepared_pdb))
                assert after.chain_sequence(peptide_chain) == candidate["peptide_sequence"]
                assert np.array_equal(after.backbone, st.backbone)
                for chain in st.unique_chains():
                    if chain != peptide_chain:
                        assert after.chain_sequence(chain) == st.chain_sequence(chain)
                metadata = {"candidate_id": unique_id, "biological_candidate_id": candidate["candidate_id"],
                    "panel_id": panel_id, "template_pdb": template_id,
                    "primary": template_id == system["primary_template"],
                    "sensitivity": template_id != system["primary_template"],
                    "is_primary_template": template_id == system["primary_template"],
                    "is_sensitivity_template": template_id != system["primary_template"],
                    "is_panel_parent": candidate["peptide_sequence"] == system["parent_peptide"],
                    "parent_biological_candidate_id": parent[0]["candidate_id"],
                    "parent_inference_candidate_id": template_id + "__" + parent[0]["candidate_id"],
                    "tcr": candidate["tcr"], "hla": candidate["hla"], "antigen": candidate["antigen"],
                    "peptide_sequence": candidate["peptide_sequence"],
                    "peptide_length": int(candidate["peptide_length"]),
                    "peptide_chain": peptide_chain, "roles": roles,
                    "panel_reference_peptide": system["parent_peptide"],
                    "native_peptide": candidate["native_peptide"],
                    "template_peptide_sequence": st.chain_sequence(peptide_chain),
                    "mutations_vs_panel_reference": candidate["mutations_vs_panel_reference"],
                    "mutations_vs_template": [change["from"] + change["author_residue_id"] + change["to"] for change in prep["mutations"]],
                    "input_mode": "fixed-backbone", "input_pdb": str(raw_pdb),
                    "template_input_pdb": str(raw_pdb), "template_input_sha256": sha256(raw_pdb),
                    "chain_sequences": {peptide_chain: candidate["peptide_sequence"]},
                    "prepared_input_pdb": str(prepared_pdb), "prepared_input_sha256": sha256(prepared_pdb),
                    "prepared_input_residues": after.num_residues,
                    "geometry_valid_count": int(after.frame_mask.sum()),
                    "all_backbone_coordinates_bit_identical_to_template": True,
                    "prepared_input_audit": str(folder / "input_preparation_audit.json"),
                    "labels_or_closest_template_used_for_preparation": False,
                    "predictions_generated": False}
                prep.update({"candidate_id": unique_id, "biological_candidate_id": candidate["candidate_id"],
                    "panel_id": panel_id, "template_pdb_id": template_id,
                    "primary": metadata["primary"], "sensitivity": metadata["sensitivity"],
                    "is_panel_parent": metadata["is_panel_parent"], "prepare_only": True,
                    "model_loaded": False, "predictions_generated": False})
                write_json(folder / "input_preparation_audit.json", prep)
                write_json(folder / "candidate_metadata.json", metadata)
                rows.append(metadata)
            print(json.dumps({"template": template_id, "prepared_candidates": len(candidates),
                "residues": st.num_residues, "peptide": st.chain_sequence(peptide_chain)}), flush=True)
    assert len(rows) == 76
    assert Counter(row["template_pdb"] for row in rows) == {"2BNQ": 14, "2BNR": 14, "1AO7": 16, "1QSE": 16, "1QSF": 16}
    assert len({row["candidate_id"] for row in rows}) == 76
    assert sum(row["is_panel_parent"] for row in rows) == 5
    manifest = {"status": "prepared_not_inferred", "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "prepare_only": True, "model_loaded": False, "predictions_generated": False,
        "candidate_count": len(rows), "biological_candidate_count": 30,
        "primary_input_count": sum(row["primary"] for row in rows),
        "sensitivity_input_count": sum(row["sensitivity"] for row in rows),
        "all_templates_applied_to_all_panel_candidates": True,
        "template_selection": "Fixed primary 2BNQ/1AO7; every predeclared sensitivity template applied to every panel candidate. Closest-template columns are ignored.",
        "recipe_vs_prepared_path": "input_pdb plus chain_sequences is the reproducible inference recipe; prepared_input_pdb is the already-generated sequence-synchronized PDB. Both have identical backbone tensors. Static all-atom sidechain completion is a separate later preparation.",
        "source_files": {str(p): sha256(p) for p in [label_path, templates_path, protocol_path, Path(__file__)]},
        "ignored_candidate_csv_fields": ["template_pdb", "template_peptide_chain", "template_peptide_sequence",
            "template_hamming", "template_selection", "template_pdb_path", "template_pdb_sha256"],
        "mutant_structure_scope": "fixed-backbone sequence substitution; mutated old side chains removed, not rebuilt or relaxed",
        "templates": template_inventory, "candidates": rows}
    write_json(OUT / "templates_preparation_inventory.json", template_inventory)
    write_json(OUT / "fixed_backbone_manifest.json", manifest)
    write_json(OUT / "fixed_backbone_preparation_audit.json", {
        "status": "passed", "manifest_sha256": sha256(OUT / "fixed_backbone_manifest.json"),
        "prepared_inputs": 76, "biological_candidates": 30, "templates": 5,
        "parents_prepared": 5, "closest_template_selection_used": False,
        "all_candidates_have_valid_complete_peptide_backbone": True,
        "all_backbone_coordinates_and_masks_bit_identical_to_respective_template": True,
        "all_nontarget_chain_sequences_unchanged": True,
        "model_loaded": False, "predictions_generated": False,
        "prepared_files": {str(Path(row["prepared_input_pdb"])): row["prepared_input_sha256"] for row in rows}})
    print(json.dumps({"status": "passed", "prepared": len(rows), "manifest": str(OUT / "fixed_backbone_manifest.json")}), flush=True)


if __name__ == "__main__":
    main()
