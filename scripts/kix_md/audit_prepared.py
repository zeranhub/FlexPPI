"""Audit downloaded prepared proteins against frozen model inputs before MD."""
import hashlib
import json
from pathlib import Path
import numpy as np

BASE = Path(__file__).resolve().parents[2] / "artifacts/kix_md_20260914"


def atoms(path):
    return {(line[21], line[22:27].strip(), line[12:16].strip()):
            np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            for line in path.read_text().splitlines() if line.startswith("ATOM")
            and line[76:78].strip() not in ("H", "D")}


def main():
    results = {}
    for state in ("apo", "myb", "mll", "both"):
        source = atoms(BASE / "inputs" / f"{state}.pdb")
        path = BASE / "remote_audit" / f"{state}_protein.pdb"
        prepared = atoms(path)
        assert set(source).issubset(prepared)
        changes = np.stack([prepared[key] - value for key, value in source.items()])
        translation = changes.mean(axis=0)
        errors = np.linalg.norm(changes - translation, axis=1)
        assert errors.max() < 0.002, "Existing heavy atoms changed beyond PDB rounding/translation"
        meta = json.loads((BASE / "remote_audit" / f"{state}_metadata.json").read_text())
        assert hashlib.sha256(path.read_bytes()).hexdigest() == meta["protein_pdb_sha256"]
        assert abs(meta["net_charge_e"]) < 1e-4
        histidines = [r for r in meta["residues"] if r["name"] == "HIS"]
        assert {int(r["author_resid"]) for r in histidines} == {592, 594, 602, 605, 651}
        assert all(r["variant"] == "HID" for r in histidines)
        assert meta["initial_periodic_heavy_gap_lower_bound_nm"] >= 2
        nonbonded = [f for f in meta["nonbonded_forces"] if f["type"] in ("NonbondedForce", "CustomNonbondedForce")]
        assert len(nonbonded) == 2
        assert all(f["getUseSwitchingFunction"] == "True" and f["getCutoffDistance"] == "1.2 nm"
                   and f["getSwitchingDistance"] == "1.0 nm" for f in nonbonded)
        results[state] = dict(passed=True, matched_source_heavy_atoms=len(source),
                             added_heavy_atoms=sorted(set(prepared) - set(source)),
                             translation_angstrom=translation.tolist(),
                             max_heavy_coordinate_deviation_after_translation_angstrom=float(errors.max()),
                             solvent_atoms=meta["total_atoms"] - meta["protein_atom_count"],
                             total_atoms=meta["total_atoms"], histidines="five HID",
                             input_matches_frozen_prediction_backbone=True)
    output = BASE / "remote_audit" / "prepared_input_audit.json"
    output.write_text(json.dumps(results, indent=2) + "\n")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
