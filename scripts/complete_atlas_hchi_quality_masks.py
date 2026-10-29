"""Complete the two quarantined ATLAS systems with explicit residue-level masks.

The 1,936 strict successes are not recomputed or altered. The entropy formula,
bins, frames, and all other geometry checks are unchanged. Only the diagnosed
LYS148 in 4g6d_B and 4g94_B is excluded from Hchi supervision.
"""
from collections import Counter
from datetime import datetime, timezone
import argparse
import json
from pathlib import Path
import shutil
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import backfill_atlas_hchi as core
import backfill_atlas_hchi_buffered as buffered

EXCLUDED = {sid: {("_", "148", "LYS")} for sid in ("4g6d_B", "4g94_B")}
REASON = "CE-NZ geometry abnormal in the source PDB and all three trajectories; not resolved by minimum-image distances"


def extract_quality_masked_replica(pdb_path, trajectory_path):
    """Core extraction with one explicitly identified residue skipped at selection."""
    import MDAnalysis as mda
    excluded = EXCLUDED[Path(pdb_path).stem]
    universe = mda.Universe(str(pdb_path), str(trajectory_path))
    quadruples, residue_info, skipped = [], [], []
    alias_count = 0
    try:
        for residue in universe.residues:
            name = str(residue.resname).upper()
            definitions = core.CHI_ATOMS.get(name, [])
            if not definitions or name not in core.THREE_TO_ONE:
                continue
            lookup = {}
            for atom in residue.atoms:
                if str(getattr(atom, "altLoc", "")).strip() not in ("", "A", "1"):
                    continue
                lookup.setdefault(str(atom.name).upper(), int(atom.index))
            if name == "ILE" and "CD1" not in lookup and "CD" in lookup:
                lookup["CD1"] = lookup["CD"]
                alias_count += 1
            if not all(atom in lookup for definition in definitions for atom in definition):
                continue
            ca = universe.atoms[lookup["CA"]]
            chain = str(getattr(ca, "chainID", "")).strip() or "_"
            insertion = str(getattr(ca, "icode", "")).strip()
            identity = (chain, str(ca.resid) + insertion, name)
            if identity in excluded:
                skipped.append(identity)
                continue
            local = []
            for definition in definitions:
                local.append(len(quadruples))
                quadruples.append([lookup[atom] for atom in definition])
            key = {"chain": chain, "residue_number": str(ca.resid), "insertion_code": insertion,
                   "amino_acid_3": name, "n_chi_angles": len(local)}
            residue_info.append((key, local))
        if len(skipped) != len(excluded) or set(skipped) != excluded:
            raise ValueError("Quality-mask residue identity is missing or ambiguous")
        if not quadruples:
            raise ValueError("No complete side-chain dihedrals")
        quad = np.asarray(quadruples, dtype=np.int64)
        frames, min_bond, max_bond = [], float("inf"), 0.0
        first_time = last_time = None
        dt_ps = float(universe.trajectory.dt)
        for ts in universe.trajectory:
            positions = ts.positions[quad]
            bonds = np.linalg.norm(np.diff(positions, axis=1), axis=-1)
            min_bond = min(min_bond, float(bonds.min()))
            max_bond = max(max_bond, float(bonds.max()))
            if min_bond < 0.5 or max_bond > 3.0:
                raise ValueError("Unmasked side-chain geometry remains invalid")
            angles = core.dihedral_degrees(positions)
            if not np.isfinite(angles).all():
                raise ValueError("Nonfinite side-chain angles")
            frames.append(np.clip(np.floor((angles + 180.0) / 120.0).astype(np.int8), 0, 2))
            if first_time is None:
                first_time = float(ts.time)
            last_time = float(ts.time)
    finally:
        universe.trajectory.close()
    states = np.stack(frames, axis=0)
    rows = []
    for key, local in residue_info:
        codes = np.zeros(len(states), dtype=np.int32)
        for power, column in enumerate(local):
            codes += states[:, column].astype(np.int32) * (3 ** power)
        rows.append(dict(key, chi_entropy=core.entropy(codes), n_frames=len(states)))
    qa = {"frames": len(states), "dt_ps": dt_ps, "start_ps": first_time, "end_ps": last_time,
          "duration_ns": (last_time - first_time) / 1000.0, "ile_cd_aliases": alias_count,
          "min_chi_bond_angstrom": min_bond, "max_chi_bond_angstrom": max_bond,
          "quality_masked_residues": [list(key) for key in sorted(excluded)],
          "quality_mask_reason": REASON}
    return rows, qa


def complete(release):
    strict_path = release / "strict_extraction_summary.json"
    if not strict_path.exists():
        summary = json.loads((release / "summary.json").read_text(encoding="utf-8"))
        if summary["done"] != 1938 or summary["complete"] != 1936 or set(summary["failures"]) != set(EXCLUDED):
            raise ValueError("Expected exactly the two diagnosed strict extraction failures")
        shutil.copyfile(release / "summary.json", strict_path)
    strict = json.loads(strict_path.read_text(encoding="utf-8"))
    diagnostic_path = release / "geometry_diagnostic.json"
    diagnostic = json.loads(diagnostic_path.read_text(encoding="utf-8"))
    if {(row["sample_id"], row["replica"]) for row in diagnostic} != {(sid, rep) for sid in EXCLUDED for rep in (1, 2, 3)}:
        raise ValueError("Need all six full-trajectory diagnostics")
    for row in diagnostic:
        if row["frames"] != 1001 or len(row["bad_pairs"]) != 1:
            raise ValueError("Unexpected geometry diagnostic scope")
        pair = row["bad_pairs"][0]
        if [(atom["resid"], atom["resname"], atom["name"]) for atom in pair["atoms"]] != [(148, "LYS", "CE"), (148, "LYS", "NZ")]:
            raise ValueError("Geometry diagnostic does not match the requested mask")
    policy = {"policy_id": "atlas_hchi_two_source_geometry_masks_v1", "reason": REASON,
              "exclusions": {sid: [list(key) for key in sorted(keys)] for sid, keys in EXCLUDED.items()},
              "approved_action": "User continued after proposal to mask only the two abnormal Hchi residues",
              "completion_code_sha256": core.sha256(Path(__file__)),
              "diagnostic_sha256": core.sha256(diagnostic_path),
              "strict_summary_sha256": core.sha256(strict_path),
              "rmsf_sasa_and_fold_membership_unchanged": True}
    core.atomic_json(release / "quality_mask_policy.json", policy)
    mappings = {"/ranze/repos/FlexPPI": str(core.ROOT),
                "/ranze/data/ATLAS_extracted": "E:/data/protein_md_datasets/ATLAS_extracted"}
    original = [core.remap_state_record(row, mappings).to_dict()
                for row in core.load_state_records(strict["source_manifest"])]
    retained_hashes = {path.name: core.sha256(path) for path in (release / "labels").glob("*.npz")
                       if path.stem not in EXCLUDED}
    if len(retained_hashes) != 1936:
        raise ValueError("Expected all 1,936 strict-success labels")
    fingerprint = dict(strict["fingerprint"], quality_mask_code_sha256=policy["completion_code_sha256"],
                       quality_mask_policy_sha256=core.sha256(release / "quality_mask_policy.json"))
    buffered._extract_replica = extract_quality_masked_replica
    for row in original:
        if row["sample_id"] not in EXCLUDED:
            continue
        report = core.process_system((row, str(release), fingerprint))
        report["quality_mask_policy"] = policy["policy_id"]
        report["record"]["metadata"]["quality_masked_residues"] = policy["exclusions"][row["sample_id"]]
        core.atomic_json(release / "systems" / (row["sample_id"] + ".json"), report)
        print(json.dumps({"sample_id": row["sample_id"], "hchi_labels": report["hchi_labels"],
                          "quality_masked_residues": policy["exclusions"][row["sample_id"]]}), flush=True)
    for name, expected in retained_hashes.items():
        if core.sha256(release / "labels" / name) != expected:
            raise ValueError("A strict-success label was altered")
    reports = [json.loads((release / "systems" / (row["sample_id"] + ".json")).read_text(encoding="utf-8")) for row in original]
    for report in reports:
        expected = fingerprint if report["sample_id"] in EXCLUDED else strict["fingerprint"]
        if report["fingerprint"] != expected or core.sha256(report["record"]["label_path"]) != report["output_sha256"]:
            raise ValueError("Release hash/provenance mismatch")
    manifest = release / "atlas_states_hchi.jsonl"
    core.write_jsonl(str(manifest), [report["record"] for report in reports])
    missing, frames = Counter(), Counter()
    for report in reports:
        missing.update(report["missing_hchi_by_resname"])
        frames.update(report["frames_per_replica"])
    summary = dict(strict, complete=1938, failures={}, manifest=str(manifest),
                   completed_utc=datetime.now(timezone.utc).isoformat(),
                   strict_success_labels_reused_unchanged=1936,
                   quality_mask_fingerprints={sid: fingerprint for sid in EXCLUDED},
                   quality_mask_policy=policy, original_strict_failures=strict["failures"],
                   rmsf_labels=sum(row["rmsf_labels"] for row in reports),
                   hchi_labels=sum(row["hchi_labels"] for row in reports),
                   missing_hchi_by_resname=dict(missing), frame_count_histogram=dict(frames))
    core.atomic_json(release / "summary.json", summary)
    core.atomic_json(release / "progress.json", {"done": 1938, "total": 1938, "complete": 1938,
                                                "failures": {}, "quality_masked_residues": 2})
    print(json.dumps({"complete": 1938, "rmsf_labels": summary["rmsf_labels"],
                      "hchi_labels": summary["hchi_labels"], "missing": dict(missing)}), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--label-release-dir", type=Path, default=core.ROOT / "artifacts/data/atlas_rmsf_hchi_v1")
    parser.add_argument("--confirm-mask-lys148", action="store_true", required=True)
    args = parser.parse_args()
    complete(args.label_release_dir.resolve())
