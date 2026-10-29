"""Validate an entire ATLAS Hchi release against original labels and structures."""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from flexppi.constants import CHI_ATOMS, METRIC_NAMES
from flexppi.data.labels import load_label_file
from flexppi.data.structures import parse_pdb


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def audit(release):
    summary = json.loads((release / "summary.json").read_text(encoding="utf-8"))
    require(summary["scope"] == "full", "Pilot release is not a full dataset")
    require(summary["requested_systems"] == summary["done"] == summary["complete"] == 1938,
            "ATLAS release must contain all 1,938 requested systems")
    require(not summary["failures"], "Extraction failures remain")
    quality_fingerprints = summary.get("quality_mask_fingerprints", {})
    quality_exclusions = summary.get("quality_mask_policy", {}).get("exclusions", {})
    if quality_fingerprints:
        expected_exclusions = {sid: [["_", "148", "LYS"]] for sid in ("4g6d_B", "4g94_B")}
        require(quality_exclusions == expected_exclusions and set(quality_fingerprints) == set(expected_exclusions),
                "Unexpected quality-mask scope")
        policy = summary["quality_mask_policy"]
        require(policy == json.loads((release / "quality_mask_policy.json").read_text(encoding="utf-8")), "Quality policy changed")
        require(policy["diagnostic_sha256"] == digest(release / "geometry_diagnostic.json"), "Diagnostic changed")
        require(policy["strict_summary_sha256"] == digest(release / "strict_extraction_summary.json"), "Strict summary changed")
        for fingerprint in quality_fingerprints.values():
            require(fingerprint == dict(summary["fingerprint"], quality_mask_code_sha256=policy["completion_code_sha256"],
                                        quality_mask_policy_sha256=digest(release / "quality_mask_policy.json")),
                    "Quality-mask extraction fingerprint mismatch")
    records = [json.loads(line) for line in Path(summary["manifest"]).read_text(encoding="utf-8").splitlines() if line]
    require(len(records) == len({row["sample_id"] for row in records}) == 1938, "Manifest coverage mismatch")
    rows, missing, frame_hist, duration_hist = [], Counter(), Counter(), Counter()
    sums = Counter()
    h = METRIC_NAMES.index("chi_entropy")
    for index, record in enumerate(records, 1):
        sid = record["sample_id"]
        report = json.loads((release / "systems" / f"{sid}.json").read_text(encoding="utf-8"))
        require(report["status"] == "complete", sid + ": incomplete report")
        require(report["fingerprint"] == quality_fingerprints.get(sid, summary["fingerprint"]), sid + ": extraction code mismatch")
        path = Path(record["label_path"])
        original_path = Path(record["metadata"]["original_label_path"])
        require(digest(path) == report["output_sha256"], sid + ": changed output file")
        require(digest(original_path) == report["source"]["old_label_sha256"], sid + ": original label changed")
        require(digest(record["structure_path"]) == report["source"]["structure_sha256"], sid + ": structure changed")
        structure = parse_pdb(record["structure_path"])
        old = load_label_file(str(original_path))
        with np.load(path, allow_pickle=False) as data:
            values, variances, mask = data["values"], data["variances"], data["mask"]
            replicas, n_chi, frames = data["hchi_replica_values"], data["hchi_chi_count"], data["hchi_replica_frames"]
            require(tuple(data["metric_names"]) == METRIC_NAMES, sid + ": metric schema mismatch")
            for key, expected in (("chain_ids", structure.chain_ids), ("author_residue_ids", structure.author_residue_ids),
                                  ("residue_names", structure.residue_names)):
                require(np.array_equal(data[key], expected), sid + ": residue identity mismatch: " + key)
        require(values.shape == variances.shape == mask.shape == old[0].shape == (structure.num_residues, 3), sid + ": array shape mismatch")
        require(mask.dtype == np.dtype(bool), sid + ": mask is not boolean")
        for col in (0, 2):
            for before, after in zip(old, (values, variances, mask)):
                require(np.array_equal(before[:, col], after[:, col], equal_nan=True), sid + ": non-Hchi labels changed")
        valid = mask[:, h]
        allowed_missing = {tuple(key) for key in quality_exclusions.get(sid, [])}
        identities = list(zip(structure.chain_ids, structure.author_residue_ids, structure.residue_names))
        for identity, eligible, count in zip(identities, valid, n_chi):
            if identity in allowed_missing:
                require(not eligible and count == 0, sid + ": quality-masked residue was not excluded")
            elif not eligible:
                require(identity[2] in ("ALA", "GLY"), sid + ": unexpected additional missing Hchi residue")
        require(allowed_missing.issubset(set(identities)), sid + ": quality-mask residue missing from structure")
        require(valid.any() and np.isfinite(values[valid, h]).all(), sid + ": invalid targets")
        require(np.isfinite(variances[valid, h]).all() and (variances[valid, h] >= 0).all(), sid + ": invalid variance")
        require(np.isnan(values[~valid, h]).all(), sid + ": masked Hchi is not missing")
        require(replicas.shape == (3, len(values)) and np.isfinite(replicas[:, valid]).all(), sid + ": invalid replica array")
        require(np.isnan(replicas[:, ~valid]).all(), sid + ": replica masks disagree")
        np.testing.assert_allclose(values[valid, h], replicas[:, valid].mean(axis=0), atol=2e-7, rtol=1e-7)
        np.testing.assert_allclose(variances[valid, h], replicas[:, valid].var(axis=0, ddof=1), atol=2e-7, rtol=1e-7)
        require((values[valid, h] >= -1e-7).all() and (values[valid, h] <= n_chi[valid] * np.log(3) + 1e-6).all(), sid + ": Hchi outside state-space bounds")
        require(np.array_equal(frames, [1001, 1001, 1001]), sid + ": unexpected frame counts")
        require(len(report["replica_qa"]) == 3, sid + ": missing replica metadata")
        for qa in report["replica_qa"]:
            require(qa.get("quality_masked_residues", []) == quality_exclusions.get(sid, []),
                    sid + ": replica quality-mask scope mismatch")
            require(qa["frames"] == 1001 and abs(qa["duration_ns"] - 100.0) < 0.001,
                    sid + ": incomplete 100 ns trajectory")
            require(abs(qa["dt_ps"] - 100.0) < 0.01 and abs(qa["start_ps"]) < 0.01,
                    sid + ": unexpected trajectory sampling")
            require(qa["min_chi_bond_angstrom"] >= 0.5 and qa["max_chi_bond_angstrom"] <= 3.0,
                    sid + ": implausible side-chain bond")
            frame_hist[qa["frames"]] += 1
            duration_hist[qa["duration_ns"]] += 1
        for name, count, eligible in zip(structure.residue_names, n_chi, valid):
            if eligible:
                require(count == len(CHI_ATOMS.get(name, [])) and count > 0, sid + ": incomplete chi definition")
            else:
                missing[name] += 1
        geometry = structure.frame_mask
        row = {"sample_id": sid, "residues": len(values), "rmsf_labels": int(mask[:, 0].sum()),
               "hchi_labels": int(valid.sum()), "hchi_geometry_valid": int((valid & geometry).sum()),
               "hchi_min": float(values[valid, h].min()), "hchi_max": float(values[valid, h].max()),
               "frames_per_replica": "1001;1001;1001", "ns_per_replica": "100;100;100",
               "original_rmsf_unchanged": True, "label_sha256": report["output_sha256"]}
        rows.append(row)
        for key in ("residues", "rmsf_labels", "hchi_labels", "hchi_geometry_valid"):
            sums[key] += row[key]
        if index % 250 == 0:
            print(f"Audited {index}/{len(records)} systems", flush=True)
    require(sums["hchi_labels"] == summary["hchi_labels"], "Summary Hchi count mismatch")
    require(sums["rmsf_labels"] == summary["rmsf_labels"], "Summary RMSF count mismatch")
    report = {"passed": True, "complete_systems": len(rows), "counts": dict(sums),
              "missing_hchi_by_residue": dict(missing), "frames_histogram": dict(frame_hist),
              "duration_ns_histogram": dict(duration_hist), "non_hchi_values_variances_masks_unchanged": True,
              "summary_sha256": digest(release / "summary.json"),
              "protocol": summary["protocol"], "all_structures_and_labels_hash_verified": True}
    report["quality_masked_residue_count"] = sum(len(keys) for keys in quality_exclusions.values())
    report["quality_masked_residues"] = quality_exclusions
    with (release / "audit_per_system.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    (release / "release_audit.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2), flush=True)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--label-release-dir", type=Path, default=ROOT / "artifacts/data/atlas_rmsf_hchi_v1")
    args = parser.parse_args()
    audit(args.label_release_dir)


if __name__ == "__main__":
    main()
