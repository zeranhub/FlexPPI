"""Reclassify completed, frozen per-record scores; never run model inference.

Biological membership must be resolved and frozen before this script may read
it for performance aggregation. The original three-stratum products remain
unchanged. Record-weighted means and family-cluster intervals are recomputed.
"""
import os
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT.parent
INPUTS = ROOT / "inputs"
OUT = ROOT / "evaluation"
CLASSIFICATION = SOURCE / "grouping_protein_complex_20261024"
GROUPS = ("protein", "protein_complex")
METRICS = ("spearman", "pearson", "z_mae", "top10_ap", "mae", "rmse", "bias")
SEED = 20261024
BOOTSTRAP = 2000


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def utc():
    return datetime.now(timezone.utc).isoformat()


def frozen_copy(source, destination):
    destination.parent.mkdir(parents=True, exist_ok=True)
    digest = sha(source)
    if destination.exists():
        assert sha(destination) == digest, f"Frozen copy differs: {destination}"
    else:
        shutil.copyfile(source, destination)
    assert sha(source) == sha(destination) == digest
    return {"source": str(source), "frozen_copy": str(destination), "sha256": digest}


def prepare_inputs():
    inventory = []
    for name in ("per_record_metrics.csv", "coverage_per_record.csv", "scoring_audit.json"):
        inventory.append(frozen_copy(SOURCE / "evaluation" / name, INPUTS / name))
    inventory.append(frozen_copy(SOURCE / "cohort/records.csv", INPUTS / "cohort_records.csv"))
    for name in ("training_exposure_protocol_sources_sha256.json", "training_exposure_membership_frozen.csv"):
        path = SOURCE / "references" / name
        if path.exists():
            inventory.append(frozen_copy(path, ROOT / "references" / name))
    manifest = INPUTS / "frozen_inputs.json"
    if manifest.exists():
        recorded = json.loads(manifest.read_text(encoding="utf-8"))
        assert all(row in inventory for row in recorded["files"])
        if recorded["files"] != inventory:
            recorded["files"] = inventory
            recorded["additional_reference_frozen_utc"] = utc()
            manifest.write_text(json.dumps(recorded, indent=2), encoding="utf-8")
    else:
        manifest.write_text(json.dumps({"created_utc": utc(), "files": inventory,
            "scope": "Completed 14-entry scores and coverage before binary biological aggregation; no inference or target changes."}, indent=2), encoding="utf-8")
    return inventory


def read_csv(path):
    return pd.read_csv(path, encoding="utf-8-sig", float_precision="round_trip")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--prepare-inputs", action="store_true")
    args = parser.parse_args()
    inventory = prepare_inputs()
    if args.prepare_inputs:
        print(json.dumps({"status": "completed_inputs_frozen", "files": len(inventory)}), flush=True)
        return

    freeze_path = CLASSIFICATION / "membership_freeze_audit.json"
    freeze = json.loads(freeze_path.read_text(encoding="utf-8"))
    assert freeze["status"] == "frozen" and freeze["unresolved_records"] == 0
    membership_path = CLASSIFICATION / "membership_frozen.csv"
    assert sha(membership_path) == freeze["membership_sha256"]
    assert sha(CLASSIFICATION / "policy_frozen.json") == freeze["policy_sha256"]
    members = read_csv(membership_path)
    assert len(members) == members.sample_id.nunique() == 3063
    assert set(members.benchmark_group) == set(GROUPS)
    assert not members.classification_status.astype(str).str.contains("unknown|unresolved").any()
    counts = {group: int((members.benchmark_group == group).sum()) for group in GROUPS}
    assert counts == freeze["group_counts"]
    evidence = []
    for name in ("membership_frozen.csv", "policy_frozen.json", "membership_freeze_audit.json"):
        evidence.append(frozen_copy(CLASSIFICATION / name, ROOT / "classification" / name))

    original = read_csv(INPUTS / "per_record_metrics.csv")
    old_audit = json.loads((INPUTS / "scoring_audit.json").read_text(encoding="utf-8"))
    methods = old_audit["complete_methods"]
    assert len(methods) == len(set(methods)) == 14
    assert set(original.method) == set(methods)
    assert len(original) == 3063 * len(methods)
    assert not original.duplicated(["method", "sample_id"]).any()
    cohort = read_csv(INPUTS / "cohort_records.csv")
    assert set(cohort.sample_id) == set(members.sample_id) == set(original.sample_id)
    reference = cohort.set_index("sample_id")
    for method, data in original.groupby("method", sort=False):
        assert len(data) == 3063 and set(data.sample_id) == set(reference.index)
        indexed = data.set_index("sample_id")
        assert indexed.n.equals(reference.loc[indexed.index, "n_valid_rmsf"].rename("n"))
        assert indexed.family.equals(reference.loc[indexed.index, "family_component_id"].rename("family"))
        for metric in METRICS:
            applicable = old_audit["methods"][method]["absolute_angstrom"] or metric not in ("mae", "rmse", "bias")
            assert bool(np.isfinite(data[metric]).all()) if applicable else bool(data[metric].isna().all()), (method, metric)

    data = original.rename(columns={"group": "legacy_group"}).merge(
        members[["sample_id", "benchmark_group", "classification_status"]],
        on="sample_id", how="left", validate="many_to_one").rename(columns={"benchmark_group": "group"})
    assert len(data) == len(original) and data.group.notna().all()
    # A merge must not change any precomputed numeric score or record order.
    for column in ("n",) + METRICS:
        assert np.array_equal(original[column].to_numpy(), data[column].to_numpy(), equal_nan=True), column
    OUT.mkdir(parents=True, exist_ok=True)
    columns = ["method", "sample_id", "group", "dataset", "family", "n"] + list(METRICS) + ["legacy_group", "classification_status"]
    data[columns].to_csv(OUT / "per_record_metrics.csv", index=False)
    metadata = cohort[["sample_id", "n_valid_rmsf", "family_component_id"]].merge(
        members[["sample_id", "benchmark_group"]], on="sample_id", validate="one_to_one")
    residues = {g: int(metadata.loc[metadata.benchmark_group == g, "n_valid_rmsf"].sum()) for g in GROUPS}
    assert sum(residues.values()) == 961321

    coverage = read_csv(INPUTS / "coverage_per_record.csv").rename(columns={"group": "legacy_group"})
    assert not coverage.duplicated(["method", "sample_id"]).any() and len(coverage) == len(data)
    assert set(coverage.method) == set(methods)
    coverage = coverage.merge(members[["sample_id", "benchmark_group"]], on="sample_id", validate="many_to_one").rename(columns={"benchmark_group": "group"})
    assert coverage.complete.astype(str).str.lower().eq("true").all()
    assert (coverage.expected == coverage.covered).all()
    coverage.to_csv(OUT / "coverage_per_record.csv", index=False)
    summary_cover = coverage.groupby(["method", "group"], sort=False).agg(
        records=("sample_id", "size"), complete_records=("complete", "sum"),
        expected_residues=("expected", "sum"), covered_residues=("covered", "sum")).reset_index()
    for row in summary_cover.itertuples():
        assert row.records == row.complete_records == counts[row.group]
        assert row.expected_residues == row.covered_residues == residues[row.group]
    summary_cover.to_csv(OUT / "coverage_summary.csv", index=False)

    rng = np.random.default_rng(SEED)
    summary, paired, bootstrap_audit = [], [], {}
    for group in GROUPS:
        group_ids = sorted(members.loc[members.benchmark_group == group, "sample_id"])
        group_data = data[data.group == group]
        family_names = sorted(group_data.family.unique())
        n_family = len(family_names)
        family_position = {family: i for i, family in enumerate(family_names)}
        family_by_record = group_data[group_data.method == "Full"].set_index("sample_id").loc[group_ids, "family"]
        codes = np.array([family_position[x] for x in family_by_record], dtype=np.int32)
        family_sizes = np.bincount(codes, minlength=n_family)
        weights = rng.multinomial(n_family, np.full(n_family, 1 / n_family), size=BOOTSTRAP) if n_family > 1 else None
        denominator = weights @ family_sizes if weights is not None else None
        bootstrap_audit[group] = {"n_records": len(group_ids), "n_families": n_family,
            "family_order_sha256": hashlib.sha256("\n".join(family_names).encode()).hexdigest(),
            "resample_weights_sha256": hashlib.sha256(weights.astype("<i4").tobytes()).hexdigest() if weights is not None else None,
            "shared_draws_for_methods_metrics_and_paired_differences": True}
        indexed = {method: group_data[group_data.method == method].set_index("sample_id").loc[group_ids] for method in methods}
        for method in methods:
            for metric in METRICS:
                values = indexed[method][metric].to_numpy(float)
                if not np.isfinite(values).any():
                    continue
                assert np.isfinite(values).all()
                sums = np.bincount(codes, weights=values, minlength=n_family)
                draws = (weights @ sums) / denominator if weights is not None else None
                lo, hi = np.quantile(draws, [.025, .975]) if draws is not None else (np.nan, np.nan)
                summary.append(dict(method=method, group=group, metric=metric, mean=float(values.mean()),
                    ci_low=float(lo), ci_high=float(hi), n_records=len(values), n_families=n_family,
                    primary_complete_cohort=True,
                    uncertainty="record_weighted_family_cluster_bootstrap" if n_family > 1 else "descriptive_single_family_no_independent_CI"))
                if method != "Full" and metric in ("spearman", "top10_ap", "z_mae", "mae"):
                    differences = indexed["Full"][metric].to_numpy(float) - values
                    dsums = np.bincount(codes, weights=differences, minlength=n_family)
                    delta_draws = (weights @ dsums) / denominator if weights is not None else None
                    dlo, dhi = np.quantile(delta_draws, [.025, .975]) if delta_draws is not None else (np.nan, np.nan)
                    paired.append(dict(group=group, comparator=method, metric=metric,
                        n_records=len(values), n_families=n_family, full_minus_comparator=float(differences.mean()),
                        ci_low=float(dlo), ci_high=float(dhi),
                        direction="negative favors Full" if metric in ("z_mae", "mae") else "positive favors Full"))
    pd.DataFrame(summary).to_csv(OUT / "group_summary.csv", index=False)
    pd.DataFrame(paired).to_csv(OUT / "paired_differences.csv", index=False)
    audit = {**old_audit, "generated_utc": utc(), "aggregation_script_sha256": sha(__file__),
        "original_scoring_audit_sha256": sha(INPUTS / "scoring_audit.json"),
        "original_per_record_metrics_sha256": sha(INPUTS / "per_record_metrics.csv"),
        "group_confirmed": True, "group_counts": counts, "group_residues": residues,
        "membership_path": str(ROOT / "classification/membership_frozen.csv"),
        "membership_sha256": sha(membership_path), "policy_sha256": freeze["policy_sha256"],
        "membership_freeze_audit_sha256": sha(freeze_path),
        "predictions_recomputed": False, "targets_or_masks_changed": False,
        "precomputed_scores_numerically_unchanged": True, "unresolved_records": 0,
        "aggregation_policy": "Same completed predictions, all 3063 records and 961321 valid residues. Biological identity groups frozen before binary performance aggregation; homomultimers and split domains of one mature protein are protein. Different proteins/independent peptide partners are protein_complex.",
        "ci_policy": "2000 shared family-cluster bootstrap draws per group; record sums divided by record counts; 95% percentile intervals; no per-record iid or unseen-family claim.",
        "bootstrap_seed": SEED, "bootstrap_draws": BOOTSTRAP, "bootstrap_audit": bootstrap_audit,
        "training_exposure_audit_sha256": sha(ROOT / "references/training_exposure_protocol_sources_sha256.json"),
        "frozen_inputs": inventory, "frozen_classification": evidence,
        "interpretation_limits": ["Mixed training exposure and model target definitions prevent an equal-training or universal-superiority claim.",
            "The 244 TCR-pMHC records remain one sequence-family dependency component within protein_complex; the new label does not create independent families.",
            "Biological annotation corrections do not rewrite frozen training exposure by original source ID; non-ATLAS exact-ID exposure checks remain limited.",
            "PLANET and RMSF-net full-cohort runs remain deferred as instructed."]}
    (OUT / "scoring_audit.json").write_text(json.dumps(audit, indent=2, ensure_ascii=False), encoding="utf-8")
    assert sha(membership_path) == freeze["membership_sha256"]
    assert all(sha(row["frozen_copy"]) == row["sha256"] for row in inventory)
    print(json.dumps({"status": "binary_aggregation_complete", "methods": len(methods), "group_counts": counts,
        "group_residues": residues, "summary_rows": len(summary), "paired_rows": len(paired),
        "membership_sha256": audit["membership_sha256"], "predictions_recomputed": False}, indent=2), flush=True)


if __name__ == "__main__":
    main()
