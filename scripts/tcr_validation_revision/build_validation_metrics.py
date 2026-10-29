"""Build direct B4.2.3 prediction-versus-MD tables without test-fitted controls.

The only evaluated prediction is deposited-5IVX input, epoch-40 fold_3.
Training-only amino-acid and global means use precisely the frozen fold_3 train
membership. The verified OOF export is used ONLY as a carrier of training labels;
no other-fold prediction is used to evaluate 5IVX or fit either baseline.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[2]
OLD = ROOT / "output/tcr_real_states_20260915"
OUT = ROOT / "output/tcr_prediction_reference_revision_20260915"
ANALYSIS = ROOT / "artifacts/evaluation_atlas_hchi_v1/analysis"
OOF = ROOT / "artifacts/evaluation_atlas_hchi_v1/oof"
SPLIT = ROOT / "output/atlas_hchi_retraining_v1_linux/payload/data/splits/state_rmsf_hchi_v2_atlas_hchi_v1"
METRICS = {"rmsf": "rmsf", "hchi": "chi_entropy"}
PREDICTORS = {"FlexPPI": "prediction", "AA-only": "aa_only_prediction", "Global-mean": "global_mean_prediction"}


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def json_lines(path):
    return [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]


def correlation(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    if len(x) < 3 or np.ptp(x) < 1e-12 or np.ptp(y) < 1e-12:
        return np.nan
    return float(spearmanr(x, y).statistic)


def regression(frame, column):
    prediction = frame[column].to_numpy(float)
    reference = frame.md_reference.to_numpy(float)
    assert len(frame) and np.isfinite([prediction, reference]).all()
    error = prediction - reference
    denominator = float(np.square(reference - reference.mean()).sum())
    return {
        "n_residues": len(frame),
        "mae": float(np.abs(error).mean()),
        "rmse": float(np.sqrt(np.square(error).mean())),
        "bias": float(error.mean()),
        "spearman": correlation(prediction, reference),
        "r2": float(1 - np.square(error).sum() / denominator) if denominator > 1e-15 else np.nan,
        "training_AA_adjusted_spearman": correlation(prediction - frame.aa_only_prediction, reference - frame.aa_only_prediction),
        "md_mean": float(reference.mean()),
        "prediction_mean": float(prediction.mean()),
    }


def build_training_priors():
    audit = json.loads((ANALYSIS / "analysis_summary.json").read_text())
    labels_path = ANALYSIS / "oof_per_residue.csv.gz"
    assert sha(labels_path) == audit["oof_per_residue_sha256"]
    manifest_path = OOF / "prediction_manifest.jsonl"
    assert sha(manifest_path) == audit["source_manifest_sha256"]
    manifest = json_lines(manifest_path)
    training = json_lines(SPLIT / "fold_3_train.jsonl")
    test = json_lines(SPLIT / "fold_3_test.jsonl")
    train_ids = {row["sample_id"] for row in training}
    test_ids = {row["sample_id"] for row in test}
    assert len(train_ids) == len(training) == 2454
    assert len(test_ids) == len(test) == 614
    assert not (train_ids & test_ids)
    assert "itrace_5ivx" in test_ids and "itrace_5ivx" not in train_ids
    assert {row["sample_id"] for row in manifest if row["fold"] != 3} == train_ids
    assert {row["sample_id"] for row in manifest if row["fold"] == 3} == test_ids
    source_hashes = manifest[0]["identity"]["source_hashes"]["3"]
    for role in ("train", "test"):
        assert sha(SPLIT / f"fold_3_{role}.jsonl") == source_hashes[role]
    print("Verified exact frozen train/test membership and original label-export hash.", flush=True)

    usecols = ["sample_id", "fold", "amino_acid", "valid_rmsf", "target_rmsf", "valid_chi_entropy", "target_chi_entropy"]
    labels = pd.read_csv(labels_path, usecols=usecols)
    assert set(labels.loc[labels.fold.ne(3), "sample_id"]) == train_ids
    # Select by audited membership, not simply by a row's numerical fold label.
    training_labels = labels.loc[labels.sample_id.isin(train_ids)].copy()
    assert training_labels.fold.ne(3).all()
    priors, global_means, counts = [], {}, {}
    for metric, old_metric in METRICS.items():
        selected = training_labels.loc[training_labels[f"valid_{old_metric}"]]
        assert np.isfinite(selected[f"target_{old_metric}"]).all()
        group = selected.groupby("amino_acid")[f"target_{old_metric}"].agg(["mean", "size"])
        global_means[metric] = float(selected[f"target_{old_metric}"].astype(float).mean())
        counts[metric] = len(selected)
        for aa, row in group.iterrows():
            priors.append({"metric": metric, "fold": 3, "amino_acid": aa,
                           "training_mean": float(row["mean"]), "training_residues": int(row["size"]),
                           "training_global_mean": global_means[metric]})
    priors = pd.DataFrame(priors)
    old_priors_path = ANALYSIS / "training_fold_aa_priors.csv"
    assert sha(old_priors_path) == audit["source_files"][old_priors_path.name]
    old_priors = pd.read_csv(old_priors_path)
    old_priors = old_priors.loc[old_priors.fold.eq(3)].copy()
    old_priors.metric = old_priors.metric.replace({"chi_entropy": "hchi"})
    compared = priors.merge(old_priors, on=["metric", "fold", "amino_acid"], suffixes=("_new", "_old"), validate="one_to_one")
    assert len(compared) == len(priors)
    # The archived per-residue CSV serializes float32 targets with decimal
    # rounding. Re-aggregation differs by <1e-9; exact membership/counts and
    # source hashes remain mandatory. This is not a relaxed scientific check.
    assert np.allclose(compared.training_mean_new, compared.training_mean_old, atol=1e-8, rtol=0)
    assert np.array_equal(compared.training_residues_new, compared.training_residues_old)
    provenance = {
        "training_samples": len(train_ids), "test_samples": len(test_ids),
        "training_valid_residues": counts, "train_sample_ids": sorted(train_ids),
        "fold3_train_manifest_sha256": sha(SPLIT / "fold_3_train.jsonl"),
        "fold3_test_manifest_sha256": sha(SPLIT / "fold_3_test.jsonl"),
        "label_export_sha256": sha(labels_path), "oof_manifest_sha256": sha(manifest_path),
        "original_training_priors_sha256": sha(old_priors_path),
        "max_difference_from_audited_prior": float(np.abs(compared.training_mean_new - compared.training_mean_old).max()),
        "definition": "Residue-pooled valid target means, fitted only on the exact frozen fold_3 training membership; no test labels or other-fold predictions used.",
    }
    return priors, global_means, provenance


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    priors, global_means, baseline_audit = build_training_priors()
    source_path = OLD / "b423_existing_holo_md_comparison.csv"
    source_audit = json.loads((OLD / "b423_holo_md_audit.json").read_text())
    checkpoint = ROOT / "artifacts/performance_first_v2_atlas_hchi_v1/fold_3/epoch_040.pt"
    assert sha(checkpoint) == source_audit["checkpoint_sha256"]
    assert sha(source_audit["cache"]) == source_audit["cache_sha256"]
    frame = pd.read_csv(source_path, dtype={"author_id": str, "topology_author_id": str})
    frame = frame.loc[frame.chain.isin(["E", "F"])].copy()
    assert len(frame) == 856 and not frame.duplicated(["chain", "author_id", "metric"]).any()
    frame["chain_role"] = frame.chain.map({"E": "TCR alpha", "F": "TCR beta"})
    frame["prediction"] = frame.predicted_deposited_fold3
    frame["unit"] = frame.metric.map({"rmsf": "Angstrom", "hchi": "nats"})
    frame["evaluation_mask"] = frame.valid & ~(frame.metric.eq("hchi") & frame.amino_acid.isin(list("AGP")))
    frame["exclusion_reason"] = np.where(~frame.valid, "original MD/input validity mask", "")
    frame.loc[frame.metric.eq("hchi") & frame.amino_acid.isin(list("AGP")), "exclusion_reason"] = "A/G lack chi angles; P near-constant label protocol: excluded from Hchi evaluation"
    frame = frame.merge(priors[["metric", "amino_acid", "training_mean"]], on=["metric", "amino_acid"], how="left", validate="many_to_one")
    frame = frame.rename(columns={"training_mean": "aa_only_prediction"})
    frame["global_mean_prediction"] = frame.metric.map(global_means)
    assert frame.loc[frame.evaluation_mask, list(PREDICTORS.values()) + ["md_reference"]].notna().all().all()
    for predictor, column in PREDICTORS.items():
        prefix = {"FlexPPI": "model", "AA-only": "aa_only", "Global-mean": "global_mean"}[predictor]
        frame[prefix + "_error"] = frame[column] - frame.md_reference
        frame[prefix + "_absolute_error"] = np.abs(frame[prefix + "_error"])
    frame["aa_adjusted_md_reference"] = frame.md_reference - frame.aa_only_prediction
    frame["aa_adjusted_prediction"] = frame.prediction - frame.aa_only_prediction
    frame["evaluation_scope"] = "Absolute holo state; combination-level heldout fold_3; homologous training proteins remain; not apo-to-holo delta validation"
    valid = frame.loc[frame.evaluation_mask].copy()
    assert valid.groupby("metric").size().to_dict() == {"hchi": 340, "rmsf": 428}
    summaries, within = [], []
    for scope, chains in (("TCR alpha beta", ["E", "F"]), ("TCR alpha", ["E"]), ("TCR beta", ["F"])):
        for metric in METRICS:
            group = valid.loc[valid.chain.isin(chains) & valid.metric.eq(metric)]
            baseline = regression(group, "aa_only_prediction")
            for predictor, column in PREDICTORS.items():
                result = regression(group, column)
                summaries.append({"scope": scope, "metric": metric, "unit": group.unit.iloc[0], "predictor": predictor,
                                  **result, "mae_reduction_vs_AA_percent": 100 * (1 - result["mae"] / baseline["mae"]),
                                  "rmse_reduction_vs_AA_percent": 100 * (1 - result["rmse"] / baseline["rmse"])})
            for aa, subgroup in group.groupby("amino_acid"):
                baseline = regression(subgroup, "aa_only_prediction")
                for predictor in ["FlexPPI", "AA-only"]:
                    result = regression(subgroup, PREDICTORS[predictor])
                    within.append({"scope": scope, "metric": metric, "amino_acid": aa, "predictor": predictor,
                                   **result, "mae_reduction_vs_AA_percent": 100 * (1 - result["mae"] / baseline["mae"]) if baseline["mae"] else np.nan})
    summary = pd.DataFrame(summaries)
    largest = valid.sort_values("model_absolute_error", ascending=False).groupby("metric", sort=False).head(20).copy()
    largest["error_rank_within_metric"] = largest.groupby("metric")["model_absolute_error"].rank(method="first", ascending=False).astype(int)
    largest = largest.sort_values(["metric", "error_rank_within_metric"])
    frame.to_csv(OUT / "validation_paired_residues.csv", index=False)
    priors.to_csv(OUT / "validation_training_only_aa_priors.csv", index=False)
    summary.to_csv(OUT / "validation_metrics.csv", index=False)
    pd.DataFrame(within).to_csv(OUT / "validation_within_aa_metrics.csv", index=False)
    largest.to_csv(OUT / "validation_largest_errors.csv", index=False)
    audit = {
        "status": "passed", "source_pair_table": str(source_path), "source_pair_table_sha256": sha(source_path),
        "prediction": "Deposited 5IVX input; frozen epoch_040.pt fold_3 only; no ensemble, refitting or calibration",
        "checkpoint_sha256": sha(checkpoint), "baseline": baseline_audit,
        "reference": source_audit,
        "coverage": {"total_TCR_residues": 428, "rmsf_valid": 428, "hchi_valid_excluding_AGP": 340,
                     "TCR_alpha_residues": 194, "TCR_beta_residues": 234},
        "statistics": "Descriptive within-one-complex residue comparisons. No independent-residue assumption, significance test, bootstrap CI, or biological replicate claim.",
        "scope": "Absolute holo-state prediction versus existing MD labels. No apo MD reference; cannot validate binding-induced delta or cellular activation.",
        "within_AA_note": "Per-type Spearman is undefined for AA-only because the type baseline is constant. Training-AA-adjusted correlation subtracts the frozen training prior, never case-specific fitted means.",
        "outputs": {name: sha(OUT / name) for name in ["validation_paired_residues.csv", "validation_training_only_aa_priors.csv", "validation_metrics.csv", "validation_within_aa_metrics.csv", "validation_largest_errors.csv"]},
    }
    (OUT / "validation_audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(summary.loc[summary.scope.eq("TCR alpha beta")].to_string(index=False), flush=True)
    print("Largest prediction errors:", flush=True)
    print(largest.loc[largest.error_rank_within_metric.le(3), ["chain_role", "author_id", "amino_acid", "metric", "md_reference", "prediction", "model_absolute_error"]].to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
