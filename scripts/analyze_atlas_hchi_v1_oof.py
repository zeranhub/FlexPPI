"""Analyze the complete, provenance-checked ATLAS-Hchi v1 OOF cache.

No model fitting or inference is performed here. Amino-acid baselines use only
the four training folds for each held-out fold. All valid residues are retained.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm

sys.path.insert(0, str(Path(__file__).resolve().parent))
from export_atlas_hchi_v1_oof_predictions import EXPECTED_COUNTS, EXPECTED_SOURCES, ROOT
from export_v2_oof_predictions import METRICS, load_cache, sha256
from flexppi.evaluation import _correlation, _spearman

RUN = ROOT / "artifacts/performance_first_v2_atlas_hchi_v1"
STATS = ["mae", "rmse", "bias", "pearson", "spearman"]


def regression(prediction, target):
    p, y = np.asarray(prediction, float), np.asarray(target, float)
    if len(y) == 0 or len(p) != len(y) or not np.isfinite([p, y]).all():
        raise ValueError("Empty, mismatched, or nonfinite regression input")
    error = p - y
    return dict(count=len(y), mae=float(np.abs(error).mean()),
                rmse=float(np.sqrt(np.mean(error ** 2))), bias=float(error.mean()),
                pearson=_correlation(p, y), spearman=_spearman(p, y),
                target_mean=float(y.mean()), target_std=float(y.std()),
                target_iqr=float(np.quantile(y, .75) - np.quantile(y, .25)),
                pred_mean=float(p.mean()))


def read_complete_cache(folder):
    summary = json.loads((folder / "run_summary.json").read_text())
    if (summary["status"] != "complete" or summary["partial"]
            or summary["processed_samples"] != 3068
            or not summary["server_reproduction_passed"]
            or summary["protocol"]["name"] != "atlas_hchi_v1_final_epoch40_oof"):
        raise ValueError("Only complete, server-verified new-run OOF predictions may be plotted")
    manifest = [json.loads(line) for line in (folder / "prediction_manifest.jsonl").read_text().splitlines()]
    if len(manifest) != 3068 or len({r["sample_id"] for r in manifest}) != 3068:
        raise ValueError("Incomplete or duplicate OOF system inventory")
    frames = []
    for i, row in enumerate(manifest):
        if row["status"] != "complete":
            raise ValueError("Incomplete cache entry")
        arrays = load_cache(Path(row["cache_path"]), row["identity"])
        meta = json.loads(str(arrays["metadata_json"]))
        if meta["protocol"]["name"] != summary["protocol"]["name"]:
            raise ValueError("Mixed old/new model caches")
        r = row["record"]
        frame = pd.DataFrame(dict(
            sample_id=row["sample_id"], dataset=r["dataset"], fold=row["fold"],
            chain=arrays["chain"], author_residue_id=arrays["author_residue_id"],
            residue_index=arrays["residue_index"], amino_acid=arrays["amino_acid"],
            biological_system=r["metadata"]["biological_system"],
            cluster_id=r["metadata"]["cluster_id"], source_group=r["source_group"],
            whole_system_length=len(arrays["chain"]), chain_count=len(set(arrays["chain"]))))
        for j, metric in enumerate(METRICS):
            for prefix, array in [("valid", "valid"), ("target", "target"),
                                  ("pred", "mean"), ("model_variance", "model_variance"),
                                  ("label_variance", "label_variance")]:
                frame[f"{prefix}_{metric}"] = arrays[array][:, j]
        frames.append(frame)
        if (i + 1) % 500 == 0:
            print(f"Read {i + 1}/3068 verified OOF caches", flush=True)
    data = pd.concat(frames, ignore_index=True)
    if data.duplicated(["sample_id", "chain", "author_residue_id"]).any():
        raise ValueError("Duplicate residue identifiers")
    observed_sources = data.groupby("dataset").sample_id.nunique().to_dict()
    if observed_sources != EXPECTED_SOURCES:
        raise ValueError("Source inventory differs from frozen run")
    for metric, count in EXPECTED_COUNTS.items():
        if int(data[f"valid_{metric}"].sum()) != count:
            raise ValueError("Invalid full residue coverage: " + metric)
    return data, manifest, summary


def attach_training_priors(data):
    rows = []
    for metric in METRICS:
        data[f"aa_prior_{metric}"] = np.nan
        for fold in range(5):
            train = data.loc[(data.fold != fold) & data[f"valid_{metric}"]]
            means = train.assign(y=train[f"target_{metric}"].astype(float)).groupby("amino_acid").y.mean()
            counts = train.groupby("amino_acid").size()
            test = (data.fold == fold) & data[f"valid_{metric}"]
            data.loc[test, f"aa_prior_{metric}"] = data.loc[test, "amino_acid"].map(means)
            for aa, value in means.items():
                rows.append(dict(metric=metric, fold=fold, amino_acid=aa,
                                 training_mean=value, training_residues=int(counts[aa])))
        if data.loc[data[f"valid_{metric}"], f"aa_prior_{metric}"].isna().any():
            raise ValueError("Held-out amino-acid type has no training prior")
    return pd.DataFrame(rows)


def system_tables(data):
    raw, baseline, residual = [], [], []
    for _, group in data.groupby("sample_id", sort=False):
        keys = ["sample_id", "dataset", "fold", "biological_system", "cluster_id",
                "whole_system_length", "chain_count"]
        base = {key: group.iloc[0][key] for key in keys}
        for metric in METRICS:
            part = group.loc[group[f"valid_{metric}"]]
            p, y, prior = [part[f"{prefix}_{metric}"] for prefix in ["pred", "target", "aa_prior"]]
            common = dict(**base, metric=metric, n_residues=len(part))
            raw.append(dict(**common, predictor="FlexPPI", **regression(p, y)))
            baseline.append(dict(**common, predictor="AA-only", **regression(prior, y)))
            if metric == "chi_entropy":
                residual.append(dict(**common, predictor="AA-corrected",
                                     **regression(p - prior, y - prior)))
    return pd.DataFrame(raw), pd.DataFrame(baseline), pd.DataFrame(residual)


def source_tables(data, raw, baseline, residual):
    rows, aa_rows = [], []
    combined = pd.concat([raw, baseline, residual], ignore_index=True)
    for source in ["all", *EXPECTED_SOURCES]:
        d = data if source == "all" else data.loc[data.dataset == source]
        s = combined if source == "all" else combined.loc[combined.dataset == source]
        for fold in [-1, 0, 1, 2, 3, 4]:
            dd = d if fold == -1 else d.loc[d.fold == fold]
            ss = s if fold == -1 else s.loc[s.fold == fold]
            for metric in METRICS:
                valid = dd.loc[dd[f"valid_{metric}"]]
                p, y, prior = [valid[f"{prefix}_{metric}"] for prefix in ["pred", "target", "aa_prior"]]
                variants = [("FlexPPI", p, y), ("AA-only", prior, y)]
                if metric == "chi_entropy":
                    variants.append(("AA-corrected", p - prior, y - prior))
                for predictor, pp, yy in variants:
                    subset = ss.loc[(ss.metric == metric) & (ss.predictor == predictor)]
                    base = dict(dataset=source, fold=fold, metric=metric, predictor=predictor,
                                n_systems=len(subset))
                    rows.append(dict(**base, scope="pooled", **regression(pp, yy)))
                    rows.append(dict(**base, scope="macro_by_system", count=len(valid),
                                     **{k: float(subset[k].mean()) for k in STATS}))
                for aa, part in valid.groupby("amino_acid"):
                    aa_rows.append(dict(dataset=source, fold=fold, metric=metric, amino_acid=aa,
                                        **regression(part[f"pred_{metric}"], part[f"target_{metric}"])))
    return pd.DataFrame(rows), pd.DataFrame(aa_rows)


def uncertainty_curves(prediction, target, model_variance, label_variance):
    p, y, mv, lv = [np.asarray(x, float) for x in [prediction, target, model_variance, label_variance]]
    if not len(y) or not np.isfinite([p, y, mv, lv]).all() or (mv <= 0).any() or (lv < 0).any():
        raise ValueError("Invalid uncertainty inputs")
    error = p - y
    order = np.argsort(mv, kind="stable")
    full_mae = float(np.abs(error).mean())
    risk, reliability, coverage = [], [], []
    for fraction in np.linspace(.1, 1, 10):
        retained = order[:max(1, int(np.ceil(fraction * len(y))))]
        mae = float(np.abs(error[retained]).mean())
        risk.append(dict(retained_fraction=float(fraction), retained_count=len(retained),
                         total_count=len(y), mae=mae, full_mae=full_mae,
                         relative_mae=mae / full_mae if full_mae else 0.0))
    for index, ids in enumerate(np.array_split(order, min(10, len(y)))):
        reliability.append(dict(bin=index + 1, count=len(ids),
                                rms_model_sd=float(np.sqrt(mv[ids].mean())),
                                rms_total_sd=float(np.sqrt((mv[ids] + lv[ids]).mean())),
                                rmse=float(np.sqrt((error[ids] ** 2).mean())),
                                mae=float(np.abs(error[ids]).mean())))
    for kind, variance in [("model_only", mv), ("model_plus_label", np.clip(mv + lv, 1e-6, 1e6))]:
        z = np.abs(error) / np.sqrt(variance)
        for level in [*.1 * np.arange(1, 10), .95, .99]:
            coverage.append(dict(variance_kind=kind, nominal_coverage=float(level),
                                 observed_coverage=float((z <= norm.ppf((1 + level) / 2)).mean()),
                                 count=len(y)))
    return risk, reliability, coverage


def uncertainty_tables(data):
    all_risk, all_rel, all_cov = [], [], []
    for (source, fold), group in data.groupby(["dataset", "fold"]):
        for metric in METRICS:
            valid = group.loc[group[f"valid_{metric}"]]
            values = [valid[f"{prefix}_{metric}"] for prefix in ["pred", "target", "model_variance", "label_variance"]]
            risk, reliability, coverage = uncertainty_curves(*values)
            base = dict(dataset=source, fold=int(fold), metric=metric)
            for target, rows in [(all_risk, risk), (all_rel, reliability), (all_cov, coverage)]:
                target.extend(dict(**base, **r) for r in rows)
    return [pd.DataFrame(rows) for rows in [all_risk, all_rel, all_cov]]


def compare_server(source, priors):
    checks = []
    for fold in range(5):
        server = json.loads((RUN / f"evaluation/fold_{fold}_test.json").read_text())
        recorded = server["metrics"]["chi_entropy"]["amino_acid_corrected"]["training_priors"]
        local_priors = priors.loc[(priors.fold == fold) & (priors.metric == "chi_entropy")].set_index("amino_acid").training_mean
        if any(abs(local_priors[aa] - value) > 1e-7 for aa, value in recorded.items()):
            raise ValueError("Training-only AA priors disagree with server")
        for dataset in ["all", *EXPECTED_SOURCES]:
            for metric in METRICS:
                block = server["metrics"][metric]
                ref = block["overall"] if dataset == "all" else block["by_dataset"][dataset]
                for scope in ["pooled", "macro_by_system"]:
                    match = source.loc[(source.fold == fold) & (source.dataset == dataset)
                                       & (source.metric == metric) & (source.scope == scope)
                                       & (source.predictor == "FlexPPI")].iloc[0]
                    for stat in ["mae", "rmse", "pearson", "spearman"]:
                        expected = ref[stat] if scope == "pooled" else ref[scope][stat]
                        diff = float(match[stat] - expected)
                        tolerance = .001 if metric == "rmsf" and stat == "rmse" else 2e-4
                        checks.append(dict(fold=fold, dataset=dataset, metric=metric, scope=scope,
                                           statistic=stat, cpu_minus_server=diff, tolerance=tolerance,
                                           strict_2e_4_passed=abs(diff) <= 2e-4, passed=abs(diff) <= tolerance))
        corrected = server["metrics"]["chi_entropy"]["amino_acid_corrected"]["residual_overall"]
        for scope in ["pooled", "macro_by_system"]:
            row = source.loc[(source.fold == fold) & (source.dataset == "all")
                             & (source.predictor == "AA-corrected") & (source.scope == scope)].iloc[0]
            for stat in ["pearson", "spearman"]:
                expected = corrected[stat] if scope == "pooled" else corrected[scope][stat]
                diff = float(row[stat] - expected)
                checks.append(dict(fold=fold, dataset="all", metric="hchi_residual", scope=scope,
                                   statistic=stat, cpu_minus_server=diff, tolerance=2e-4,
                                   strict_2e_4_passed=abs(diff) <= 2e-4, passed=abs(diff) <= 2e-4))
    checks = pd.DataFrame(checks)
    if not checks.passed.all():
        raise ValueError("Analysis disagrees with server metrics")
    return checks


def run(folder, output):
    output.mkdir(parents=True, exist_ok=True)
    data, manifest, summary = read_complete_cache(folder)
    priors = attach_training_priors(data)
    raw, baseline, residual = system_tables(data)
    print("System and training-only baseline metrics computed", flush=True)
    sources, aa = source_tables(data, raw, baseline, residual)
    checks = compare_server(sources, priors)
    risk, reliability, coverage = uncertainty_tables(data)
    tables = {
        "per_system_metrics": raw, "aa_only_per_system_metrics": baseline,
        "hchi_residual_per_system_metrics": residual, "source_fold_metrics": sources,
        "training_fold_aa_priors": priors, "within_aa_metrics": aa,
        "uncertainty_risk": risk, "uncertainty_reliability": reliability,
        "uncertainty_coverage": coverage, "cpu_server_reproduction": checks,
    }
    for name, frame in tables.items():
        frame.to_csv(output / (name + ".csv"), index=False, encoding="utf-8")
    data.to_csv(output / "oof_per_residue.csv.gz", index=False, compression="gzip")
    report = dict(status="complete", model_run=RUN.name, n_systems=3068,
                  rows=len(data), source_systems=EXPECTED_SOURCES,
                  valid_residues=EXPECTED_COUNTS,
                  invalid_metric_rows={m: int((~data[f"valid_{m}"]).sum()) for m in METRICS},
                  undefined_correlations={name: int(frame.spearman.isna().sum()) for name, frame in
                                          [("raw", raw), ("aa_only", baseline), ("residual", residual)]},
                  test_selection="Fixed epoch 40; one held-out fold model per system; no ensemble.",
                  baseline="Residue-pooled AA-type means fitted on four training folds only.",
                  residual="Subtract training-fold AA mean from both prediction and target; MAE unchanged.",
                  uncertainty="Model variance is the single model's variance head, not five-fold ensemble variance.",
                  coverage_warning="Gaussian intervals are uncalibrated diagnostics, not validated confidence intervals. model_plus_label uses held-out label variance and is not deployment uncertainty.",
                  statistics="Descriptive only. Five-fold mean +/- sample SD where shown; no significance testing or independence claim for residues/folds.",
                  server_max_absolute_difference=float(checks.cpu_minus_server.abs().max()),
                  strict_2e_4_exceedances=int((~checks.strict_2e_4_passed).sum()),
                  platform_tolerance="RMSF RMSE 0.001 Angstrom, other metrics 0.0002; initial strict audit retained.",
                  source_manifest_sha256=sha256(folder / "prediction_manifest.jsonl"),
                  oof_per_residue_sha256=sha256(output / "oof_per_residue.csv.gz"),
                  analysis_script_sha256=sha256(Path(__file__)),
                  source_files={name + ".csv": sha256(output / (name + ".csv")) for name in tables})
    (output / "analysis_summary.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=ROOT / "artifacts/evaluation_atlas_hchi_v1/oof")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "artifacts/evaluation_atlas_hchi_v1/analysis")
    args = parser.parse_args()
    run(args.input_dir.resolve(), args.output_dir.resolve())
