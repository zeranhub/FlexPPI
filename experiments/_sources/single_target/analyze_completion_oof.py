"""Matched target comparison and held-out uncertainty diagnostics (CPU only).

Optional two-way calibration uses only the other component-disjoint half of
the *same model's held-out fold*. It never calibrates an evaluated residue on
its own error, nor uses the other four models' potentially contaminated OOF
errors. Gaussian intervals are descriptive and are not confidence guarantees.
"""
from __future__ import annotations

import argparse
import hashlib
from pathlib import Path
from statistics import NormalDist

import numpy as np
import pandas as pd

from completion_common import ARM_METRICS, correlation, fold_summary, jsonl, read_cache, sha256, statistics, write_json

LEVELS = (0.5, 0.68, 0.8, 0.9, 0.95)
RETAINED = tuple(np.arange(0.1, 1.01, 0.1))


def component_side(component, seed):
    return int(hashlib.sha256((str(seed) + ":" + component).encode()).hexdigest(), 16) % 2


def collect(evaluation, variants):
    rows = [dict(row, evaluation=str(evaluation)) for row in jsonl(Path(evaluation) / "prediction_manifest.jsonl") if row["variant"] in variants]
    keys = [(row["variant"], row["sample_id"]) for row in rows]
    if len(set(keys)) != len(keys):
        raise ValueError("Repeated sample/variant prediction")
    return rows


def save_macro(frame, out, prefix, keys, values):
    if frame.empty:
        return
    folds, summary = fold_summary(frame, values, keys)
    folds.to_csv(out / (prefix + "_source_fold_macro.csv"), index=False)
    summary.to_csv(out / (prefix + "_five_fold_mean_sample_sd.csv"), index=False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dual-evaluation", type=Path, required=True)
    parser.add_argument("--single-evaluation", type=Path, help="May be omitted for baseline-only uncertainty diagnostics")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--group-audit", type=Path, help="Frozen record_overlap_audit.csv; required for calibration fitting")
    parser.add_argument("--fit-group-calibration", action="store_true")
    parser.add_argument("--split-seed", type=int, default=1729)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    manifests = collect(args.dual_evaluation, ["full"])
    if args.single_evaluation:
        manifests += collect(args.single_evaluation, ["rmsf_only", "chi_entropy_only"])
    if args.fit_group_calibration and args.group_audit is None:
        raise ValueError("Calibration requires a frozen component audit")
    counts = {}
    for row in manifests:
        counts[row["variant"]] = counts.get(row["variant"], 0) + 1
    if not counts or any(count != 3063 for count in counts.values()):
        raise ValueError("Each analyzed arm must contain all 3063 frozen OOF records: " + str(counts))
    index = {(row["variant"], row["sample_id"]): row for row in manifests}
    per_record, coverage, risk, reliability, diagnostic = [], [], [], [], []
    calibration_parameters, calibrated_coverage, calibrated_errors, assignments = [], [], [], []
    groups = None
    if args.group_audit:
        audit = pd.read_csv(args.group_audit)
        audit = audit[audit.role == "training"].copy()
        if audit.sample_id.duplicated().any():
            raise ValueError("Repeated sample in sequence component audit")
        groups = audit.set_index("sample_id")["component_id"].to_dict()
        for row in manifests:
            component = groups[row["sample_id"]]
            row["component"] = component
            row["calibration_half"] = component_side(component, args.split_seed)
        dual_rows = [row for row in manifests if row["variant"] == "full"]
        assignments = [dict(sample_id=row["sample_id"], dataset=row["dataset"], fold=row["fold"],
                            component_id=row["component"], half=row["calibration_half"]) for row in dual_rows]
        # Write the label-blind split before examining residuals or fitting any parameter.
        pd.DataFrame(assignments).to_csv(args.output / "calibration_group_assignment.csv", index=False)
    fit_sufficient = []
    for row in manifests:
        for metric in ARM_METRICS[row["variant"]]:
            cache = read_cache(row["evaluation"], row, metric)
            mask = cache["valid"]
            if not mask.any():
                continue
            target = cache["target"][mask]
            prediction = cache["prediction"][mask]
            model_variance = np.clip(cache["model_variance"][mask], 1e-6, 1e6)
            label_variance = np.maximum(np.nan_to_num(cache["label_variance"][mask], nan=0.0), 0.0)
            error = prediction - target
            identity = dict(variant=row["variant"], fold=int(row["fold"]), dataset=row["dataset"], metric=metric, sample_id=row["sample_id"])
            per_record.append(identity | statistics(target, prediction))
            for kind, variance in [("model", model_variance), ("model_plus_label", np.clip(model_variance + label_variance, 1e-6, 1e6))]:
                nll = 0.5 * (error ** 2 / variance + np.log(variance) + np.log(2 * np.pi))
                diagnostic.append(identity | dict(variance_kind=kind, count=len(error), gaussian_nll=float(nll.mean()),
                     error_variance_spearman=correlation(abs(error), variance, True), mean_variance=float(variance.mean()),
                     squared_error=float(np.mean(error ** 2)), mean_interval_sd=float(np.sqrt(variance).mean())))
                standard_error = abs(error) / np.sqrt(variance)
                for level in LEVELS:
                    bound = NormalDist().inv_cdf((1 + level) / 2)
                    coverage.append(identity | dict(variance_kind=kind, nominal_coverage=level,
                                                     observed_coverage=float(np.mean(standard_error <= bound)),
                                                     mean_interval_width=float(2 * bound * np.sqrt(variance).mean()), count=len(error)))
                order = np.argsort(variance, kind="stable")
                for fraction in RETAINED:
                    retained = order[:max(1, int(np.ceil(len(order) * fraction)))]
                    risk.append(identity | dict(variance_kind=kind, retained_fraction=float(round(fraction, 2)),
                              retained_mae=float(abs(error[retained]).mean()), retained_rmse=float(np.sqrt(np.mean(error[retained] ** 2))),
                              count=len(retained)))
                for bin_index, selected in enumerate(np.array_split(order, 10)):
                    if len(selected):
                        reliability.append(identity | dict(variance_kind=kind, bin_index=bin_index + 1,
                                  mean_variance=float(variance[selected].mean()), squared_error=float(np.mean(error[selected] ** 2)), count=len(selected)))
            if args.fit_group_calibration:
                fit_sufficient.append(identity | dict(component=row["component"], half=row["calibration_half"],
                    standardized_square=float(np.mean(error ** 2 / model_variance))))
    per = pd.DataFrame(per_record)
    per.to_csv(args.output / "per_record_metrics.csv", index=False)
    save_macro(per, args.output, "metrics", ["variant", "metric"], ["mae", "rmse", "bias", "spearman", "pearson"])
    tables = [("uncertainty_coverage", coverage, ["variant", "metric", "variance_kind", "nominal_coverage"], ["observed_coverage", "mean_interval_width"]),
              ("uncertainty_risk", risk, ["variant", "metric", "variance_kind", "retained_fraction"], ["retained_mae", "retained_rmse"]),
              ("uncertainty_reliability", reliability, ["variant", "metric", "variance_kind", "bin_index"], ["mean_variance", "squared_error"]),
              ("uncertainty_diagnostics", diagnostic, ["variant", "metric", "variance_kind"], ["gaussian_nll", "error_variance_spearman", "mean_variance", "squared_error", "mean_interval_sd"])]
    for prefix, values, keys, numbers in tables:
        table = pd.DataFrame(values)
        table.to_csv(args.output / (prefix + "_per_record.csv"), index=False)
        save_macro(table, args.output, prefix, keys, numbers)
    if args.single_evaluation:
        paired = []
        for row in manifests:
            arm = row["variant"]
            if arm == "full":
                continue
            dual = index["full", row["sample_id"]]
            assert row["fold"] == dual["fold"] and row["dataset"] == dual["dataset"]
            assert row["input_sha256"] == dual["input_sha256"] and row["label_sha256"] == dual["label_sha256"]
            metric = ARM_METRICS[arm][0]
            a, b = read_cache(row["evaluation"], row, metric), read_cache(dual["evaluation"], dual, metric)
            assert np.array_equal(a["valid"], b["valid"])
            mask = a["valid"] & b["valid"]
            assert np.array_equal(a["target"][mask], b["target"][mask])
            if not mask.any():
                continue
            single_stats = statistics(a["target"][mask], a["prediction"][mask])
            dual_stats = statistics(b["target"][mask], b["prediction"][mask])
            numbers = {"single_" + key: single_stats[key] for key in ["mae", "rmse", "spearman", "pearson"]}
            numbers.update({"dual_" + key: dual_stats[key] for key in ["mae", "rmse", "spearman", "pearson"]})
            numbers.update({"single_minus_dual_" + key: single_stats[key] - dual_stats[key] for key in ["mae", "rmse", "spearman", "pearson"]})
            paired.append(dict(variant=arm, metric=metric, fold=row["fold"], sample_id=row["sample_id"], dataset=row["dataset"], count=int(mask.sum()), **numbers))
        paired = pd.DataFrame(paired)
        paired.to_csv(args.output / "single_vs_dual_per_record.csv", index=False)
        save_macro(paired, args.output, "single_vs_dual", ["variant", "metric"], [column for column in paired if column.startswith(("single_", "dual_"))])
    if args.fit_group_calibration:
        sufficient = pd.DataFrame(fit_sufficient)
        multipliers = {}
        for (variant, fold, metric), data in sufficient.groupby(["variant", "fold", "metric"]):
            for test_half in [0, 1]:
                calibration = data[data.half != test_half]
                test = data[data.half == test_half]
                if calibration.empty or test.empty:
                    raise ValueError("An entire component half lacks calibration/test labels")
                assert not set(calibration.component) & set(test.component)
                alpha = max(float(calibration.standardized_square.mean()), 1e-12)
                multipliers[variant, int(fold), metric, test_half] = alpha
                calibration_parameters.append(dict(variant=variant, fold=int(fold), metric=metric, test_half=test_half,
                    variance_multiplier=alpha, calibration_records=len(calibration), test_records=len(test),
                    calibration_components=int(calibration.component.nunique()), test_components=int(test.component.nunique()),
                    fitted_on="opposite component half of same held-out fold; equal record weight"))
        for row in manifests:
            for metric in ARM_METRICS[row["variant"]]:
                cache = read_cache(row["evaluation"], row, metric)
                mask = cache["valid"]
                if not mask.any():
                    continue
                error = cache["prediction"][mask] - cache["target"][mask]
                base_variance = np.clip(cache["model_variance"][mask], 1e-6, 1e6)
                alpha = multipliers[row["variant"], int(row["fold"]), metric, row["calibration_half"]]
                identity = dict(variant=row["variant"], fold=int(row["fold"]), dataset=row["dataset"], metric=metric, sample_id=row["sample_id"])
                for kind, variance in [("uncalibrated_model", base_variance), ("group_crossfit_scaled_model", base_variance * alpha)]:
                    calibrated_errors.append(identity | dict(variance_kind=kind, gaussian_nll=float(np.mean(0.5 * (error ** 2 / variance + np.log(variance) + np.log(2 * np.pi))))))
                    for level in LEVELS:
                        bound = NormalDist().inv_cdf((1 + level) / 2)
                        calibrated_coverage.append(identity | dict(variance_kind=kind, nominal_coverage=level,
                               observed_coverage=float(np.mean(abs(error) <= bound * np.sqrt(variance))),
                               mean_interval_width=float(2 * bound * np.sqrt(variance).mean())))
        pd.DataFrame(calibration_parameters).to_csv(args.output / "group_calibration_parameters.csv", index=False)
        table = pd.DataFrame(calibrated_coverage)
        table.to_csv(args.output / "group_calibration_coverage_per_record.csv", index=False)
        save_macro(table, args.output, "group_calibration_coverage", ["variant", "metric", "variance_kind", "nominal_coverage"], ["observed_coverage", "mean_interval_width"])
        table = pd.DataFrame(calibrated_errors)
        table.to_csv(args.output / "group_calibration_nll_per_record.csv", index=False)
        save_macro(table, args.output, "group_calibration_nll", ["variant", "metric", "variance_kind"], ["gaussian_nll"])
    inputs = dict(dual_manifest_sha256=sha256(args.dual_evaluation / "prediction_manifest.jsonl"))
    if args.single_evaluation:
        inputs["single_manifest_sha256"] = sha256(args.single_evaluation / "prediction_manifest.jsonl")
    if args.group_audit:
        inputs["group_audit_sha256"] = sha256(args.group_audit)
    write_json(args.output / "analysis_audit.json", dict(status="complete", records_per_arm=counts, inputs=inputs,
       code_sha256=sha256(__file__), group_calibration_fitted=args.fit_group_calibration, split_seed=args.split_seed,
       calibration_rule="Two complementary component-disjoint halves within each model's true held-out fold; no use of evaluated residue error during fitting",
       aggregation="Within record first, equal record macro per source and fold; five-fold mean and sample SD",
       interval_interpretation="Unbounded Gaussian diagnostics for frozen reference label means; not thermodynamic uncertainty, experimental confidence intervals, or conformal guarantees",
       label_variance_interpretation="Across available replicas/windows; single observation zero is a placeholder; model_plus_label is a sensitivity diagnostic, not variance of the mean"))
    print(args.output, flush=True)


if __name__ == "__main__":
    main()
