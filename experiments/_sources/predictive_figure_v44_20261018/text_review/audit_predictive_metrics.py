"""Independent direct-threshold and pairwise spot checks of v44 metrics.

Does not import or call the metric worker's implementation. All source inputs
are read-only. Only this agent's JSON audit is written.
"""
from pathlib import Path
import hashlib
import json
import sys
from bisect import bisect_left

ROOT = Path("E:/Project/FlexPPI/experiments/_sources/predictive_figure_v44_20261018")
EXP = ROOT.parent
CORE = Path("E:/Project/FlexPPI/experiments/core_corrected_pipeline")
sys.path.insert(0, str(EXP / "completed_results_20261012/plot_runtime"))
import numpy as np
import pandas as pd


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def evaluate_direct(labels, scores):
    """Enumerate actual thresholds and compare all positive/negative pairs."""
    npos = int(labels.sum())
    nneg = len(labels) - npos
    thresholds = [np.inf] + sorted(set(map(float, scores)), reverse=True)
    recall, precision, fpr = [0.0], [1.0], [0.0]
    ap, last_recall = 0.0, 0.0
    for threshold in thresholds[1:]:
        selected = scores >= threshold
        hits = int(np.count_nonzero(labels & selected))
        rr, pp = hits / npos, hits / int(selected.sum())
        ap += (rr - last_recall) * pp
        last_recall = rr
        recall.append(rr)
        precision.append(pp)
        fpr.append((int(selected.sum()) - hits) / nneg)
    positive, negative = scores[labels], scores[~labels]
    pairdiff = positive[:, None] - negative[None, :]
    auc = (np.count_nonzero(pairdiff > 0) + 0.5 * np.count_nonzero(pairdiff == 0)) / pairdiff.size
    macro = []
    for desired in np.linspace(0.0, 1.0, 101):
        if desired == 0:
            macro.append(1.0)
        else:
            first = next(i for i, actual in enumerate(recall) if actual >= desired)
            macro.append(precision[first])
    return ap, auc, {"threshold": np.array(thresholds), "recall": np.array(recall), "precision": np.array(precision), "fpr": np.array(fpr), "tpr": np.array(recall)}, np.array(macro)


def main():
    analysis = ROOT / "analysis"
    summary = json.loads((analysis / "summary.json").read_text(encoding="utf-8"))
    data = pd.read_csv(analysis / "per_record_metrics.csv")
    codes = pd.read_csv(analysis / "record_codes.csv")
    types = pd.read_csv(analysis / "within_type_per_record_metrics.csv")
    macro = pd.read_csv(analysis / "macro_pr_curve.csv")
    typed_macro = pd.read_csv(analysis / "within_type_record_macro_metrics.csv")
    endpoints = pd.read_csv(analysis / "all_selected_endpoint.csv")
    schema = json.loads((CORE / "static_corrected/training_schema.json").read_text(encoding="utf-8"))["records"]
    static = np.load(CORE / "static_corrected/static_sequence_complex_oof.npz")["prediction"]
    aa_only = np.load(CORE / "static_corrected/aa_only_oof.npz")["prediction"]
    rng = np.random.default_rng(20261018)
    chosen = set(map(int, rng.choice(len(codes), 6, replace=False)))
    for target in ["rmsf", "chi_entropy"]:
        gg = data[(data.target == target) & (data.method == "Full") & data.eligible]
        for col, which in [("ap", "min"), ("ap", "max"), ("prevalence", "max"), ("n_positions", "max")]:
            idx = gg[col].idxmin() if which == "min" else gg[col].idxmax()
            chosen.add(int(gg.loc[idx, "record_code"]))
    chosen = sorted(chosen)
    raw = {}
    segments = {}
    for t in ["rmsf", "chi_entropy"]:
        for m in ["Full", "Intra", "Static", "AA"]:
            with np.load(analysis / f"raw_curves_{t}_{m}.npz") as loaded:
                arr = {k: loaded[k].copy() for k in loaded.files}
            raw[(t, m)] = arr
            start = np.r_[0, np.flatnonzero(np.diff(arr["record_code"]) != 0) + 1]
            stop = np.r_[start[1:], len(arr["record_code"])]
            segments[(t, m)] = {int(arr["record_code"][s]): slice(int(s), int(e)) for s, e in zip(start, stop)}
    checked = []
    typed_checked = []
    errors = []
    max_abs_ap = max_abs_auc = max_abs_curve = max_abs_mae = max_abs_type_ap = max_abs_type_auc = 0.0
    for code in chosen:
        record = schema[code]
        sid, start, stop = record["sample_id"], record["start"], record["stop"]
        with np.load(CORE / "neural_exports/evaluation/full" / (sid + ".npz")) as z:
            full = {k: z[k].copy() for k in z.files}
        with np.load(CORE / "neural_exports/evaluation/intrachain" / (sid + ".npz")) as z:
            intra = {k: z[k].copy() for k in z.files}
        preds = {"Full": full["mean"], "Intra": intra["mean"], "Static": static[start:stop], "AA": aa_only[start:stop]}
        for j, target in enumerate(["rmsf", "chi_entropy"]):
            valid = full["valid"][:, j] & intra["valid"][:, j] & np.isfinite(full["target"][:, j])
            for p in preds.values():
                valid &= np.isfinite(p[:, j])
            y = full["target"][valid, j].astype(float)
            q90 = float(np.quantile(y, 0.9, method="linear"))
            labels = y >= q90
            for method, pred in preds.items():
                score = pred[valid, j].astype(float)
                row = data[(data.record_code == code) & (data.target == target) & (data.method == method)].iloc[0]
                assert int(row.n_positions) == len(y) and int(row.n_positive) == int(labels.sum())
                assert abs(row.threshold_q90 - q90) < 1e-12
                mae = float(np.mean(np.abs(score - y)))
                max_abs_mae = max(max_abs_mae, abs(mae - row.mae))
                if labels.all() or not labels.any():
                    assert not row.eligible and np.isnan(row.ap) and np.isnan(row.auroc)
                    continue
                ap, auc, curve, sampled = evaluate_direct(labels, score)
                max_abs_ap = max(max_abs_ap, abs(ap - row.ap))
                max_abs_auc = max(max_abs_auc, abs(auc - row.auroc))
                rr = raw[(target, method)]
                keep = segments[(target, method)][code]
                for key, values in curve.items():
                    got = rr[key][keep]
                    assert got.shape == values.shape, (sid, target, method, key)
                    if key == "threshold":
                        assert np.isinf(got[0]) and np.allclose(got[1:], values[1:], rtol=0, atol=0)
                    else:
                        max_abs_curve = max(max_abs_curve, float(np.max(np.abs(got - values))))
                checked.append({"record_code": code, "sample_id": sid, "target": target, "method": method, "n_positions": len(y), "n_positive": int(labels.sum()), "n_score_ties": int(len(score) - len(np.unique(score))), "ap": ap, "auroc": auc})
            if target == "chi_entropy":
                aa_names = full["amino_acid"][valid].astype(str)
                for aa_name in sorted(set(aa_names)):
                    sel_type = aa_names == aa_name
                    yy = y[sel_type]
                    tq = float(np.quantile(yy, 0.9, method="linear"))
                    labels_type = yy >= tq
                    eligible_type = len(yy) >= 4 and 0 < int(labels_type.sum()) < len(yy)
                    for method, pred in preds.items():
                        row_type = types[(types.record_code == code) & (types.amino_acid == aa_name) & (types.method == method)].iloc[0]
                        assert bool(row_type.eligible) == eligible_type
                        assert int(row_type.n_positions) == len(yy) and int(row_type.n_positive) == int(labels_type.sum())
                        assert abs(row_type.threshold_q90 - tq) < 1e-12
                        if eligible_type:
                            sc_type = pred[valid, j][sel_type].astype(float)
                            ap_type, auc_type, _, _ = evaluate_direct(labels_type, sc_type)
                            max_abs_type_ap = max(max_abs_type_ap, abs(ap_type - row_type.ap))
                            max_abs_type_auc = max(max_abs_type_auc, abs(auc_type - row_type.auroc))
                            typed_checked.append({"record_code": code, "sample_id": sid, "amino_acid": aa_name, "method": method, "n_positions": len(yy), "n_positive": int(labels_type.sum()), "ap": ap_type, "auroc": auc_type})
    checks = {
        "3063_record_metadata_complete": len(codes) == 3063 and codes.record_code.is_unique,
        "24504_target_method_record_rows_unique": len(data) == 3063 * 2 * 4 and not data[["record_code", "target", "method"]].duplicated().any(),
        "direct_threshold_AP_spotchecks": max_abs_ap < 1e-12,
        "direct_positive_negative_pair_AUC_spotchecks": max_abs_auc < 1e-12,
        "all_saved_threshold_points_for_spotchecks": max_abs_curve < 1e-12,
        "source_raw_unit_MAE_spotchecks": max_abs_mae < 1e-12,
        "within_type_raw_source_AP_spotchecks": max_abs_type_ap < 1e-12,
        "within_type_raw_source_pair_AUC_spotchecks": max_abs_type_auc < 1e-12,
        "score_ties_present_in_spotchecks": any(r["n_score_ties"] > 0 for r in checked),
    }
    aa = types[(types.method == "AA") & types.eligible]
    checks["all_within_type_AA_AP_equals_prevalence"] = bool(np.allclose(aa.ap, aa.prevalence, rtol=0, atol=1e-12))
    checks["all_within_type_AA_AUC_equals_half"] = bool(np.allclose(aa.auroc, 0.5, rtol=0, atol=1e-12))
    checks["all_ineligible_type_reasons_explicit"] = bool(types.loc[~types.eligible, "undefined_reason"].notna().all())
    checks["all_invalid_AP_AUC_retained_as_NaN"] = bool(data.loc[~data.eligible, ["ap", "auroc"]].isna().all().all() and types.loc[~types.eligible, ["ap", "auroc"]].isna().all().all())
    checks["actual_q90_prevalence_at_least_ten_percent"] = bool((data.loc[data.eligible, "prevalence"] >= 0.1 - 1e-12).all())
    checks["all_18_Hchi_record_types_including_absent_types_retained"] = bool(len(types) == 3063 * 18 * 4 and types.groupby(["record_code", "method"]).size().eq(18).all())
    for target in ["rmsf", "chi_entropy"]:
        for method in ["Full", "Intra", "Static", "AA"]:
            g = data[(data.target == target) & (data.method == method)]
            ss = summary["metrics"][target][method]
            checks[f"summary_from_all_records_{target}_{method}"] = bool(abs(g.ap.mean() - ss["mean_ap"]) < 1e-12 and abs(g.auroc.mean() - ss["mean_auroc"]) < 1e-12 and int(g.eligible.sum()) == ss["n_records"] and int(g.n_positions.sum()) == ss["n_positions"])
            rr = raw[(target, method)]
            grids = []
            for code, group in g[g.eligible].iterrows():
                cc = int(group.record_code)
                keep = segments[(target, method)][cc]
                recall, precision = rr["recall"][keep], rr["precision"][keep]
                assert recall[-1] == 1.0 and abs(precision[-1] - group.prevalence) < 1e-12
                # Direct first-complete-threshold lookup, without interpolation/envelope.
                recall_list = recall.tolist()
                grids.append([1.0 if r == 0 else precision[bisect_left(recall_list, r)] for r in np.linspace(0, 1, 101)])
            mean_grid = np.mean(grids, axis=0)
            gg = macro[(macro.target == target) & (macro.method == method)].sort_values("recall")
            checks[f"macroPR_all_records_from_raw_thresholds_{target}_{method}"] = bool(np.allclose(mean_grid, gg.precision, rtol=0, atol=1e-12))
            checks[f"PR_no_skill_actual_prevalence_{target}_{method}"] = bool(np.allclose(gg.no_skill_precision, g.loc[g.eligible, "prevalence"].mean(), rtol=0, atol=1e-12))
            endpoint = endpoints[(endpoints.target == target) & (endpoints.method == method)].iloc[0]
            checks[f"full_selection_endpoint_{target}_{method}"] = bool(endpoint.recall == 1.0 and abs(endpoint.precision - g.loc[g.eligible, "prevalence"].mean()) < 1e-12 and endpoint.n_records == g.eligible.sum())
    for method in ["Full", "Intra", "Static", "AA"]:
        group = types[(types.method == method) & types.eligible]
        equal_type = group.groupby("record_code")[["ap", "auroc"]].mean()
        ss = summary["within_type"][method]
        checks[f"within_type_equal_type_then_record_{method}"] = bool(abs(equal_type.ap.mean() - ss["mean_ap"]) < 1e-12 and abs(equal_type.auroc.mean() - ss["mean_auroc"]) < 1e-12 and len(equal_type) == ss["n_records"] and len(group) == ss["n_eligible_record_types"])
    for target in ["rmsf", "chi_entropy"]:
        with np.load(analysis / f"prediction_reference_{target}.npz") as z:
            ref, pred, code = z["reference"].astype(float), z["prediction"].astype(float), z["record_code"]
            counts = np.bincount(code, minlength=3063)
            weight = 1.0 / (3063 * counts[code])
            observed_hist = np.histogram2d(ref, pred, bins=(z["x_edges"], z["y_edges"]), weights=weight)[0]
            checks[f"density_all_record_equal_mass_{target}"] = bool(len(ref) == summary["metrics"][target]["Full"]["n_positions"] and np.allclose(observed_hist, z["histogram_mass"], rtol=0, atol=1e-12) and abs(weight.sum() - 1) < 1e-12)
            checks[f"density_no_value_clipping_{target}"] = bool(ref.min() >= z["x_edges"][0] and ref.max() <= z["x_edges"][-1] and pred.min() >= z["y_edges"][0] and pred.max() <= z["y_edges"][-1])
            macro_mae = float(np.mean(np.bincount(code, weights=np.abs(ref - pred), minlength=3063) / counts))
            checks[f"density_original_unit_macroMAE_{target}"] = abs(macro_mae - summary["metrics"][target]["Full"]["macro_mae"]) < 1e-12
    checks = {k: bool(v) for k, v in checks.items()}
    report = {
        "status": "pass" if all(checks.values()) else "fail",
        "checks": checks,
        "summary_sha256": sha(analysis / "summary.json"),
        "protocol_sha256": sha(analysis / "PRESPECIFIED_PROTOCOL.md"),
        "endpoint_source_sha256": sha(analysis / "all_selected_endpoint.csv"),
        "metric_worker_code_sha256": sha(analysis / "compute_metrics.py"),
        "independent_audit_code_sha256": sha(__file__),
        "spotchecked_records": chosen,
        "spotchecked_target_method_comparisons": checked,
        "spotchecked_within_type_comparisons": typed_checked,
        "maximum_absolute_differences": {"AP": max_abs_ap, "AUROC": max_abs_auc, "raw_PR_ROC_points": max_abs_curve, "MAE": max_abs_mae, "within_type_AP": max_abs_type_ap, "within_type_AUROC": max_abs_type_auc},
        "coverage": summary["coverage"],
        "within_type_coverage": {m: summary["within_type"][m]["n_eligible_record_types"] for m in summary["within_type"]},
        "interpretation": [
            "AP/AUROC evaluate a derived q90 ranking task; they are not classification-training accuracies.",
            "Whole-record Hchi ranking includes amino-acid-type prior; within-type AP is a separate auxiliary task.",
            "Reference threshold ties can raise positive prevalence above ten percent.",
            "Mean record AP is not the finite-grid area of the macro PR curve.",
            "Density and trend use original predictions and do not calibrate the model.",
            "Record bootstrap intervals do not represent independently resampled biological families.",
        ],
    }
    out = ROOT / "text_review/independent_metrics_audit.json"
    out.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({"status": report["status"], "n_spotchecked_records": len(chosen), "n_comparisons": len(checked), "maximum_absolute_differences": report["maximum_absolute_differences"], "failed_checks": [k for k, v in checks.items() if not v]}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
