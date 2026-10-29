"""Read-only investigation of the protein-group gap; no fitting or inference.

Every protein is included in summary tables. Worst-case residue profiles are
explicitly diagnostic selections and never replace the frozen benchmark cohort.
"""
import os
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[name] = "1"
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[1]
NEW = ROOT / "protein_complex_benchmark"
METHODS = ("Full", "Intra", "backflip-2.1 (chains)", "backflip-2.1 (assembly)", "DynaProt-M", "backflip-1.0 (chains)", "backflip-1.0-seq (chains)")


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def metrics(y, p):
    return dict(n=len(y), rho=float(np.corrcoef(pd.Series(y).rank().to_numpy(), pd.Series(p).rank().to_numpy())[0, 1]) if len(y) >= 4 and np.ptp(y) > 0 and np.ptp(p) > 0 else np.nan,
        mae=float(np.mean(abs(y - p))), bias=float(np.mean(p - y)), target_mean=float(y.mean()), prediction_mean=float(p.mean()),
        target_std=float(y.std()), prediction_std=float(p.std()),
        std_ratio=float(p.std() / y.std()) if y.std() else np.nan,
        target_p95=float(np.quantile(y, .95)), prediction_p95=float(np.quantile(p, .95)),
        min_prediction=float(p.min()), max_prediction=float(p.max()))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    paths = [NEW / "classification/membership_frozen.csv", NEW / "evaluation/per_record_metrics.csv",
        ROOT / "cohort/records.csv", ROOT / "evaluation/aligned_predictions.npz",
        ROOT.parent / "baselines/cohort/truth_residues.csv.gz"]
    source_hashes = {str(path): sha(path) for path in paths}
    members = pd.read_csv(paths[0], encoding="utf-8-sig")
    per = pd.read_csv(paths[1], float_precision="round_trip")
    cohort = pd.read_csv(paths[2], encoding="utf-8-sig")
    protein_ids = set(members.loc[members.benchmark_group == "protein", "sample_id"])
    assert len(protein_ids) == 2203
    record = cohort[cohort.sample_id.isin(protein_ids)].copy()
    wide = per[per.group == "protein"].pivot(index="sample_id", columns="method", values=["spearman", "mae", "top10_ap", "bias"])
    wide.columns = [f"{metric}::{method}" for metric, method in wide.columns]
    record = record.merge(wide, on="sample_id", validate="one_to_one")
    record["delta_rho_full_minus_bf2"] = record["spearman::Full"] - record["spearman::backflip-2.1 (chains)"]
    record["delta_mae_full_minus_bf2"] = record["mae::Full"] - record["mae::backflip-2.1 (chains)"]
    record["length_bin"] = pd.cut(record.n_residues, [0, 100, 200, 400, 1024], labels=["1-100", "101-200", "201-400", "401-1024"])
    record.sort_values("delta_rho_full_minus_bf2").to_csv(OUT / "all_2203_proteins_ranked.csv", index=False)
    summaries = []
    for stratum in ("dataset", "original_fold", "length_bin", "n_chains"):
        for key, frame in record.groupby(stratum, observed=True):
            summaries.append(dict(stratifier=stratum, level=str(key), records=len(frame),
                full_rho=frame["spearman::Full"].mean(), bf2_rho=frame["spearman::backflip-2.1 (chains)"].mean(),
                delta_rho=frame.delta_rho_full_minus_bf2.mean(), full_wins=int((frame.delta_rho_full_minus_bf2 > 0).sum()),
                full_mae=frame["mae::Full"].mean(), bf2_mae=frame["mae::backflip-2.1 (chains)"].mean()))
    pd.DataFrame(summaries).to_csv(OUT / "gap_by_record_attributes.csv", index=False)
    fam = record.groupby("family_component_id").agg(records=("sample_id", "size"), full_rho=("spearman::Full", "mean"),
        bf2_rho=("spearman::backflip-2.1 (chains)", "mean"), delta_rho=("delta_rho_full_minus_bf2", "mean"))
    fam.to_csv(OUT / "gap_by_family.csv")

    truth = pd.read_csv(paths[4], keep_default_na=False, float_precision="round_trip", dtype={"chain": str, "author_residue_id": str})
    assert len(truth) == 961321
    assert sha(paths[4]) == json.loads((NEW / "evaluation/scoring_audit.json").read_text(encoding="utf-8"))["truth_sha256"]
    with np.load(paths[3]) as arrays:
        predictions = {name: arrays[name] for name in METHODS}
    assert all(len(p) == len(truth) for p in predictions.values())
    selected = truth.sample_id.isin(protein_ids) & truth.valid_rmsf.astype(bool)
    assert int(selected.sum()) == 500964
    pooled = {name: metrics(truth.loc[selected, "target_rmsf"].to_numpy(), p[selected]) for name, p in predictions.items()}
    records_by_id = record.set_index("sample_id")
    profile_rows, regions, top_profiles = [], [], []
    worst_ids = set(record.nsmallest(20, "delta_rho_full_minus_bf2").sample_id)
    for sample, rows in truth[selected].groupby("sample_id", sort=False):
        y = rows.target_rmsf.to_numpy(float)
        indices = rows.index.to_numpy()
        # Rank each actual input chain; ten residues at each end are a diagnostic
        # region, not a new scoring mask in the benchmark. Chains<=20 have no core.
        edge = np.zeros(len(rows), dtype=bool)
        for chain in rows.chain.unique():
            positions = np.flatnonzero(rows.chain.to_numpy() == chain)
            edge[positions[:10]] = True
            edge[positions[-10:]] = True
        for name, p in predictions.items():
            values = p[indices]
            desc = metrics(y, values)
            saved = per[(per.sample_id == sample) & (per.method == name)].iloc[0]
            assert abs(desc["rho"] - saved.spearman) < 1e-12 and abs(desc["mae"] - saved.mae) < 1e-12
            profile_rows.append(dict(sample_id=sample, dataset=records_by_id.loc[sample, "dataset"], method=name, **desc))
            for region, mask in (("ten_residues_per_chain_end", edge), ("remaining_chain_core", ~edge)):
                if mask.any():
                    regions.append(dict(sample_id=sample, dataset=records_by_id.loc[sample, "dataset"], method=name,
                        region=region, **metrics(y[mask], values[mask])))
        if sample in worst_ids:
            profile = rows[["sample_id", "residue_index", "amino_acid", "chain", "author_residue_id", "target_rmsf"]].copy()
            for name in METHODS:
                profile[name] = predictions[name][indices]
            profile["diagnostic_chain_end"] = edge
            top_profiles.append(profile)
    profile = pd.DataFrame(profile_rows)
    profile.to_csv(OUT / "all_protein_prediction_profiles.csv", index=False)
    region = pd.DataFrame(regions)
    region.to_csv(OUT / "all_protein_region_metrics.csv", index=False)
    pd.concat(top_profiles, ignore_index=True).to_csv(OUT / "worst20_residue_profiles_diagnostic_only.csv", index=False)
    rows = []
    for keys, frame in profile.groupby(["dataset", "method"], sort=False):
        rows.append(dict(dataset=keys[0], method=keys[1], records=len(frame),
            **{f"mean_{column}": float(frame[column].mean()) for column in ("rho", "mae", "bias", "target_mean", "prediction_mean", "target_std", "prediction_std", "std_ratio", "target_p95", "prediction_p95")},
            median_std_ratio=float(frame.std_ratio.median()), std_ratio_below_half=int((frame.std_ratio < .5).sum()),
            negative_rho=int((frame.rho < 0).sum())))
    pd.DataFrame(rows).to_csv(OUT / "prediction_profile_by_dataset.csv", index=False)
    region.groupby(["dataset", "method", "region"], sort=False).agg(records=("sample_id", "size"), rho=("rho", "mean"),
        mae=("mae", "mean"), bias=("bias", "mean"), std_ratio=("std_ratio", "median")).to_csv(OUT / "region_summary.csv")
    delta = record.delta_rho_full_minus_bf2
    trim = []
    for n_drop in (0, 10, 20, 100):
        remaining = record.sort_values("delta_rho_full_minus_bf2").iloc[n_drop:]
        trim.append(dict(diagnostic_only=True, worst_gap_records_removed=n_drop, n_records=len(remaining),
            full_rho=float(remaining["spearman::Full"].mean()), bf2_rho=float(remaining["spearman::backflip-2.1 (chains)"].mean()),
            delta_rho=float(remaining.delta_rho_full_minus_bf2.mean())))
    pd.DataFrame(trim).to_csv(OUT / "outlier_influence_diagnostic_only.csv", index=False)
    audit = dict(created_utc=datetime.now(timezone.utc).isoformat(), script_sha256=sha(__file__), source_sha256=source_hashes,
        records=2203, residues=500964, full_rho_wins=int((delta > 0).sum()), bf2_rho_wins=int((delta < 0).sum()),
        delta_rho_quantiles={str(q): float(delta.quantile(q)) for q in (0, .1, .25, .5, .75, .9, 1)},
        family_equal_descriptive_rho=dict(full=float(fam.full_rho.mean()), bf2=float(fam.bf2_rho.mean())),
        pooled_residue_metrics_diagnostic_only=pooled, no_inference=True, no_fitting=True, no_benchmark_change=True,
        source_files_stable=all(sha(path) == digest for path, digest in source_hashes.items()),
        interpretation="Descriptive investigation only. Pooled-residue correlations and removing worst cases are not replacements for the frozen per-record benchmark. Prediction variance shrinkage does not establish a training cause. No model-tuned masks or calibration.")
    assert audit["source_files_stable"]
    (OUT / "profile_audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(json.dumps({"records": 2203, "full_wins": audit["full_rho_wins"], "bf2_wins": audit["bf2_rho_wins"],
        "by_attributes": summaries, "source_files_stable": True}, indent=2), flush=True)


if __name__ == "__main__":
    main()
