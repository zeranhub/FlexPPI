"""Build the external H_chi benchmark for FlexPPI on ITrace trajectories.

Ground truth: per-residue side-chain rotamer entropy (H_chi) computed from the
ITrace MD trajectories with the FlexPPI label convention (see
``extract_itrace_hchi.py``), averaged over the three replicas per complex.
Prediction: FlexPPI ``chi_entropy`` on the complete-complex input.

Metric protocol mirrors the RMSF benchmark: per-chain Spearman/Pearson on each
subset (>=4 residues), equal-chain-weight per-complex mean, equal-complex
macro mean with complex bootstrap 95% CI, and the same cohort filters
(all_245, external_clean without 3KXF, replica-stable external clean).
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, Iterable, Optional

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr


ROOT = Path(r"E:\Project\FlexPPI")
OUTPUT_DIR = ROOT / "output" / "itrace_hchi_benchmark"
RAW_DIR = OUTPUT_DIR / "raw"
METRICS_DIR = OUTPUT_DIR / "metrics"
PAPER_DIR = OUTPUT_DIR / "paper"
CONTAMINATED = {"3kxf"}
BOOTSTRAP_REPETITIONS = 2000
SEED = 2026

SUBSETS = (
    "overall",
    "interface",
    "non_interface",
    "tcr",
    "pmhc",
    "mhc",
    "b2m",
    "peptide",
    "tcr_alpha",
    "tcr_beta",
    "cdr_all",
    "cdr1",
    "cdr2",
    "cdr3",
    "tcr_framework",
)


def normalize_key(frame: pd.DataFrame) -> pd.DataFrame:
    frame = frame.copy()
    frame["pdb_id"] = frame["pdb_id"].astype(str).str.lower()
    frame["chain"] = frame["chain"].astype(str)
    frame["residue_number"] = frame["residue_number"].astype(str)
    frame["insertion_code"] = frame["insertion_code"].fillna("").astype(str).str.strip()
    frame["amino_acid_3"] = frame["amino_acid_3"].astype(str).str.upper()
    return frame


def subset_mask(frame: pd.DataFrame, subset: str) -> np.ndarray:
    if subset == "overall":
        return np.ones(len(frame), dtype=bool)
    if subset == "interface":
        return frame["is_tcr_pmhc_interface_8A"].fillna(False).to_numpy(bool)
    if subset == "non_interface":
        return ~frame["is_tcr_pmhc_interface_8A"].fillna(False).to_numpy(bool)
    if subset == "tcr":
        return frame["chain_role"].isin(["tcr_alpha", "tcr_beta"]).to_numpy()
    if subset == "pmhc":
        return frame["chain_role"].isin(["mhc", "b2m", "peptide"]).to_numpy()
    if subset in {"mhc", "b2m", "peptide", "tcr_alpha", "tcr_beta"}:
        return (frame["chain_role"] == subset).to_numpy()
    if subset == "cdr_all":
        return frame["region_group"].isin(["cdr1", "cdr2", "cdr3"]).to_numpy()
    if subset in {"cdr1", "cdr2", "cdr3", "tcr_framework"}:
        return (frame["region_group"] == subset).to_numpy()
    raise ValueError(f"Unknown subset: {subset}")


def chain_metrics(target: np.ndarray, predicted: np.ndarray) -> Optional[Dict[str, float]]:
    valid = np.isfinite(target) & np.isfinite(predicted)
    target = target[valid]
    predicted = predicted[valid]
    if len(target) < 4 or np.all(target == target[0]) or np.all(predicted == predicted[0]):
        return None
    return {
        "spearman": float(spearmanr(target, predicted).statistic),
        "pearson": float(pearsonr(target, predicted).statistic),
        "residue_count": int(len(target)),
    }


def bootstrap_mean_ci(values: np.ndarray, repetitions: int = BOOTSTRAP_REPETITIONS, seed: int = SEED) -> tuple:
    values = np.asarray(values, dtype=float)
    rng = np.random.default_rng(seed)
    samples = rng.choice(values, size=(repetitions, values.size), replace=True)
    means = samples.mean(axis=1)
    return float(np.quantile(means, 0.025)), float(np.quantile(means, 0.975))


def main() -> None:
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    METRICS_DIR.mkdir(parents=True, exist_ok=True)
    PAPER_DIR.mkdir(parents=True, exist_ok=True)

    backbone = normalize_key(
        pd.read_csv(
            ROOT / "output" / "itrace_benchmark" / "raw" / "common_residue_predictions.csv",
            usecols=[
                "pdb_id", "chain", "residue_number", "insertion_code", "amino_acid_3",
                "region", "chain_role", "is_tcr_pmhc_interface_8A", "region_group",
            ],
            low_memory=False,
        )
    )
    replica = normalize_key(
        pd.read_csv(RAW_DIR / "chi_entropy_per_replica.csv", low_memory=False)
    )
    predictions = normalize_key(
        pd.read_csv(
            ROOT / "output" / "itrace_residue_biology" / "raw" / "all_condition_predictions.csv",
            usecols=[
                "pdb_id", "condition", "chain", "residue_number", "insertion_code",
                "amino_acid_3", "chi_entropy", "chi_entropy_std",
            ],
            low_memory=False,
        )
    )
    predictions = predictions.loc[predictions["condition"].eq("full")].drop(columns=["condition"])

    key = ["pdb_id", "chain", "residue_number", "insertion_code", "amino_acid_3"]
    consensus = (
        replica.groupby(key, dropna=False)["chi_entropy"]
        .agg(["mean", "std", "count"])
        .reset_index()
        .rename(
            columns={
                "mean": "target_chi_entropy",
                "std": "target_chi_entropy_std",
                "count": "target_chi_replica_count",
            }
        )
    )
    common = (
        backbone.merge(consensus, on=key, how="left")
        .merge(predictions, on=key, how="left")
        .rename(
            columns={
                "chi_entropy": "flexppi_chi_entropy",
                "chi_entropy_std": "flexppi_chi_entropy_std",
            }
        )
    )
    common.to_csv(RAW_DIR / "hchi_consensus_per_residue.csv", index=False)

    # ---- replica reliability (per-chain pairwise Spearman, complex macro) ----
    reliability_rows = []
    replica_id = ["pdb_id", "chain", "residue_number", "insertion_code", "amino_acid_3"]
    for (pdb_id, chain), group in replica.groupby(["pdb_id", "chain"], sort=True):
        piv = group.pivot_table(
            index=replica_id,
            columns="traj_id",
            values="chi_entropy",
            aggfunc="first",
        )
        rep_columns = sorted(piv.columns)
        pairs = []
        for i in range(len(rep_columns)):
            for j in range(i + 1, len(rep_columns)):
                a = piv[rep_columns[i]].to_numpy(float)
                b = piv[rep_columns[j]].to_numpy(float)
                mask = np.isfinite(a) & np.isfinite(b)
                a, b = a[mask], b[mask]
                if len(a) < 4 or np.all(a == a[0]) or np.all(b == b[0]):
                    continue
                pairs.append(float(spearmanr(a, b).statistic))
        if pairs:
            reliability_rows.append(
                {
                    "pdb_id": pdb_id,
                    "chain": chain,
                    "pair_count": len(pairs),
                    "mean_pairwise_spearman": float(np.mean(pairs)),
                    "median_pairwise_spearman": float(np.median(pairs)),
                }
            )
    reliability = pd.DataFrame(reliability_rows)
    complex_reliability = (
        reliability.groupby("pdb_id", as_index=False)["mean_pairwise_spearman"]
        .mean()
        .rename(columns={"mean_pairwise_spearman": "complex_mean_pairwise_spearman"})
    )
    reliability.to_csv(METRICS_DIR / "replica_reliability.csv", index=False)

    # ---- per-chain metrics ----
    chain_rows = []
    for subset in SUBSETS:
        selected = common.loc[subset_mask(common, subset)]
        for (pdb_id, chain), group in selected.groupby(["pdb_id", "chain"], sort=True):
            metrics = chain_metrics(
                group["target_chi_entropy"].to_numpy(float),
                group["flexppi_chi_entropy"].to_numpy(float),
            )
            if metrics is None:
                continue
            chain_rows.append(
                {
                    "pdb_id": pdb_id,
                    "chain": chain,
                    "subset": subset,
                    **metrics,
                }
            )
    per_chain = pd.DataFrame(chain_rows)
    per_chain.to_csv(METRICS_DIR / "per_chain_metrics.csv", index=False)

    per_complex = (
        per_chain.groupby(["pdb_id", "subset"], as_index=False)[["spearman", "pearson"]]
        .mean()
        .rename(columns={"spearman": "flexppi__spearman", "pearson": "flexppi__pearson"})
    )
    residue_sums = per_chain.groupby(["pdb_id", "subset"], as_index=False)["residue_count"].sum()
    per_complex = per_complex.merge(residue_sums, on=["pdb_id", "subset"], how="left")
    per_complex = per_complex.sort_values(["subset", "pdb_id"])
    per_complex.to_csv(METRICS_DIR / "per_complex_metrics.csv", index=False)

    # ---- cohorts and summary ----
    all_ids = set(per_complex.loc[per_complex["subset"].eq("overall"), "pdb_id"].astype(str))
    clean_ids = all_ids - {value.lower() for value in CONTAMINATED}
    stable_ids = set(
        complex_reliability.loc[
            complex_reliability["complex_mean_pairwise_spearman"] >= 0.5, "pdb_id"
        ].astype(str)
    )
    cohorts = {
        "all_245": all_ids,
        "external_clean": clean_ids,
        "replica_stable_external_clean": clean_ids & stable_ids,
    }
    summary_rows = []
    for cohort_name, identifiers in cohorts.items():
        for subset in SUBSETS:
            values = per_complex.loc[
                per_complex["subset"].eq(subset) & per_complex["pdb_id"].isin(identifiers),
                "flexppi__spearman",
            ].dropna()
            if values.empty:
                continue
            low, high = bootstrap_mean_ci(values.to_numpy(float))
            summary_rows.append(
                {
                    "cohort": cohort_name,
                    "subset": subset,
                    "method": "flexppi",
                    "metric": "spearman",
                    "macro_mean": float(values.mean()),
                    "macro_median": float(values.median()),
                    "ci95_low": low,
                    "ci95_high": high,
                    "complex_count": int(len(values)),
                }
            )
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(METRICS_DIR / "summary_metrics.csv", index=False)

    # ---- context: correlation with per-complex RMSF Spearman ----
    rmsf_complex = pd.read_csv(
        ROOT / "output" / "itrace_benchmark" / "metrics" / "per_complex_metrics.csv",
        usecols=["pdb_id", "label_definition", "subset", "flexppi__spearman"],
    )
    rmsf_overall = rmsf_complex.loc[
        rmsf_complex["subset"].eq("overall")
        & rmsf_complex["label_definition"].eq("chain_fitted_internal_rmsf")
    ].copy()
    rmsf_overall = rmsf_overall.rename(columns={"flexppi__spearman": "rmsf_spearman"})
    hchi_overall = per_complex.loc[per_complex["subset"].eq("overall")].copy()
    joined = hchi_overall.merge(rmsf_overall[["pdb_id", "rmsf_spearman"]], on="pdb_id", how="inner")
    joined = joined.loc[joined["pdb_id"].astype(str).isin(clean_ids)]
    joined = joined.dropna(subset=["flexppi__spearman", "rmsf_spearman"])
    cross = {
        "n_complexes": int(len(joined)),
        "pearson_rmsf_vs_hchi_spearman": float(
            pearsonr(joined["rmsf_spearman"], joined["flexppi__spearman"]).statistic
        ),
        "spearman_rmsf_vs_hchi_spearman": float(
            spearmanr(joined["rmsf_spearman"], joined["flexppi__spearman"]).statistic
        ),
    }
    joined.to_csv(METRICS_DIR / "per_complex_rmsf_vs_hchi.csv", index=False)

    # ---- diagnostics: distribution + amino-acid-type prior ----
    external = common.loc[
        common["pdb_id"].isin(clean_ids)
        & common["target_chi_entropy"].notna()
        & common["flexppi_chi_entropy"].notna()
    ].copy()
    overall_values = per_complex.loc[
        per_complex["subset"].eq("overall") & per_complex["pdb_id"].isin(clean_ids),
        "flexppi__spearman",
    ].dropna()
    aa_prior_rows = []
    aa_means = external.groupby("amino_acid_3")["target_chi_entropy"].mean()
    for (pdb_id, chain), group in external.groupby(["pdb_id", "chain"]):
        prior = group["amino_acid_3"].map(aa_means).to_numpy(float)
        target = group["target_chi_entropy"].to_numpy(float)
        valid = np.isfinite(prior) & np.isfinite(target)
        if valid.sum() < 4:
            continue
        prior = prior[valid]
        target = target[valid]
        if np.all(prior == prior[0]) or np.all(target == target[0]):
            continue
        aa_prior_rows.append(
            {
                "pdb_id": pdb_id,
                "chain": chain,
                "aa_type_prior_spearman": float(spearmanr(prior, target).statistic),
            }
        )
    aa_prior = pd.DataFrame(aa_prior_rows)
    aa_complex = (
        aa_prior.groupby("pdb_id", as_index=False)["aa_type_prior_spearman"]
        .mean()
        .rename(columns={"aa_type_prior_spearman": "complex_aa_prior_spearman"})
    )
    aa_complex.to_csv(METRICS_DIR / "per_complex_aa_type_prior.csv", index=False)

    # ---- residue coverage summary ----
    coverage = {
        "backbone_residues": int(len(common)),
        "residues_with_hchi_target": int(common["target_chi_entropy"].notna().sum()),
        "residues_with_target_and_prediction": int(
            (common["target_chi_entropy"].notna() & common["flexppi_chi_entropy"].notna()).sum()
        ),
        "per_replica_chi_rows": int(len(replica)),
        "residue_keys_with_any_replica": int(
            replica.groupby(["pdb_id", "chain", "residue_number", "insertion_code", "amino_acid_3"]).ngroups
        ),
    }
    diagnostics = {
        "external_clean_overall_per_complex_spearman": {
            "n_complexes": int(len(overall_values)),
            "mean": float(overall_values.mean()),
            "median": float(overall_values.median()),
            "min": float(overall_values.min()),
            "max": float(overall_values.max()),
            "q25": float(overall_values.quantile(0.25)),
            "q75": float(overall_values.quantile(0.75)),
        },
        "amino_acid_type_prior_diagnostic": {
            "method": (
                "baseline prediction = dataset-wide mean H_chi per amino-acid type; "
                "per-chain Spearman vs MD, equal-chain-weight per complex"
            ),
            "n_complexes": int(len(aa_complex)),
            "mean_spearman": float(aa_complex["complex_aa_prior_spearman"].mean()),
        },
        "replica_reliability_external_clean": {
            "n_chains": int(
                reliability.loc[
                    reliability["pdb_id"].isin(clean_ids)
                ].shape[0]
            ),
            "mean_chain_pairwise_spearman": float(
                reliability.loc[
                    reliability["pdb_id"].isin(clean_ids), "mean_pairwise_spearman"
                ].mean()
            ),
        },
    }
    report = {
        "protocol": (
            "H_chi per residue per replica from ITrace trajectories using the FlexPPI label "
            "convention (CHI_ATOMS, 3-bin rotamers, base-3 code, Shannon entropy over frames); "
            "three replicas averaged per residue; per-chain Spearman on each subset (>=4 residues); "
            "equal-chain-weight per-complex mean; equal-complex macro mean with complex bootstrap "
            "95% CI (2000 resamples, seed 2026); 3KXF excluded from external_clean."
        ),
        "coverage": coverage,
        "diagnostics": diagnostics,
        "rmsf_vs_hchi_correlation": cross,
    }
    report_path = OUTPUT_DIR / "itrace_hchi_benchmark_report.json"
    report_path.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")

    print(json.dumps(report, ensure_ascii=False, indent=2))
    print()
    print("== external_clean summary (FlexPPI H_chi Spearman vs MD) ==")
    table = summary.loc[summary["cohort"].eq("external_clean")]
    for _, row in table.iterrows():
        print(
            f"{row['subset']:<16} n={int(row['complex_count']):>3}  "
            f"mean={row['macro_mean']:.3f}  CI=[{row['ci95_low']:.3f},{row['ci95_high']:.3f}]"
        )


if __name__ == "__main__":
    main()
