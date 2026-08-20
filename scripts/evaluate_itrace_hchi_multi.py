"""Multi-method external H_chi benchmark on ITrace.

Tier 1 (direct): FlexPPI's dedicated H_chi head vs MD rotamer-entropy truth.
Tier 2 (cross-descriptor proxy): each static flexibility score (FlexPPI RMSF,
FastProtFlex, BackFlip-structure/seq, ProDy-ANM/GNM, RMSF-net-PDB) is evaluated
as a ranking/detection proxy for side-chain entropy; these methods were NOT
trained on H_chi, so the comparison is explicitly cross-descriptor.
Controls: amino-acid-type prior (dataset means), MD RMSF<->H_chi target
orthogonality, and FlexPPI input-condition ablation (full / pMHC / TCR /
isolated chains).

Metric protocol mirrors the RMSF benchmark: per-chain metrics on each subset
(>=4 residues), equal-chain-weight per-complex mean, equal-complex macro mean
with complex bootstrap 95% CI, external_clean cohort (3KXF excluded).
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, Optional

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr


ROOT = Path(r"E:\Project\FlexPPI")
OUTPUT_DIR = ROOT / "output" / "itrace_hchi_multi_benchmark"
METRICS_DIR = OUTPUT_DIR / "metrics"
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

# method label -> score column in the merged frame
METHOD_COLUMNS: Dict[str, str] = {
    "FlexPPI-Hchi": "flexppi_chi_entropy",
    "FlexPPI-RMSF": "flexppi_pretrain_rmsf_angstrom",
    "FastProtFlex": "fastprotflex_score",
    "BackFlip-structure": "backflip_structure_rmsf_angstrom",
    "BackFlip-seq": "backflip_seq_rmsf_angstrom",
    "ProDy-ANM": "anm_sq_fluct",
    "ProDy-GNM": "gnm_sq_fluct",
    "RMSF-net-PDB": "rmsf_net_pdb_score",
}

KEY = ["pdb_id", "chain", "residue_number", "insertion_code", "amino_acid_3"]


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


def auroc_ap_pk(target: np.ndarray, score: np.ndarray) -> Optional[Dict[str, float]]:
    """Top-20% detection metrics on one chain-subset."""

    valid = np.isfinite(target) & np.isfinite(score)
    target = target[valid]
    score = score[valid]
    if len(target) < 8:
        return None
    positive = target >= np.quantile(target, 0.8)
    n_pos = int(positive.sum())
    n_neg = len(positive) - n_pos
    if n_pos == 0 or n_neg == 0:
        return None
    order = np.argsort(-score, kind="mergesort")
    ranked_positive = positive[order]
    # AUROC via Mann-Whitney U on ranks
    ranks = np.arange(len(score), 0, -1)
    pos_ranks = ranks[ranked_positive].sum()
    auroc = (pos_ranks - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg)
    # average precision over the ranking
    cumulative = np.cumsum(ranked_positive)
    precisions = cumulative[ranked_positive] / np.arange(1, len(score) + 1)[ranked_positive]
    ap = float(precisions.mean()) if len(precisions) else float("nan")
    # precision at K (K = number of positives)
    top_k = ranked_positive[:n_pos]
    pk = float(top_k.mean()) if n_pos else float("nan")
    return {
        "auroc": float(auroc),
        "average_precision": ap,
        "precision_at_k": pk,
    }


def chain_metrics(target: np.ndarray, score: np.ndarray) -> Optional[Dict[str, float]]:
    valid = np.isfinite(target) & np.isfinite(score)
    target = target[valid]
    score = score[valid]
    if len(target) < 4 or np.all(target == target[0]) or np.all(score == score[0]):
        return None
    out = {
        "spearman": float(spearmanr(target, score).statistic),
        "pearson": float(pearsonr(target, score).statistic),
        "residue_count": int(len(target)),
    }
    detection = auroc_ap_pk(target, score)
    if detection:
        out.update(detection)
    return out


def bootstrap_mean_ci(values: np.ndarray, repetitions: int = BOOTSTRAP_REPETITIONS, seed: int = SEED) -> tuple:
    values = np.asarray(values, dtype=float)
    rng = np.random.default_rng(seed)
    samples = rng.choice(values, size=(repetitions, values.size), replace=True)
    means = samples.mean(axis=1)
    return float(np.quantile(means, 0.025)), float(np.quantile(means, 0.975))


def main() -> None:
    METRICS_DIR.mkdir(parents=True, exist_ok=True)

    consensus = normalize_key(
        pd.read_csv(
            ROOT / "output" / "itrace_hchi_benchmark" / "raw" / "hchi_consensus_per_residue.csv",
            low_memory=False,
        )
    )
    baseline_scores = normalize_key(
        pd.read_csv(
            ROOT / "artifacts" / "evaluation" / "additional_baselines" / "results"
            / "itrace_common_residue_predictions.csv",
            low_memory=False,
        )
    )
    conditions = normalize_key(
        pd.read_csv(
            ROOT / "output" / "itrace_residue_biology" / "raw" / "all_condition_predictions.csv",
            usecols=[
                "pdb_id", "condition", "chain", "residue_number", "insertion_code",
                "amino_acid_3", "chi_entropy",
            ],
            low_memory=False,
        )
    )
    condition_columns = {}
    for condition in ("full", "pmhc", "tcr", "isolated_A", "isolated_B", "isolated_C",
                      "isolated_D", "isolated_E"):
        subset = conditions.loc[conditions["condition"].eq(condition), KEY + ["chi_entropy"]]
        subset = subset.rename(columns={"chi_entropy": f"flexppi_chi_{condition}"})
        condition_columns[f"flexppi_chi_{condition}"] = subset

    frame = consensus.merge(
        baseline_scores[
            KEY
            + ["target_chainfit_rmsf_angstrom"]
            + [col for col in METHOD_COLUMNS.values() if col != "flexppi_chi_entropy"]
        ],
        on=KEY,
        how="left",
    )
    for column, subset in condition_columns.items():
        frame = frame.merge(subset, on=KEY, how="left")

    # amino-acid-type prior: dataset-wide mean target H_chi per residue type
    aa_means = frame.loc[frame["target_chi_entropy"].notna()].groupby("amino_acid_3")["target_chi_entropy"].mean()
    frame["aa_type_prior"] = frame["amino_acid_3"].map(aa_means)
    method_columns = dict(METHOD_COLUMNS)
    method_columns["AA-type prior"] = "aa_type_prior"

    # ---- per-chain metrics ----
    chain_rows = []
    for subset in SUBSETS:
        selected = frame.loc[subset_mask(frame, subset)]
        for (pdb_id, chain), group in selected.groupby(["pdb_id", "chain"], sort=True):
            target = group["target_chi_entropy"].to_numpy(float)
            for method, column in method_columns.items():
                metrics = chain_metrics(target, group[column].to_numpy(float))
                if metrics is None:
                    continue
                chain_rows.append(
                    {
                        "pdb_id": pdb_id,
                        "chain": chain,
                        "subset": subset,
                        "method": method,
                        **metrics,
                    }
                )
    per_chain = pd.DataFrame(chain_rows)
    per_chain.to_csv(METRICS_DIR / "per_chain_metrics.csv", index=False)

    metric_columns = ["spearman", "pearson", "auroc", "average_precision", "precision_at_k"]
    per_complex = (
        per_chain.groupby(["pdb_id", "subset", "method"], as_index=False)[metric_columns]
        .mean()
        .rename(
            columns={
                "spearman": "flexppi__spearman",
                "pearson": "flexppi__pearson",
                "auroc": "flexppi__top20_auroc",
                "average_precision": "flexppi__top20_average_precision",
                "precision_at_k": "flexppi__top20_precision_at_k",
            }
        )
    )
    residue_sums = per_chain.groupby(["pdb_id", "subset", "method"], as_index=False)["residue_count"].sum()
    per_complex = per_complex.merge(residue_sums, on=["pdb_id", "subset", "method"], how="left")
    per_complex = per_complex.sort_values(["subset", "method", "pdb_id"])
    per_complex.to_csv(METRICS_DIR / "per_complex_metrics.csv", index=False)

    # ---- cohorts and summary ----
    all_ids = set(per_complex.loc[per_complex["subset"].eq("overall"), "pdb_id"].astype(str))
    clean_ids = all_ids - {value.lower() for value in CONTAMINATED}
    cohorts = {"all_245": all_ids, "external_clean": clean_ids}
    summary_rows = []
    for cohort_name, identifiers in cohorts.items():
        for subset in SUBSETS:
            for method in method_columns:
                for metric in metric_columns:
                    key_metric = {
                        "spearman": "flexppi__spearman",
                        "pearson": "flexppi__pearson",
                        "auroc": "flexppi__top20_auroc",
                        "average_precision": "flexppi__top20_average_precision",
                        "precision_at_k": "flexppi__top20_precision_at_k",
                    }[metric]
                    values = per_complex.loc[
                        per_complex["subset"].eq(subset)
                        & per_complex["method"].eq(method)
                        & per_complex["pdb_id"].isin(identifiers),
                        key_metric,
                    ].dropna()
                    if values.empty:
                        continue
                    low, high = bootstrap_mean_ci(values.to_numpy(float))
                    summary_rows.append(
                        {
                            "cohort": cohort_name,
                            "subset": subset,
                            "method": method,
                            "metric": metric,
                            "macro_mean": float(values.mean()),
                            "macro_median": float(values.median()),
                            "ci95_low": low,
                            "ci95_high": high,
                            "complex_count": int(len(values)),
                        }
                    )
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(METRICS_DIR / "summary_metrics.csv", index=False)

    # ---- FlexPPI condition ablation, matched on chain role ----
    # For each chain role, compare the conditions that actually contain that
    # chain: full, the group input (pmhc or tcr), and the isolated chain.
    role_groups = {
        "mhc": "pmhc",
        "b2m": "pmhc",
        "peptide": "pmhc",
        "tcr_alpha": "tcr",
        "tcr_beta": "tcr",
    }
    chain_letter = (
        frame.loc[subset_mask(frame, "overall")]
        .groupby(["pdb_id", "chain_role"])["chain"]
        .first()
        .reset_index()
        .rename(columns={"chain": "chain_letter"})
    )
    condition_rows = []
    condition_labels = {
        "flexppi_chi_full": "full",
        "flexppi_chi_pmhc": "pmhc",
        "flexppi_chi_tcr": "tcr",
        "flexppi_chi_isolated_A": "isolated_A",
        "flexppi_chi_isolated_B": "isolated_B",
        "flexppi_chi_isolated_C": "isolated_C",
        "flexppi_chi_isolated_D": "isolated_D",
        "flexppi_chi_isolated_E": "isolated_E",
    }
    condition_frame = frame.merge(chain_letter, on=["pdb_id", "chain_role"], how="left")
    for role, group_condition in role_groups.items():
        role_frame = condition_frame.loc[condition_frame["chain_role"].eq(role)]
        for base_condition in ("full", group_condition):
            column = f"flexppi_chi_{base_condition}"
            for (pdb_id, chain), group in role_frame.groupby(["pdb_id", "chain"], sort=True):
                metrics = chain_metrics(
                    group["target_chi_entropy"].to_numpy(float),
                    group[column].to_numpy(float),
                )
                if metrics is None:
                    continue
                condition_rows.append(
                    {
                        "pdb_id": pdb_id,
                        "chain_role": role,
                        "condition": base_condition,
                        "spearman": metrics["spearman"],
                        "residue_count": metrics["residue_count"],
                    }
                )
        # isolated-chain condition: role -> chain letter -> isolated_<letter>
        for (pdb_id, chain_letter), group in role_frame.groupby(["pdb_id", "chain_letter"], sort=True):
            column = f"flexppi_chi_isolated_{chain_letter}"
            if column not in condition_labels:
                continue
            metrics = chain_metrics(
                group["target_chi_entropy"].to_numpy(float),
                group[column].to_numpy(float),
            )
            if metrics is None:
                continue
            condition_rows.append(
                {
                    "pdb_id": pdb_id,
                    "chain_role": role,
                    "condition": column.replace("flexppi_chi_", ""),
                    "spearman": metrics["spearman"],
                    "residue_count": metrics["residue_count"],
                }
            )
    condition_per_chain = pd.DataFrame(condition_rows)
    condition_per_chain.to_csv(METRICS_DIR / "per_chain_condition_hchi.csv", index=False)
    condition_summary = []
    for (role, condition), group in condition_per_chain.groupby(["chain_role", "condition"]):
        values = group.loc[group["pdb_id"].isin(clean_ids), "spearman"].dropna()
        if values.empty:
            continue
        low, high = bootstrap_mean_ci(values.to_numpy(float))
        condition_summary.append(
            {
                "chain_role": role,
                "condition": condition,
                "n_complexes": int(len(values)),
                "mean_spearman": float(values.mean()),
                "ci95_low": low,
                "ci95_high": high,
            }
        )
    condition_summary = pd.DataFrame(condition_summary)
    condition_summary.to_csv(METRICS_DIR / "condition_summary.csv", index=False)

    # ---- target orthogonality: MD RMSF vs MD H_chi ----
    ortho_rows = []
    for (pdb_id, chain), group in frame.loc[subset_mask(frame, "overall")].groupby(
        ["pdb_id", "chain"], sort=True
    ):
        valid = group["target_chi_entropy"].notna() & group["target_chainfit_rmsf_angstrom"].notna()
        sub = group.loc[valid]
        if len(sub) >= 4:
            t = sub["target_chainfit_rmsf_angstrom"].to_numpy(float)
            h = sub["target_chi_entropy"].to_numpy(float)
            if not (np.all(t == t[0]) or np.all(h == h[0])):
                ortho_rows.append(
                    {
                        "pdb_id": pdb_id,
                        "chain": chain,
                        "spearman_md_rmsf_vs_md_hchi": float(spearmanr(t, h).statistic),
                        "residue_count": int(len(sub)),
                    }
                )
    ortho = pd.DataFrame(ortho_rows)
    ortho_complex = (
        ortho.groupby("pdb_id", as_index=False)["spearman_md_rmsf_vs_md_hchi"]
        .mean()
        .rename(columns={"spearman_md_rmsf_vs_md_hchi": "complex_orthogonality_spearman"})
    )
    ortho_complex = ortho_complex.loc[ortho_complex["pdb_id"].isin(clean_ids)]
    ortho_complex.to_csv(METRICS_DIR / "target_orthogonality.csv", index=False)

    # ---- report ----
    ext_overall = summary.loc[
        summary["cohort"].eq("external_clean") & summary["subset"].eq("overall")
    ]
    headline = {}
    for method in method_columns:
        row = ext_overall.loc[ext_overall["method"].eq(method) & ext_overall["metric"].eq("spearman")]
        if not row.empty:
            headline[method] = {
                "spearman": float(row["macro_mean"].iloc[0]),
                "ci95": [float(row["ci95_low"].iloc[0]), float(row["ci95_high"].iloc[0])],
            }
    report = {
        "protocol": (
            "Tier 1: FlexPPI H_chi head vs MD rotamer-entropy truth (direct). "
            "Tier 2: static flexibility scores as cross-descriptor proxies for side-chain "
            "entropy (not trained on H_chi). Per-chain metrics per subset (>=4 residues), "
            "equal-chain-weight per-complex mean, equal-complex macro mean with complex "
            "bootstrap 95% CI (2000, seed 2026), external_clean excludes 3KXF. "
            "Top-20% detection: positives = target >= 80th percentile within chain-subset."
        ),
        "external_clean_overall_spearman": headline,
        "flexppi_condition_ablation": condition_summary.to_dict("records"),
        "md_target_orthogonality": {
            "n_chains": int(len(ortho)),
            "mean_chain_spearman_md_rmsf_vs_md_hchi": float(ortho["spearman_md_rmsf_vs_md_hchi"].mean()),
        },
    }
    (OUTPUT_DIR / "itrace_hchi_multi_benchmark_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    print(json.dumps(report, ensure_ascii=False, indent=2))
    print()
    print("== external_clean overall: Spearman / top-20% AP / P@K ==")
    for method in method_columns:
        row = ext_overall.loc[ext_overall["method"].eq(method)]
        if row.empty:
            continue
        def val(metric):
            r = row.loc[row["metric"].eq(metric)]
            return f"{r['macro_mean'].iloc[0]:.3f}" if not r.empty else "  -"
        print(f"{method:<18} spearman={val('spearman')}  auroc={val('auroc')}  "
              f"ap={val('average_precision')}  pk={val('precision_at_k')}")
    print()
    print("== FlexPPI condition ablation, matched by chain role (external_clean) ==")
    for _, row in condition_summary.iterrows():
        print(f"{row['chain_role']:<10} {row['condition']:<10} n={int(row['n_complexes']):>3}  "
              f"mean={row['mean_spearman']:.3f}  CI=[{row['ci95_low']:.3f},{row['ci95_high']:.3f}]")


if __name__ == "__main__":
    main()
