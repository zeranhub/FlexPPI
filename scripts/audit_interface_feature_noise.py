"""Post-hoc noise, redundancy, and sampling-sensitivity audit.

This script consumes the trajectory descriptors produced by
``analyze_dynarepo_interface_dynamics.py``.  It does not import, train, or
modify FlexPPI.  Correlations are first computed within an apo/holo pair and
then averaged at PDB level so that long proteins do not dominate the result.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


DYNAMIC_FEATURES = [
    "interface_persistence_7a",
    "interface_persistence_8a",
    "interface_persistence_9a",
    "distance_mean",
    "distance_std",
    "contact_count_mean",
    "contact_count_std",
    "edge_entropy_mean",
    "edge_entropy_sum",
    "switching_edge_fraction",
    "stable_edge_fraction",
    "partner_identity_entropy",
    "initial_contact_retention",
]

STATIC_CONTROLS = [
    "static_distance",
    "static_distance_exp",
    "static_burial_proxy",
    "static_partner_count",
]

GLOBAL_FEATURES = [
    "state_entropy_025_norm",
    "state_entropy_035_norm",
    "state_entropy_045_norm",
    "contact_turnover",
    "consensus_dispersion",
    "total_contacts_cv",
    "pose_translation_rms",
    "pose_rotation_rms_deg",
    "pose_rmsd_rms",
]


def safe_spearman(x: pd.Series | np.ndarray, y: pd.Series | np.ndarray) -> float:
    x_array = np.asarray(x, dtype=np.float64)
    y_array = np.asarray(y, dtype=np.float64)
    valid = np.isfinite(x_array) & np.isfinite(y_array)
    if valid.sum() < 6 or np.nanstd(x_array[valid]) < 1e-12 or np.nanstd(y_array[valid]) < 1e-12:
        return float("nan")
    return float(spearmanr(x_array[valid], y_array[valid]).statistic)


def pdb_bootstrap_mean(
    table: pd.DataFrame,
    value_column: str,
    repetitions: int = 5000,
    seed: int = 17,
) -> tuple[float, float, float, int]:
    by_pdb = table.groupby("pdb_id", as_index=False)[value_column].mean(numeric_only=True)
    values = by_pdb[value_column].to_numpy(dtype=np.float64)
    values = values[np.isfinite(values)]
    if not len(values):
        return float("nan"), float("nan"), float("nan"), 0
    rng = np.random.default_rng(seed)
    draws = rng.choice(values, size=(repetitions, len(values)), replace=True).mean(axis=1)
    return (
        float(values.mean()),
        float(np.quantile(draws, 0.025)),
        float(np.quantile(draws, 0.975)),
        int(len(values)),
    )


def static_redundancy(residue: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for feature in DYNAMIC_FEATURES:
        for control in STATIC_CONTROLS:
            pair_rows = []
            for (pdb_id, pair_id), group in residue.groupby(["pdb_id", "pair_id"], sort=False):
                pair_rows.append({
                    "pdb_id": pdb_id,
                    "pair_id": pair_id,
                    "correlation": safe_spearman(group[feature], group[control]),
                })
            pair_table = pd.DataFrame(pair_rows)
            mean, low, high, count = pdb_bootstrap_mean(pair_table, "correlation")
            rows.append({
                "feature": feature,
                "static_control": control,
                "macro_spearman": mean,
                "ci_low": low,
                "ci_high": high,
                "pdbs": count,
            })
    return pd.DataFrame(rows)


def dynamic_redundancy(residue: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for index, first in enumerate(DYNAMIC_FEATURES):
        for second in DYNAMIC_FEATURES[index + 1 :]:
            pair_rows = []
            for (pdb_id, pair_id), group in residue.groupby(["pdb_id", "pair_id"], sort=False):
                pair_rows.append({
                    "pdb_id": pdb_id,
                    "pair_id": pair_id,
                    "correlation": safe_spearman(group[first], group[second]),
                })
            pair_table = pd.DataFrame(pair_rows)
            mean, low, high, count = pdb_bootstrap_mean(pair_table, "correlation")
            rows.append({
                "feature_1": first,
                "feature_2": second,
                "macro_spearman": mean,
                "ci_low": low,
                "ci_high": high,
                "pdbs": count,
            })
    return pd.DataFrame(rows)


def residue_stride_sensitivity(base: pd.DataFrame, dense: pd.DataFrame) -> pd.DataFrame:
    keys = ["pair_id", "pdb_id", "target_local_index"]
    merged = base[keys + DYNAMIC_FEATURES].merge(
        dense[keys + DYNAMIC_FEATURES],
        on=keys,
        suffixes=("_1ns", "_05ns"),
        validate="one_to_one",
    )
    rows: list[dict[str, object]] = []
    for feature in DYNAMIC_FEATURES:
        pair_rows = []
        for (pdb_id, pair_id), group in merged.groupby(["pdb_id", "pair_id"], sort=False):
            first = group[feature + "_1ns"].to_numpy(dtype=np.float64)
            second = group[feature + "_05ns"].to_numpy(dtype=np.float64)
            scale = max(float(np.nanstd(first)), 1e-8)
            pair_rows.append({
                "pdb_id": pdb_id,
                "pair_id": pair_id,
                "spearman": safe_spearman(first, second),
                "normalized_mae": float(np.nanmean(np.abs(first - second)) / scale),
            })
        pair_table = pd.DataFrame(pair_rows)
        spearman = pdb_bootstrap_mean(pair_table, "spearman")
        error = pdb_bootstrap_mean(pair_table, "normalized_mae")
        rows.append({
            "feature": feature,
            "macro_spearman": spearman[0],
            "spearman_ci_low": spearman[1],
            "spearman_ci_high": spearman[2],
            "normalized_mae": error[0],
            "normalized_mae_ci_low": error[1],
            "normalized_mae_ci_high": error[2],
            "pdbs": spearman[3],
        })
    return pd.DataFrame(rows)


def complex_stride_sensitivity(base: pd.DataFrame, dense: pd.DataFrame) -> pd.DataFrame:
    keys = ["pair_id", "pdb_id"]
    merged = base[keys + GLOBAL_FEATURES].merge(
        dense[keys + GLOBAL_FEATURES],
        on=keys,
        suffixes=("_1ns", "_05ns"),
        validate="one_to_one",
    )
    rows = []
    for feature in GLOBAL_FEATURES:
        first = merged[feature + "_1ns"].to_numpy(dtype=np.float64)
        second = merged[feature + "_05ns"].to_numpy(dtype=np.float64)
        denominator = np.maximum(np.maximum(np.abs(first), np.abs(second)), 1e-8)
        rows.append({
            "feature": feature,
            "spearman_across_pdbs": safe_spearman(first, second),
            "median_symmetric_relative_difference": float(np.median(np.abs(first - second) / denominator)),
            "pdbs": len(merged),
        })
    return pd.DataFrame(rows)


def cutoff_sensitivity(
    base_residue: pd.DataFrame,
    base_complex: pd.DataFrame,
    variants: dict[str, tuple[pd.DataFrame, pd.DataFrame]],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    residue_rows: list[dict[str, object]] = []
    complex_rows: list[dict[str, object]] = []
    residue_keys = ["pair_id", "pdb_id", "target_local_index"]
    complex_keys = ["pair_id", "pdb_id"]
    for cutoff, (variant_residue, variant_complex) in variants.items():
        merged_residue = base_residue[residue_keys + DYNAMIC_FEATURES].merge(
            variant_residue[residue_keys + DYNAMIC_FEATURES],
            on=residue_keys,
            suffixes=("_8a", "_variant"),
            validate="one_to_one",
        )
        for feature in DYNAMIC_FEATURES:
            pair_rows = []
            for (pdb_id, pair_id), group in merged_residue.groupby(["pdb_id", "pair_id"], sort=False):
                first = group[feature + "_8a"].to_numpy(dtype=np.float64)
                second = group[feature + "_variant"].to_numpy(dtype=np.float64)
                scale = max(float(np.nanstd(first)), 1e-8)
                pair_rows.append({
                    "pdb_id": pdb_id,
                    "pair_id": pair_id,
                    "spearman": safe_spearman(first, second),
                    "normalized_mae": float(np.nanmean(np.abs(first - second)) / scale),
                })
            pair_table = pd.DataFrame(pair_rows)
            correlation = pdb_bootstrap_mean(pair_table, "spearman")
            error = pdb_bootstrap_mean(pair_table, "normalized_mae")
            residue_rows.append({
                "contact_cutoff": cutoff,
                "feature": feature,
                "macro_spearman_vs_8a": correlation[0],
                "spearman_ci_low": correlation[1],
                "spearman_ci_high": correlation[2],
                "normalized_mae": error[0],
                "pdbs": correlation[3],
            })

        merged_complex = base_complex[complex_keys + GLOBAL_FEATURES].merge(
            variant_complex[complex_keys + GLOBAL_FEATURES],
            on=complex_keys,
            suffixes=("_8a", "_variant"),
            validate="one_to_one",
        )
        for feature in GLOBAL_FEATURES:
            first = merged_complex[feature + "_8a"].to_numpy(dtype=np.float64)
            second = merged_complex[feature + "_variant"].to_numpy(dtype=np.float64)
            denominator = np.maximum(np.maximum(np.abs(first), np.abs(second)), 1e-8)
            complex_rows.append({
                "contact_cutoff": cutoff,
                "feature": feature,
                "spearman_across_pdbs_vs_8a": safe_spearman(first, second),
                "median_symmetric_relative_difference": float(np.median(np.abs(first - second) / denominator)),
                "pdbs": len(merged_complex),
            })
    return pd.DataFrame(residue_rows), pd.DataFrame(complex_rows)


def label_stability(
    residue: pd.DataFrame,
    block_residue: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    keys = ["pdb_id", "pair_id", "target_local_index"]
    pivot = block_residue.pivot_table(
        index=keys,
        columns="block",
        values="delta_rmsf_block",
        aggfunc="first",
    )
    values = pivot.to_numpy(dtype=np.float64)
    count = np.sum(np.isfinite(values), axis=1)
    mean = np.nanmean(values, axis=1)
    standard_deviation = np.nanstd(values, axis=1, ddof=1)
    standard_error = standard_deviation / np.sqrt(np.maximum(count, 1))
    # Student-t 97.5% quantile for five blocks (four degrees of freedom).
    ci_half_width = 2.7764451051977987 * standard_error
    positive_fraction = np.nanmean(values > 0, axis=1)
    negative_fraction = np.nanmean(values < 0, axis=1)
    sign_agreement = np.maximum(positive_fraction, negative_fraction)
    all_same_sign = (positive_fraction == 1.0) | (negative_fraction == 1.0)
    ci_excludes_zero = (mean - ci_half_width > 0) | (mean + ci_half_width < 0)

    stability = pivot.reset_index()[keys].copy()
    stability["block_mean_delta_rmsf"] = mean
    stability["block_sd_delta_rmsf"] = standard_deviation
    stability["block_se_delta_rmsf"] = standard_error
    stability["ci95_low"] = mean - ci_half_width
    stability["ci95_high"] = mean + ci_half_width
    stability["sign_agreement"] = sign_agreement
    stability["all_five_blocks_same_sign"] = all_same_sign
    stability["ci95_excludes_zero"] = ci_excludes_zero
    annotations = residue[keys + ["static_distance", "delta_rmsf", "delta_rmsf_z"]]
    stability = stability.merge(annotations, on=keys, how="left", validate="one_to_one")
    stability["region"] = np.where(stability["static_distance"] <= 8.0, "interface_8a", "non_interface")

    per_pdb_rows = []
    for (pdb_id, region), group in stability.groupby(["pdb_id", "region"], sort=False):
        per_pdb_rows.append({
            "pdb_id": pdb_id,
            "region": region,
            "residues": len(group),
            "fraction_all_five_blocks_same_sign": float(group["all_five_blocks_same_sign"].mean()),
            "fraction_ci95_excludes_zero": float(group["ci95_excludes_zero"].mean()),
            "fraction_sign_agreement_ge_0_8": float((group["sign_agreement"] >= 0.8).mean()),
            "median_abs_delta_rmsf": float(np.median(np.abs(group["delta_rmsf"]))),
            "median_block_sd": float(np.median(group["block_sd_delta_rmsf"])),
        })
    per_pdb = pd.DataFrame(per_pdb_rows)

    summary_rows = []
    for region in ["all", "interface_8a", "non_interface"]:
        region_table = per_pdb if region == "all" else per_pdb[per_pdb["region"] == region]
        # For the all-residue summary, average the raw residue statistics within
        # each PDB rather than giving equal weight to interface/non-interface.
        if region == "all":
            raw_rows = []
            for pdb_id, group in stability.groupby("pdb_id", sort=False):
                raw_rows.append({
                    "pdb_id": pdb_id,
                    "fraction_all_five_blocks_same_sign": float(group["all_five_blocks_same_sign"].mean()),
                    "fraction_ci95_excludes_zero": float(group["ci95_excludes_zero"].mean()),
                    "fraction_sign_agreement_ge_0_8": float((group["sign_agreement"] >= 0.8).mean()),
                    "median_abs_delta_rmsf": float(np.median(np.abs(group["delta_rmsf"]))),
                    "median_block_sd": float(np.median(group["block_sd_delta_rmsf"])),
                })
            region_table = pd.DataFrame(raw_rows)
        row: dict[str, object] = {"region": region, "pdbs": int(region_table["pdb_id"].nunique())}
        for metric in [
            "fraction_all_five_blocks_same_sign",
            "fraction_ci95_excludes_zero",
            "fraction_sign_agreement_ge_0_8",
            "median_abs_delta_rmsf",
            "median_block_sd",
        ]:
            mean, low, high, _ = pdb_bootstrap_mean(region_table, metric)
            row[metric] = mean
            row[metric + "_ci_low"] = low
            row[metric + "_ci_high"] = high
        summary_rows.append(row)
    return stability, per_pdb, pd.DataFrame(summary_rows)


def binding_response_landscape(
    stability: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    annotated = stability.copy()
    annotated["robust_softening"] = annotated["ci95_low"] > 0
    annotated["robust_stiffening"] = annotated["ci95_high"] < 0
    annotated["distance_bin"] = pd.cut(
        annotated["static_distance"],
        bins=[-np.inf, 7.0, 9.0, 12.0, np.inf],
        labels=["<=7A", "7-9A", "9-12A", ">12A"],
        right=True,
    ).astype(str)

    per_pdb_rows = []
    for (pdb_id, region), group in annotated.groupby(["pdb_id", "region"], sort=False):
        per_pdb_rows.append({
            "pdb_id": pdb_id,
            "region": region,
            "residues": len(group),
            "mean_delta_rmsf": float(group["delta_rmsf"].mean()),
            "median_delta_rmsf": float(group["delta_rmsf"].median()),
            "fraction_stiffening": float((group["delta_rmsf"] < 0).mean()),
            "fraction_softening": float((group["delta_rmsf"] > 0).mean()),
            "fraction_robust_stiffening": float(group["robust_stiffening"].mean()),
            "fraction_robust_softening": float(group["robust_softening"].mean()),
            "has_both_robust_directions": bool(group["robust_stiffening"].any() and group["robust_softening"].any()),
        })
    per_pdb = pd.DataFrame(per_pdb_rows)

    summary_rows = []
    for region in ["interface_8a", "non_interface"]:
        region_table = per_pdb[per_pdb["region"] == region]
        row: dict[str, object] = {"region": region, "pdbs": int(region_table["pdb_id"].nunique())}
        for metric in [
            "mean_delta_rmsf",
            "median_delta_rmsf",
            "fraction_stiffening",
            "fraction_softening",
            "fraction_robust_stiffening",
            "fraction_robust_softening",
        ]:
            mean, low, high, _ = pdb_bootstrap_mean(region_table, metric)
            row[metric] = mean
            row[metric + "_ci_low"] = low
            row[metric + "_ci_high"] = high
        row["fraction_pdbs_with_both_robust_directions"] = float(region_table["has_both_robust_directions"].mean())
        summary_rows.append(row)

    distance_rows = []
    distance_order = ["<=7A", "7-9A", "9-12A", ">12A"]
    for distance_bin in distance_order:
        subset = annotated[annotated["distance_bin"] == distance_bin]
        pdb_rows = []
        for pdb_id, group in subset.groupby("pdb_id", sort=False):
            pdb_rows.append({
                "pdb_id": pdb_id,
                "residues": len(group),
                "mean_delta_rmsf": float(group["delta_rmsf"].mean()),
                "mean_abs_delta_rmsf": float(np.abs(group["delta_rmsf"]).mean()),
                "fraction_stiffening": float((group["delta_rmsf"] < 0).mean()),
                "fraction_robust_stiffening": float(group["robust_stiffening"].mean()),
                "fraction_robust_softening": float(group["robust_softening"].mean()),
            })
        pdb_table = pd.DataFrame(pdb_rows)
        row = {"distance_bin": distance_bin, "pdbs": int(pdb_table["pdb_id"].nunique()), "residues": len(subset)}
        for metric in [
            "mean_delta_rmsf",
            "mean_abs_delta_rmsf",
            "fraction_stiffening",
            "fraction_robust_stiffening",
            "fraction_robust_softening",
        ]:
            mean, low, high, _ = pdb_bootstrap_mean(pdb_table, metric)
            row[metric] = mean
            row[metric + "_ci_low"] = low
            row[metric + "_ci_high"] = high
        distance_rows.append(row)
    return per_pdb, pd.DataFrame(summary_rows), pd.DataFrame(distance_rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-dir", type=Path, default=Path("artifacts/dynarepo_interface_dynamics"))
    parser.add_argument("--dense-dir", type=Path, default=Path("artifacts/dynarepo_interface_dynamics_stride05"))
    parser.add_argument("--cutoff07-dir", type=Path, default=Path("artifacts/dynarepo_interface_dynamics_cutoff07"))
    parser.add_argument("--cutoff09-dir", type=Path, default=Path("artifacts/dynarepo_interface_dynamics_cutoff09"))
    args = parser.parse_args()

    residue = pd.read_csv(args.base_dir / "residue_features.csv")
    dense_residue = pd.read_csv(args.dense_dir / "residue_features.csv")
    block_residue = pd.read_csv(args.base_dir / "block_residue_features.csv")
    complex_table = pd.read_csv(args.base_dir / "complex_features.csv")
    dense_complex = pd.read_csv(args.dense_dir / "complex_features.csv")
    cutoff07_residue = pd.read_csv(args.cutoff07_dir / "residue_features.csv")
    cutoff07_complex = pd.read_csv(args.cutoff07_dir / "complex_features.csv")
    cutoff09_residue = pd.read_csv(args.cutoff09_dir / "residue_features.csv")
    cutoff09_complex = pd.read_csv(args.cutoff09_dir / "complex_features.csv")

    static_redundancy(residue).to_csv(args.base_dir / "static_redundancy.csv", index=False)
    dynamic_redundancy(residue).to_csv(args.base_dir / "dynamic_feature_redundancy.csv", index=False)
    residue_stride_sensitivity(residue, dense_residue).to_csv(
        args.base_dir / "stride_sensitivity_residue.csv", index=False
    )
    complex_stride_sensitivity(complex_table, dense_complex).to_csv(
        args.base_dir / "stride_sensitivity_complex.csv", index=False
    )
    cutoff_residue, cutoff_complex = cutoff_sensitivity(
        residue,
        complex_table,
        {
            "7a": (cutoff07_residue, cutoff07_complex),
            "9a": (cutoff09_residue, cutoff09_complex),
        },
    )
    cutoff_residue.to_csv(args.base_dir / "cutoff_sensitivity_residue.csv", index=False)
    cutoff_complex.to_csv(args.base_dir / "cutoff_sensitivity_complex.csv", index=False)
    stability, stability_per_pdb, stability_summary = label_stability(residue, block_residue)
    stability.to_csv(args.base_dir / "label_stability_residue.csv", index=False)
    stability_per_pdb.to_csv(args.base_dir / "label_stability_per_pdb.csv", index=False)
    stability_summary.to_csv(args.base_dir / "label_stability_summary.csv", index=False)
    response_per_pdb, response_summary, distance_response = binding_response_landscape(stability)
    response_per_pdb.to_csv(args.base_dir / "binding_response_per_pdb.csv", index=False)
    response_summary.to_csv(args.base_dir / "binding_response_summary.csv", index=False)
    distance_response.to_csv(args.base_dir / "distance_bin_response.csv", index=False)


if __name__ == "__main__":
    main()
