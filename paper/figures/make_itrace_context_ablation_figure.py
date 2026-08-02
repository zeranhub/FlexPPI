# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a) paired performance scatter → MarginalDensity/plot_VariableCorrelation.py → param inherit (plain paired scatter; marginal densities are not scientifically required)
# (b) metric-effect forest → BarComparison/plot_comparison_Trajectory.py → param inherit (confidence intervals replace bars)
# (c) chain-role effect forest → BarComparison/plot_comparison_Trajectory.py → param inherit (signed effects require a zero reference)
# (d) grouped RMSF-scale comparison → BarComparison/plot_comparison_Trajectory.py → param inherit (point intervals preserve the grouped comparison visual system)
# RULE: "native run" = load pre-rendered PNG via Image.open().ax.imshow().
#       "param inherit" = drawing function below that copies Class A/B/C values.
#       If a panel says "native run" and you write a drawing function, you broke the contract.

# Academic Figure Skill Typography Baseline — COPY VERBATIM, place at TOP of script
import matplotlib as mpl
mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans"],
    "font.size": 8,
    "axes.titlesize": 8,
    "axes.labelsize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 8,
    "figure.titlesize": 9,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.6,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "legend.frameon": False,
})

# Academic Figure Skill Nature/Cell/Science Color Palette -- COPY VERBATIM
CATEGORICAL = ["#2166AC", "#B2182B", "#1B7837", "#F1A340", "#762A83", "#666666"]
CATEGORICAL_EXTENDED = [
    "#2166AC", "#B2182B", "#1B7837", "#F1A340", "#762A83", "#666666",
    "#4393C3", "#D6604D", "#5AAE61", "#B35806", "#9970AB", "#999999",
]
DIVERGING   = ["#2166AC", "#F7F7F7", "#B2182B"]
SEQUENTIAL  = ["#F7FBFF", "#6BAED6", "#08306B"]
ACCENT_RED  = "#B2182B"
GREY        = "#999999"
BLACK       = "#222222"

# Academic Figure Skill Export Baseline — COPY VERBATIM
mpl.rcParams.update({
    "pdf.fonttype": 42,         # TrueType font embedding
    "svg.fonttype": "none",     # editable text in SVG
    "savefig.bbox": "tight",    # trim whitespace
    "savefig.dpi": 300,
})

def save_cns_figure(fig, filename):
    """Standard Academic Figure Skill export: vector PDF + 300dpi PNG preview."""
    fig.savefig(f"{filename}.pdf", bbox_inches="tight", dpi=300)
    fig.savefig(f"{filename}.png", bbox_inches="tight", dpi=300)


import json
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Sequence

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
FIGURE_DIR = Path(__file__).resolve().parent
DATA_DIR = PROJECT_ROOT / "paper" / "data" / "itrace_context"
CONTEXT_DIR = PROJECT_ROOT / "output" / "itrace_zero_shot_context"
OUTPUT_STEM = FIGURE_DIR / "fig_itrace_context_ablation"
STATISTICS_PATH = FIGURE_DIR / "fig_itrace_context_ablation_statistics.json"

MM_TO_INCH = 1.0 / 25.4
FULL_COLOR = CATEGORICAL[1]
ISOLATED_COLOR = CATEGORICAL[0]
TARGET_COLOR = BLACK
POSITIVE_COLOR = CATEGORICAL[1]
NEGATIVE_COLOR = CATEGORICAL[0]
LIGHT_GREY = "#D0D0D0"

METRIC_ORDER = [
    "spearman",
    "pearson",
    "z_mae",
    "z_rmse",
    "top20_auroc",
    "top20_average_precision",
    "top20_precision_at_k",
]
METRIC_LABELS = {
    "spearman": "Spearman ρ",
    "pearson": "Pearson r",
    "z_mae": "z-MAE",
    "z_rmse": "z-RMSE",
    "top20_auroc": "Top-20% AUROC",
    "top20_average_precision": "Top-20% AP",
    "top20_precision_at_k": "Top-20% P@K",
}
ROLE_ORDER = ["mhc", "b2m", "peptide", "tcr_alpha", "tcr_beta"]
ROLE_LABELS = {
    "mhc": "MHC",
    "b2m": "β2m",
    "peptide": "Peptide",
    "tcr_alpha": "TCRα",
    "tcr_beta": "TCRβ",
}
METHOD_ORDER = ["target", "full_complex", "isolated_chain"]
METHOD_LABELS = {
    "target": "ITrace MD",
    "full_complex": "Full complex",
    "isolated_chain": "Isolated chain",
}
METHOD_COLORS = {
    "target": TARGET_COLOR,
    "full_complex": FULL_COLOR,
    "isolated_chain": ISOLATED_COLOR,
}
METHOD_MARKERS = {
    "target": "D",
    "full_complex": "o",
    "isolated_chain": "s",
}


def load_json(path: Path) -> Dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def bootstrap_mean_ci(
    values: Sequence[float],
    repetitions: int = 2000,
    seed: int = 7,
) -> List[float]:
    array = np.asarray(values, dtype=float)
    array = array[np.isfinite(array)]
    if len(array) < 2:
        return [float("nan"), float("nan")]
    generator = np.random.default_rng(seed)
    means = generator.choice(
        array,
        size=(repetitions, len(array)),
        replace=True,
    ).mean(axis=1)
    return [float(np.quantile(means, 0.025)), float(np.quantile(means, 0.975))]


def benjamini_hochberg(p_values: Iterable[float]) -> np.ndarray:
    values = np.asarray(list(p_values), dtype=float)
    order = np.argsort(values)
    ranked = values[order]
    adjusted = ranked * len(values) / np.arange(1, len(values) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    output = np.empty_like(adjusted)
    output[order] = np.clip(adjusted, 0.0, 1.0)
    return output


def validate_inputs(
    per_complex: pd.DataFrame,
    per_chain: pd.DataFrame,
    residues: pd.DataFrame,
    excluded_ids: set[str],
) -> Dict[str, Any]:
    clean_complex = per_complex[~per_complex["pdb_id"].isin(excluded_ids)]
    clean_chain = per_chain[~per_chain["pdb_id"].isin(excluded_ids)]
    clean_residue = residues[~residues["pdb_id"].isin(excluded_ids)]
    required_residue = [
        "target_chainfit_rmsf_angstrom",
        "full_flexppi_rmsf_angstrom",
        "isolated_flexppi_rmsf_angstrom",
    ]
    checks = {
        "external_clean_complexes": int(clean_complex["pdb_id"].nunique()),
        "external_clean_chains": int(len(clean_chain)),
        "external_clean_residues": int(len(clean_residue)),
        "per_complex_unique": bool(not clean_complex["pdb_id"].duplicated().any()),
        "per_chain_unique": bool(
            not clean_chain.duplicated(["pdb_id", "chain"]).any()
        ),
        "missing_required_residue_values": int(
            clean_residue[required_residue].isna().sum().sum()
        ),
        "chain_roles": sorted(clean_chain["chain_role"].astype(str).unique()),
    }
    expected = {
        "external_clean_complexes": 244,
        "external_clean_chains": 1220,
        "external_clean_residues": 199668,
        "per_complex_unique": True,
        "per_chain_unique": True,
        "missing_required_residue_values": 0,
        "chain_roles": sorted(ROLE_ORDER),
    }
    if checks != expected:
        raise ValueError(
            "ITrace context source-data validation failed:\n"
            + json.dumps({"observed": checks, "expected": expected}, indent=2)
        )
    return checks


def prepare_source_data() -> Dict[str, Any]:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    report = load_json(CONTEXT_DIR / "itrace_context_report.json")
    per_complex = pd.read_csv(
        CONTEXT_DIR / "metrics" / "per_complex_metrics.csv"
    )
    per_chain = pd.read_csv(
        CONTEXT_DIR / "metrics" / "per_chain_metrics.csv"
    )
    residues = pd.read_csv(
        CONTEXT_DIR / "raw" / "common_context_predictions.csv",
        dtype={"residue_number": str, "insertion_code": str},
    )
    for frame in (per_complex, per_chain, residues):
        frame["pdb_id"] = frame["pdb_id"].astype(str).str.lower()
    excluded = {
        str(value).lower()
        for value in report["cohort"]["excluded_exact_training_overlap"]
    }
    validation = validate_inputs(per_complex, per_chain, residues, excluded)
    clean_complex = per_complex[~per_complex["pdb_id"].isin(excluded)].copy()
    clean_chain = per_chain[~per_chain["pdb_id"].isin(excluded)].copy()
    clean_residue = residues[~residues["pdb_id"].isin(excluded)].copy()

    panel_a = clean_complex[
        ["pdb_id", "full_spearman", "isolated_spearman"]
    ].copy()
    panel_a["full_advantage_spearman"] = (
        panel_a["full_spearman"] - panel_a["isolated_spearman"]
    )
    panel_a["winner"] = np.where(
        panel_a["full_advantage_spearman"] > 0,
        "full_complex",
        np.where(
            panel_a["full_advantage_spearman"] < 0,
            "isolated_chain",
            "tie",
        ),
    )
    panel_a.to_csv(DATA_DIR / "panel_a_per_complex_spearman.csv", index=False)

    result = report["external_clean_results"]
    panel_b_rows: List[Dict[str, Any]] = []
    for metric in METRIC_ORDER:
        paired = result["paired_full_vs_isolated"][metric]
        panel_b_rows.append(
            {
                "metric": metric,
                "metric_label": METRIC_LABELS[metric],
                "full_complex_macro_mean": result["full_complex"][metric][
                    "macro_mean"
                ],
                "isolated_chain_macro_mean": result["isolated_chain"][metric][
                    "macro_mean"
                ],
                "mean_advantage": paired["mean_advantage"],
                "ci_low": paired["complex_bootstrap_95_ci"][0],
                "ci_high": paired["complex_bootstrap_95_ci"][1],
                "win_rate": paired["win_rate"],
                "wilcoxon_two_sided_p": paired["wilcoxon_two_sided_p"],
                "n_complexes": paired["paired_complex_count"],
            }
        )
    panel_b = pd.DataFrame(panel_b_rows)
    panel_b["bh_fdr_q"] = benjamini_hochberg(
        panel_b["wilcoxon_two_sided_p"]
    )
    panel_b.to_csv(DATA_DIR / "panel_b_metric_advantages.csv", index=False)

    role_results = report["sensitivity_analysis"]["by_chain_role"]
    panel_c_rows: List[Dict[str, Any]] = []
    for role in ROLE_ORDER:
        role_result = role_results[role]
        paired = role_result["paired_full_vs_isolated"]["spearman"]
        panel_c_rows.append(
            {
                "chain_role": role,
                "chain_role_label": ROLE_LABELS[role],
                "full_complex_spearman": role_result["full_complex"]["spearman"][
                    "macro_mean"
                ],
                "isolated_chain_spearman": role_result["isolated_chain"][
                    "spearman"
                ]["macro_mean"],
                "mean_advantage": paired["mean_advantage"],
                "ci_low": paired["complex_bootstrap_95_ci"][0],
                "ci_high": paired["complex_bootstrap_95_ci"][1],
                "wilcoxon_two_sided_p": paired["wilcoxon_two_sided_p"],
                "n_complexes": paired["paired_complex_count"],
            }
        )
    panel_c = pd.DataFrame(panel_c_rows)
    panel_c["bh_fdr_q"] = benjamini_hochberg(
        panel_c["wilcoxon_two_sided_p"]
    )
    panel_c.to_csv(DATA_DIR / "panel_c_role_spearman.csv", index=False)

    method_columns = {
        "target": "target_chainfit_rmsf_angstrom",
        "full_complex": "full_flexppi_rmsf_angstrom",
        "isolated_chain": "isolated_flexppi_rmsf_angstrom",
    }
    panel_d_parts = []
    for method, column in method_columns.items():
        part = (
            clean_residue.groupby(["pdb_id", "chain_role"], sort=True)[column]
            .median()
            .rename("complex_role_median_rmsf_angstrom")
            .reset_index()
        )
        part["method"] = method
        panel_d_parts.append(part)
    panel_d = pd.concat(panel_d_parts, ignore_index=True)
    panel_d = panel_d[
        panel_d["chain_role"].astype(str).isin(ROLE_ORDER)
    ].copy()
    panel_d.to_csv(DATA_DIR / "panel_d_role_rmsf_scale.csv", index=False)

    scale_rows = []
    for (role, method), group in panel_d.groupby(
        ["chain_role", "method"],
        sort=False,
    ):
        values = group["complex_role_median_rmsf_angstrom"].to_numpy(dtype=float)
        interval = bootstrap_mean_ci(values)
        scale_rows.append(
            {
                "chain_role": role,
                "chain_role_label": ROLE_LABELS[str(role)],
                "method": method,
                "method_label": METHOD_LABELS[str(method)],
                "macro_mean_of_complex_medians": float(np.mean(values)),
                "macro_median_of_complex_medians": float(np.median(values)),
                "ci_low": interval[0],
                "ci_high": interval[1],
                "n_complexes": int(len(values)),
            }
        )
    panel_d_summary = pd.DataFrame(scale_rows)
    panel_d_summary.to_csv(
        DATA_DIR / "panel_d_role_rmsf_scale_summary.csv",
        index=False,
    )

    clean_complex.to_csv(DATA_DIR / "per_complex_metrics_external_clean.csv", index=False)
    clean_chain.to_csv(DATA_DIR / "per_chain_metrics_external_clean.csv", index=False)
    return {
        "report": report,
        "validation": validation,
        "panel_a": panel_a,
        "panel_b": panel_b,
        "panel_c": panel_c,
        "panel_d": panel_d,
        "panel_d_summary": panel_d_summary,
    }


def add_panel_label(axis: Any, label: str) -> None:
    axis.text(
        -0.15,
        1.06,
        label,
        transform=axis.transAxes,
        fontsize=9,
        fontweight="bold",
        va="bottom",
        ha="left",
    )


def style_axis(axis: Any) -> None:
    axis.spines["left"].set_linewidth(0.6)
    axis.spines["bottom"].set_linewidth(0.6)
    axis.tick_params(length=3, width=0.6)


def asymmetric_error(value: float, low: float, high: float) -> np.ndarray:
    return np.asarray([[value - low], [high - value]], dtype=float)


def panel_a_scatter(axis: Any, data: pd.DataFrame, report: Mapping[str, Any]) -> None:
    positive = data["full_advantage_spearman"] > 0
    axis.scatter(
        data.loc[~positive, "isolated_spearman"],
        data.loc[~positive, "full_spearman"],
        s=15,
        color=ISOLATED_COLOR,
        marker="s",
        alpha=0.60,
        edgecolors="white",
        linewidths=0.25,
        rasterized=True,
        label="Isolated chain better",
        zorder=2,
    )
    axis.scatter(
        data.loc[positive, "isolated_spearman"],
        data.loc[positive, "full_spearman"],
        s=15,
        color=FULL_COLOR,
        marker="o",
        alpha=0.60,
        edgecolors="white",
        linewidths=0.25,
        rasterized=True,
        label="Full complex better",
        zorder=3,
    )
    lower = float(
        min(data["isolated_spearman"].min(), data["full_spearman"].min())
    )
    upper = float(
        max(data["isolated_spearman"].max(), data["full_spearman"].max())
    )
    margin = 0.05 * (upper - lower)
    limits = (lower - margin, upper + margin)
    axis.plot(
        limits,
        limits,
        color=GREY,
        linestyle="--",
        linewidth=0.7,
        alpha=0.8,
        zorder=1,
    )
    paired = report["external_clean_results"]["paired_full_vs_isolated"][
        "spearman"
    ]
    annotation = (
        f"Δρ = {paired['mean_advantage']:+.3f} "
        f"[{paired['complex_bootstrap_95_ci'][0]:.3f}, "
        f"{paired['complex_bootstrap_95_ci'][1]:.3f}]; "
        f"win = {100 * paired['win_rate']:.1f}%"
    )
    axis.text(
        0.04,
        0.06,
        annotation,
        transform=axis.transAxes,
        va="bottom",
        ha="left",
        fontsize=6.5,
        linespacing=1.3,
        bbox={
            "facecolor": "white",
            "edgecolor": "none",
            "alpha": 0.88,
            "pad": 1.5,
        },
        zorder=4,
    )
    axis.set_xlim(limits)
    axis.set_ylim(limits)
    axis.set_aspect("equal", adjustable="box")
    axis.set_xlabel("Isolated-chain Spearman ρ")
    axis.set_ylabel("Full-complex Spearman ρ")
    axis.set_title("Paired complex-level performance", loc="left")
    style_axis(axis)


def effect_forest(
    axis: Any,
    data: pd.DataFrame,
    label_column: str,
    title: str,
    x_label: str,
    color_by_sign: bool,
) -> None:
    ordered = data.reset_index(drop=True)
    y = np.arange(len(ordered))[::-1]
    axis.axvline(
        0.0,
        color=GREY,
        linestyle="--",
        linewidth=0.7,
        alpha=0.8,
        zorder=0,
    )
    for position, (_, row) in zip(y, ordered.iterrows()):
        value = float(row["mean_advantage"])
        color = (
            POSITIVE_COLOR
            if not color_by_sign or value >= 0
            else NEGATIVE_COLOR
        )
        marker = "o" if value >= 0 else "s"
        axis.errorbar(
            value,
            position,
            xerr=asymmetric_error(
                value,
                float(row["ci_low"]),
                float(row["ci_high"]),
            ),
            fmt=marker,
            color=color,
            markerfacecolor=color,
            markeredgecolor="white",
            markeredgewidth=0.35,
            markersize=5.2,
            elinewidth=1.2,
            capsize=2.2,
            capthick=0.8,
            zorder=2,
        )
    axis.set_yticks(y)
    axis.set_yticklabels(ordered[label_column])
    axis.set_xlabel(x_label)
    axis.set_title(title, loc="left")
    axis.margins(x=0.12, y=0.12)
    style_axis(axis)


def panel_d_scale(axis: Any, summary: pd.DataFrame) -> None:
    base = np.arange(len(ROLE_ORDER), dtype=float)
    offsets = {"target": -0.18, "full_complex": 0.0, "isolated_chain": 0.18}
    for method in METHOD_ORDER:
        ordered = (
            summary[summary["method"] == method]
            .set_index("chain_role")
            .loc[ROLE_ORDER]
            .reset_index()
        )
        values = ordered["macro_mean_of_complex_medians"].to_numpy(dtype=float)
        lows = ordered["ci_low"].to_numpy(dtype=float)
        highs = ordered["ci_high"].to_numpy(dtype=float)
        axis.errorbar(
            base + offsets[method],
            values,
            yerr=np.vstack([values - lows, highs - values]),
            fmt=METHOD_MARKERS[method],
            color=METHOD_COLORS[method],
            markerfacecolor=(
                "white" if method == "target" else METHOD_COLORS[method]
            ),
            markeredgecolor=METHOD_COLORS[method],
            markeredgewidth=0.8,
            markersize=4.8,
            elinewidth=1.0,
            capsize=2.0,
            capthick=0.8,
            label=METHOD_LABELS[method],
            zorder=3,
        )
    axis.set_yscale("log")
    axis.set_ylim(0.25, 16.0)
    axis.set_xticks(base)
    axis.set_xticklabels([ROLE_LABELS[role] for role in ROLE_ORDER])
    axis.set_ylabel("Complex-level median RMSF (Å, log scale)")
    axis.set_title("Absolute-scale stress test", loc="left")
    axis.grid(
        axis="y",
        which="major",
        color="#E0E0E0",
        linewidth=0.35,
        alpha=0.6,
        zorder=0,
    )
    axis.legend(
        loc="upper left",
        bbox_to_anchor=(0.0, 1.01),
        ncol=3,
        fontsize=6.5,
        handletextpad=0.3,
        columnspacing=0.8,
        borderaxespad=0.0,
    )
    style_axis(axis)


def make_figure(data: Mapping[str, Any]) -> None:
    figure = plt.figure(
        figsize=(183 * MM_TO_INCH, 132 * MM_TO_INCH),
        constrained_layout=False,
    )
    grid = figure.add_gridspec(
        2,
        2,
        width_ratios=[1.08, 0.92],
        height_ratios=[1.0, 1.0],
        left=0.085,
        right=0.985,
        bottom=0.10,
        top=0.95,
        wspace=0.38,
        hspace=0.47,
    )
    axes = [
        figure.add_subplot(grid[0, 0]),
        figure.add_subplot(grid[0, 1]),
        figure.add_subplot(grid[1, 0]),
        figure.add_subplot(grid[1, 1]),
    ]
    panel_a_scatter(axes[0], data["panel_a"], data["report"])
    effect_forest(
        axes[1],
        data["panel_b"],
        "metric_label",
        "Cross-metric consistency",
        "Full-complex advantage",
        color_by_sign=False,
    )
    effect_forest(
        axes[2],
        data["panel_c"],
        "chain_role_label",
        "Chain-role heterogeneity",
        "Full-complex Spearman advantage",
        color_by_sign=True,
    )
    panel_d_scale(axes[3], data["panel_d_summary"])
    for label, axis in zip("abcd", axes):
        add_panel_label(axis, label)
    save_cns_figure(figure, str(OUTPUT_STEM))
    figure.savefig(
        f"{OUTPUT_STEM}.svg",
        bbox_inches="tight",
        dpi=300,
    )
    plt.close(figure)


def write_statistics(data: Mapping[str, Any]) -> None:
    report = data["report"]
    main = report["external_clean_results"]
    non_peptide = report["sensitivity_analysis"]["exclude_peptide_chain"]
    statistics = {
        "figure": "ITrace frozen-checkpoint full-complex context ablation",
        "core_conclusion": (
            "The same frozen FlexPPI checkpoint benefits from full-complex "
            "context overall, but the gain is chain-role dependent and the "
            "isolated peptide condition exhibits severe length extrapolation."
        ),
        "cohort": {
            "independent_unit": "PDB complex",
            "n_complexes": 244,
            "n_chains": 1220,
            "n_residues": 199668,
            "excluded_exact_training_pdb": ["3kxf"],
        },
        "inference": {
            "parameter_updates": 0,
            "checkpoint_sha256": report["model"]["checkpoint_sha256"],
            "full_input": "bound-state five-chain complex",
            "ablation_input": (
                "one chain extracted from the same bound-state topology with "
                "retained atom coordinates unchanged"
            ),
            "target_used_during_inference": False,
        },
        "statistics": {
            "center": (
                "metric computed within chain, equal-chain mean within complex, "
                "then equal-complex macro mean"
            ),
            "interval": "2,000 complex-level bootstrap resamples, percentile 95% CI",
            "paired_test": "two-sided Wilcoxon signed-rank test across complexes",
            "multiple_comparison": (
                "Benjamini-Hochberg FDR within the seven panel-b metrics and "
                "within the five panel-c chain roles"
            ),
        },
        "main_results": main,
        "non_peptide_sensitivity": non_peptide,
        "panel_b_source": data["panel_b"].to_dict("records"),
        "panel_c_source": data["panel_c"].to_dict("records"),
        "panel_d_summary": data["panel_d_summary"].to_dict("records"),
        "data_validation": data["validation"],
        "source_files": {
            "context_report": str(
                CONTEXT_DIR / "itrace_context_report.json"
            ),
            "per_complex_metrics": str(
                CONTEXT_DIR / "metrics" / "per_complex_metrics.csv"
            ),
            "per_chain_metrics": str(
                CONTEXT_DIR / "metrics" / "per_chain_metrics.csv"
            ),
            "residue_predictions": str(
                CONTEXT_DIR / "raw" / "common_context_predictions.csv"
            ),
        },
    }
    STATISTICS_PATH.write_text(
        json.dumps(statistics, indent=2, ensure_ascii=False, allow_nan=False),
        encoding="utf-8",
    )


def main() -> None:
    data = prepare_source_data()
    make_figure(data)
    write_statistics(data)
    print(
        json.dumps(
            {
                "figure_pdf": str(OUTPUT_STEM.with_suffix(".pdf")),
                "figure_png": str(OUTPUT_STEM.with_suffix(".png")),
                "figure_svg": str(OUTPUT_STEM.with_suffix(".svg")),
                "statistics": str(STATISTICS_PATH),
                "data_validation": data["validation"],
            },
            ensure_ascii=False,
        )
    )


if __name__ == "__main__":
    main()
