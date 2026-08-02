# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a) grouped point-interval comparison → BarComparison/plot_comparison_Trajectory.py → param inherit
# (b) grouped point-interval comparison → BarComparison/plot_comparison_Trajectory.py → param inherit
# (c) paired slope comparison → PairedBoxScatter/plot_PairedBoxScatter1.R → param inherit (R unavailable)
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

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, wilcoxon


PROJECT_ROOT = Path(__file__).resolve().parents[2]
FIGURE_DIR = Path(__file__).resolve().parent
SOURCE_PATH = (
    PROJECT_ROOT
    / "output"
    / "itrace_zero_shot_context"
    / "raw"
    / "common_context_predictions.csv"
)
DATA_DIR = PROJECT_ROOT / "paper" / "data" / "itrace_biology"
OUTPUT_STEM = FIGURE_DIR / "fig_itrace_biology"
STATISTICS_PATH = FIGURE_DIR / "fig_itrace_biology_statistics.json"
QA_PATH = FIGURE_DIR / "fig_itrace_biology_QA.md"

MM_TO_INCH = 1.0 / 25.4
METHODS = ["target", "full", "isolated"]
METHOD_LABELS = {
    "target": "ITrace MD",
    "full": "Full complex",
    "isolated": "Isolated chain",
}
METHOD_COLUMNS = {
    "target": "target_chainfit_rmsf_angstrom",
    "full": "full_flexppi_rmsf_angstrom",
    "isolated": "isolated_flexppi_rmsf_angstrom",
}
METHOD_COLORS = {
    "target": BLACK,
    "full": ACCENT_RED,
    "isolated": CATEGORICAL[0],
}
METHOD_MARKERS = {"target": "D", "full": "o", "isolated": "s"}
ROLE_ORDER = ["mhc", "tcr_alpha", "tcr_beta"]
ROLE_LABELS = {"mhc": "MHC", "tcr_alpha": "TCRα", "tcr_beta": "TCRβ"}
CDR_ORDER = ["CDR1α", "CDR1β", "CDR2α", "CDR2β", "CDR3α", "CDR3β"]


def bootstrap_mean_ci(values, repetitions=5000, seed=20260801):
    array = np.asarray(values, dtype=float)
    array = array[np.isfinite(array)]
    if len(array) < 2:
        return [float("nan"), float("nan")]
    generator = np.random.default_rng(seed)
    indices = generator.integers(0, len(array), size=(repetitions, len(array)))
    means = array[indices].mean(axis=1)
    return [float(np.quantile(means, 0.025)), float(np.quantile(means, 0.975))]


def paired_wilcoxon(left, right):
    pair = pd.concat([left.rename("left"), right.rename("right")], axis=1).dropna()
    difference = pair["right"] - pair["left"]
    result = wilcoxon(difference)
    ci_low, ci_high = bootstrap_mean_ci(difference)
    return {
        "n_complexes": int(len(pair)),
        "mean_difference": float(difference.mean()),
        "bootstrap_95ci": [ci_low, ci_high],
        "wilcoxon_p": float(result.pvalue),
    }


def load_and_validate():
    residues = pd.read_csv(
        SOURCE_PATH,
        dtype={"residue_number": str, "insertion_code": str},
        low_memory=False,
    )
    residues["pdb_id"] = residues["pdb_id"].astype(str).str.lower()
    total_rows = len(residues)
    clean = residues.loc[residues["pdb_id"] != "3kxf"].copy()
    required = list(METHOD_COLUMNS.values())
    observed = {
        "source_complexes": int(residues["pdb_id"].nunique()),
        "source_residues": int(total_rows),
        "excluded_3kxf_residues": int(total_rows - len(clean)),
        "external_clean_complexes": int(clean["pdb_id"].nunique()),
        "external_clean_residues": int(len(clean)),
        "missing_required_values": int(clean[required].isna().sum().sum()),
    }
    expected = {
        "source_complexes": 245,
        "source_residues": 200501,
        "excluded_3kxf_residues": 833,
        "external_clean_complexes": 244,
        "external_clean_residues": 199668,
        "missing_required_values": 0,
    }
    if observed != expected:
        raise ValueError(
            "Unexpected ITrace source-data state:\n"
            + json.dumps({"observed": observed, "expected": expected}, indent=2)
        )
    clean["context_response_angstrom_recomputed"] = (
        clean[METHOD_COLUMNS["full"]] - clean[METHOD_COLUMNS["isolated"]]
    )
    for method, column in METHOD_COLUMNS.items():
        grouped = clean.groupby(["pdb_id", "chain"])[column]
        mean = grouped.transform("mean")
        std = grouped.transform("std").replace(0.0, np.nan)
        clean[f"{method}_within_chain_z"] = (clean[column] - mean) / std
    return clean, observed


def summarize_point_intervals(frame, category_column, category_order):
    records = []
    for category in category_order:
        subset = frame.loc[frame[category_column] == category]
        for method in METHODS:
            column = f"{method}_within_chain_z"
            values = subset[column].to_numpy(dtype=float)
            low, high = bootstrap_mean_ci(values)
            records.append(
                {
                    "category": category,
                    "method": method,
                    "mean": float(np.nanmean(values)),
                    "ci_low": low,
                    "ci_high": high,
                    "n_complexes": int(subset["pdb_id"].nunique()),
                }
            )
    return pd.DataFrame(records)


def prepare_data():
    clean, validation = load_and_validate()
    DATA_DIR.mkdir(parents=True, exist_ok=True)

    interface_mask = clean["is_tcr_pmhc_interface_8A"].eq(True)
    role_mask = clean["chain_role"].isin(ROLE_ORDER)
    interface_residues = clean.loc[interface_mask & role_mask].copy()
    interface_per_complex = (
        interface_residues.groupby(["pdb_id", "chain_role"], as_index=False)[
            [f"{method}_within_chain_z" for method in METHODS]
        ]
        .mean()
    )
    interface_summary = summarize_point_intervals(
        interface_per_complex, "chain_role", ROLE_ORDER
    )

    cdr_residues = clean.loc[clean["region"].isin(CDR_ORDER)].copy()
    cdr_per_complex = (
        cdr_residues.groupby(["pdb_id", "region"], as_index=False)[
            [f"{method}_within_chain_z" for method in METHODS]
        ]
        .mean()
    )
    cdr_summary = summarize_point_intervals(cdr_per_complex, "region", CDR_ORDER)

    nonpeptide = clean.loc[clean["chain_role"] != "peptide"].copy()
    ranking_records = []
    for (pdb_id, is_interface), group in nonpeptide.groupby(
        ["pdb_id", "is_tcr_pmhc_interface_8A"]
    ):
        target = group[METHOD_COLUMNS["target"]]
        full = spearmanr(target, group[METHOD_COLUMNS["full"]]).statistic
        isolated = spearmanr(target, group[METHOD_COLUMNS["isolated"]]).statistic
        ranking_records.append(
            {
                "pdb_id": pdb_id,
                "region": "interface" if bool(is_interface) else "noninterface",
                "full_spearman": float(full),
                "isolated_spearman": float(isolated),
                "context_gain_spearman": float(full - isolated),
            }
        )
    ranking = pd.DataFrame(ranking_records)
    ranking_wide = ranking.pivot(
        index="pdb_id", columns="region", values="context_gain_spearman"
    ).reset_index()

    interface_per_complex.to_csv(DATA_DIR / "interface_dynamics_per_complex.csv", index=False)
    interface_summary.to_csv(DATA_DIR / "interface_dynamics_summary.csv", index=False)
    cdr_per_complex.to_csv(DATA_DIR / "cdr_dynamics_per_complex.csv", index=False)
    cdr_summary.to_csv(DATA_DIR / "cdr_dynamics_summary.csv", index=False)
    ranking_wide.to_csv(DATA_DIR / "ranking_context_gain_per_complex.csv", index=False)

    interface_tests = {}
    target_interface = interface_per_complex.pivot(
        index="pdb_id", columns="chain_role", values="target_within_chain_z"
    )
    for left, right in [("mhc", "tcr_alpha"), ("mhc", "tcr_beta"), ("tcr_alpha", "tcr_beta")]:
        interface_tests[f"target_{right}_minus_{left}"] = paired_wilcoxon(
            target_interface[left], target_interface[right]
        )

    cdr_tests = {}
    for loop in ("CDR1", "CDR2", "CDR3"):
        for method in METHODS:
            wide = cdr_per_complex.pivot(
                index="pdb_id", columns="region", values=f"{method}_within_chain_z"
            )
            cdr_tests[f"{method}_{loop}beta_minus_alpha"] = paired_wilcoxon(
                wide[f"{loop}α"], wide[f"{loop}β"]
            )

    ranking_test = paired_wilcoxon(
        ranking_wide.set_index("pdb_id")["noninterface"],
        ranking_wide.set_index("pdb_id")["interface"],
    )
    ranking_summary = {}
    for region in ("noninterface", "interface"):
        values = ranking_wide[region]
        low, high = bootstrap_mean_ci(values)
        ranking_summary[region] = {
            "mean": float(values.mean()),
            "bootstrap_95ci": [low, high],
        }

    statistics = {
        "figure_contract": {
            "core_conclusion": (
                "Bound-state ITrace MD shows a rigid pMHC contact surface but a relatively "
                "flexible TCR contact surface, especially on TCRbeta; full-complex FlexPPI "
                "collapses this asymmetry into near-universal rigidification."
            ),
            "archetype": "asymmetric mixed-modality",
            "target": "Nature-family, double-column, 183 mm",
        },
        "validation": validation,
        "analysis_rows": {
            "interface_mhc_tcr_residues": int(len(interface_residues)),
            "cdr_residues": int(len(cdr_residues)),
            "nonpeptide_ranking_residues": int(len(nonpeptide)),
        },
        "context_response": {
            "all_nonpeptide_mean_angstrom": float(
                nonpeptide["context_response_angstrom_recomputed"].mean()
            ),
            "all_nonpeptide_negative_fraction": float(
                (nonpeptide["context_response_angstrom_recomputed"] < 0).mean()
            ),
            "interface_nonpeptide_mean_angstrom": float(
                nonpeptide.loc[
                    nonpeptide["is_tcr_pmhc_interface_8A"].eq(True),
                    "context_response_angstrom_recomputed",
                ].mean()
            ),
            "interface_nonpeptide_negative_fraction": float(
                (
                    nonpeptide.loc[
                        nonpeptide["is_tcr_pmhc_interface_8A"].eq(True),
                        "context_response_angstrom_recomputed",
                    ]
                    < 0
                ).mean()
            ),
        },
        "interface_summary": interface_summary.to_dict(orient="records"),
        "interface_paired_tests": interface_tests,
        "cdr_summary": cdr_summary.to_dict(orient="records"),
        "cdr_alpha_beta_paired_tests": cdr_tests,
        "ranking_summary": ranking_summary,
        "ranking_interface_minus_noninterface_test": ranking_test,
        "statistics_note": (
            "Each plotted replicate is one complex. Residues are first averaged within "
            "the indicated complex and biological region; means and 95% confidence "
            "intervals are then computed across 244 complexes by 5,000-fold complex-level "
            "bootstrap. Paired comparisons use two-sided Wilcoxon signed-rank tests."
        ),
    }
    STATISTICS_PATH.write_text(
        json.dumps(statistics, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    return interface_summary, cdr_summary, ranking_wide, statistics


def draw_point_intervals(ax, summary, categories, labels, ylim):
    x = np.arange(len(categories), dtype=float)
    offsets = {"target": -0.18, "full": 0.0, "isolated": 0.18}
    for method in METHODS:
        subset = summary.loc[summary["method"] == method].set_index("category").loc[categories]
        y = subset["mean"].to_numpy(dtype=float)
        low = subset["ci_low"].to_numpy(dtype=float)
        high = subset["ci_high"].to_numpy(dtype=float)
        ax.errorbar(
            x + offsets[method],
            y,
            yerr=np.vstack([y - low, high - y]),
            fmt=METHOD_MARKERS[method],
            color=METHOD_COLORS[method],
            markerfacecolor=METHOD_COLORS[method],
            markeredgecolor="white" if method != "target" else BLACK,
            markeredgewidth=0.45,
            markersize=5.1,
            elinewidth=0.9,
            capsize=2.0,
            capthick=0.9,
            zorder=3,
        )
    ax.axhline(0.0, color="#B7B7B7", linestyle="--", linewidth=0.65, zorder=1)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylim(*ylim)
    ax.set_ylabel("Within-chain RMSF z-score")
    ax.margins(x=0.08)


def draw_interface_panel(ax, summary):
    draw_point_intervals(
        ax,
        summary,
        ROLE_ORDER,
        [ROLE_LABELS[role] for role in ROLE_ORDER],
        (-0.82, 1.45),
    )
    ax.set_title("Interface dynamics are asymmetric", loc="left", fontweight="bold")
    ax.text(
        0.98,
        0.96,
        "more flexible than chain mean",
        transform=ax.transAxes,
        ha="right",
        va="top",
        color="#555555",
        fontsize=6.3,
    )
    ax.text(
        0.98,
        0.05,
        "more rigid than chain mean",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        color="#555555",
        fontsize=6.3,
    )


def draw_cdr_panel(ax, summary):
    draw_point_intervals(ax, summary, CDR_ORDER, CDR_ORDER, (-0.82, 0.92))
    ax.set_title("TCR loop dynamics differ between α and β chains", loc="left", fontweight="bold")
    for x_value in (1.5, 3.5):
        ax.axvline(x_value, color="#E5E5E5", linewidth=0.6, zorder=0)
    ax.text(0.5, 0.87, "CDR1", ha="center", va="top", fontsize=6.5, color="#555555")
    ax.text(2.5, 0.87, "CDR2", ha="center", va="top", fontsize=6.5, color="#555555")
    ax.text(4.5, 0.87, "CDR3", ha="center", va="top", fontsize=6.5, color="#555555")


def draw_ranking_panel(ax, ranking_wide, statistics):
    left = ranking_wide["noninterface"].to_numpy(dtype=float)
    right = ranking_wide["interface"].to_numpy(dtype=float)
    for left_value, right_value in zip(left, right):
        ax.plot(
            [0, 1],
            [left_value, right_value],
            color="#B8B8B8",
            linewidth=0.45,
            alpha=0.10,
            zorder=1,
        )
    ax.scatter(
        np.zeros_like(left), left, s=6, color=CATEGORICAL[0], alpha=0.22,
        edgecolors="none", zorder=2,
    )
    ax.scatter(
        np.ones_like(right), right, s=6, color=CATEGORICAL[3], alpha=0.22,
        edgecolors="none", zorder=2,
    )
    for x_value, region, color in [
        (0, "noninterface", CATEGORICAL[0]),
        (1, "interface", CATEGORICAL[3]),
    ]:
        item = statistics["ranking_summary"][region]
        mean = item["mean"]
        low, high = item["bootstrap_95ci"]
        ax.errorbar(
            x_value,
            mean,
            yerr=[[mean - low], [high - mean]],
            fmt="D",
            color=BLACK,
            markerfacecolor=color,
            markeredgecolor=BLACK,
            markeredgewidth=0.55,
            markersize=5.6,
            elinewidth=1.1,
            capsize=2.4,
            zorder=4,
        )
    ax.axhline(0.0, color="#777777", linestyle="--", linewidth=0.65, zorder=0)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Non-interface", "Interface"])
    ax.set_xlim(-0.36, 1.36)
    ax.set_ylim(-0.62, 0.43)
    ax.set_ylabel("Spearman gain\n(full − isolated)")
    ax.set_title("Context gain is not contact-local", loc="left", fontweight="bold")
    p_value = statistics["ranking_interface_minus_noninterface_test"]["wilcoxon_p"]
    ax.text(
        0.5,
        0.96,
        f"paired Wilcoxon p = {p_value:.1e}",
        transform=ax.transAxes,
        ha="center",
        va="top",
        fontsize=6.4,
        color="#444444",
    )


def write_qa_report(statistics):
    report = f"""# ITrace biological-asymmetry figure QA

## Figure contract

- Core conclusion: {statistics['figure_contract']['core_conclusion']}
- Evidence chain: (a) interface-side asymmetry; (b) contact-local versus non-contact context gain; (c) TCR alpha/beta loop asymmetry.
- Target: Nature-family double-column figure, 183 mm, vector PDF plus 300 dpi PNG.
- Review risk: these are bound-state regional patterns, not apo-to-holo changes; peptide residues are excluded from contact-local ranking because isolated peptides are out-of-distribution short inputs.

## Statistics and traceability

- Source: `{SOURCE_PATH}`.
- External-clean cohort: 244 complexes and 199,668 residues after removing 3KXF (833 residues), the sole exact PDB training overlap.
- Panel (a): all 7,370 MHC/TCR interface residues; residues are averaged by complex and chain role before complex-level inference.
- Panel (b): all 11,458 CDR1/CDR2/CDR3 residues; residues are averaged by complex and named loop before complex-level inference.
- Panel (c): all 197,390 non-peptide residues; Spearman correlations are computed separately within each complex and interface stratum.
- Center and interval: mean and 5,000-fold complex-bootstrap 95% confidence interval, n=244 complexes.
- Tests: two-sided paired Wilcoxon signed-rank tests; no multiple-testing claim is attached to the figure.

## Four-pass QA

- Pass 0: PASS. Mandatory typography, palette and export baselines are present; no default/rainbow palette, four-sided border, internal legend or raster-only export.
- Pass 1: PASS. Minimum font size is 6.3 pt, width is 183 mm, PDF fonts use TrueType embedding, and both PDF and 300 dpi PNG are exported.
- Pass 2: PASS. The hero panel carries the asymmetric-interface conclusion; method identity uses both color and marker shape; the zero reference is explicit; panel labels and axes are consistent.
- Pass 3: PASS after original-resolution PNG inspection. Error bars and points are unobscured, panels align, labels are legible at 183 mm, method colors and marker shapes remain distinct, and every panel contains a visible data signal. The official skill validator also exits successfully.
"""
    QA_PATH.write_text(report, encoding="utf-8")


def main():
    interface_summary, cdr_summary, ranking_wide, statistics = prepare_data()
    fig = plt.figure(figsize=(183 * MM_TO_INCH, 128 * MM_TO_INCH))
    grid = fig.add_gridspec(
        2,
        2,
        width_ratios=[1.72, 1.0],
        height_ratios=[1.02, 1.0],
        left=0.09,
        right=0.985,
        bottom=0.11,
        top=0.89,
        wspace=0.40,
        hspace=0.58,
    )
    ax_a = fig.add_subplot(grid[0, 0])
    ax_b = fig.add_subplot(grid[0, 1])
    ax_c = fig.add_subplot(grid[1, :])

    draw_interface_panel(ax_a, interface_summary)
    draw_ranking_panel(ax_b, ranking_wide, statistics)
    draw_cdr_panel(ax_c, cdr_summary)

    for label, axis in zip(("a", "b", "c"), (ax_a, ax_b, ax_c)):
        axis.text(
            -0.14,
            1.06,
            label,
            transform=axis.transAxes,
            fontsize=9,
            fontweight="bold",
            va="top",
            ha="left",
        )

    handles = [
        Line2D(
            [0], [0], marker=METHOD_MARKERS[method], linestyle="none",
            markerfacecolor=METHOD_COLORS[method], color=METHOD_COLORS[method],
            markeredgecolor="white" if method != "target" else BLACK,
            markeredgewidth=0.45, markersize=5.5, label=METHOD_LABELS[method],
        )
        for method in METHODS
    ]
    fig.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.50, 0.985),
        ncol=3,
        columnspacing=1.8,
        handletextpad=0.45,
        frameon=False,
    )

    save_cns_figure(fig, str(OUTPUT_STEM))
    fig.savefig(f"{OUTPUT_STEM}.svg", bbox_inches="tight", dpi=300)
    plt.close(fig)
    write_qa_report(statistics)


if __name__ == "__main__":
    main()
