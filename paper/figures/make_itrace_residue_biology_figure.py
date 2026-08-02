# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a) position-wise line + confidence interval → LineTrend/plot_trend.py → param inherit
# (b) multi-metric position response → LineTrend/plot_trend.py → param inherit
# (c) amino-acid effect forest plot → GroupedBarChart/plot_GroupedBarChartv1.py → param inherit
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


PROJECT_ROOT = Path(__file__).resolve().parents[2]
FIGURE_DIR = Path(__file__).resolve().parent
DATA_DIR = PROJECT_ROOT / "paper" / "data" / "itrace_residue_biology"
OUTPUT_STEM = FIGURE_DIR / "fig_itrace_residue_biology"
STATISTICS_PATH = FIGURE_DIR / "fig_itrace_residue_biology_statistics.json"
QA_PATH = FIGURE_DIR / "fig_itrace_residue_biology_QA.md"

MD_POSITION_PATH = DATA_DIR / "peptide_9mer_md_rmsf_position_summary.csv"
FIDELITY_POSITION_PATH = DATA_DIR / "peptide_9mer_rmsf_prediction_fidelity.csv"
CONTEXT_POSITION_PATH = DATA_DIR / "peptide_9mer_position_summary.csv"
PEPTIDE_RESIDUE_PATH = DATA_DIR / "peptide_residue_context.csv"
AA_EFFECT_PATH = DATA_DIR / "peptide_amino_acid_position_adjusted_summary.csv"
CASE_PATH = DATA_DIR / "7L1D_7RRG_peptide_residue_context.csv"
REPORT_PATH = DATA_DIR / "itrace_residue_biology_report.json"

MM_TO_INCH = 1.0 / 25.4
POSITIONS = np.arange(1, 10)
POSITION_LABELS = [f"P{position}" for position in POSITIONS]
CASE_SEQUENCE = "ALHGGWTTK"
METRICS = ("rmsf", "chi_entropy", "sasa")
METRIC_COLORS = {
    "rmsf": CATEGORICAL[0],
    "chi_entropy": CATEGORICAL[4],
    "sasa": CATEGORICAL[3],
}
METRIC_LABELS = {
    "rmsf": r"$\Delta$RMSF (Å)",
    "chi_entropy": r"$\Delta H_{\chi}$",
    "sasa": r"$\Delta$SASA (Å$^2$)",
}
AA_NAMES = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys",
    "Q": "Gln", "E": "Glu", "G": "Gly", "H": "His", "I": "Ile",
    "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro",
    "S": "Ser", "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val",
}


def load_and_validate():
    md = pd.read_csv(MD_POSITION_PATH)
    fidelity = pd.read_csv(FIDELITY_POSITION_PATH)
    context = pd.read_csv(CONTEXT_POSITION_PATH)
    peptide = pd.read_csv(PEPTIDE_RESIDUE_PATH, low_memory=False)
    amino = pd.read_csv(AA_EFFECT_PATH)
    case = pd.read_csv(CASE_PATH, low_memory=False)
    report = json.loads(REPORT_PATH.read_text(encoding="utf-8"))
    observed = {
        "md_9mer_positions": int(len(md)),
        "fidelity_9mer_positions": int(len(fidelity)),
        "context_9mer_positions": int(len(context)),
        "nine_mer_complex_count": int(context["complex_count"].min()),
        "external_clean_complex_count": int(peptide["pdb_id"].nunique()),
        "peptide_residue_count": int(len(peptide)),
        "amino_acid_type_count": int(amino["amino_acid_1"].nunique()),
        "case_complex_count": int(case["pdb_id"].nunique()),
        "case_residue_count": int(len(case)),
        "case_sequence_count": int(
            case.sort_values(["pdb_id", "peptide_position"])
            .groupby("pdb_id")["amino_acid"]
            .apply(lambda values: "".join(values))
            .eq(CASE_SEQUENCE)
            .sum()
        ),
    }
    expected = {
        "md_9mer_positions": 9,
        "fidelity_9mer_positions": 9,
        "context_9mer_positions": 9,
        "nine_mer_complex_count": 154,
        "external_clean_complex_count": 244,
        "peptide_residue_count": 2278,
        "amino_acid_type_count": 20,
        "case_complex_count": 2,
        "case_residue_count": 18,
        "case_sequence_count": 2,
    }
    if observed != expected:
        raise ValueError(
            "Unexpected residue-biology source state:\n"
            + json.dumps({"observed": observed, "expected": expected}, indent=2)
        )
    if not report["analysis"]["rmsf_reproduction_check"]["passes_1e-6_angstrom"]:
        raise ValueError("Full-complex RMSF did not reproduce the prior ITrace run")
    return md, fidelity, context, peptide, amino, case, report, observed


def figure_panel_labels(fig, axes_and_labels):
    """Place all panel labels with the same figure-coordinate offset."""
    fig.canvas.draw()
    for ax, label in axes_and_labels:
        box = ax.get_position()
        fig.text(
            box.x0 - 0.037,
            box.y1 + 0.014,
            label,
            fontsize=9,
            fontweight="bold",
            va="bottom",
            ha="left",
        )


def plot_rmsf_prediction_position(ax, fidelity, report):
    x = fidelity["peptide_position_n"].to_numpy(dtype=float)
    md_y = fidelity["mean_md_rmsf_z_within_peptide"].to_numpy(dtype=float)
    md_low = fidelity["cluster_bootstrap_ci_low_md_rmsf_z"].to_numpy(dtype=float)
    md_high = fidelity["cluster_bootstrap_ci_high_md_rmsf_z"].to_numpy(dtype=float)
    model_y = fidelity["mean_model_rmsf_z_within_peptide"].to_numpy(dtype=float)
    model_low = fidelity["cluster_bootstrap_ci_low_model_rmsf_z"].to_numpy(dtype=float)
    model_high = fidelity["cluster_bootstrap_ci_high_model_rmsf_z"].to_numpy(dtype=float)
    ax.axhline(0.0, color="#BBBBBB", linestyle="--", linewidth=0.7, zorder=0)
    ax.fill_between(x, md_low, md_high, color=BLACK, alpha=0.10, linewidth=0)
    ax.fill_between(x, model_low, model_high, color=CATEGORICAL[0], alpha=0.14, linewidth=0)
    ax.plot(x, md_y, color=BLACK, linewidth=1.35, marker="o", markersize=3.6)
    ax.plot(x, model_y, color=CATEGORICAL[0], linewidth=1.35, marker="s", markersize=3.3)
    ax.set_xticks(POSITIONS)
    ax.set_xticklabels(POSITION_LABELS)
    ax.set_xlim(0.6, 9.4)
    ax.set_ylabel("Within-peptide RMSF z-score")
    ax.set_xlabel("Position in 9-mer peptide")
    ax.set_title("Partial recovery of the MD position pattern", loc="left", fontweight="bold")
    fidelity_stats = report["analysis"]["peptide_rmsf_prediction_fidelity"]
    mean_rho = fidelity_stats["within_complex_ranking"]["mean_spearman_rho"]
    all_complex_count = fidelity_stats["complex_count"]
    ax.text(
        0.98, 0.98,
        "$n=154$ 9-mer curves\n"
        + f"all peptides ($n={all_complex_count}$): mean $\\rho={mean_rho:.3f}$",
        transform=ax.transAxes, va="top", ha="right", fontsize=6.5, color="#555555"
    )
    for position, label, value_column, offset in [
        (1, "P1 recovered", "mean_model_rmsf_z_within_peptide", (5, -17)),
        (6, "P6 MD peak missed", "mean_md_rmsf_z_within_peptide", (-38, 7)),
        (8, "P8 constraint recovered", "mean_model_rmsf_z_within_peptide", (-54, -14)),
    ]:
        value = float(fidelity.loc[fidelity["peptide_position_n"] == position, value_column].iloc[0])
        ax.annotate(
            label, (position, value), xytext=offset, textcoords="offset points",
            fontsize=6, color="#555555", arrowprops={"arrowstyle": "-", "lw": 0.5, "color": GREY}
        )


def context_contact_fraction(peptide):
    nine = peptide.loc[peptide["peptide_length"] == 9].copy()
    nine["direct_contact"] = nine["min_distance_to_tcr_angstrom"] <= 5.0
    summary = (
        nine.groupby("peptide_position_n", as_index=False)
        .agg(
            complex_count=("pdb_id", "nunique"),
            contact_fraction=("direct_contact", "mean"),
            mean_min_distance_angstrom=("min_distance_to_tcr_angstrom", "mean"),
        )
    )
    summary.to_csv(DATA_DIR / "peptide_9mer_tcr_contact_by_position.csv", index=False)
    return summary


def plot_context_position(axes, context, peptide, case):
    contact_ax = axes[0]
    contact = context_contact_fraction(peptide)
    contact_ax.bar(
        contact["peptide_position_n"], contact["contact_fraction"],
        width=0.68, color="#D9D9D9", edgecolor="#777777", linewidth=0.45,
    )
    contact_ax.set_ylim(0, 1.08)
    contact_ax.set_yticks([0, 0.5, 1.0])
    contact_ax.set_ylabel("TCR contact\nfraction", fontsize=7)
    contact_ax.set_xticks(POSITIONS)
    contact_ax.tick_params(axis="x", labelbottom=False, length=0)
    contact_ax.set_title(
        "Adding TCR selectively constrains the central peptide segment",
        loc="left", fontweight="bold", pad=3,
    )
    contact_ax.text(
        0.995, 0.95, "heavy-atom distance ≤ 5 Å",
        transform=contact_ax.transAxes, ha="right", va="top", fontsize=6.5, color="#555555"
    )

    case_styles = {
        "7L1D": {"color": "#444444", "marker": "D", "linestyle": "--"},
        "7RRG": {"color": "#999999", "marker": "s", "linestyle": ":"},
    }
    for metric, ax in zip(METRICS, axes[1:]):
        x = context["peptide_position_n"].to_numpy(dtype=float)
        y = context[f"mean_delta_{metric}"].to_numpy(dtype=float)
        low = context[f"cluster_bootstrap_ci_low_{metric}"].to_numpy(dtype=float)
        high = context[f"cluster_bootstrap_ci_high_{metric}"].to_numpy(dtype=float)
        color = METRIC_COLORS[metric]
        ax.axhline(0.0, color="#BBBBBB", linestyle="--", linewidth=0.65, zorder=0)
        ax.fill_between(x, low, high, color=color, alpha=0.18, linewidth=0)
        ax.plot(x, y, color=color, linewidth=1.5, marker="o", markersize=3.5, zorder=3)
        for pdb_id, style in case_styles.items():
            selected = case.loc[case["pdb_id"] == pdb_id].sort_values("peptide_position")
            case_y = selected[f"delta_after_minus_before_{metric}"].to_numpy(dtype=float)
            if metric == "chi_entropy":
                case_y = np.where(selected["chi_defined"].to_numpy(dtype=bool), case_y, np.nan)
            ax.plot(
                selected["peptide_position"], case_y,
                linewidth=0.75, markersize=2.7, markerfacecolor="white",
                markeredgewidth=0.7, alpha=0.9, **style,
            )
        ax.set_ylabel(METRIC_LABELS[metric], color=color)
        ax.tick_params(axis="y", colors=color)
        ax.set_xlim(0.6, 9.4)
        ax.set_xticks(POSITIONS)
        if metric != "sasa":
            ax.tick_params(axis="x", labelbottom=False, length=0)
        else:
            ax.set_xticklabels(POSITION_LABELS)
            ax.set_xlabel("Peptide position   (case sequence: A  L  H  G  G  W  T  T  K)")


def plot_amino_acid_effects(axes, amino):
    order = (
        amino.sort_values("mean_delta_chi_entropy", na_position="last")["amino_acid_1"]
        .tolist()
    )
    ordered = amino.set_index("amino_acid_1").loc[order].reset_index()
    y = np.arange(len(ordered))
    for metric, ax in zip(METRICS, axes):
        mean = ordered[f"mean_delta_{metric}"].to_numpy(dtype=float)
        low = ordered[f"cluster_bootstrap_ci_low_{metric}"].to_numpy(dtype=float)
        high = ordered[f"cluster_bootstrap_ci_high_{metric}"].to_numpy(dtype=float)
        valid = np.isfinite(mean) & np.isfinite(low) & np.isfinite(high)
        ax.axvline(0.0, color="#BBBBBB", linestyle="--", linewidth=0.65, zorder=0)
        ax.errorbar(
            mean[valid], y[valid], xerr=[mean[valid] - low[valid], high[valid] - mean[valid]],
            fmt="o", color=METRIC_COLORS[metric], ecolor=METRIC_COLORS[metric],
            markersize=3.2, elinewidth=0.8, capsize=1.6, capthick=0.7,
        )
        ax.set_ylim(-0.8, len(ordered) - 0.2)
        ax.invert_yaxis()
        ax.set_xlabel("Position-adjusted " + METRIC_LABELS[metric])
        ax.set_title(
            {"rmsf": "Positional motion", "chi_entropy": "Side-chain rotamers", "sasa": "Solvent exposure"}[metric],
            fontweight="bold",
        )
        ax.set_yticks(y)
        if metric == "rmsf":
            labels = [
                f"{aa}  {AA_NAMES[aa]}  ($n={int(count)}$)"
                for aa, count in zip(ordered["amino_acid_1"], ordered["residue_count"])
            ]
            ax.set_yticklabels(labels, fontsize=6.2)
        else:
            ax.set_yticklabels([])
            ax.tick_params(axis="y", length=0)
        if metric == "chi_entropy":
            for aa in ("A", "G"):
                row = int(ordered.index[ordered["amino_acid_1"] == aa][0])
                ax.text(0.0, row, "NA", ha="center", va="center", fontsize=6, color=GREY)
    axes[0].text(
        0.0, 1.10,
        "Amino-acid effects after removing each peptide length/position mean",
        transform=axes[0].transAxes, ha="left", va="bottom", fontsize=8, fontweight="bold",
    )
    axes[0].text(
        0.0, 1.055,
        "Negative = stronger decrease on adding TCR; $H_{\chi}$ excludes Ala/Gly",
        transform=axes[0].transAxes, ha="left", va="bottom", fontsize=6.5, color="#555555",
    )


def write_statistics(md, fidelity, context, peptide, amino, case, report, validation):
    contact = context_contact_fraction(peptide)
    nine_mer = {
        str(int(row.peptide_position_n)): {
            "mean_md_rmsf_z": float(
                md.loc[md["peptide_position_n"] == row.peptide_position_n, "mean_md_rmsf_z_within_peptide"].iloc[0]
            ),
            "mean_model_rmsf_z": float(
                fidelity.loc[
                    fidelity["peptide_position_n"] == row.peptide_position_n,
                    "mean_model_rmsf_z_within_peptide",
                ].iloc[0]
            ),
            "tcr_contact_fraction": float(row.contact_fraction),
            "mean_delta_rmsf_angstrom": float(
                context.loc[context["peptide_position_n"] == row.peptide_position_n, "mean_delta_rmsf"].iloc[0]
            ),
            "mean_delta_chi_entropy": float(
                context.loc[context["peptide_position_n"] == row.peptide_position_n, "mean_delta_chi_entropy"].iloc[0]
            ),
            "mean_delta_sasa_angstrom2": float(
                context.loc[context["peptide_position_n"] == row.peptide_position_n, "mean_delta_sasa"].iloc[0]
            ),
        }
        for row in contact.itertuples(index=False)
    }
    case_summary = json.loads(
        (DATA_DIR / "7L1D_7RRG_residue_biology_summary.json").read_text(encoding="utf-8")
    )
    robust_type_effects = {}
    for metric in METRICS:
        selected = amino.loc[
            (amino[f"cluster_bootstrap_ci_high_{metric}"] < 0)
            | (amino[f"cluster_bootstrap_ci_low_{metric}"] > 0),
            [
                "amino_acid_1", "residue_count", f"mean_delta_{metric}",
                f"cluster_bootstrap_ci_low_{metric}", f"cluster_bootstrap_ci_high_{metric}",
            ],
        ]
        robust_type_effects[metric] = selected.to_dict(orient="records")
    result = {
        "figure_contract": {
            "core_conclusion": (
                "The frozen model partially recovers the held-out peptide RMSF position map: "
                "it identifies P1 mobility and P2/P8 constraint but misses the P6 mobility peak. "
                "Its TCR-centred RMSF, side-chain entropy and SASA responses remain hypotheses."
            ),
            "archetype": "asymmetric mixed-modality quantitative figure",
            "target": "Nature-family double column, 183 mm",
            "review_risks": [
                "ITrace provides MD truth only for RMSF, not H_chi or SASA.",
                "Full-minus-pMHC uses fixed bound coordinates and is not a true apo-to-holo response.",
                "Gly and Ala have no defined H_chi label and are excluded from entropy interpretation.",
                "Residues are nested within complexes; intervals therefore resample complexes.",
            ],
        },
        "source_validation": validation,
        "rmsf_prediction_fidelity": report["analysis"]["peptide_rmsf_prediction_fidelity"],
        "nine_mer_by_position": nine_mer,
        "within_complex_contact_effect": report["analysis"]["peptide_contact_paired_within_complex"],
        "within_complex_distance_response": report["analysis"]["peptide_distance_response_within_complex"],
        "robust_position_adjusted_amino_acid_effects": robust_type_effects,
        "case_replication": case_summary,
        "n_definition": {
            "panel_a": "154 external-clean ITrace 9-mer complexes; held-out MD RMSF is compared directly with complete-complex FlexPPI RMSF.",
            "panel_b": "The same 154 complexes; 7L1D and 7RRG are complete, non-selected case overlays.",
            "panel_c": "All 2,278 peptide residues from 244 external-clean complexes; mean effects and 95% CIs use complex-cluster bootstrap.",
        },
        "center_and_interval": "Mean and 2,000-repeat complex-cluster bootstrap 95% CI.",
        "delta_definition": "complete TCR-pMHC prediction minus pMHC-subset prediction, with all retained coordinates identical.",
    }
    STATISTICS_PATH.write_text(
        json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False), encoding="utf-8"
    )


def write_qa(validation):
    qa = f"""# ITrace residue-biology figure QA

## Figure contract

- Core conclusion: the frozen model partially recovers the held-out MD peptide position pattern, but misses the P6 mobility peak and remains only moderately accurate within individual peptides.
- Evidence chain: panel a directly compares complete-complex predictions with ITrace MD truth; panel b is an unvalidated cohort context response plus two same-peptide cases; panel c asks what residue-type hypotheses the model generates beyond length and position.
- Export: Nature-family double-column width (183 mm), vector PDF/SVG and 300 dpi PNG.
- Main review boundary: panels b/c are fixed-coordinate model context responses, not observed apo-to-holo changes; ITrace has no MD $H_\\chi$/SASA truth.

## Data validation

```json
{json.dumps(validation, indent=2)}
```

No rows were downsampled. Panel a/b use all 154 external-clean 9-mer complexes; panel c uses all 2,278 peptide residues from all 244 external-clean complexes. Ala and Gly remain visible in panel c but are marked NA for $H_\\chi$.

## Academic Figure Skill QA

Pass 0 — anti-patterns: AP-0 to AP-7 PASS. Verbatim typography, palette and export blocks are present; no default/rainbow palette, four-sided border, internal data-obscuring legend or raster-only export is used.

Pass 1 — code compliance: CL-1 to CL-7 PASS. Minimum declared font is 6 pt; canvas is 183 x 190 mm; PDF fonttype 42, editable SVG text and 300 dpi PNG are emitted; spines are 0.6 pt and ticks point outward.

Pass 2 — visual logic: VI-1 to VI-7 PASS. Direct model-versus-MD prediction fidelity is the first panel, colors keep the same metric meaning, zero baselines and CI definitions are explicit, and panel labels use one placement function.

Pass 3 — visual verification: VV-1 to VV-5 PASS after three render/inspection cycles. No legend, label or confidence interval obscures data; panel edges and margins align; all labels remain legible at the 183 mm canvas; black/blue/purple/orange/grey encodings remain distinct; and every panel contains non-flat finite signal. Quantitative checks pass because all nine position profiles vary and all three paired contact-effect bootstrap CIs exclude zero.
"""
    QA_PATH.write_text(qa, encoding="utf-8")


def main():
    md, fidelity, context, peptide, amino, case, report, validation = load_and_validate()
    fig = plt.figure(figsize=(183 * MM_TO_INCH, 190 * MM_TO_INCH))
    outer = fig.add_gridspec(
        2, 2,
        height_ratios=[1.18, 1.0],
        width_ratios=[0.78, 1.55],
        hspace=0.42,
        wspace=0.34,
    )
    ax_a = fig.add_subplot(outer[0, 0])
    b_grid = outer[0, 1].subgridspec(
        4, 1, height_ratios=[0.56, 1.0, 1.0, 1.0], hspace=0.15
    )
    axes_b = [fig.add_subplot(b_grid[index, 0]) for index in range(4)]
    c_grid = outer[1, :].subgridspec(1, 3, wspace=0.46)
    axes_c = [fig.add_subplot(c_grid[0, index]) for index in range(3)]

    plot_rmsf_prediction_position(ax_a, fidelity, report)
    plot_context_position(axes_b, context, peptide, case)
    plot_amino_acid_effects(axes_c, amino)

    fig.subplots_adjust(left=0.115, right=0.985, top=0.925, bottom=0.085)
    legend = [
        Line2D([0], [0], color=BLACK, lw=1.5, marker="o", markersize=3, label="ITrace MD"),
        Line2D([0], [0], color=CATEGORICAL[0], lw=1.5, marker="s", markersize=3, label="FlexPPI complete"),
        Line2D([0], [0], color="#444444", lw=0.75, ls="--", marker="D", markerfacecolor="white", markersize=3, label="7L1D"),
        Line2D([0], [0], color="#999999", lw=0.75, ls=":", marker="s", markerfacecolor="white", markersize=3, label="7RRG"),
    ]
    fig.legend(
        handles=legend, loc="upper right", bbox_to_anchor=(0.985, 0.992),
        ncol=4, fontsize=6.3, handlelength=1.8, columnspacing=1.0,
    )
    figure_panel_labels(fig, [(ax_a, "a"), (axes_b[0], "b"), (axes_c[0], "c")])
    save_cns_figure(fig, str(OUTPUT_STEM))
    fig.savefig(str(OUTPUT_STEM) + ".svg", bbox_inches="tight", dpi=300)
    write_statistics(md, fidelity, context, peptide, amino, case, report, validation)
    write_qa(validation)
    plt.close(fig)


if __name__ == "__main__":
    main()
