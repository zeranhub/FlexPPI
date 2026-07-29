# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a) paired checkpoint dumbbell → PairedBoxScatter/plot_PairedBoxScatter1.R → param inherit (R unavailable; aggregate estimates are not sample-level boxplot data)
# (b) paired model comparison → PairedBoxScatter/plot_PairedBoxScatter1.R → param inherit (R unavailable; aggregate estimates with bootstrap intervals)
# (c) response-method comparison → BarComparison/plot_comparison_Trajectory.py → param inherit (source structure has two targets and bootstrap intervals rather than repeated bars)
# (d) chain-level distribution → Violin/plot_Violin.py → param inherit (one continuous distribution rather than a multi-group significance plot)
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
import warnings

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import spearmanr


PROJECT_ROOT = Path(__file__).resolve().parents[2]
FIGURE_DIR = Path(__file__).resolve().parent
OUTPUT_STEM = FIGURE_DIR / "fig_completed_experiments"

PRETRAIN_STATE = PROJECT_ROOT / "artifacts/evaluation/state_validation_pretrain.json"
FINETUNE_STATE = PROJECT_ROOT / "artifacts/evaluation/state_validation_finetune.json"
PAIR_REPORT = PROJECT_ROOT / "artifacts/evaluation/pair_validation.json"
ATLAS_COMPARISON = (
    PROJECT_ROOT
    / "artifacts/evaluation/fastprotflex_comparison/atlas_comparison_report.json"
)
CONTEXT_COMPARISON = (
    PROJECT_ROOT
    / "artifacts/evaluation/fastprotflex_comparison/chain_context/"
    "chain_context_comparison_report.json"
)
TCR3D_REPORT = PROJECT_ROOT / "artifacts/evaluation/tcr3d_full/tcr3d_report.json"
TCR3D_RESIDUES = (
    PROJECT_ROOT / "artifacts/evaluation/tcr3d_full/tcr3d_residue_predictions.csv"
)

BLUE = CATEGORICAL[0]
ORANGE = CATEGORICAL[3]
PURPLE = CATEGORICAL[4]
DARK_GREY = CATEGORICAL[5]
LIGHT_GREY = "#C7C7C7"
PALE_BLUE = "#DCEAF7"
PALE_ORANGE = "#F8E1C2"


def load_json(path):
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def ci_array(value):
    if isinstance(value, str):
        return np.asarray([float(part) for part in value.split()], dtype=float)
    return np.asarray(value, dtype=float)


def metric_estimate(metric_block, metric_name):
    block = metric_block[metric_name]["overall"]["macro_by_system"]
    return {
        "estimate": float(block["spearman"]),
        "ci": ci_array(block["bootstrap_95ci"]["spearman"]),
        "n": int(block["correlation_system_count"]),
    }


def method_spearman(method_block):
    block = method_block["spearman"]
    return {
        "estimate": float(block["macro_mean"]),
        "ci": ci_array(block["protein_bootstrap_95_ci"]),
        "n": int(block["protein_count"]),
    }


def macro_mae(method_block):
    overall = method_block["overall"]
    block = overall["macro_by_system"]
    return {
        "estimate": float(block["mae"]),
        "ci": ci_array(block["bootstrap_95ci"]["mae"]),
        "n": int(overall["system_count"]),
    }


def compute_chain_spearman(csv_path):
    required = ["sample_id", "predicted_rmsf", "bfactor_rmsf_proxy"]
    frame = pd.read_csv(csv_path, usecols=required)
    input_rows = int(len(frame))
    finite = frame.replace([np.inf, -np.inf], np.nan).dropna(
        subset=["predicted_rmsf", "bfactor_rmsf_proxy"]
    )
    finite_rows = int(len(finite))
    records = []
    excluded_constant_or_short = 0
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for sample_id, group in finite.groupby("sample_id", sort=True):
            if (
                len(group) < 3
                or group["predicted_rmsf"].nunique() < 2
                or group["bfactor_rmsf_proxy"].nunique() < 2
            ):
                excluded_constant_or_short += 1
                continue
            rho = spearmanr(
                group["predicted_rmsf"].to_numpy(dtype=float),
                group["bfactor_rmsf_proxy"].to_numpy(dtype=float),
            ).statistic
            if np.isfinite(rho):
                records.append(
                    {
                        "sample_id": sample_id,
                        "spearman": float(rho),
                        "residue_count": int(len(group)),
                    }
                )
            else:
                excluded_constant_or_short += 1
    result = pd.DataFrame.from_records(records)
    return result, {
        "input_residue_rows": input_rows,
        "finite_residue_rows": finite_rows,
        "dropped_nonfinite_residue_rows": input_rows - finite_rows,
        "correlatable_chains": int(len(result)),
        "excluded_constant_or_short_chains": int(excluded_constant_or_short),
    }


def asymmetric_error(estimate, ci):
    return np.asarray(
        [[estimate - float(ci[0])], [float(ci[1]) - estimate]], dtype=float
    )


def add_panel_label(ax, label):
    ax.text(
        -0.14,
        1.08,
        label,
        transform=ax.transAxes,
        fontsize=9,
        fontweight="bold",
        va="top",
        ha="left",
    )


def style_axis(ax):
    ax.spines["left"].set_linewidth(0.6)
    ax.spines["bottom"].set_linewidth(0.6)
    ax.tick_params(length=3, width=0.6)


def panel_a(ax, state_data):
    metrics = ["RMSF", r"Side-chain entropy ($H_{\chi}$)"]
    y_positions = np.asarray([1.0, 0.0])
    for y, metric in zip(y_positions, metrics):
        pre = state_data[metric]["Pretrain"]
        fine = state_data[metric]["Fine-tune"]
        ax.plot(
            [pre["estimate"], fine["estimate"]],
            [y, y],
            color=LIGHT_GREY,
            linewidth=1.2,
            zorder=1,
        )
        ax.errorbar(
            pre["estimate"],
            y,
            xerr=asymmetric_error(pre["estimate"], pre["ci"]),
            fmt="o",
            markersize=5.2,
            markerfacecolor=BLUE,
            markeredgecolor="white",
            markeredgewidth=0.5,
            ecolor=BLUE,
            elinewidth=0.9,
            capsize=2,
            zorder=3,
        )
        ax.errorbar(
            fine["estimate"],
            y,
            xerr=asymmetric_error(fine["estimate"], fine["ci"]),
            fmt="s",
            markersize=4.8,
            markerfacecolor=ORANGE,
            markeredgecolor="white",
            markeredgewidth=0.5,
            ecolor=ORANGE,
            elinewidth=0.9,
            capsize=2,
            zorder=3,
        )
        delta = fine["estimate"] - pre["estimate"]
        x_text = min(pre["estimate"], fine["estimate"]) - 0.012
        ax.text(
            x_text,
            y + 0.16,
            rf"$\Delta={delta:+.3f}$",
            color=ACCENT_RED if delta < -0.02 else DARK_GREY,
            fontsize=6.5,
            ha="left",
            va="center",
        )
    ax.set_yticks(y_positions, metrics)
    ax.set_xlim(0.70, 0.96)
    ax.set_ylim(-0.55, 1.52)
    ax.set_xlabel(r"Protein-equal Spearman $\rho$")
    ax.set_title("Absolute-state retention", loc="left", fontweight="bold")
    handles = [
        Line2D(
            [0], [0], marker="o", linestyle="none", color=BLUE,
            markerfacecolor=BLUE, markeredgecolor="white", markersize=5,
            label="Pretrain",
        ),
        Line2D(
            [0], [0], marker="s", linestyle="none", color=ORANGE,
            markerfacecolor=ORANGE, markeredgecolor="white", markersize=5,
            label="Fine-tune",
        ),
    ]
    ax.legend(
        handles=handles,
        loc="upper right",
        ncol=2,
        handletextpad=0.35,
        columnspacing=1.0,
        borderaxespad=0.35,
    )
    add_panel_label(ax, "a")
    style_axis(ax)


def panel_b(ax, comparison_data):
    datasets = [
        "ATLAS\nsingle-chain",
        "DynaRepo\nsingle-chain",
        "DynaRepo\ncomplexes",
    ]
    y_base = np.asarray([2.0, 1.0, 0.0])
    offset = 0.11
    for y, dataset in zip(y_base, datasets):
        fast = comparison_data[dataset]["FastProtFlex"]
        flex = comparison_data[dataset]["FlexPPI"]
        ax.plot(
            [fast["estimate"], flex["estimate"]],
            [y - offset, y + offset],
            color=LIGHT_GREY,
            linewidth=0.9,
            zorder=1,
        )
        ax.errorbar(
            fast["estimate"],
            y - offset,
            xerr=asymmetric_error(fast["estimate"], fast["ci"]),
            fmt="s",
            markersize=4.5,
            markerfacecolor=DARK_GREY,
            markeredgecolor="white",
            markeredgewidth=0.5,
            ecolor=DARK_GREY,
            elinewidth=0.8,
            capsize=2,
            zorder=3,
        )
        ax.errorbar(
            flex["estimate"],
            y + offset,
            xerr=asymmetric_error(flex["estimate"], flex["ci"]),
            fmt="o",
            markersize=5.0,
            markerfacecolor=BLUE,
            markeredgecolor="white",
            markeredgewidth=0.5,
            ecolor=BLUE,
            elinewidth=0.8,
            capsize=2,
            zorder=3,
        )
        ax.text(
            max(fast["ci"][1], flex["ci"][1]) + 0.006,
            y,
            rf"$+{flex['estimate'] - fast['estimate']:.3f}$",
            fontsize=6.2,
            color=BLUE,
            ha="left",
            va="center",
        )
    ax.set_yticks(y_base, datasets)
    ax.set_xlim(0.70, 0.895)
    ax.set_ylim(-0.55, 2.55)
    ax.set_xlabel(r"Protein/complex-equal Spearman $\rho$")
    ax.set_title("Static RMSF ranking", loc="left", fontweight="bold")
    handles = [
        Line2D(
            [0], [0], marker="o", linestyle="none", color=BLUE,
            markerfacecolor=BLUE, markeredgecolor="white", markersize=5,
            label="FlexPPI",
        ),
        Line2D(
            [0], [0], marker="s", linestyle="none", color=DARK_GREY,
            markerfacecolor=DARK_GREY, markeredgecolor="white", markersize=5,
            label="FastProtFlex",
        ),
    ]
    ax.legend(
        handles=handles,
        loc="upper left",
        ncol=2,
        handletextpad=0.35,
        columnspacing=1.0,
        borderaxespad=0.35,
    )
    add_panel_label(ax, "b")
    style_axis(ax)


def response_facet(ax, title, values, show_ylabels):
    methods = ["Response head", "State subtraction", "Zero change"]
    colors = [ORANGE, PURPLE, DARK_GREY]
    markers = ["o", "D", "s"]
    y_positions = np.asarray([2.0, 1.0, 0.0])
    baseline = values["Zero change"]["estimate"]
    ax.axvline(
        baseline,
        color=GREY,
        linewidth=0.7,
        linestyle=(0, (3, 2)),
        alpha=0.8,
        zorder=0,
    )
    for y, method, color, marker in zip(y_positions, methods, colors, markers):
        item = values[method]
        ax.errorbar(
            item["estimate"],
            y,
            xerr=asymmetric_error(item["estimate"], item["ci"]),
            fmt=marker,
            markersize=4.5,
            markerfacecolor=color,
            markeredgecolor="white",
            markeredgewidth=0.45,
            ecolor=color,
            elinewidth=0.85,
            capsize=2,
            zorder=3,
        )
    ax.set_yticks(y_positions)
    ax.set_yticklabels(methods if show_ylabels else [])
    ax.set_ylim(-0.55, 2.55)
    ax.set_title(title, loc="center", fontsize=7.5, fontweight="bold")
    ax.set_xlabel("MAE")
    style_axis(ax)


def panel_c(ax_left, ax_right, response_data):
    response_facet(
        ax_left,
        r"$\Delta$RMSF ($\AA$)",
        response_data["RMSF"],
        show_ylabels=True,
    )
    response_facet(
        ax_right,
        r"$\Delta H_{\chi}$",
        response_data[r"$H_{\chi}$"],
        show_ylabels=False,
    )
    ax_left.set_xlim(0.10, 0.39)
    ax_right.set_xlim(0.075, 0.185)
    ax_left.text(
        0.0,
        1.22,
        "Paired response error",
        transform=ax_left.transAxes,
        fontsize=8,
        fontweight="bold",
        ha="left",
        va="bottom",
    )
    add_panel_label(ax_left, "c")


def panel_d(ax, chain_values, report_summary):
    values = chain_values["spearman"].to_numpy(dtype=float)
    violin = ax.violinplot(
        values,
        positions=[0.0],
        vert=False,
        widths=0.72,
        showmeans=False,
        showmedians=False,
        showextrema=False,
        points=200,
        bw_method=0.22,
    )
    for body in violin["bodies"]:
        body.set_facecolor(PALE_ORANGE)
        body.set_edgecolor(ORANGE)
        body.set_linewidth(0.8)
        body.set_alpha(1.0)
    q1, median, q3 = np.quantile(values, [0.25, 0.5, 0.75])
    ax.plot([q1, q3], [0.0, 0.0], color=DARK_GREY, linewidth=5.0, solid_capstyle="butt")
    ax.plot([median, median], [-0.08, 0.08], color="white", linewidth=1.2)
    mean = float(report_summary["spearman"])
    ci = ci_array(report_summary["ci"])
    ax.errorbar(
        mean,
        0.0,
        xerr=asymmetric_error(mean, ci),
        fmt="o",
        markerfacecolor=ORANGE,
        markeredgecolor="white",
        markeredgewidth=0.5,
        markersize=5.0,
        ecolor=ORANGE,
        elinewidth=1.0,
        capsize=2,
        zorder=4,
    )
    ax.axvline(
        0.0,
        color=GREY,
        linewidth=0.7,
        linestyle=(0, (3, 2)),
        alpha=0.8,
        zorder=0,
    )
    ax.text(
        0.98,
        0.93,
        f"n = {len(values):,} chains\nmean = {mean:.3f}\n95% CI [{ci[0]:.3f}, {ci[1]:.3f}]",
        transform=ax.transAxes,
        fontsize=6.5,
        ha="right",
        va="top",
        color=BLACK,
    )
    ax.set_xlim(-1.0, 1.0)
    ax.set_ylim(-0.62, 0.62)
    ax.set_yticks([])
    ax.set_xlabel(r"Within-chain Spearman $\rho$")
    ax.set_title("TCR3D external B-factor proxy", loc="left", fontweight="bold")
    add_panel_label(ax, "d")
    style_axis(ax)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", left=False)


def build_data():
    pretrain = load_json(PRETRAIN_STATE)
    finetune = load_json(FINETUNE_STATE)
    pair_report = load_json(PAIR_REPORT)
    atlas_report = load_json(ATLAS_COMPARISON)
    context_report = load_json(CONTEXT_COMPARISON)
    tcr_report = load_json(TCR3D_REPORT)

    state_data = {
        "RMSF": {
            "Pretrain": metric_estimate(pretrain["metrics"], "rmsf"),
            "Fine-tune": metric_estimate(finetune["metrics"], "rmsf"),
        },
        r"Side-chain entropy ($H_{\chi}$)": {
            "Pretrain": metric_estimate(pretrain["metrics"], "chi_entropy"),
            "Fine-tune": metric_estimate(finetune["metrics"], "chi_entropy"),
        },
    }

    dynarepo_single = context_report["single_chain"]["dynarepo_single_protein_chain"][
        "methods"
    ]
    complexes = context_report["protein_complex"]["equal_complex_methods"]
    comparison_data = {
        "ATLAS\nsingle-chain": {
            "FastProtFlex": method_spearman(atlas_report["methods"]["fastprotflex"]),
            "FlexPPI": method_spearman(atlas_report["methods"]["flexppi_pretrain"]),
        },
        "DynaRepo\nsingle-chain": {
            "FastProtFlex": method_spearman(dynarepo_single["fastprotflex"]),
            "FlexPPI": method_spearman(dynarepo_single["flexppi_pretrain"]),
        },
        "DynaRepo\ncomplexes": {
            "FastProtFlex": method_spearman(complexes["fastprotflex_chain"]),
            "FlexPPI": method_spearman(complexes["flexppi_pretrain_complex"]),
        },
    }

    method_key = {
        "Response head": "response_head",
        "State subtraction": "state_subtraction",
        "Zero change": "zero_change",
    }
    response_data = {}
    for display_metric, key in [("RMSF", "rmsf"), (r"$H_{\chi}$", "chi_entropy")]:
        response_data[display_metric] = {}
        methods = pair_report["metrics"][key]["method_comparison"]
        for display_method, raw_method in method_key.items():
            response_data[display_metric][display_method] = macro_mae(methods[raw_method])

    chain_values, chain_audit = compute_chain_spearman(TCR3D_RESIDUES)
    tcr_macro = tcr_report["bfactor_rank_proxy"]["within_chain_zscore"][
        "macro_by_system"
    ]
    tcr_summary = {
        "spearman": float(tcr_macro["spearman"]),
        "ci": ci_array(tcr_macro["bootstrap_95ci"]["spearman"]),
        "n": int(tcr_macro["correlation_system_count"]),
    }
    calculated_mean = float(chain_values["spearman"].mean())
    if not np.isclose(calculated_mean, tcr_summary["spearman"], atol=1e-12):
        raise ValueError(
            "Recomputed TCR3D chain mean does not match the archived report: "
            f"{calculated_mean} vs {tcr_summary['spearman']}"
        )
    if len(chain_values) != tcr_summary["n"]:
        raise ValueError(
            "Recomputed TCR3D chain count does not match the archived report: "
            f"{len(chain_values)} vs {tcr_summary['n']}"
        )
    return (
        state_data,
        comparison_data,
        response_data,
        chain_values,
        chain_audit,
        tcr_summary,
    )


def validate_signal(state_data, comparison_data, response_data, chain_values):
    state_values = np.asarray(
        [
            item["estimate"]
            for metric in state_data.values()
            for item in metric.values()
        ]
    )
    comparison_deltas = np.asarray(
        [
            methods["FlexPPI"]["estimate"] - methods["FastProtFlex"]["estimate"]
            for methods in comparison_data.values()
        ]
    )
    response_ranges = np.asarray(
        [
            max(item["estimate"] for item in methods.values())
            - min(item["estimate"] for item in methods.values())
            for methods in response_data.values()
        ]
    )
    chain_std = float(chain_values["spearman"].std(ddof=1))
    checks = {
        "panel_a_value_range": float(np.ptp(state_values)),
        "panel_b_min_FlexPPI_minus_FastProtFlex": float(comparison_deltas.min()),
        "panel_c_min_method_range": float(response_ranges.min()),
        "panel_d_chain_spearman_sd": chain_std,
    }
    if checks["panel_a_value_range"] <= 0.05 * float(state_values.max()):
        raise ValueError("Panel (a) lacks sufficient visual signal.")
    if np.any(comparison_deltas <= 0):
        raise ValueError("Panel (b) does not have a consistent positive method difference.")
    if np.any(response_ranges <= 0.05 * np.asarray(
        [max(item["estimate"] for item in methods.values()) for methods in response_data.values()]
    )):
        raise ValueError("Panel (c) lacks sufficient between-method range.")
    if chain_std <= 0.05:
        raise ValueError("Panel (d) distribution is nearly degenerate.")
    return checks


def write_source_data(
    state_data,
    comparison_data,
    response_data,
    chain_values,
    chain_audit,
    tcr_summary,
    signal_checks,
):
    records = []
    for metric, methods in state_data.items():
        for method, item in methods.items():
            records.append(
                {
                    "panel": "a",
                    "sample_id": "",
                    "group": metric,
                    "method": method,
                    "value": item["estimate"],
                    "ci_low": item["ci"][0],
                    "ci_high": item["ci"][1],
                    "n": item["n"],
                    "value_definition": "protein-equal macro Spearman",
                }
            )
    for dataset, methods in comparison_data.items():
        for method, item in methods.items():
            records.append(
                {
                    "panel": "b",
                    "sample_id": "",
                    "group": dataset.replace("\n", " "),
                    "method": method,
                    "value": item["estimate"],
                    "ci_low": item["ci"][0],
                    "ci_high": item["ci"][1],
                    "n": item["n"],
                    "value_definition": "protein/complex-equal macro Spearman",
                }
            )
    for metric, methods in response_data.items():
        for method, item in methods.items():
            records.append(
                {
                    "panel": "c",
                    "sample_id": "",
                    "group": metric,
                    "method": method,
                    "value": item["estimate"],
                    "ci_low": item["ci"][0],
                    "ci_high": item["ci"][1],
                    "n": item["n"],
                    "value_definition": "system-equal macro MAE",
                }
            )
    for row in chain_values.itertuples(index=False):
        records.append(
            {
                "panel": "d",
                "sample_id": row.sample_id,
                "group": "TCR3D",
                "method": "FlexPPI finetune",
                "value": row.spearman,
                "ci_low": np.nan,
                "ci_high": np.nan,
                "n": row.residue_count,
                "value_definition": "within-chain Spearman",
            }
        )
    pd.DataFrame.from_records(records).to_csv(
        FIGURE_DIR / "fig_completed_experiments_source_data.csv", index=False
    )
    statistics = {
        "figure_contract": {
            "core_conclusion": (
                "Absolute-state pretraining learns residue-level flexibility and "
                "outperforms FastProtFlex on shared static tasks, whereas paired "
                "response remains weak and TCR3D provides only moderate proxy evidence."
            ),
            "archetype": "quantitative_grid",
            "target": "Nature-family standard, double column, 183 mm",
        },
        "input_files": [
            str(PRETRAIN_STATE.relative_to(PROJECT_ROOT)),
            str(FINETUNE_STATE.relative_to(PROJECT_ROOT)),
            str(PAIR_REPORT.relative_to(PROJECT_ROOT)),
            str(ATLAS_COMPARISON.relative_to(PROJECT_ROOT)),
            str(CONTEXT_COMPARISON.relative_to(PROJECT_ROOT)),
            str(TCR3D_REPORT.relative_to(PROJECT_ROOT)),
            str(TCR3D_RESIDUES.relative_to(PROJECT_ROOT)),
        ],
        "tcr3d_data_audit": chain_audit,
        "tcr3d_reported_summary": {
            "mean_spearman": tcr_summary["spearman"],
            "bootstrap_95_ci": tcr_summary["ci"].tolist(),
            "correlatable_chains": tcr_summary["n"],
        },
        "signal_checks": signal_checks,
        "statistics": {
            "panels_a_b_c": (
                "Point estimates are equal-system means; intervals are archived "
                "system-level nonparametric bootstrap 95% confidence intervals."
            ),
            "panel_d": (
                "The violin contains all correlatable chains. The orange point is "
                "the equal-chain mean and the whisker is the archived chain bootstrap "
                "95% confidence interval."
            ),
            "hypothesis_tests": "No p-values are asserted in this figure.",
            "multiple_comparison_correction": "Not applicable.",
        },
    }
    with (FIGURE_DIR / "fig_completed_experiments_statistics.json").open(
        "w", encoding="utf-8"
    ) as handle:
        json.dump(statistics, handle, ensure_ascii=False, indent=2)


def main():
    (
        state_data,
        comparison_data,
        response_data,
        chain_values,
        chain_audit,
        tcr_summary,
    ) = build_data()
    signal_checks = validate_signal(
        state_data, comparison_data, response_data, chain_values
    )

    mm_to_inch = 1.0 / 25.4
    fig = plt.figure(figsize=(183 * mm_to_inch, 122 * mm_to_inch))
    outer = fig.add_gridspec(
        2,
        2,
        left=0.09,
        right=0.985,
        bottom=0.11,
        top=0.91,
        wspace=0.44,
        hspace=0.53,
    )
    ax_a = fig.add_subplot(outer[0, 0])
    ax_b = fig.add_subplot(outer[0, 1])
    inner_c = outer[1, 0].subgridspec(1, 2, wspace=0.20)
    ax_c_left = fig.add_subplot(inner_c[0, 0])
    ax_c_right = fig.add_subplot(inner_c[0, 1])
    ax_d = fig.add_subplot(outer[1, 1])

    panel_a(ax_a, state_data)
    panel_b(ax_b, comparison_data)
    panel_c(ax_c_left, ax_c_right, response_data)
    panel_d(ax_d, chain_values, tcr_summary)

    save_cns_figure(fig, str(OUTPUT_STEM))
    plt.close(fig)
    write_source_data(
        state_data,
        comparison_data,
        response_data,
        chain_values,
        chain_audit,
        tcr_summary,
        signal_checks,
    )
    print(
        json.dumps(
            {
                "pdf": str(OUTPUT_STEM.with_suffix(".pdf")),
                "png": str(OUTPUT_STEM.with_suffix(".png")),
                "source_rows": 4 + 6 + 6 + len(chain_values),
                "tcr3d_chain_count": len(chain_values),
                "signal_checks": signal_checks,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
