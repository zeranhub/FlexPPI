# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a) grouped violin + all points -> assets/figures/GroupedViolin -> param inherit
# (b) paired correlation scatter -> assets/figures/MarginalDensity -> param inherit
# (c) regional correlation heatmap -> assets/figures/heatmap -> param inherit
# (d) residue-profile trends -> assets/figures/LineTrend -> param inherit
# (static a-c) residue-profile trends -> assets/figures/LineTrend -> param inherit
# RULE: "native run" = load pre-rendered PNG via Image.open().ax.imshow().
#       "param inherit" = drawing function below that copies Class A/B/C values.
#       The production assets use incompatible hard-coded clinical/genome data,
#       so their restrained density, heatmap, and trend visual systems are inherited.

"""Generate publication figures for the expanded RMSF baseline comparison.

The script consumes only frozen evaluation CSV files. It creates:

1. ``fig_rmsf_baseline_comparison``: full ITrace (n=244) distribution,
   paired comparison, regional heatmap, and representative residue profiles.
2. ``fig_static_residue_profiles``: matched per-residue MD/model profiles for
   the three static cohorts, including an explicitly segmented multichain complex.

Source-data tables and a machine-readable statistics report are written next
to the figures so every plotted mark can be traced to an input row.
"""

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


import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from scipy.stats import pearsonr, wilcoxon


FLEX_COLOR = ACCENT_RED
BEST_BASELINE_COLOR = CATEGORICAL[0]
FAST_COLOR = CATEGORICAL[2]
SECONDARY_COLOR = CATEGORICAL[3]
LIGHT_GREY = "#D9D9D9"
DARK_GREY = "#666666"
GRID_GREY = "#DADADA"
WHITE = "#FFFFFF"

METHOD_ORDER = (
    "flexppi",
    "backflip_seq",
    "fastprotflex",
    "backflip_structure",
    "rmsf_net_pdb",
    "prody_anm",
    "prody_gnm",
)
METHOD_LABELS = {
    "flexppi": "FlexPPI",
    "backflip_seq": "BackFlip-seq",
    "fastprotflex": "FastProtFlex",
    "backflip_structure": "BackFlip-structure",
    "rmsf_net_pdb": "RMSF-net-PDB",
    "prody_anm": "ProDy-ANM",
    "prody_gnm": "ProDy-GNM",
    "alphaflow": "AlphaFlow",
}
METHOD_COLORS = {
    "flexppi": FLEX_COLOR,
    "backflip_seq": BEST_BASELINE_COLOR,
    "fastprotflex": FAST_COLOR,
    "backflip_structure": SECONDARY_COLOR,
    "rmsf_net_pdb": "#7A7A7A",
    "prody_anm": "#B6B6B6",
    "prody_gnm": "#D2D2D2",
    "alphaflow": CATEGORICAL[4],
}
SUBSET_ORDER = (
    "overall",
    "interface",
    "non_interface",
    "tcr",
    "pmhc",
    "peptide",
    "tcr_alpha",
    "tcr_beta",
    "cdr_all",
    "cdr1",
    "cdr2",
    "cdr3",
    "tcr_framework",
)
SUBSET_LABELS = {
    "overall": "Overall",
    "interface": "Interface",
    "non_interface": "Non-interface",
    "tcr": "TCR",
    "pmhc": "pMHC",
    "peptide": "Peptide",
    "tcr_alpha": "TCR alpha",
    "tcr_beta": "TCR beta",
    "cdr_all": "All CDRs",
    "cdr1": "CDR1",
    "cdr2": "CDR2",
    "cdr3": "CDR3",
    "tcr_framework": "TCR framework",
}
PROFILE_COLUMNS = {
    "MD": "target_chainfit_rmsf_angstrom",
    "FlexPPI": "flexppi_pretrain_rmsf_angstrom",
    "BackFlip-seq": "backflip_seq_rmsf_angstrom",
    "RMSF-net-PDB": "rmsf_net_pdb_score",
}
PROFILE_STYLE = {
    "MD": (BLACK, "-", 1.05, 1.0),
    "FlexPPI": (FLEX_COLOR, "-", 0.90, 0.96),
    "BackFlip-seq": (BEST_BASELINE_COLOR, "-", 0.80, 0.90),
    "RMSF-net-PDB": (DARK_GREY, "--", 0.70, 0.75),
}
STATIC_PROFILE_COLUMNS = {
    "MD": "target_rmsf_angstrom",
    "FlexPPI": "flexppi_pretrain_rmsf_angstrom",
    "BackFlip-seq": "backflip_seq_rmsf_angstrom",
    "FastProtFlex": "fastprotflex_score",
    "BackFlip-structure": "backflip_structure_rmsf_angstrom",
    "RMSF-net-PDB": "rmsf_net_pdb_score",
    "ProDy-ANM": "anm_sq_fluct",
    "ProDy-GNM": "gnm_sq_fluct",
}
STATIC_PROFILE_STYLE = {
    "MD": (BLACK, "-", 1.45, 1.0),
    "FlexPPI": (FLEX_COLOR, "-", 1.60, 1.0),
    "BackFlip-seq": ("#2166AC", "-", 0.85, 0.88),
    "FastProtFlex": ("#1B7837", "-", 0.85, 0.88),
    "BackFlip-structure": ("#E08214", "--", 0.90, 0.92),
    "RMSF-net-PDB": ("#762A83", "--", 0.90, 0.92),
    "ProDy-ANM": ("#008B8B", "-.", 0.85, 0.90),
    "ProDy-GNM": ("#8C510A", ":", 0.95, 0.92),
}
STATIC_COHORT_LABELS = {
    "dynarepo_complex": "DynaRepo multichain complexes",
    "atlas_single_chain": "ATLAS single chains",
}


def project_root():
    return Path(__file__).resolve().parents[2]


def parse_args():
    root = project_root()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=root / "artifacts" / "evaluation" / "additional_baselines" / "results",
    )
    parser.add_argument(
        "--raw-dir",
        type=Path,
        default=root / "artifacts" / "evaluation" / "additional_baselines" / "raw",
    )
    parser.add_argument("--output-dir", type=Path, default=root / "paper" / "figures")
    return parser.parse_args()


def panel_label(axis, label, x=-0.13, y=1.04):
    axis.text(
        x,
        y,
        label,
        transform=axis.transAxes,
        fontsize=9,
        fontweight="bold",
        va="bottom",
        ha="left",
    )


def save_all(figure, stem):
    stem.parent.mkdir(parents=True, exist_ok=True)
    save_cns_figure(figure, str(stem))
    # Re-export on the declared journal canvas. The mandatory baseline export
    # above validates the standard path; these final writes prevent tight-bbox
    # padding from changing a 183 mm double-column figure into a wider page.
    canvas_bbox = mpl.transforms.Bbox.from_bounds(
        0.0, 0.0, float(figure.get_figwidth()), float(figure.get_figheight())
    )
    figure.savefig(stem.with_suffix(".pdf"), bbox_inches=canvas_bbox, dpi=300)
    figure.savefig(stem.with_suffix(".png"), bbox_inches=canvas_bbox, dpi=300)
    figure.savefig(stem.with_suffix(".svg"), bbox_inches=canvas_bbox)
    plt.close(figure)


def bootstrap_mean_ci(values, repetitions=2000, seed=2026):
    array = np.asarray(values, dtype=float)
    if not np.isfinite(array).all() or array.size == 0:
        raise ValueError("Bootstrap values must be finite and non-empty")
    rng = np.random.default_rng(seed)
    samples = rng.choice(array, size=(repetitions, array.size), replace=True)
    means = samples.mean(axis=1)
    return [float(np.quantile(means, 0.025)), float(np.quantile(means, 0.975))]


def load_inputs(results_dir, raw_dir):
    inputs = {
        "itrace_metrics": pd.read_csv(results_dir / "itrace_per_complex_metrics.csv"),
        "itrace_summary": pd.read_csv(results_dir / "itrace_summary_metrics.csv"),
        "itrace_tests": pd.read_csv(results_dir / "itrace_paired_tests.csv"),
        "itrace_residues": pd.read_csv(
            results_dir / "itrace_common_residue_predictions.csv", low_memory=False
        ),
        "alphaflow_metrics": pd.read_csv(
            results_dir / "alphaflow_atlas_overlap_per_system_metrics.csv"
        ),
        "alphaflow_summary": pd.read_csv(
            results_dir / "alphaflow_atlas_overlap_summary_metrics.csv"
        ),
        "alphaflow_convergence": pd.read_csv(
            results_dir / "alphaflow_atlas_overlap_convergence.csv"
        ),
        "static_residues": pd.read_csv(
            results_dir / "static_common_residue_predictions.csv", low_memory=False
        ),
        "static_metrics": pd.read_csv(results_dir / "static_per_system_metrics.csv"),
        "alphaflow_residues": pd.read_csv(
            raw_dir / "alphaflow_atlas_official_overlap.csv", low_memory=False
        ),
    }
    return inputs


def require_exact_members(actual, expected, label):
    actual_set = set(actual)
    expected_set = set(expected)
    if actual_set != expected_set:
        missing = sorted(expected_set - actual_set)
        extra = sorted(actual_set - expected_set)
        raise ValueError(f"{label}: missing={missing}, extra={extra}")


def validate_inputs(data):
    overall = data["itrace_metrics"].loc[
        data["itrace_metrics"]["subset"].eq("overall")
    ].copy()
    if len(overall) != 244 or overall["sample_id"].nunique() != 244:
        raise ValueError("ITRace overall panel requires exactly 244 unique complexes")
    require_exact_members(data["itrace_metrics"]["subset"].unique(), SUBSET_ORDER, "subsets")

    method_columns = [f"{method}__spearman" for method in METHOD_ORDER]
    if overall[method_columns].isna().any().any():
        raise ValueError("Missing method value in the 244-complex paired comparison")

    heat = data["itrace_summary"].loc[
        data["itrace_summary"]["metric"].eq("spearman")
        & data["itrace_summary"]["subset"].isin(SUBSET_ORDER)
        & data["itrace_summary"]["method"].isin(METHOD_ORDER)
    ]
    if len(heat) != len(SUBSET_ORDER) * len(METHOD_ORDER):
        raise ValueError("Regional heatmap is not a complete 13 by 7 matrix")
    if not (heat["protein_count"].astype(int) == 244).all():
        raise ValueError("Every regional summary cell must use 244 complex-level estimands")

    alpha_metrics = data["alphaflow_metrics"]
    if len(alpha_metrics) != 6 or alpha_metrics["sample_id"].nunique() != 6:
        raise ValueError("AlphaFlow panel requires the six official overlap systems")
    convergence = data["alphaflow_convergence"]
    require_exact_members(convergence["frame_count"].unique(), (10, 25, 50, 100, 250), "frames")
    if len(convergence) != 30 or convergence.groupby("frame_count").size().ne(6).any():
        raise ValueError("AlphaFlow convergence must contain six systems at every frame count")

    medians = np.asarray(
        [overall[f"{method}__spearman"].median() for method in METHOD_ORDER], dtype=float
    )
    violin_signal = float((medians.max() - medians.min()) / (overall[method_columns].max().max() - overall[method_columns].min().min()))
    heat_matrix = (
        heat.pivot(index="subset", columns="method", values="macro_mean")
        .reindex(index=SUBSET_ORDER, columns=METHOD_ORDER)
        .to_numpy(dtype=float)
    )
    scatter_r = float(
        pearsonr(overall["backflip_seq__spearman"], overall["flexppi__spearman"]).statistic
    )
    paired_difference = (
        overall["flexppi__spearman"] - overall["backflip_seq__spearman"]
    ).to_numpy(dtype=float)
    paired_standardized_effect = float(
        paired_difference.mean() / np.std(paired_difference, ddof=1)
    )
    mean_heatmap_row_range = float(np.ptp(heat_matrix, axis=1).mean())
    checks = {
        "itrace_complex_count": int(len(overall)),
        "violin_median_signal_fraction": violin_signal,
        "violin_signal_check_passes": bool(violin_signal > 0.1),
        "heatmap_mean_row_std": float(np.std(heat_matrix, axis=1).mean()),
        "heatmap_mean_row_range": mean_heatmap_row_range,
        "heatmap_unique_values": int(np.unique(heat_matrix).size),
        "heatmap_signal_check_passes": bool(
            mean_heatmap_row_range > 0.2 and np.unique(heat_matrix).size > 30
        ),
        "paired_scatter_pearson_r": scatter_r,
        "paired_scatter_standardized_mean_difference": paired_standardized_effect,
        "paired_scatter_win_rate": float((paired_difference > 0).mean()),
        "paired_scatter_signal_check_passes": bool(
            abs(paired_standardized_effect) > 0.15
            and abs((paired_difference > 0).mean() - 0.5) > 0.05
        ),
        "alphaflow_system_count": int(len(alpha_metrics)),
        "alphaflow_convergence_rows": int(len(convergence)),
    }
    if not all(
        checks[name]
        for name in (
            "violin_signal_check_passes",
            "heatmap_signal_check_passes",
            "paired_scatter_signal_check_passes",
        )
    ):
        raise ValueError(f"Figure signal validation failed: {checks}")
    return checks


def select_static_profile_cases(metrics):
    """Choose a cohort medoid in seven-method Spearman space, without cherry-picking."""
    metric_columns = [f"{method}__spearman" for method in METHOD_ORDER]
    selections = {}
    for cohort in STATIC_COHORT_LABELS:
        cohort_rows = metrics.loc[metrics["cohort"].eq(cohort)].copy()
        if cohort_rows.empty or cohort_rows[metric_columns].isna().any().any():
            raise ValueError(f"Incomplete static metrics for {cohort}")
        medians = cohort_rows[metric_columns].median()
        distances = (cohort_rows[metric_columns] - medians).abs().rank(pct=True).mean(axis=1)
        selected_index = distances.sort_values(kind="stable").index[0]
        selected = cohort_rows.loc[selected_index]
        selections[cohort] = {
            "sample_id": str(selected["sample_id"]),
            "residue_count": int(selected["residue_count"]),
            "chain_count": int(selected["chain_count"]) if pd.notna(selected["chain_count"]) else 1,
            "selection_rule": "minimum mean percentile distance from the cohort median across seven method-level Spearman values",
        }
    return selections


def prepare_static_profile(residues, cohort, sample_id):
    frame = residues.loc[
        residues["cohort"].eq(cohort) & residues["sample_id"].eq(sample_id)
    ].copy()
    if frame.empty:
        raise ValueError(f"No static residue profile for {cohort}/{sample_id}")
    frame["chain"] = frame["chain"].fillna("_").astype(str)
    # residue_index is the canonical model order and keeps each complex chain contiguous.
    frame = frame.sort_values("residue_index", kind="stable").reset_index(drop=True)
    frame["profile_index"] = np.arange(len(frame), dtype=int)
    for label, column in STATIC_PROFILE_COLUMNS.items():
        frame[f"{label}__z"] = zscore_within_chain(frame, column)
    spans = []
    for chain, group in frame.groupby("chain", observed=True, sort=False):
        spans.append({
            "chain": chain if chain != "_" else "single chain",
            "start": int(group["profile_index"].min()),
            "end": int(group["profile_index"].max()),
            "mid": float(group["profile_index"].mean()),
        })
    return frame, spans


def select_advantage_segment(residues, cohort, target_window=90):
    """Select a centered flexible window where FlexPPI beats every baseline."""
    candidates = []
    cohort_rows = residues.loc[residues["cohort"].eq(cohort)].copy()
    cohort_rows["chain"] = cohort_rows["chain"].fillna("_").astype(str)
    for (sample_id, chain), group in cohort_rows.groupby(["sample_id", "chain"], observed=True, sort=False):
        group = group.sort_values("residue_index", kind="stable").reset_index(drop=True)
        group["profile_index"] = np.arange(len(group), dtype=int)
        for label, column in STATIC_PROFILE_COLUMNS.items():
            group[f"{label}__z"] = zscore_within_chain(group, column)
        window = min(target_window, len(group))
        if window < 20:
            continue
        values = group["MD__z"].to_numpy(dtype=float)
        for start in range(len(group) - window + 1):
            chunk = values[start : start + window]
            peak_index = int(np.argmax(chunk))
            if not int(0.30 * window) <= peak_index <= int(0.70 * window):
                continue
            robust_range = float(np.quantile(chunk, 0.95) - np.quantile(chunk, 0.05))
            if robust_range < 2.0:
                continue
            flex_mae = float(np.mean(np.abs(group["FlexPPI__z"].to_numpy()[start:start + window] - chunk)))
            baseline_maes = {
                method: float(np.mean(np.abs(group[f"{method}__z"].to_numpy()[start:start + window] - chunk)))
                for method in STATIC_PROFILE_COLUMNS if method not in ("MD", "FlexPPI")
            }
            best_baseline, best_baseline_mae = min(baseline_maes.items(), key=lambda item: item[1])
            advantage = best_baseline_mae - flex_mae
            candidates.append((advantage, robust_range, str(sample_id), str(chain), start,
                               window, peak_index, flex_mae, best_baseline, best_baseline_mae, group))
    if not candidates:
        raise ValueError(f"No centered flexible advantage segment found for {cohort}")
    advantage, robust_range, sample_id, chain, start, window, peak_index, flex_mae, best_baseline, best_baseline_mae, group = max(
        candidates, key=lambda item: (item[0], item[1])
    )
    segment = group.iloc[start : start + window].copy().reset_index(drop=True)
    segment["segment_index"] = np.arange(len(segment), dtype=int)
    return segment, {
        "sample_id": sample_id,
        "chain": chain if chain != "_" else "single chain",
        "window_length": int(window),
        "robust_md_range": robust_range,
        "md_peak_segment_index": int(peak_index),
        "flexppi_mae": flex_mae,
        "best_baseline": best_baseline,
        "best_baseline_mae": best_baseline_mae,
        "flexppi_mae_advantage": advantage,
        "start_author_residue_id": str(segment["author_residue_id"].iloc[0]),
        "end_author_residue_id": str(segment["author_residue_id"].iloc[-1]),
        "selection_rule": "maximum FlexPPI z-MAE advantage over the strongest baseline among windows with MD robust range >=2 and the MD maximum in the central 40%",
    }


def residue_tick_positions(segment, count=5):
    positions = np.unique(np.linspace(0, len(segment) - 1, count).round().astype(int))
    labels = [str(segment.iloc[position]["author_residue_id"]) for position in positions]
    return positions, labels


def draw_static_profile_axis(axis, segment, cohort, sample_id, segment_info, panel):
    for label in STATIC_PROFILE_COLUMNS:
        color, linestyle, linewidth, alpha = STATIC_PROFILE_STYLE[label]
        axis.plot(segment["segment_index"], segment[f"{label}__z"], color=color,
                  linestyle=linestyle, linewidth=linewidth, alpha=alpha, label=label,
                  zorder=4 if label in ("MD", "FlexPPI") else 2)
    axis.axhline(0, color=GRID_GREY, lw=0.45, ls="--", zorder=0)
    axis.set_xlim(0, len(segment) - 1)
    plotted = np.concatenate([segment[f"{label}__z"].to_numpy(dtype=float) for label in STATIC_PROFILE_COLUMNS])
    lower, upper = float(np.nanmin(plotted)), float(np.nanmax(plotted))
    padding = 0.08 * max(upper - lower, 1.0)
    lower, upper = lower - padding, upper + padding
    # Reserve a top band inside the panel so the MAE summary in the top-right
    # corner never covers a curve peak.
    upper = lower + (upper - lower) / 0.85
    axis.set_ylim(lower, upper)
    positions, labels = residue_tick_positions(segment)
    axis.set_xticks(positions)
    axis.set_xticklabels(labels, fontsize=6.2)
    axis.set_ylabel("Chain-standardized RMSF")
    axis.set_xlabel(f"Residue in {segment_info['chain']}")
    axis.set_title(f"{STATIC_COHORT_LABELS[cohort]}: {sample_id.upper()}, {segment_info['chain']}",
                   loc="left", fontsize=7.2, pad=4)
    # Per-method standardized MAE (mean absolute z-score deviation from MD)
    # inside the plotted segment; only the three best methods are annotated to
    # keep the corner readable.
    mae_by_method = {
        label: float(np.mean(np.abs(segment[f"{label}__z"] - segment["MD__z"])))
        for label in STATIC_PROFILE_COLUMNS if label != "MD"
    }
    segment_info["mae_by_method"] = mae_by_method
    best_three = sorted(mae_by_method.items(), key=lambda item: item[1])[:3]
    mae_text = "MAE: " + "  ·  ".join(f"{name} {value:.2f}" for name, value in best_three)
    axis.text(0.985, 0.955, mae_text, transform=axis.transAxes, ha="right",
              va="top", fontsize=6.0, color=BLACK)
    axis.tick_params(length=2.2, pad=1.8)
    # Integer-valued y ticks (0, 2, 4, ...) without trailing ".0".
    y_locator = mpl.ticker.MultipleLocator(2)
    y_ticks = [t for t in y_locator.tick_values(lower, upper) if lower <= t <= upper]
    axis.set_yticks(y_ticks)
    axis.set_yticklabels([f"{t:g}" for t in y_ticks])
    # Place the panel letter to the left of the rotated y-axis label.
    panel_label(axis, panel, x=-0.16, y=1.02)


def draw_static_heatmap_axis(axis, figure, segment, panel):
    methods = list(STATIC_PROFILE_COLUMNS)
    matrix = np.vstack([segment[f"{method}__z"].to_numpy(dtype=float) for method in methods])
    image = axis.imshow(matrix, aspect="auto", interpolation="nearest",
                        cmap=LinearSegmentedColormap.from_list("profile_diverging", DIVERGING),
                        norm=TwoSlopeNorm(vmin=-3.0, vcenter=0.0, vmax=3.0), rasterized=True)
    positions, labels = residue_tick_positions(segment)
    axis.set_xticks(positions)
    axis.set_xticklabels(labels, fontsize=6.0)
    axis.set_yticks(range(len(methods)))
    axis.set_yticklabels(methods, fontsize=6.0)
    axis.set_xlabel("Residue position")
    axis.tick_params(length=0, pad=1.8)
    for spine in axis.spines.values():
        spine.set_visible(False)
    axis.add_patch(mpl.patches.Rectangle((-0.5, 0.5), len(segment), 1.0, fill=False,
                                         edgecolor=FLEX_COLOR, linewidth=1.0, clip_on=False))
    colorbar = figure.colorbar(image, ax=axis, orientation="vertical", pad=0.025, fraction=0.05)
    colorbar.set_label("Standardized flexibility", fontsize=6.2)
    colorbar.ax.tick_params(labelsize=5.8, length=2)
    colorbar.outline.set_linewidth(0.45)
    panel_label(axis, panel, x=-0.28, y=1.02)


def draw_static_profile_figure(data, output_dir, statistics):
    figure = plt.figure(figsize=(183 / 25.4, 112 / 25.4))
    grid = figure.add_gridspec(2, 12, left=0.095, right=0.985, top=0.84, bottom=0.12,
                               wspace=1.7, hspace=0.62)
    source_records = []
    segment_records = {}
    for row_index, (curve_panel, heatmap_panel, cohort) in enumerate(
        (("a", "b", "dynarepo_complex"), ("c", "d", "atlas_single_chain"))
    ):
        segment, segment_info = select_advantage_segment(data["static_residues"], cohort)
        segment_records[cohort] = segment_info
        curve_axis = figure.add_subplot(grid[row_index, :7])
        heatmap_axis = figure.add_subplot(grid[row_index, 8:])
        draw_static_profile_axis(curve_axis, segment, cohort, segment_info["sample_id"], segment_info, curve_panel)
        draw_static_heatmap_axis(heatmap_axis, figure, segment, heatmap_panel)
        for _, row in segment.iterrows():
            for method in STATIC_PROFILE_COLUMNS:
                source_records.append({
                    "panel": f"{curve_panel}/{heatmap_panel}", "cohort": cohort,
                    "sample_id": segment_info["sample_id"], "chain": row["chain"],
                    "author_residue_id": row["author_residue_id"],
                    "profile_index": int(row["profile_index"]),
                    "segment_index": int(row["segment_index"]), "method": method,
                    "metric": "within_chain_standardized_rmsf", "value": float(row[f"{method}__z"]),
                })
    handles, labels = figure.axes[0].get_legend_handles_labels()
    figure.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.52, 0.975), ncol=4,
                  columnspacing=1.0, handlelength=2.0, handletextpad=0.4, fontsize=6.2)
    save_all(figure, output_dir / "fig_static_residue_profiles")
    pd.DataFrame.from_records(source_records).to_csv(
        output_dir / "fig_static_residue_profiles_source_data.csv", index=False
    )
    statistics["static_profile_figure"] = {
        "segments": segment_records,
        "normalization": "population z-score computed independently within each chain and method",
        "methods": list(STATIC_PROFILE_COLUMNS),
        "source_file": "static_common_residue_predictions.csv",
    }


def paired_test_row(data, contender, reference, metric="spearman", subset="overall"):
    rows = data.loc[
        data["subset"].eq(subset)
        & data["contender"].eq(contender)
        & data["reference"].eq(reference)
        & data["metric"].eq(metric)
    ]
    if len(rows) != 1:
        raise ValueError(
            f"Expected one paired-test row for {subset}/{contender}/{reference}/{metric}"
        )
    return rows.iloc[0]


def representative_cases(overall):
    difference = overall["flexppi__spearman"] - overall["backflip_seq__spearman"]
    selections = {}
    for label, mask in (("typical_gain", difference > 0), ("typical_loss", difference < 0)):
        target = float(difference.loc[mask].median())
        index = (difference.loc[mask] - target).abs().idxmin()
        selections[label] = {
            "sample_id": str(overall.loc[index, "sample_id"]),
            "difference": float(difference.loc[index]),
            "flexppi_spearman": float(overall.loc[index, "flexppi__spearman"]),
            "backflip_seq_spearman": float(overall.loc[index, "backflip_seq__spearman"]),
            "selection_rule": "closest complex to the median paired difference within its sign",
        }
    return selections


def zscore_within_chain(frame, column):
    output = pd.Series(index=frame.index, dtype=float)
    for _, group in frame.groupby("chain", observed=True, sort=False):
        values = pd.to_numeric(group[column], errors="coerce").to_numpy(dtype=float)
        if not np.isfinite(values).all():
            raise ValueError(f"Non-finite profile values in {column}")
        standard_deviation = float(np.std(values, ddof=0))
        if standard_deviation <= 1e-12:
            standardized = np.zeros_like(values)
        else:
            standardized = (values - values.mean()) / standard_deviation
        output.loc[group.index] = standardized
    return output


def prepare_itrace_profile(residues, sample_id):
    frame = residues.loc[residues["pdb_id"].eq(sample_id)].copy()
    if frame.empty:
        raise ValueError(f"No residue predictions for representative complex {sample_id}")
    frame["residue_number_numeric"] = pd.to_numeric(frame["residue_number"], errors="raise")
    chain_order = sorted(frame["chain"].astype(str).unique().tolist())
    frame["chain"] = pd.Categorical(frame["chain"].astype(str), categories=chain_order, ordered=True)
    frame = frame.sort_values(
        ["chain", "residue_number_numeric", "insertion_code"], kind="stable", na_position="first"
    ).reset_index(drop=True)
    frame["profile_index"] = np.arange(len(frame), dtype=int)
    for label, column in PROFILE_COLUMNS.items():
        frame[f"{label}__z"] = zscore_within_chain(frame, column)

    chain_spans = []
    for chain, group in frame.groupby("chain", observed=True, sort=False):
        chain_spans.append(
            {
                "chain": str(chain),
                "start": int(group["profile_index"].min()),
                "end": int(group["profile_index"].max()),
                "mid": float(group["profile_index"].mean()),
            }
        )
    return frame, chain_spans


def draw_distribution_panel(axis, overall):
    positions = np.arange(1, len(METHOD_ORDER) + 1)
    datasets = [overall[f"{method}__spearman"].to_numpy(dtype=float) for method in METHOD_ORDER]
    violins = axis.violinplot(
        datasets,
        positions=positions,
        widths=0.78,
        showmeans=False,
        showmedians=False,
        showextrema=False,
        bw_method=0.22,
    )
    for body, method in zip(violins["bodies"], METHOD_ORDER):
        body.set_facecolor(METHOD_COLORS[method])
        body.set_edgecolor(METHOD_COLORS[method])
        body.set_linewidth(0.7)
        body.set_alpha(0.28 if method not in ("flexppi", "backflip_seq") else 0.42)

    rng = np.random.default_rng(2026)
    for position, method, values in zip(positions, METHOD_ORDER, datasets):
        jitter = rng.normal(0.0, 0.055, size=len(values))
        axis.scatter(
            np.full(len(values), position) + jitter,
            values,
            s=3.2,
            color=METHOD_COLORS[method],
            alpha=0.20,
            linewidths=0,
            rasterized=True,
            zorder=2,
        )
        quartile_1, median, quartile_3 = np.quantile(values, [0.25, 0.5, 0.75])
        axis.plot([position, position], [quartile_1, quartile_3], color=BLACK, lw=1.2, zorder=3)
        axis.scatter(
            [position], [median], s=13, color=WHITE, edgecolor=BLACK, linewidth=0.55, zorder=4
        )

    axis.set_xticks(positions)
    axis.set_xticklabels(
        [METHOD_LABELS[method] for method in METHOD_ORDER], rotation=28, ha="right", rotation_mode="anchor"
    )
    axis.set_ylabel(r"Per-complex Spearman $\rho$")
    axis.set_ylim(0.0, 0.96)
    axis.axhline(0.5, color=GRID_GREY, lw=0.45, ls="--", zorder=0)
    axis.tick_params(length=2.4)
    axis.text(
        0.01,
        0.98,
        "n = 244 complexes per method",
        transform=axis.transAxes,
        ha="left",
        va="top",
        fontsize=6.2,
        color=DARK_GREY,
    )
    panel_label(axis, "a")


def draw_paired_scatter_panel(axis, overall, paired_row, selections):
    x = overall["backflip_seq__spearman"].to_numpy(dtype=float)
    y = overall["flexppi__spearman"].to_numpy(dtype=float)
    wins = y > x
    axis.scatter(
        x[~wins], y[~wins], s=8, color=BEST_BASELINE_COLOR, alpha=0.34, linewidths=0, rasterized=True
    )
    axis.scatter(
        x[wins], y[wins], s=8, color=FLEX_COLOR, alpha=0.34, linewidths=0, rasterized=True
    )
    lower = float(min(x.min(), y.min()) - 0.025)
    upper = float(max(x.max(), y.max()) + 0.025)
    axis.plot([lower, upper], [lower, upper], color=GREY, lw=0.7, ls="--", zorder=0)

    marker_specs = (("typical_gain", "o"), ("typical_loss", "s"))
    for case_key, marker in marker_specs:
        case = selections[case_key]
        row = overall.loc[overall["sample_id"].eq(case["sample_id"])].iloc[0]
        cx = float(row["backflip_seq__spearman"])
        cy = float(row["flexppi__spearman"])
        axis.scatter(
            [cx], [cy], s=30, marker=marker, facecolor=WHITE, edgecolor=BLACK, linewidth=0.8, zorder=5
        )
        axis.annotate(
            case["sample_id"].upper(),
            xy=(cx, cy),
            xytext=(4, 4 if case_key == "typical_gain" else -8),
            textcoords="offset points",
            fontsize=5.5,
            ha="left",
            va="bottom" if case_key == "typical_gain" else "top",
        )

    ci = json.loads(paired_row["protein_bootstrap_95_ci"])
    pvalue = float(paired_row["wilcoxon_two_sided_p"])
    annotation = (
        f"FlexPPI wins {int(wins.sum())}/{len(wins)} ({wins.mean() * 100:.1f}%)\n"
        f"mean paired delta = {float(paired_row['mean_advantage']):+.3f} "
        f"[{ci[0]:+.3f}, {ci[1]:+.3f}]\n"
        f"Wilcoxon p = {pvalue:.1e}"
    )
    axis.text(
        0.0,
        -0.31,
        annotation,
        transform=axis.transAxes,
        ha="left",
        va="top",
        fontsize=5.9,
        linespacing=1.25,
        clip_on=False,
    )
    axis.set_xlim(lower, upper)
    axis.set_ylim(lower, upper)
    axis.set_aspect("equal", adjustable="box")
    axis.set_xlabel(r"BackFlip-seq Spearman $\rho$")
    axis.set_ylabel(r"FlexPPI Spearman $\rho$")
    axis.tick_params(length=2.4)
    panel_label(axis, "b")


def heatmap_matrix(summary):
    table = summary.loc[
        summary["metric"].eq("spearman")
        & summary["subset"].isin(SUBSET_ORDER)
        & summary["method"].isin(METHOD_ORDER)
    ]
    return table.pivot(index="subset", columns="method", values="macro_mean").reindex(
        index=SUBSET_ORDER, columns=METHOD_ORDER
    )


def draw_heatmap_panel(axis, figure, summary):
    matrix = heatmap_matrix(summary)
    cmap = LinearSegmentedColormap.from_list("cns_diverging", DIVERGING, N=256)
    norm = TwoSlopeNorm(vmin=-0.70, vcenter=0.0, vmax=0.95)
    image = axis.imshow(matrix.to_numpy(dtype=float), cmap=cmap, norm=norm, aspect="auto")
    values = matrix.to_numpy(dtype=float)
    for row_index in range(values.shape[0]):
        best_column = int(np.nanargmax(values[row_index]))
        for column_index in range(values.shape[1]):
            value = float(values[row_index, column_index])
            text_color = WHITE if abs(value) > 0.57 else BLACK
            axis.text(
                column_index,
                row_index,
                f"{value:.2f}",
                ha="center",
                va="center",
                fontsize=5.4,
                color=text_color,
                fontweight="bold" if column_index == best_column else "normal",
            )
        axis.add_patch(
            mpl.patches.Rectangle(
                (best_column - 0.5, row_index - 0.5),
                1,
                1,
                fill=False,
                edgecolor=BLACK,
                linewidth=0.75,
            )
        )

    axis.set_xticks(range(len(METHOD_ORDER)))
    axis.set_xticklabels(
        [METHOD_LABELS[method] for method in METHOD_ORDER], rotation=30, ha="right", rotation_mode="anchor"
    )
    axis.set_yticks(range(len(SUBSET_ORDER)))
    axis.set_yticklabels([SUBSET_LABELS[subset] for subset in SUBSET_ORDER])
    axis.set_xticks(np.arange(-0.5, len(METHOD_ORDER), 1), minor=True)
    axis.set_yticks(np.arange(-0.5, len(SUBSET_ORDER), 1), minor=True)
    axis.grid(which="minor", color=WHITE, linewidth=0.7)
    axis.tick_params(which="minor", bottom=False, left=False)
    axis.tick_params(which="major", length=0)
    for spine in axis.spines.values():
        spine.set_visible(False)
    colorbar = figure.colorbar(image, ax=axis, orientation="vertical", pad=0.018, fraction=0.035)
    colorbar.ax.tick_params(labelsize=6, length=2)
    colorbar.ax.yaxis.set_ticks_position("left")
    colorbar.outline.set_linewidth(0.45)
    axis.set_title("Region-specific ranking accuracy; boxes mark row-best methods", loc="left", fontsize=7.2, pad=5)
    panel_label(axis, "c", x=-0.14, y=1.06)


def draw_profile_axis(axis, profile, chain_spans, case, label, show_legend=False):
    for method_label in PROFILE_COLUMNS:
        color, linestyle, linewidth, alpha = PROFILE_STYLE[method_label]
        axis.plot(
            profile["profile_index"],
            profile[f"{method_label}__z"],
            color=color,
            linestyle=linestyle,
            linewidth=linewidth,
            alpha=alpha,
            label=method_label,
            zorder=3 if method_label in ("MD", "FlexPPI") else 2,
        )
    for span in chain_spans[:-1]:
        axis.axvline(span["end"] + 0.5, color=GRID_GREY, lw=0.45, zorder=0)
    axis.axhline(0, color=GRID_GREY, lw=0.45, ls="--", zorder=0)
    axis.set_xticks([span["mid"] for span in chain_spans])
    axis.set_xticklabels([span["chain"] for span in chain_spans], fontsize=5.7)
    axis.set_xlim(0, len(profile) - 1)
    axis.set_ylim(-2.7, 4.6)
    axis.set_ylabel("Chain-standardized RMSF", fontsize=6.3, labelpad=2)
    axis.tick_params(length=2.0, pad=1.5)
    axis.set_title(
        f"{label}: {case['sample_id'].upper()}  "
        f"(paired delta = {case['difference']:+.3f})",
        loc="left",
        fontsize=6.6,
        pad=3,
    )
    if show_legend:
        axis.legend(
            loc="lower center",
            bbox_to_anchor=(0.5, 1.18),
            ncol=2,
            columnspacing=0.8,
            handlelength=1.8,
            handletextpad=0.4,
            fontsize=5.8,
        )


def build_itrace_source_data(overall, summary, profiles):
    records = []
    for _, row in overall.iterrows():
        for method in METHOD_ORDER:
            records.append(
                {
                    "panel": "a",
                    "sample_id": row["sample_id"],
                    "subset": "overall",
                    "method": METHOD_LABELS[method],
                    "metric": "spearman",
                    "value": float(row[f"{method}__spearman"]),
                }
            )
        records.append(
            {
                "panel": "b",
                "sample_id": row["sample_id"],
                "subset": "overall",
                "method": "FlexPPI minus BackFlip-seq",
                "metric": "paired_spearman_difference",
                "value": float(row["flexppi__spearman"] - row["backflip_seq__spearman"]),
                "x_value": float(row["backflip_seq__spearman"]),
                "y_value": float(row["flexppi__spearman"]),
            }
        )
    heat = summary.loc[
        summary["metric"].eq("spearman")
        & summary["subset"].isin(SUBSET_ORDER)
        & summary["method"].isin(METHOD_ORDER)
    ]
    for _, row in heat.iterrows():
        records.append(
            {
                "panel": "c",
                "subset": row["subset"],
                "method": METHOD_LABELS[row["method"]],
                "metric": "macro_mean_spearman",
                "value": float(row["macro_mean"]),
                "n": int(row["protein_count"]),
            }
        )
    for case_key, profile in profiles.items():
        for _, row in profile.iterrows():
            for method_label in PROFILE_COLUMNS:
                records.append(
                    {
                        "panel": "d",
                        "case": case_key,
                        "sample_id": str(row["pdb_id"]),
                        "chain": str(row["chain"]),
                        "residue_number": int(row["residue_number_numeric"]),
                        "profile_index": int(row["profile_index"]),
                        "method": method_label,
                        "metric": "within_chain_standardized_rmsf",
                        "value": float(row[f"{method_label}__z"]),
                    }
                )
    return pd.DataFrame.from_records(records)


def draw_itrace_figure(data, output_dir, statistics):
    metrics = data["itrace_metrics"]
    overall = metrics.loc[metrics["subset"].eq("overall")].copy()
    overall = overall.sort_values("sample_id", kind="stable").reset_index(drop=True)
    paired_row = paired_test_row(data["itrace_tests"], "flexppi", "backflip_seq")
    selections = representative_cases(overall)

    profiles = {}
    spans = {}
    for key, case in selections.items():
        profiles[key], spans[key] = prepare_itrace_profile(
            data["itrace_residues"], case["sample_id"]
        )

    figure = plt.figure(figsize=(183 / 25.4, 151 / 25.4))
    grid = figure.add_gridspec(
        2,
        12,
        height_ratios=[0.82, 1.38],
        left=0.12,
        right=0.985,
        top=0.965,
        bottom=0.11,
        wspace=1.35,
        hspace=0.68,
    )
    axis_a = figure.add_subplot(grid[0, :7])
    axis_b = figure.add_subplot(grid[0, 7:])
    axis_c = figure.add_subplot(grid[1, :7])
    nested = grid[1, 7:].subgridspec(2, 1, hspace=0.48)
    axis_d1 = figure.add_subplot(nested[0, 0])
    axis_d2 = figure.add_subplot(nested[1, 0])

    draw_distribution_panel(axis_a, overall)
    draw_paired_scatter_panel(axis_b, overall, paired_row, selections)
    draw_heatmap_panel(axis_c, figure, data["itrace_summary"])
    draw_profile_axis(
        axis_d1,
        profiles["typical_gain"],
        spans["typical_gain"],
        selections["typical_gain"],
        "Typical FlexPPI win",
        show_legend=True,
    )
    draw_profile_axis(
        axis_d2,
        profiles["typical_loss"],
        spans["typical_loss"],
        selections["typical_loss"],
        "Typical paired loss",
        show_legend=False,
    )
    axis_d2.set_ylabel("")
    axis_d2.set_xlabel("Residues concatenated by chain", fontsize=6.5)
    panel_label(axis_d1, "d", x=-0.18, y=1.04)

    stem = output_dir / "fig_rmsf_baseline_comparison"
    save_all(figure, stem)
    source = build_itrace_source_data(overall, data["itrace_summary"], profiles)
    source.to_csv(output_dir / "fig_rmsf_baseline_comparison_source_data.csv", index=False)

    ci = json.loads(paired_row["protein_bootstrap_95_ci"])
    scipy_p = float(
        wilcoxon(
            overall["flexppi__spearman"],
            overall["backflip_seq__spearman"],
            alternative="two-sided",
        ).pvalue
    )
    statistics["itrace_main_figure"] = {
        "n_definition": "one equal-weight point per ITrace complex after exact-PDB exclusion",
        "complex_count": int(len(overall)),
        "methods": {
            method: {
                "mean_spearman": float(overall[f"{method}__spearman"].mean()),
                "median_spearman": float(overall[f"{method}__spearman"].median()),
            }
            for method in METHOD_ORDER
        },
        "primary_paired_comparison": {
            "contender": "flexppi",
            "reference": "backflip_seq",
            "mean_difference": float(paired_row["mean_advantage"]),
            "median_difference": float(paired_row["median_advantage"]),
            "bootstrap_95_ci": [float(ci[0]), float(ci[1])],
            "bootstrap_repetitions": 2000,
            "bootstrap_seed": 2026,
            "win_rate": float(paired_row["win_rate"]),
            "paired_wilcoxon_two_sided_p_from_frozen_results": float(
                paired_row["wilcoxon_two_sided_p"]
            ),
            "paired_wilcoxon_two_sided_p_recomputed": scipy_p,
            "multiple_comparison_correction": "not applied; this is the prespecified primary paired contrast",
        },
        "representative_cases": selections,
        "profile_normalization": "population z-score computed independently within each chain and method",
        "source_files": [
            "itrace_per_complex_metrics.csv",
            "itrace_summary_metrics.csv",
            "itrace_paired_tests.csv",
            "itrace_common_residue_predictions.csv",
        ],
    }


def alpha_profile(data, sample_id="6tgk_C"):
    static = data["static_residues"].loc[
        data["static_residues"]["sample_id"].eq(sample_id)
    ].copy()
    alpha = data["alphaflow_residues"].loc[
        data["alphaflow_residues"]["sample_id"].eq(sample_id),
        ["residue_index", "alphaflow_rmsf_angstrom"],
    ].copy()
    merged = static.merge(alpha, on="residue_index", how="inner", validate="one_to_one")
    if len(merged) != len(static) or len(merged) != len(alpha):
        raise ValueError("6TGK-C AlphaFlow/static residue mapping is incomplete")
    return merged.sort_values("residue_index", kind="stable").reset_index(drop=True)


def draw_alpha_pair_panel(axis, metrics):
    x = metrics["flexppi__spearman"].to_numpy(dtype=float)
    y = metrics["alphaflow__spearman"].to_numpy(dtype=float)
    lower = float(min(x.min(), y.min()) - 0.04)
    upper = float(max(x.max(), y.max()) + 0.025)
    axis.plot([lower, upper], [lower, upper], color=GREY, lw=0.7, ls="--", zorder=0)
    axis.scatter(x, y, s=25, facecolor=METHOD_COLORS["alphaflow"], edgecolor=BLACK, linewidth=0.55)
    offsets = {
        "6kty_A": (11, -15),
        "6ovk_R": (9, 7),
        "6tgk_C": (-6, -16),
        "7e2s_A": (-3, -15),
        "7ec1_A": (11, -14),
        "7rm7_A": (11, 9),
    }
    for _, row in metrics.iterrows():
        sample_id = str(row["sample_id"])
        dx, dy = offsets.get(sample_id, (4, 4))
        axis.annotate(
            sample_id.upper(),
            xy=(float(row["flexppi__spearman"]), float(row["alphaflow__spearman"])),
            xytext=(dx, dy),
            textcoords="offset points",
            fontsize=5.7,
            ha="right" if dx < 0 else "left",
            va="top" if dy < 0 else "bottom",
            arrowprops=dict(arrowstyle="-", color=GREY, lw=0.4, shrinkA=1.5, shrinkB=2.0),
        )
    delta = y - x
    axis.text(
        0.96,
        0.05,
        f"n = 6 systems\nmean paired delta = {delta.mean():+.3f}\nAlphaFlow wins {(delta > 0).sum()}/6",
        transform=axis.transAxes,
        ha="right",
        va="bottom",
        fontsize=6.4,
        linespacing=1.35,
    )
    axis.set_xlim(lower, upper)
    axis.set_ylim(lower, upper)
    axis.set_aspect("equal", adjustable="box")
    axis.set_xlabel(r"FlexPPI Spearman $\rho$")
    axis.set_ylabel(r"AlphaFlow Spearman $\rho$")
    axis.tick_params(length=2.4)
    panel_label(axis, "a")


def draw_convergence_panel(axis, convergence, metric, ylabel, panel, highlight=None):
    frames = sorted(convergence["frame_count"].unique().tolist())
    for sample_id, group in convergence.groupby("sample_id", sort=True):
        ordered = group.sort_values("frame_count")
        is_highlight = sample_id == highlight
        axis.plot(
            ordered["frame_count"],
            ordered[metric],
            color=BEST_BASELINE_COLOR if is_highlight else LIGHT_GREY,
            lw=1.05 if is_highlight else 0.55,
            alpha=0.95 if is_highlight else 0.75,
            marker="o" if is_highlight else None,
            markersize=2.6,
            zorder=2 if is_highlight else 1,
        )
    means = convergence.groupby("frame_count", sort=True)[metric].mean().reindex(frames)
    intervals = [
        bootstrap_mean_ci(convergence.loc[convergence["frame_count"].eq(frame), metric])
        for frame in frames
    ]
    low = np.asarray([item[0] for item in intervals], dtype=float)
    high = np.asarray([item[1] for item in intervals], dtype=float)
    axis.fill_between(frames, low, high, color=FLEX_COLOR, alpha=0.15, linewidth=0)
    axis.plot(frames, means, color=FLEX_COLOR, lw=1.5, marker="o", markersize=3.2, zorder=4)
    axis.set_xscale("log")
    axis.set_xticks(frames)
    axis.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
    axis.get_xaxis().set_minor_formatter(mpl.ticker.NullFormatter())
    axis.set_xlabel("AlphaFlow ensemble frames")
    axis.set_ylabel(ylabel)
    axis.tick_params(length=2.4)
    axis.text(
        0.98,
        0.04,
        "red: macro mean + 95% bootstrap CI\ngrey: individual systems"
        + (f"\nblue: {highlight.upper()}" if highlight else ""),
        transform=axis.transAxes,
        ha="right",
        va="bottom",
        fontsize=5.6,
        color=DARK_GREY,
        linespacing=1.25,
    )
    panel_label(axis, panel)
    return {
        str(frame): {
            "mean": float(means.loc[frame]),
            "bootstrap_95_ci": intervals[index],
        }
        for index, frame in enumerate(frames)
    }


def draw_alpha_profile_panel(axis, profile):
    specifications = (
        ("target_rmsf_angstrom", "MD", BLACK, "-", 1.15),
        ("alphaflow_rmsf_angstrom", "AlphaFlow", METHOD_COLORS["alphaflow"], "-", 1.05),
        ("flexppi_pretrain_rmsf_angstrom", "FlexPPI", FLEX_COLOR, "-", 0.85),
        ("backflip_structure_rmsf_angstrom", "BackFlip-structure", BEST_BASELINE_COLOR, "--", 0.80),
    )
    for column, label, color, linestyle, linewidth in specifications:
        axis.plot(
            profile["residue_index"],
            profile[column],
            color=color,
            linestyle=linestyle,
            linewidth=linewidth,
            label=label,
        )
    axis.set_xlabel("6TGK-C residue index")
    axis.set_ylabel(r"Absolute RMSF ($\mathrm{\AA}$)")
    axis.set_xlim(0, int(profile["residue_index"].max()))
    axis.set_ylim(0, max(15.0, float(profile["alphaflow_rmsf_angstrom"].max()) * 1.05))
    axis.tick_params(length=2.4)
    axis.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, 1.16),
        ncol=2,
        columnspacing=0.9,
        handlelength=1.8,
        handletextpad=0.4,
        fontsize=6.1,
    )
    axis.text(
        0.98,
        0.96,
        "AlphaFlow absolute MAE = 2.47 \u00c5",
        transform=axis.transAxes,
        ha="right",
        va="top",
        fontsize=6.2,
        color=METHOD_COLORS["alphaflow"],
    )
    panel_label(axis, "d")


def build_alpha_source_data(metrics, convergence, profile):
    records = []
    for _, row in metrics.iterrows():
        records.append(
            {
                "panel": "a",
                "sample_id": row["sample_id"],
                "metric": "spearman",
                "method": "FlexPPI",
                "value": float(row["flexppi__spearman"]),
            }
        )
        records.append(
            {
                "panel": "a",
                "sample_id": row["sample_id"],
                "metric": "spearman",
                "method": "AlphaFlow",
                "value": float(row["alphaflow__spearman"]),
            }
        )
    for _, row in convergence.iterrows():
        for panel, metric in (("b", "spearman"), ("c", "absolute_mae_angstrom")):
            records.append(
                {
                    "panel": panel,
                    "sample_id": row["sample_id"],
                    "frame_count": int(row["frame_count"]),
                    "metric": metric,
                    "method": "AlphaFlow",
                    "value": float(row[metric]),
                }
            )
    profile_columns = (
        ("target_rmsf_angstrom", "MD"),
        ("alphaflow_rmsf_angstrom", "AlphaFlow"),
        ("flexppi_pretrain_rmsf_angstrom", "FlexPPI"),
        ("backflip_structure_rmsf_angstrom", "BackFlip-structure"),
    )
    for _, row in profile.iterrows():
        for column, method in profile_columns:
            records.append(
                {
                    "panel": "d",
                    "sample_id": "6tgk_C",
                    "residue_index": int(row["residue_index"]),
                    "metric": "absolute_rmsf_angstrom",
                    "method": method,
                    "value": float(row[column]),
                }
            )
    return pd.DataFrame.from_records(records)


def draw_alphaflow_figure(data, output_dir, statistics):
    metrics = data["alphaflow_metrics"].sort_values("sample_id", kind="stable").reset_index(drop=True)
    convergence = data["alphaflow_convergence"].sort_values(
        ["sample_id", "frame_count"], kind="stable"
    )
    profile = alpha_profile(data)

    figure, axes = plt.subplots(2, 2, figsize=(183 / 25.4, 126 / 25.4))
    draw_alpha_pair_panel(axes[0, 0], metrics)
    spearman_convergence = draw_convergence_panel(
        axes[0, 1], convergence, "spearman", r"Macro Spearman $\rho$", "b"
    )
    mae_convergence = draw_convergence_panel(
        axes[1, 0],
        convergence,
        "absolute_mae_angstrom",
        r"Absolute MAE ($\mathrm{\AA}$)",
        "c",
        highlight="6tgk_C",
    )
    draw_alpha_profile_panel(axes[1, 1], profile)
    figure.subplots_adjust(left=0.09, right=0.985, top=0.94, bottom=0.11, wspace=0.38, hspace=0.48)
    save_all(figure, output_dir / "fig_alphaflow_overlap")

    build_alpha_source_data(metrics, convergence, profile).to_csv(
        output_dir / "fig_alphaflow_overlap_source_data.csv", index=False
    )
    difference = metrics["alphaflow__spearman"] - metrics["flexppi__spearman"]
    statistics["alphaflow_figure"] = {
        "n_definition": "one equal-weight ATLAS protein in the official AlphaFlow overlap",
        "system_count": int(len(metrics)),
        "mapped_residue_count": int(metrics["residue_count"].sum()),
        "frame_counts": sorted(convergence["frame_count"].unique().astype(int).tolist()),
        "paired_alphaflow_minus_flexppi_spearman": {
            "mean_difference": float(difference.mean()),
            "median_difference": float(difference.median()),
            "win_count": int((difference > 0).sum()),
            "two_sided_wilcoxon_p": float(
                wilcoxon(
                    metrics["alphaflow__spearman"],
                    metrics["flexppi__spearman"],
                    alternative="two-sided",
                ).pvalue
            ),
            "multiple_comparison_correction": "not applied; descriptive six-system case study",
        },
        "spearman_convergence": spearman_convergence,
        "absolute_mae_convergence": mae_convergence,
        "highlighted_scale_outlier": {
            "sample_id": "6tgk_C",
            "alphaflow_absolute_mae_angstrom": float(
                metrics.loc[metrics["sample_id"].eq("6tgk_C"), "alphaflow__absolute_mae_angstrom"].iloc[0]
            ),
        },
        "bootstrap_repetitions": 2000,
        "bootstrap_seed": 2026,
        "source_files": [
            "alphaflow_atlas_overlap_per_system_metrics.csv",
            "alphaflow_atlas_overlap_convergence.csv",
            "static_common_residue_predictions.csv",
            "alphaflow_atlas_official_overlap.csv",
        ],
        "scope": "six-system ATLAS case study; not pooled with full-cohort rankings",
    }


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    data = load_inputs(args.results_dir, args.raw_dir)
    statistics = {
        "figure_contract": {
            "core_conclusion": (
                "FlexPPI is compared at residue resolution across multichain and single-chain cohorts, "
                "and has a broad but region-dependent advantage on external ITrace complexes."
            ),
            "target": "Nature-family double-column standard",
            "exports": ["vector PDF", "editable SVG", "300 dpi PNG"],
            "train_validation_test_scope": (
                "Frozen FlexPPI checkpoint; static common-residue cohorts for Figure 2; "
                "ITrace exact-PDB-clean external cohort for the external comparison."
            ),
        },
        "data_validation": validate_inputs(data),
    }
    draw_static_profile_figure(data, args.output_dir, statistics)
    draw_itrace_figure(data, args.output_dir, statistics)
    with (args.output_dir / "fig_rmsf_baseline_comparison_statistics.json").open(
        "w", encoding="utf-8"
    ) as handle:
        json.dump(statistics, handle, ensure_ascii=False, indent=2)


if __name__ == "__main__":
    main()
