# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (01a-f) All-residue density -> MarginalDensity/plot_VariableCorrelation.py -> param inherit
# (02a-d) System distributions -> Violin/plot_Violin.py -> param inherit
# (03a-d) Error stratification -> MarginalDensity/plot_VariableCorrelation.py -> param inherit
# (04a-d) Baseline and AA correction -> LineTrend/plot_sweep.py -> param inherit
# (05a-d) Uncertainty diagnostics -> LineTrend/plot_sweep.py -> param inherit
# User-requested white/boxed-header/muted reference styling overrides asset colors.
# Assets assume genome filters/log fits, independent significance tests or embedded
# optimization data; these are incompatible with the OOF statistical protocol.
# Preserve clean spines, group fills, marker/line emphasis; do not copy those tests.
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


import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LogNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.ticker import FixedLocator, FuncFormatter, MaxNLocator

ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = ROOT / "artifacts/evaluation_atlas_hchi_v1/analysis"
OUTPUT = ROOT / "output/flexppi_atlas_hchi_v1_extended_figures"
SOURCES = ["atlas", "dynarepo", "itrace"]
LABELS = dict(atlas="ATLAS", dynarepo="DynaRepo", itrace="ITrace")
COLORS = dict(atlas="#4E617F", dynarepo="#76A99B", itrace="#AAA56B")
MARKERS = dict(atlas="o", dynarepo="s", itrace="^")
LINES = dict(atlas="-", dynarepo="--", itrace=":")
NAMES = dict(rmsf="RMSF", chi_entropy="Hχ")
UNITS = dict(rmsf="Å", chi_entropy="nats")
ACCENT = "#D44D4C"
DENSITY_CMAP = LinearSegmentedColormap.from_list("ink_density", ["#F1F4F7", "#8195AE", "#273E5B"])
AA_ORDER = list("CSTV" + "DFHINPWY" + "EQM" + "KR" + "L")
# Alphabetical order has no implied physical continuum; no connecting lines.
AA_ORDER = sorted(set(AA_ORDER))
AUDIT = {"panels": [], "figures": [], "data_downsampling": False,
         "style": "User-supplied reference overrides asset palette and boxed facet titles.",
         "statistical_tests": "None; descriptive figures only."}

# Keep all axis/annotation text readable at the final 183 mm page width.
mpl.rcParams.update({"xtick.labelsize": 8, "ytick.labelsize": 8,
                     "mathtext.fontset": "custom", "mathtext.rm": "Arial",
                     "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold"})


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def compose_figure(rows, cols, height=140, right=.965, top=.885, bottom=.15):
    """Symmetric quantitative_grid at actual journal double-column dimensions."""
    fig, axes = plt.subplots(rows, cols, figsize=(183 / 25.4, height / 25.4), squeeze=False)
    fig.subplots_adjust(left=.095, right=right, bottom=bottom, top=top,
                        wspace=.40 if cols == 2 else .40, hspace=.58)
    return fig, axes


def facet(ax, label, title):
    ax.add_patch(Rectangle((0, 1.025), 1, .115, transform=ax.transAxes,
                           fc="white", ec=BLACK, lw=.7, clip_on=False))
    ax.text(.5, 1.082, title, transform=ax.transAxes, ha="center", va="center", fontsize=8)
    ax.text(-.18, 1.15, label, transform=ax.transAxes, ha="left", va="top",
            fontsize=10, fontweight="bold")
    ax.tick_params(length=3, pad=2)
    ax.yaxis.set_major_locator(MaxNLocator(4))


def source_legend(fig, y=.965):
    handles = [Line2D([], [], color=COLORS[s], marker=MARKERS[s], ms=4,
                      lw=1.2, linestyle=LINES[s], label=LABELS[s]) for s in SOURCES]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(.54, y),
               ncol=3, columnspacing=1.8, handlelength=2.0, frameon=False)


def footer(fig, text):
    fig.text(.095, .026, text, ha="left", va="bottom", fontsize=8, linespacing=1.4, color="#444444")


def finish(fig, stem, out, description):
    """Fixed-canvas vector exports preserve exactly 183 mm rather than recropping."""
    fig.canvas.draw()
    with mpl.rc_context({"savefig.bbox": None}):
        for ext in ["pdf", "svg", "png"]:
            fig.savefig(out / f"{stem}.{ext}", dpi=600 if ext == "png" else 300,
                        bbox_inches=None, facecolor="white")
    widths = fig.get_size_inches() * 25.4
    texts = [t for t in fig.findobj(mpl.text.Text) if t.get_visible() and t.get_text()]
    AUDIT["figures"].append(dict(stem=stem, description=description,
                                  width_mm=float(widths[0]), height_mm=float(widths[1]),
                                  minimum_text_size_pt=min(t.get_fontsize() for t in texts),
                                  formats=["pdf", "svg", "png"], raster_dpi=600))
    plt.close(fig)
    print("Rendered " + stem, flush=True)


def read_data(folder):
    summary = json.loads((folder / "analysis_summary.json").read_text())
    if summary["status"] != "complete" or summary["n_systems"] != 3068:
        raise ValueError("Full validated new-run analysis is required")
    frames = {}
    for name in ["per_system_metrics", "aa_only_per_system_metrics", "hchi_residual_per_system_metrics",
                 "source_fold_metrics", "within_aa_metrics", "uncertainty_coverage", "uncertainty_risk"]:
        path = folder / (name + ".csv")
        if sha256(path) != summary["source_files"][path.name]:
            raise ValueError("Analysis table changed: " + name)
        frames[name] = pd.read_csv(path)
    if sha256(folder / "oof_per_residue.csv.gz") != summary["oof_per_residue_sha256"]:
        raise ValueError("Per-residue source data changed after validation")
    data = pd.read_csv(folder / "oof_per_residue.csv.gz", dtype={
        "sample_id": str, "chain": str, "author_residue_id": str,
        "dataset": "category", "amino_acid": str}, keep_default_na=False)
    for metric, count in summary["valid_residues"].items():
        if data[f"valid_{metric}"].dtype != bool or int(data[f"valid_{metric}"].sum()) != count:
            raise ValueError("Invalid plotting mask or incomplete data")
    AUDIT["analysis_summary"] = str((folder / "analysis_summary.json").resolve())
    AUDIT["analysis_summary_sha256"] = sha256(folder / "analysis_summary.json")
    AUDIT["source_counts"] = summary["source_systems"]
    AUDIT["valid_residues"] = summary["valid_residues"]
    return data, frames


def density_figure(data, tables, out):
    fig, axes = compose_figure(2, 3, height=143, right=.85, top=.87, bottom=.15)
    hexes, max_count = [], 1
    pooled = tables["source_fold_metrics"]
    for row, metric in enumerate(["rmsf", "chi_entropy"]):
        valid_all = data.loc[data[f"valid_{metric}"]]
        p_all, y_all = valid_all[f"pred_{metric}"].to_numpy(float), valid_all[f"target_{metric}"].to_numpy(float)
        transform = np.arcsinh if metric == "rmsf" else lambda x: x
        lower = float(min(0., transform(p_all).min(), transform(y_all).min()))
        upper = float(max(transform(p_all).max(), transform(y_all).max()) * 1.015)
        for col, source in enumerate(SOURCES):
            ax = axes[row, col]
            selected = valid_all.loc[valid_all.dataset == source]
            p, y = selected[f"pred_{metric}"].to_numpy(float), selected[f"target_{metric}"].to_numpy(float)
            h = ax.hexbin(transform(y), transform(p), gridsize=58, mincnt=1,
                          extent=(lower, upper, lower, upper), cmap=DENSITY_CMAP,
                          linewidths=0, rasterized=False)
            if abs(float(h.get_array().sum()) - len(y)) > .01:
                raise ValueError("Density binning lost data points")
            hexes.append(h)
            max_count = max(max_count, float(h.get_array().max()))
            ax.plot([lower, upper], [lower, upper], color="#666666", lw=.7, ls=(0, (3, 3)), zorder=3)
            ax.set_xlim(lower, upper)
            ax.set_ylim(lower, upper)
            ax.set_aspect("equal", adjustable="box")
            facet(ax, "abcdef"[row * 3 + col], f"{LABELS[source]} · {NAMES[metric]}")
            if metric == "rmsf":
                raw_ticks = np.array([-10, -2, 0, 1, 3, 10, 30, 100.])
                ticks = [(transform(v), f"{v:g}") for v in raw_ticks if lower <= transform(v) <= upper]
                ax.set_xticks([x for x, _ in ticks], [s for _, s in ticks])
                ax.set_yticks([x for x, _ in ticks], [s for _, s in ticks])
            else:
                ax.xaxis.set_major_locator(MaxNLocator(4))
                ax.yaxis.set_major_locator(MaxNLocator(4))
            stat = pooled.loc[(pooled.dataset == source) & (pooled.metric == metric)
                              & (pooled.fold == -1) & (pooled.predictor == "FlexPPI")
                              & (pooled.scope == "pooled")].iloc[0]
            ax.text(.045, .955, f"n = {len(y):,}\nρ = {stat.spearman:.3f}", transform=ax.transAxes,
                    fontsize=8, ha="left", va="top", linespacing=1.45)
            ax.set_xlabel(f"MD value ({UNITS[metric]})")
            ax.set_ylabel(f"Prediction ({UNITS[metric]})" if col == 0 else "")
            AUDIT["panels"].append(dict(figure="01", panel="abcdef"[row * 3 + col],
                                          n=len(y), dataset=source, task=metric,
                                          n_definition="valid held-out residues", aggregation="hexagonal counts",
                                          center="none", interval="none", downsampled=0,
                                          axis_transform="asinh, linear width 1 Å" if metric == "rmsf" else "linear"))
    norm = LogNorm(vmin=1, vmax=max_count)
    for h in hexes:
        h.set_norm(norm)
    cax = fig.add_axes([.877, .305, .016, .36])
    cb = fig.colorbar(hexes[0], cax=cax)
    ticks = [10 ** k for k in range(8) if 10 ** k <= max_count]
    cb.set_ticks(ticks, labels=[f"{x:,}" for x in ticks])
    cb.ax.minorticks_off()
    cb.set_label("Residues per hexagon", fontsize=8)
    cb.outline.set_visible(False)
    # Matplotlib otherwise rasterizes long colorbars despite vector panel data.
    cb.solids.set_rasterized(False)
    cb.solids.set_edgecolor("face")  # Avoid white segment seams in PDF viewers.
    footer(fig, "All valid OOF residues; dashed line: prediction = MD value.\nRMSF axes use asinh spacing; tick labels remain in Å. No outliers removed.")
    finish(fig, "fig01_prediction_density", out, "All-residue predicted versus MD-derived values")


def distribution_figure(tables, out):
    systems = tables["per_system_metrics"]
    fig, axes = compose_figure(2, 2, height=147)
    for index, (metric, stat) in enumerate([("rmsf", "spearman"), ("chi_entropy", "spearman"),
                                          ("rmsf", "mae"), ("chi_entropy", "mae")]):
        ax = axes.flat[index]
        facet(ax, "abcd"[index], f"{NAMES[metric]} · {'Spearman' if stat == 'spearman' else 'MAE'}")
        plot_rows = systems.loc[systems.metric == metric]
        use_log = stat == "mae"
        ranges = []
        counts = []
        for x, source in enumerate(SOURCES):
            part = plot_rows.loc[plot_rows.dataset == source].sort_values("sample_id")
            raw_values = part[stat].to_numpy(float)
            if not np.isfinite(raw_values).all() or (use_log and (raw_values <= 0).any()):
                raise ValueError("Unexpected undefined statistic in full system distribution")
            values = np.log10(raw_values) if use_log else raw_values
            counts.append(len(values))
            ranges.extend(values)
            violin = ax.violinplot([values], positions=[x], widths=.64,
                                   showmeans=False, showmedians=False, showextrema=False)
            for body in violin["bodies"]:
                body.set_facecolor(COLORS[source]); body.set_edgecolor(COLORS[source]); body.set_alpha(.43)
                path = body.get_paths()[0]
                path.vertices[:, 0] = np.minimum(path.vertices[:, 0], x)
            jitter = np.random.default_rng(20 + index * 3 + x).uniform(.04, .27, len(values))
            ax.scatter(x + jitter, values, c=COLORS[source], s=2.6, alpha=.24,
                       marker=MARKERS[source], linewidths=0, rasterized=False)
            quartiles = np.quantile(raw_values, [.25, .5, .75])
            q25, median, q75 = np.log10(quartiles) if use_log else quartiles
            ax.plot([x, x], [q25, q75], lw=2.5, color=BLACK, solid_capstyle="round")
            ax.plot(x, median, "o", ms=3.5, mfc="white", mec=BLACK, mew=.6)
        ax.set_xlim(-.5, 2.5)
        ax.set_xticks(range(3), [f"{LABELS[s]}\nn = {n:,}" for s, n in zip(SOURCES, counts)])
        if use_log:
            low, high = min(ranges), max(ranges)
            candidates = ([.001, .003, .01, .03, .1, .3, 1., 3., 10., 30.]
                          if metric == "rmsf" else [.03, .05, .1, .2, .3, .5, 1.])
            ticks = [v for v in candidates if low - .1 <= np.log10(v) <= high + .1]
            ax.set_yticks(np.log10(ticks), [f"{v:g}" for v in ticks])
            ax.set_ylabel(f"MAE ({UNITS[metric]}; log scale)")
            ax.set_ylim(low - .08, high + .10)
        else:
            ax.set_ylabel("Per-system Spearman ρ")
            ax.set_ylim(max(-1.05, min(ranges) - .05), 1.045)
            candidates = np.arange(-1, 1.001, .25 if metric == "rmsf" else .05)
            ax.set_yticks([x for x in candidates if ax.get_ylim()[0] <= x <= 1.001])
        AUDIT["panels"].append(dict(figure="02", panel="abcd"[index], task=metric, statistic=stat,
                                      n=sum(counts), n_definition="one point per whole system",
                                      center="median", interval="interquartile range; violin is KDE",
                                      downsampled=0, kde_space="log10(MAE)" if use_log else "raw correlation"))
    footer(fig, "One point per held-out system; white dot: median; black bar: interquartile range.\nAll 3,068 systems shown in each task. MAE axes are logarithmic.")
    finish(fig, "fig02_system_distributions", out, "Performance distributions beyond five-fold averages")


def strata_figure(tables, out):
    systems = tables["per_system_metrics"]
    fig, axes = compose_figure(2, 2, height=145, top=.86)
    source_legend(fig)
    for index, (metric, xcol) in enumerate([("rmsf", "whole_system_length"), ("chi_entropy", "whole_system_length"),
                                          ("rmsf", "target_iqr"), ("chi_entropy", "target_iqr")]):
        ax = axes.flat[index]
        facet(ax, "abcd"[index], f"{NAMES[metric]} · {'length' if xcol == 'whole_system_length' else 'label variation'}")
        rows = systems.loc[systems.metric == metric]
        if (rows.mae <= 0).any() or not np.isfinite(rows[[xcol, "mae"]]).all().all():
            raise ValueError("Cannot show every system on error axes")
        for source in SOURCES:
            sub = rows.loc[rows.dataset == source]
            ax.scatter(sub[xcol], sub.mae, s=5, c=COLORS[source], marker=MARKERS[source],
                       alpha=.42, edgecolors="none", rasterized=False)
        ax.set_yscale("log")
        ax.yaxis.set_major_locator(FixedLocator([.01, .1, 1, 10] if metric == "rmsf"
                                              else [.03, .05, .1, .2, .3, .5]))
        ax.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:g}"))
        ax.minorticks_off()
        ax.set_ylim(rows.mae.min() / 1.18, rows.mae.max() * 1.18)
        ax.set_ylabel(f"MAE ({UNITS[metric]}; log scale)")
        if xcol == "whole_system_length":
            ax.set_xscale("log")
            ax.set_xticks([50, 100, 300, 1000], ["50", "100", "300", "1,000"])
            ax.set_xlim(rows[xcol].min() / 1.05, rows[xcol].max() * 1.05)
            ax.minorticks_off()
            ax.set_xlabel("Whole-system length (residues; log)")
        else:
            ax.set_xscale("function", functions=(np.arcsinh, np.sinh))
            maximum = rows[xcol].max()
            ticks = [x for x in [0, .5, 1, 2, 5, 10, 20] if x <= maximum * 1.05]
            ax.set_xticks(ticks, [f"{x:g}" for x in ticks])
            ax.set_xlim(0, maximum * 1.05)
            ax.set_xlabel(f"Target IQR ({UNITS[metric]}; asinh)")
        AUDIT["panels"].append(dict(figure="03", panel="abcd"[index], n=len(rows), task=metric,
                                      n_definition="whole systems", x=xcol, y="MAE", center="none",
                                      interval="none", correlation_claim="descriptive; no causal inference", downsampled=0))
    footer(fig, "All 3,068 systems per panel; IQR = the middle 50% range of residue labels.\nLength includes all chains. Source groups have different sequence-sharing protocols.")
    finish(fig, "fig03_error_stratification", out, "Error dependence on system length and target variation")


def paired_fold_panel(ax, table, metric, variants, categories):
    for j, source in enumerate(SOURCES):
        positions = np.array([0, 1]) + (j - 1) * .035
        values = []
        for variant in variants:
            part = table.loc[(table.dataset == source) & (table.fold >= 0) & (table.metric == metric)
                             & (table.predictor == variant) & (table.scope == "macro_by_system")].sort_values("fold")
            if len(part) != 5:
                raise ValueError("Missing fold statistic")
            values.append(part.spearman.to_numpy())
        values = np.asarray(values)
        for fold in range(5):
            ax.plot(positions, values[:, fold], c=COLORS[source], lw=.55, alpha=.25)
            ax.scatter(positions, values[:, fold], c=COLORS[source], s=5, alpha=.5, marker=MARKERS[source], linewidth=0)
        ax.errorbar(positions, values.mean(1), yerr=values.std(1, ddof=1), color=COLORS[source],
                    marker=MARKERS[source], ms=4.5, lw=1.3, ls=LINES[source], capsize=2, capthick=.8)
    ax.set_xlim(-.23, 1.23)
    ax.set_xticks([0, 1], categories)
    ax.set_ylabel("Macro Spearman ρ")


def amino_acid_figure(tables, out):
    fig, axes = compose_figure(2, 2, height=149, top=.86)
    source_legend(fig)
    table = tables["source_fold_metrics"]
    for index, metric in enumerate(["rmsf", "chi_entropy"]):
        ax = axes.flat[index]
        facet(ax, "ab"[index], f"{NAMES[metric]} · training-only baseline")
        paired_fold_panel(ax, table, metric, ["AA-only", "FlexPPI"], ["AA-only", "FlexPPI"])
        ax.set_ylim((.15, 1.02) if metric == "rmsf" else (.84, 1.01))
    ax = axes[1, 0]
    facet(ax, "c", "Hχ · beyond amino-acid type")
    paired_fold_panel(ax, table, "chi_entropy", ["FlexPPI", "AA-corrected"], ["Raw Hχ", "AA-adjusted Hχ"])
    ax.set_ylim(.40, 1.02)
    ax = axes[1, 1]
    facet(ax, "d", "Hχ · within each amino-acid type")
    aa = tables["within_aa_metrics"]
    selected = aa.loc[(aa.metric == "chi_entropy") & (aa.fold >= 0) & aa.dataset.isin(SOURCES)]
    order = sorted(selected.amino_acid.unique())
    for j, source in enumerate(SOURCES):
        group = selected.loc[selected.dataset == source].groupby("amino_acid").spearman.agg(["mean", "std"]).reindex(order)
        ax.errorbar(np.arange(len(order)) + (j - 1) * .2, group["mean"], yerr=group["std"],
                    color=COLORS[source], marker=MARKERS[source], ms=2.6, lw=.65,
                    ls="none", capsize=0, alpha=.95)
    ax.set_xlim(-.8, len(order) - .2)
    ax.set_xticks(range(len(order)), order)
    ax.set_xlabel("Amino-acid one-letter code")
    ax.set_ylabel("Within-type Spearman ρ")
    low = selected.spearman.min()
    ax.set_ylim(min(0., low - .08), 1.02)
    ax.axhline(0, color="#BDBDBD", ls=(0, (3, 3)), lw=.6, zorder=0)
    ax.text(.04, .17, "P: near-constant MD labels", transform=ax.transAxes,
            fontsize=8, color="#555555", va="bottom")
    for index in range(4):
        AUDIT["panels"].append(dict(figure="04", panel="abcd"[index],
                                      n_definition="5 held-out folds; all systems/residues in each stratum",
                                      center="equal-fold mean", interval="sample SD (ddof=1)",
                                      baseline="training-fold AA means only", downsampled=0))
    footer(fig, "Means ± SD across 5 folds; faint points/lines in a–c show individual folds.\nAA-adjusted: subtract training-fold AA means from prediction and target; no refitting.")
    finish(fig, "fig04_amino_acid_controls", out, "Separating amino-acid priors from environment-specific signal")


def uncertainty_figure(tables, out):
    fig, axes = compose_figure(2, 2, height=144, top=.86)
    source_legend(fig)
    cover, risk = tables["uncertainty_coverage"], tables["uncertainty_risk"]
    for column, metric in enumerate(["rmsf", "chi_entropy"]):
        ax = axes[0, column]
        facet(ax, "ab"[column], f"{NAMES[metric]} · Gaussian interval check")
        for source in SOURCES:
            rows = cover.loc[(cover.metric == metric) & (cover.dataset == source) & (cover.variance_kind == "model_only")]
            group = rows.groupby("nominal_coverage").observed_coverage.agg(["mean", "std"])
            x = group.index.to_numpy(float); mean = group["mean"].to_numpy(); sd = group["std"].to_numpy()
            ax.plot(x, mean, c=COLORS[source], marker=MARKERS[source], ms=2.5, ls=LINES[source], lw=1.2)
            ax.fill_between(x, np.clip(mean - sd, 0, 1), np.clip(mean + sd, 0, 1), color=COLORS[source], alpha=.13)
        ax.plot([0, 1], [0, 1], color="#666666", lw=.7, ls=(0, (3, 3)))
        ax.set_xlim(0, 1); ax.set_ylim(0, 1)
        ax.set_xticks([0, .25, .5, .75, 1])
        ax.set_yticks([0, .25, .5, .75, 1])
        ax.set_xlabel("Nominal Gaussian coverage")
        ax.set_ylabel("Observed residue coverage")
        ax = axes[1, column]
        facet(ax, "cd"[column], f"{NAMES[metric]} · uncertainty-based selection")
        maxima = [1]
        minima = [1]
        for source in SOURCES:
            rows = risk.loc[(risk.metric == metric) & (risk.dataset == source)]
            group = rows.groupby("retained_fraction").relative_mae.agg(["mean", "std"])
            x = group.index.to_numpy(float); mean = group["mean"].to_numpy(); sd = group["std"].to_numpy()
            ax.plot(x, mean, c=COLORS[source], marker=MARKERS[source], ms=2.5, ls=LINES[source], lw=1.2)
            ax.fill_between(x, mean - sd, mean + sd, color=COLORS[source], alpha=.13)
            maxima.extend(mean + sd); minima.extend(mean - sd)
        ax.axhline(1, color="#666666", lw=.7, ls=(0, (3, 3)))
        ax.set_xlim(.08, 1.02)
        ax.set_ylim(max(0, min(minima) - .07), max(maxima) + .10)
        ax.set_xticks([.1, .25, .5, .75, 1], ["10", "25", "50", "75", "100"])
        ax.set_xlabel("Retained residues (lowest model SD first; %)")
        ax.set_ylabel("MAE / all-residue MAE")
    footer(fig, "Mean ± SD across 5 folds. Uses the single-model variance head, not an ensemble.\nGaussian intervals are diagnostic only; no post-hoc calibration or label variance used.")
    for index in range(4):
        AUDIT["panels"].append(dict(figure="05", panel="abcd"[index],
                                      n_definition="5 held-out folds, all valid residues",
                                      center="equal-fold mean", interval="sample SD (ddof=1)",
                                      variance="model-only; no label variance or five-model ensemble",
                                      downsampled=0, calibration_fitted=False))
    finish(fig, "fig05_uncertainty_diagnostics", out, "Reliability and usefulness of predicted uncertainty")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis-dir", type=Path, default=ANALYSIS)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    data, tables = read_data(args.analysis_dir)
    density_figure(data, tables, args.output_dir)
    distribution_figure(tables, args.output_dir)
    strata_figure(tables, args.output_dir)
    amino_acid_figure(tables, args.output_dir)
    uncertainty_figure(tables, args.output_dir)
    AUDIT["plot_script_sha256"] = sha256(Path(__file__))
    (args.output_dir / "statistics_and_provenance.json").write_text(json.dumps(AUDIT, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
