# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a) ranked lollipop/dumbbell with CI -> cross-type inherit (BarCategorical) -> param inherit
# (b) joint scatter + marginals -> cross-type inherit (MarginalDensity) -> param inherit
# (c) ECDF per-complex -> cross-type inherit (LineTrend) -> param inherit
# RULE: all panels use drawing functions below copying Class A/B/C values.

# Academic Figure Skill Typography Baseline -- COPY VERBATIM, place at TOP of script
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

# Academic Figure Skill Export Baseline -- COPY VERBATIM
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

import os
import ast
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

OUT = r"E:\Project\FlexPPI\output\flexppi_performance_figures"
os.makedirs(OUT, exist_ok=True)

# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------
SUMMARY_CSV = r"E:\Project\FlexPPI\artifacts\evaluation\additional_baselines\results\itrace_summary_metrics.csv"
PER_COMPLEX_CSV = r"E:\Project\FlexPPI\artifacts\evaluation\additional_baselines\results\itrace_per_complex_metrics.csv"
RESIDUE_CSV = r"E:\Project\FlexPPI\artifacts\evaluation\additional_baselines\results\itrace_common_residue_predictions.csv"

summ = pd.read_csv(SUMMARY_CSV)
summ = summ[(summ.subset == "overall") & (summ.metric == "spearman")].copy()
summ["ci_lo"] = summ.protein_bootstrap_95_ci.map(lambda s: ast.literal_eval(s)[0])
summ["ci_hi"] = summ.protein_bootstrap_95_ci.map(lambda s: ast.literal_eval(s)[1])

method_meta = [
    ("flexppi", "FlexPPI (ours)", ACCENT_RED, 1.0),
    ("backflip_seq", "BackFlip-seq", "#4393C3", 0.9),
    ("fastprotflex", "FastProtFlex", "#4393C3", 0.9),
    ("backflip_structure", "BackFlip-structure", "#4393C3", 0.9),
    ("rmsf_net_pdb", "RMSF-net-PDB", "#4393C3", 0.9),
    ("prody_anm", "ProDy ANM", "#999999", 0.8),
    ("prody_gnm", "ProDy GNM", "#999999", 0.8),
]
rows = []
for key, lab, col, _ in method_meta:
    r = summ[summ.method == key].iloc[0]
    rows.append((lab, r.macro_mean, r.ci_lo, r.ci_hi, col))
methods = pd.DataFrame(rows, columns=["label", "value", "lo", "hi", "color"]).sort_values("value")

per = pd.read_csv(PER_COMPLEX_CSV)
per = per[per.subset == "overall"].copy()
ecdf_cols = {key: f"{key}__spearman" for key, _, _, _ in method_meta}

win_fpf = (per["flexppi__spearman"] > per["fastprotflex__spearman"]).mean()
win_bfs = (per["flexppi__spearman"] > per["backflip_seq__spearman"]).mean()

res = pd.read_csv(RESIDUE_CSV, usecols=["target_chainfit_rmsf_angstrom", "flexppi_pretrain_rmsf_angstrom"])
res = res.dropna()
n_res = len(res)
tgt = res["target_chainfit_rmsf_angstrom"].to_numpy()
prd = res["flexppi_pretrain_rmsf_angstrom"].to_numpy()
pooled_pearson, _ = stats.pearsonr(tgt, prd)

RMSF_MEAN = summ[summ.method == "flexppi"].iloc[0].macro_mean
RMSF_LO = summ[summ.method == "flexppi"].iloc[0].ci_lo
RMSF_HI = summ[summ.method == "flexppi"].iloc[0].ci_hi

# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------
mm = 1 / 25.4
fig = plt.figure(figsize=(183 * mm, 114 * mm))

SEQ_CMAP = LinearSegmentedColormap.from_list("seq", SEQUENTIAL)

# ---- (a) hero: ranked 7-method Spearman ------------------------------------
axa = fig.add_axes([0.055, 0.44, 0.40, 0.50])
y = np.arange(len(methods))[::-1]
for yi, m in zip(y, methods.itertuples()):
    axa.hlines(yi, 0.30, m.value, color=m.color, linewidth=1.0, zorder=2, alpha=0.9)
    axa.errorbar(m.value, yi, xerr=[[m.value - m.lo], [m.hi - m.value]],
                 fmt="o", ms=4.2 if m.label != "FlexPPI (ours)" else 6.2,
                 color=m.color, ecolor="#555555", elinewidth=0.7, capsize=2.2,
                 zorder=3, linewidth=0.5)
    axa.text(m.hi + 0.022, yi, f"{m.value:.3f}", va="center", ha="left",
             fontsize=7, color=BLACK, fontweight="bold" if m.label == "FlexPPI (ours)" else "normal")
axa.set_yticks(y)
axa.set_yticklabels(methods.label.values, fontsize=7)
axa.set_xlim(0.30, 0.87)
axa.set_xticks(np.arange(0.30, 0.90, 0.10))
axa.set_xlabel("Spearman $\\rho$ vs MD RMSF (244 complexes)", fontsize=7)
axa.set_title("Seven static-structure methods", fontsize=8, loc="left", pad=4)
axa.text(-0.11, 1.08, "a", transform=axa.transAxes, fontsize=9, fontweight="bold", va="top")

# ---- (b) joint scatter MD vs FlexPPI RMSF ------------------------------------
bx = fig.add_axes([0.57, 0.50, 0.24, 0.42])
bx_top = fig.add_axes([0.57, 0.945, 0.24, 0.035], sharex=bx)
bx_right = fig.add_axes([0.82, 0.50, 0.06, 0.42], sharey=bx)
bx_top.axis("off")
bx_right.axis("off")

lo, hi = 0.0, 10.0
hb = bx.hexbin(tgt, prd, gridsize=48, bins="log", cmap=SEQ_CMAP,
               mincnt=1, linewidths=0.0, extent=(lo, hi, lo, hi))
bx.plot([lo, hi], [lo, hi], color=GREY, linestyle="--", linewidth=0.6, alpha=0.6, zorder=2)
bx.set_xlim(lo, hi)
bx.set_ylim(lo, hi)

bin_edges = np.linspace(lo, hi, 50)
bx_top.hist(tgt, bins=bin_edges, color="#9ECAE1", edgecolor="white", linewidth=0.1)
bx_right.hist(prd, bins=bin_edges, color="#9ECAE1", orientation="horizontal",
              edgecolor="white", linewidth=0.1)
bx_top.set_ylim(0, bx_top.get_ylim()[1] * 1.4)
bx_right.set_xlim(0, bx_right.get_xlim()[1] * 1.4)
bx.set_xlabel("MD RMSF (\\AA)", fontsize=7)
bx.set_ylabel("FlexPPI RMSF (\\AA)", fontsize=7)
bx.text(-0.17, 1.08, "b", transform=bx.transAxes, fontsize=9, fontweight="bold", va="top")

# colorbar + statistics BELOW panel (b), outside the plot area
cax = fig.add_axes([0.57, 0.46, 0.30, 0.014])
cb = fig.colorbar(hb, cax=cax, orientation="horizontal")
cb.set_label("log$_{10}$(residue count)", fontsize=6.5)
cb.ax.tick_params(labelsize=6)
c_lo, c_hi = hb.get_clim()
c_ticks = [10.0 ** k for k in range(int(np.ceil(np.log10(c_lo))),
                                    int(np.floor(np.log10(c_hi))) + 1)]
cb.set_ticks(c_ticks)
c_labels = []
for v in c_ticks:
    k = int(round(np.log10(v)))
    c_labels.append("1" if k == 0 else ("10" if k == 1 else
                    ("100" if k == 2 else f"$10^{{{k}}}$")))
cb.set_ticklabels(c_labels)

fig.text(0.57, 0.428, "per-chain Spearman (complex mean):", fontsize=6.5, color=BLACK)
fig.text(0.57, 0.406, f"$\\rho$ = {RMSF_MEAN:.3f}   [{RMSF_LO:.3f}, {RMSF_HI:.3f}]",
         fontsize=7, color=ACCENT_RED, fontweight="bold")
fig.text(0.57, 0.384, f"pooled Pearson = {pooled_pearson:.3f}", fontsize=6.5, color=BLACK)
fig.text(0.57, 0.362, f"n = {n_res:,} residues, 244 complexes", fontsize=6.5, color=BLACK)

# ---- (c) ECDF of per-complex Spearman -----------------------------------------
axc = fig.add_axes([0.055, 0.06, 0.90, 0.26])
order = methods.label.values[::-1]
color_map = dict(zip(methods.label, methods.color))
for lab in order:
    key = [k for k, l, _, _ in method_meta if l == lab][0]
    col = ecdf_cols[key]
    vals = np.sort(per[col].dropna().values)
    ecdf = np.arange(1, len(vals) + 1) / len(vals)
    lw = 1.6 if lab == "FlexPPI (ours)" else 0.8
    alpha = 1.0 if lab == "FlexPPI (ours)" else 0.75
    axc.plot(vals, ecdf, color=color_map[lab], linewidth=lw, alpha=alpha,
             label=lab, zorder=3 if lab == "FlexPPI (ours)" else 2)
axc.axvline(0.5, color="#E0E0E0", linewidth=0.4, linestyle=":", zorder=1)
axc.text(0.76, 0.06,
         f"FlexPPI wins per complex vs FastProtFlex in {win_fpf:.0%},\n"
         f"vs BackFlip-seq in {win_bfs:.0%}",
         fontsize=6.5, color=BLACK, ha="left", va="bottom")
axc.set_xlim(0.0, 1.0)
axc.set_ylim(0.0, 1.04)
axc.set_xlabel("Per-complex Spearman $\\rho$ vs MD RMSF", fontsize=7)
axc.set_ylabel("Fraction of complexes", fontsize=7)
axc.legend(loc="upper left", fontsize=6, ncol=2, handlelength=1.5, frameon=False)
axc.text(-0.055, 1.22, "c", transform=axc.transAxes, fontsize=9, fontweight="bold", va="top")

save_cns_figure(fig, os.path.join(OUT, "fig_flexppi_rmsf_benchmark"))
print("saved", os.path.join(OUT, "fig_flexppi_rmsf_benchmark"))
