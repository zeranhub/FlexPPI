# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a) joint scatter + marginals -> cross-type inherit (MarginalDensity) -> param inherit
# (b) ranked horizontal bar with CI -> cross-type inherit (BarCategorical) -> param inherit
# (c) histogram + KDE + reference lines -> cross-type inherit (KernelDensity) -> param inherit
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
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D

OUT = r"E:\Project\FlexPPI\output\flexppi_performance_figures"
os.makedirs(OUT, exist_ok=True)

# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------
RESIDUE_CSV = r"E:\Project\FlexPPI\output\itrace_hchi_benchmark\raw\hchi_consensus_per_residue.csv"
SUMMARY_CSV = r"E:\Project\FlexPPI\output\itrace_hchi_benchmark\metrics\summary_metrics.csv"
REPLICA_CSV = r"E:\Project\FlexPPI\output\itrace_hchi_benchmark\metrics\replica_reliability.csv"

res = pd.read_csv(RESIDUE_CSV, low_memory=False)
res = res[["target_chi_entropy", "flexppi_chi_entropy"]].dropna()
n_res = len(res)
target = res["target_chi_entropy"].to_numpy()
pred = res["flexppi_chi_entropy"].to_numpy()

pooled_pearson, _ = stats.pearsonr(target, pred)
print(f"Hchi residues={n_res} pooled pearson={pooled_pearson:.4f}")

summ = pd.read_csv(SUMMARY_CSV)
summ = summ[(summ.cohort == "external_clean") & (summ.method == "flexppi") & (summ.metric == "spearman")]
overall = summ[summ.subset == "overall"].iloc[0]
HCHI_MEAN = overall.macro_mean
HCHI_LO, HCHI_HI = overall.ci95_low, overall.ci95_high

repl = pd.read_csv(REPLICA_CSV)
repl = repl[repl.pdb_id.str.lower() != "3kxf"]
rep_mean = repl.mean_pairwise_spearman.mean()

AA_PRIOR = 0.8818686970496346

regions = [
    ("non_interface", "Non-interface"),
    ("mhc", "MHC"),
    ("tcr_alpha", "TCR alpha"),
    ("tcr_framework", "TCR framework"),
    ("tcr_beta", "TCR beta"),
    ("overall", "Overall"),
    ("cdr_all", "All CDR"),
    ("cdr2", "CDR2"),
    ("cdr3", "CDR3"),
    ("interface", "Interface"),
    ("peptide", "Peptide"),
    ("cdr1", "CDR1"),
]
rows = []
for key, lab in regions:
    r = summ[summ.subset == key].iloc[0]
    rows.append((lab, r.macro_mean, r.ci95_low, r.ci95_high, r.complex_count))
region_df = pd.DataFrame(rows, columns=["label", "value", "lo", "hi", "n"]).sort_values("value")

# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------
mm = 1 / 25.4
fig = plt.figure(figsize=(183 * mm, 124 * mm))

SEQ_CMAP = LinearSegmentedColormap.from_list("seq", SEQUENTIAL)

# ---- (a) hero: joint scatter MD vs FlexPPI H_chi --------------------------
ax = fig.add_axes([0.055, 0.47, 0.34, 0.45])
ax_top = fig.add_axes([0.055, 0.945, 0.34, 0.035], sharex=ax)
ax_right = fig.add_axes([0.405, 0.47, 0.045, 0.45], sharey=ax)
ax_top.axis("off")
ax_right.axis("off")

lo = min(target.min(), pred.min())
hi = max(target.max(), pred.max())
hb = ax.hexbin(target, pred, gridsize=55, bins="log", cmap=SEQ_CMAP,
               mincnt=1, linewidths=0.0, extent=(lo, hi, lo, hi))
ax.plot([lo, hi], [lo, hi], color=GREY, linestyle="--", linewidth=0.6, alpha=0.6, zorder=2)
ax.set_xlim(lo, hi)
ax.set_ylim(lo, hi)

bin_edges = np.linspace(lo, hi, 60)
ax_top.hist(target, bins=bin_edges, color="#9ECAE1", edgecolor="white", linewidth=0.2)
ax_right.hist(pred, bins=bin_edges, color="#9ECAE1", orientation="horizontal",
              edgecolor="white", linewidth=0.2)
ax_top.set_ylim(0, ax_top.get_ylim()[1] * 1.35)
ax_right.set_xlim(0, ax_right.get_xlim()[1] * 1.35)

ax.set_xlabel("MD $H_\\chi$ (rotamer entropy, nats)")
ax.set_ylabel("FlexPPI $H_\\chi$ prediction")
ax.text(-0.14, 1.05, "a", transform=ax.transAxes, fontsize=9, fontweight="bold", va="top")

# colorbar + statistics BELOW panel (a), outside the plot area
cax = fig.add_axes([0.055, 0.425, 0.34, 0.014])
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

fig.text(0.055, 0.390, "per-chain Spearman (complex mean):", fontsize=6.5, color=BLACK)
fig.text(0.055, 0.367, f"$\\rho$ = {HCHI_MEAN:.3f}   [{HCHI_LO:.3f}, {HCHI_HI:.3f}]",
         fontsize=7, color=ACCENT_RED, fontweight="bold")
fig.text(0.055, 0.345, f"pooled Pearson = {pooled_pearson:.3f}", fontsize=6.5, color=BLACK)
fig.text(0.055, 0.323, f"n = {n_res:,} residues, 244 complexes", fontsize=6.5, color=BLACK)

# ---- (b) region-ranked Spearman bars --------------------------------------
axb = fig.add_axes([0.50, 0.32, 0.44, 0.63])
y = np.arange(len(region_df))
err = np.vstack([region_df.value.values - region_df.lo.values,
                 region_df.hi.values - region_df.value.values])
colors = [ACCENT_RED if lab == "Overall" else "#4393C3" for lab in region_df.label]
axb.barh(y, region_df.value.values, xerr=err, color=colors,
         height=0.62, error_kw=dict(ecolor="#555555", lw=0.7, capsize=2.2), zorder=3)
axb.axvline(AA_PRIOR, color=GREY, linestyle="--", linewidth=0.8, zorder=2)
axb.text(AA_PRIOR + 0.002, len(region_df) - 0.20, "AA-type prior = 0.882",
         color="#444444", fontsize=6, va="bottom")
for yi, r in zip(y, region_df.itertuples()):
    axb.text(r.value + err[1][yi] + 0.016, yi, f"{r.value:.3f}",
             va="center", ha="left", fontsize=6.5, color=BLACK)
axb.set_yticks(y)
axb.set_yticklabels(region_df.label.values, fontsize=6.5)
axb.set_xlim(0.70, 1.12)
axb.set_xlabel("Spearman $\\rho$ (MD $H_\\chi$)", fontsize=7)
axb.text(-0.05, 1.04, "b", transform=axb.transAxes, fontsize=9, fontweight="bold", va="top")

# ---- (c) replica reliability, full-width bottom panel -----------------------
axc = fig.add_axes([0.055, 0.05, 0.90, 0.22])
vals = repl.mean_pairwise_spearman.to_numpy()
axc.hist(vals, bins=np.linspace(0.4, 1.0, 25), color="#BDBDBD",
         edgecolor="white", linewidth=0.3, density=True, zorder=3)
kde = stats.gaussian_kde(vals)
xx = np.linspace(0.4, 1.0, 300)
axc.plot(xx, kde(xx), color="#2166AC", linewidth=1.0, zorder=4)
axc.axvline(rep_mean, color=BLACK, linestyle="-", linewidth=1.0, zorder=4)
axc.axvline(HCHI_MEAN, color=ACCENT_RED, linestyle="--", linewidth=1.0, zorder=4)
axc.set_xlim(0.40, 1.0)
axc.set_xlabel("Pairwise replica Spearman $\\rho$", fontsize=7)
axc.set_ylabel("Density", fontsize=7)
axc.text(-0.035, 1.16, "c", transform=axc.transAxes, fontsize=9, fontweight="bold", va="top")
axc.legend(handles=[
    Line2D([0], [0], color=BLACK, lw=1.0, label=f"MD replicas mean {rep_mean:.3f}"),
    Line2D([0], [0], color=ACCENT_RED, lw=1.0, ls="--", label=f"FlexPPI {HCHI_MEAN:.3f}"),
], loc="upper left", fontsize=6.5, frameon=False, handlelength=1.5)

save_cns_figure(fig, os.path.join(OUT, "fig_flexppi_hchi_performance"))
print("saved", os.path.join(OUT, "fig_flexppi_hchi_performance"))
