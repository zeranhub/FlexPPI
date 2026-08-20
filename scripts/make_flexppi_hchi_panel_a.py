# Standalone panel (a) of the H_chi performance figure:
# MD-vs-prediction joint density (hexbin) + marginals + colorbar.
# No statistics text inside or below the panel; axis tick formats are
# identical on both axes (plain integers 0-4).

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

OUT = r"E:\Project\FlexPPI\output\flexppi_performance_figures"
os.makedirs(OUT, exist_ok=True)

RESIDUE_CSV = r"E:\Project\FlexPPI\output\itrace_hchi_benchmark\raw\hchi_consensus_per_residue.csv"

res = pd.read_csv(RESIDUE_CSV, low_memory=False)
res = res[["target_chi_entropy", "flexppi_chi_entropy"]].dropna()
target = res["target_chi_entropy"].to_numpy()
pred = res["flexppi_chi_entropy"].to_numpy()

SEQ_CMAP = LinearSegmentedColormap.from_list("seq", SEQUENTIAL)
mm = 1 / 25.4
fig = plt.figure(figsize=(89 * mm, 86 * mm))

# main axes + marginals
ax = fig.add_axes([0.13, 0.18, 0.60, 0.68])
ax_top = fig.add_axes([0.13, 0.88, 0.60, 0.055], sharex=ax)
ax_right = fig.add_axes([0.745, 0.18, 0.06, 0.68], sharey=ax)
ax_top.axis("off")
ax_right.axis("off")

lo, hi = 0.0, 4.05
hb = ax.hexbin(target, pred, gridsize=55, bins="log", cmap=SEQ_CMAP,
               mincnt=1, linewidths=0.0, extent=(lo, hi, lo, hi))
ax.plot([lo, hi], [lo, hi], color=GREY, linestyle="--", linewidth=0.6, alpha=0.6, zorder=2)
ax.set_xlim(lo, hi)
ax.set_ylim(lo, hi)

# identical integer ticks on both axes
ticks = [0, 1, 2, 3, 4]
ax.set_xticks(ticks)
ax.set_yticks(ticks)
ax.set_xticklabels([f"{t}" for t in ticks])
ax.set_yticklabels([f"{t}" for t in ticks])

bin_edges = np.linspace(lo, hi, 60)
ax_top.hist(target, bins=bin_edges, color="#9ECAE1", edgecolor="white", linewidth=0.2)
ax_right.hist(pred, bins=bin_edges, color="#9ECAE1", orientation="horizontal",
              edgecolor="white", linewidth=0.2)
ax_top.set_ylim(0, ax_top.get_ylim()[1] * 1.35)
ax_right.set_xlim(0, ax_right.get_xlim()[1] * 1.35)

ax.set_xlabel("MD $H_\\chi$ (rotamer entropy, nats)")
ax.set_ylabel("FlexPPI $H_\\chi$ prediction")

save_cns_figure(fig, os.path.join(OUT, "fig_flexppi_hchi_panel_a"))
print("saved", os.path.join(OUT, "fig_flexppi_hchi_panel_a"))
