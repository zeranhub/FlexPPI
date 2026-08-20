# Standalone RMSF joint-scatter panel (b of the RMSF benchmark figure) with
# the same treatment as the H_chi panel (a): no colorbar, no statistics text,
# identical integer tick formatting on both axes.
# Also renders the combined H_chi + RMSF two-panel joint-scatter figure.

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
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

OUT = r"E:\Project\FlexPPI\output\flexppi_performance_figures"
os.makedirs(OUT, exist_ok=True)

SEQ_CMAP = LinearSegmentedColormap.from_list("seq", SEQUENTIAL)

# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------
h = pd.read_csv(r"E:\Project\FlexPPI\output\itrace_hchi_benchmark\raw\hchi_consensus_per_residue.csv",
                low_memory=False)
h = h[["target_chi_entropy", "flexppi_chi_entropy"]].dropna()
h_target = h.target_chi_entropy.to_numpy()
h_pred = h.flexppi_chi_entropy.to_numpy()

r = pd.read_csv(r"E:\Project\FlexPPI\artifacts\evaluation\additional_baselines\results\itrace_common_residue_predictions.csv",
                usecols=["target_chainfit_rmsf_angstrom", "flexppi_pretrain_rmsf_angstrom"])
r = r.dropna()
r_target = r.target_chainfit_rmsf_angstrom.to_numpy()
r_pred = r.flexppi_pretrain_rmsf_angstrom.to_numpy()

# ---------------------------------------------------------------------------
# Drawing helper: joint density + marginals, no colorbar, no stats text
# ---------------------------------------------------------------------------
def draw_joint(ax, ax_top, ax_right, target, pred, lo, hi, ticks,
               xlabel, ylabel):
    hb = ax.hexbin(target, pred, gridsize=55, bins="log", cmap=SEQ_CMAP,
                   mincnt=1, linewidths=0.0, extent=(lo, hi, lo, hi))
    ax.plot([lo, hi], [lo, hi], color=GREY, linestyle="--", linewidth=0.6,
            alpha=0.6, zorder=2)
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    ax.set_xticklabels([f"{t}" for t in ticks])
    ax.set_yticklabels([f"{t}" for t in ticks])

    bin_edges = np.linspace(lo, hi, 60)
    ax_top.hist(target, bins=bin_edges, color="#9ECAE1", edgecolor="white",
                linewidth=0.2)
    ax_right.hist(pred, bins=bin_edges, color="#9ECAE1", orientation="horizontal",
                  edgecolor="white", linewidth=0.2)
    ax_top.set_ylim(0, ax_top.get_ylim()[1] * 1.35)
    ax_right.set_xlim(0, ax_right.get_xlim()[1] * 1.35)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    return hb

# ---------------------------------------------------------------------------
# 1) Combined H_chi + RMSF figure (183 mm)
# ---------------------------------------------------------------------------
mm = 1 / 25.4
fig = plt.figure(figsize=(183 * mm, 100 * mm))

ax_a = fig.add_axes([0.06, 0.18, 0.36, 0.68])
ax_a_top = fig.add_axes([0.06, 0.88, 0.36, 0.055], sharex=ax_a)
ax_a_right = fig.add_axes([0.43, 0.18, 0.045, 0.68], sharey=ax_a)
ax_a_top.axis("off")
ax_a_right.axis("off")
draw_joint(ax_a, ax_a_top, ax_a_right, h_target, h_pred, 0.0, 4.05,
           [0, 1, 2, 3, 4],
           "MD $H_\\chi$ (rotamer entropy, nats)",
           "FlexPPI $H_\\chi$ prediction")
ax_a.text(-0.14, 1.05, "a", transform=ax_a.transAxes, fontsize=9,
          fontweight="bold", va="top")

ax_b = fig.add_axes([0.54, 0.18, 0.36, 0.68])
ax_b_top = fig.add_axes([0.54, 0.88, 0.36, 0.055], sharex=ax_b)
ax_b_right = fig.add_axes([0.91, 0.18, 0.045, 0.68], sharey=ax_b)
ax_b_top.axis("off")
ax_b_right.axis("off")
draw_joint(ax_b, ax_b_top, ax_b_right, r_target, r_pred, 0.0, 10.0,
           [0, 2, 4, 6, 8, 10],
           "MD RMSF (\\AA)",
           "FlexPPI RMSF (\\AA)")
ax_b.text(-0.14, 1.05, "b", transform=ax_b.transAxes, fontsize=9,
          fontweight="bold", va="top")

save_cns_figure(fig, os.path.join(OUT, "fig_flexppi_hchi_rmsf_joint"))
print("saved", os.path.join(OUT, "fig_flexppi_hchi_rmsf_joint"))

# ---------------------------------------------------------------------------
# 2) Standalone RMSF panel b (89 mm, same geometry as the H_chi standalone)
# ---------------------------------------------------------------------------
fig2 = plt.figure(figsize=(89 * mm, 86 * mm))
ax2 = fig2.add_axes([0.13, 0.18, 0.60, 0.68])
ax2_top = fig2.add_axes([0.13, 0.88, 0.60, 0.055], sharex=ax2)
ax2_right = fig2.add_axes([0.745, 0.18, 0.06, 0.68], sharey=ax2)
ax2_top.axis("off")
ax2_right.axis("off")
draw_joint(ax2, ax2_top, ax2_right, r_target, r_pred, 0.0, 10.0,
           [0, 2, 4, 6, 8, 10],
           "MD RMSF (\\AA)",
           "FlexPPI RMSF (\\AA)")
save_cns_figure(fig2, os.path.join(OUT, "fig_flexppi_rmsf_panel_b"))
print("saved", os.path.join(OUT, "fig_flexppi_rmsf_panel_b"))
