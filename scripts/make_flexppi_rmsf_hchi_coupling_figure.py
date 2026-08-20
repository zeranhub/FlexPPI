# Experiment: does the model reproduce the MD-truth residue-level coupling
# between RMSF and H_chi? Two-panel version:
# (a) MD truth: residue RMSF vs residue H_chi
# (b) predictions (full input): residue RMSF vs H_chi
# No colorbars; identical integer tick formatting on both axes of both panels;
# the Spearman annotations share one color.

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

KEYS = ["pdb_id", "chain", "residue_number", "insertion_code"]
ALL = r"E:\Project\FlexPPI\output\itrace_residue_biology\raw\all_condition_predictions.csv"
HCHI = r"E:\Project\FlexPPI\output\itrace_hchi_benchmark\raw\hchi_consensus_per_residue.csv"
CC = r"E:\Project\FlexPPI\output\itrace_zero_shot_context\raw\common_context_predictions.csv"

df = pd.read_csv(ALL, usecols=KEYS + ["condition", "rmsf", "chi_entropy"],
                 low_memory=False)
h = pd.read_csv(HCHI, usecols=KEYS + ["target_chi_entropy"], low_memory=False)
cc = pd.read_csv(CC, usecols=KEYS + ["target_chainfit_rmsf_angstrom"], low_memory=False)
df = df.merge(h, on=KEYS, how="left").merge(cc, on=KEYS, how="left")
full = df[(df.condition == "full") & df.target_chainfit_rmsf_angstrom.notna()
          & df.target_chi_entropy.notna() & df.rmsf.notna() & df.chi_entropy.notna()]

rho_md, _ = stats.spearmanr(full.target_chainfit_rmsf_angstrom, full.target_chi_entropy)
rho_pd, _ = stats.spearmanr(full.rmsf, full.chi_entropy)
print(f"MD rho={rho_md:.3f} | pred rho={rho_pd:.3f} | n={len(full)}")

# ---------------------------------------------------------------------------
# Figure (two panels)
# ---------------------------------------------------------------------------
mm = 1 / 25.4
SEQ_CMAP = LinearSegmentedColormap.from_list("seq", SEQUENTIAL)
fig = plt.figure(figsize=(183 * mm, 88 * mm))

def draw_joint(xpos, xdata, ydata, xlabel, ylabel, rho_text):
    ax = fig.add_axes([xpos, 0.20, 0.36, 0.66])
    ax_top = fig.add_axes([xpos, 0.885, 0.36, 0.06], sharex=ax)
    ax_right = fig.add_axes([xpos + 0.37, 0.20, 0.045, 0.66], sharey=ax)
    ax_top.axis("off"); ax_right.axis("off")

    ax.hexbin(xdata, ydata, gridsize=45, bins="log", cmap=SEQ_CMAP,
              mincnt=1, linewidths=0.0, extent=(0, 10, 0, 4))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 4)
    # identical integer tick formatting on both axes
    ax.set_xticks([0, 2, 4, 6, 8, 10])
    ax.set_yticks([0, 1, 2, 3, 4])
    ax.set_xticklabels(["0", "2", "4", "6", "8", "10"])
    ax.set_yticklabels(["0", "1", "2", "3", "4"])

    ax_top.hist(xdata, bins=np.linspace(0, 10, 50),
                color="#9ECAE1", edgecolor="white", linewidth=0.2)
    ax_right.hist(ydata, bins=np.linspace(0, 4, 40),
                  color="#9ECAE1", orientation="horizontal",
                  edgecolor="white", linewidth=0.2)
    ax_top.set_ylim(0, ax_top.get_ylim()[1] * 1.4)
    ax_right.set_xlim(0, ax_right.get_xlim()[1] * 1.4)

    ax.set_xlabel(xlabel, fontsize=7)
    ax.set_ylabel(ylabel, fontsize=7)
    ax.text(0.97, 0.035, rho_text, transform=ax.transAxes, fontsize=6.5,
            color=BLACK, ha="right", va="bottom")
    return ax

axa = draw_joint(0.055, full.target_chainfit_rmsf_angstrom,
                 full.target_chi_entropy,
                 "MD RMSF (\\AA)", "MD $H_\\chi$ (nats)",
                 f"Spearman $\\rho$ = {rho_md:.3f}")
axa.text(-0.13, 1.15, "a", transform=axa.transAxes, fontsize=9,
         fontweight="bold", va="top")

axb = draw_joint(0.52, full.rmsf, full.chi_entropy,
                 "Predicted RMSF (\\AA)", "Predicted $H_\\chi$ (nats)",
                 f"Spearman $\\rho$ = {rho_pd:.3f}")
axb.text(-0.13, 1.15, "b", transform=axb.transAxes, fontsize=9,
         fontweight="bold", va="top")

save_cns_figure(fig, os.path.join(OUT, "fig_flexppi_rmsf_hchi_coupling"))
print("saved", os.path.join(OUT, "fig_flexppi_rmsf_hchi_coupling"))
