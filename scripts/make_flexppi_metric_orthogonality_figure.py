# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a) per-complex orthogonality scatter + marginals -> cross-type inherit (MarginalDensity) -> param inherit
# (b) within-chain MD descriptor correlation histogram + KDE -> cross-type inherit (KernelDensity) -> param inherit
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

OUT = r"E:\Project\FlexPPI\output\flexppi_performance_figures"
os.makedirs(OUT, exist_ok=True)

# per-complex RMSF vs H_chi performance (chainfit RMSF label definition)
v = pd.read_csv(r"E:\Project\FlexPPI\output\itrace_hchi_benchmark\metrics\per_complex_rmsf_vs_hchi.csv")
v = v[v.subset == "overall"].groupby("pdb_id").head(1)
x = v.rmsf_spearman.to_numpy()          # per-complex FlexPPI RMSF Spearman
y = v.flexppi__spearman.to_numpy()      # per-complex FlexPPI H_chi Spearman
pr, pp = stats.pearsonr(x, y)
sr, sp = stats.spearmanr(x, y)

# within-chain MD RMSF <-> H_chi (both are MD ground-truth descriptors)
h = pd.read_csv(r"E:\Project\FlexPPI\output\itrace_hchi_benchmark\raw\hchi_consensus_per_residue.csv",
                usecols=["pdb_id", "chain", "residue_number", "insertion_code", "target_chi_entropy"])
r = pd.read_csv(r"E:\Project\FlexPPI\artifacts\evaluation\additional_baselines\results\itrace_common_residue_predictions.csv",
                usecols=["pdb_id", "chain", "residue_number", "insertion_code", "target_chainfit_rmsf_angstrom"])
m = h.merge(r, on=["pdb_id", "chain", "residue_number", "insertion_code"], how="inner")
m = m[m.target_chi_entropy.notna() & m.target_chainfit_rmsf_angstrom.notna()]
within = []
for _, d in m.groupby(["pdb_id", "chain"]):
    if len(d) >= 4:
        spc, _ = stats.spearmanr(d.target_chi_entropy, d.target_chainfit_rmsf_angstrom)
        within.append(spc)
within = np.array(within)

# ---------------------------------------------------------------------------
# Figure (single column)
# ---------------------------------------------------------------------------
mm = 1 / 25.4
fig = plt.figure(figsize=(89 * mm, 136 * mm))

# ---- (a) per-complex performance orthogonality ------------------------------
ax = fig.add_axes([0.15, 0.63, 0.55, 0.30])
ax_top = fig.add_axes([0.15, 0.945, 0.55, 0.03], sharex=ax)
ax_right = fig.add_axes([0.715, 0.63, 0.06, 0.30], sharey=ax)
ax_top.axis("off")
ax_right.axis("off")

ax.scatter(x, y, s=8, color="#4393C3", alpha=0.65, edgecolors="none", zorder=3)
ax.set_xlim(0.2, 1.0)
ax.set_ylim(0.65, 1.0)
ax.set_xlabel("FlexPPI RMSF Spearman $\\rho$", fontsize=7)
ax.set_ylabel("FlexPPI $H_\\chi$ Spearman $\\rho$", fontsize=7)

binx = np.linspace(0.2, 1.0, 24)
biny = np.linspace(0.65, 1.0, 24)
ax_top.hist(x, bins=binx, color="#9ECAE1", edgecolor="white", linewidth=0.2)
ax_right.hist(y, bins=biny, color="#9ECAE1", orientation="horizontal",
              edgecolor="white", linewidth=0.2)
ax_top.set_ylim(0, ax_top.get_ylim()[1] * 1.4)
ax_right.set_xlim(0, ax_right.get_xlim()[1] * 1.4)
ax.text(-0.22, 1.08, "a", transform=ax.transAxes, fontsize=9, fontweight="bold", va="top")

# statistics BELOW panel (a), outside the plot area
fig.text(0.15, 0.585, f"Pearson r = {pr:.3f}   (p = {pp:.2g})",
         fontsize=7, color=BLACK)
fig.text(0.15, 0.563, f"Spearman $\\rho$ = {sr:.3f}   (p = {sp:.2g})",
         fontsize=7, color=BLACK)
fig.text(0.15, 0.541, f"n = {len(x)} complexes", fontsize=7, color=BLACK)

# ---- (b) within-chain MD descriptor correlation -----------------------------
axb = fig.add_axes([0.15, 0.14, 0.62, 0.33])
axb.hist(within, bins=np.linspace(-1.0, 1.0, 41), color="#BDBDBD",
         edgecolor="white", linewidth=0.3, density=True, zorder=3)
kde = stats.gaussian_kde(within)
xx = np.linspace(-1.0, 1.0, 300)
axb.plot(xx, kde(xx), color="#2166AC", linewidth=1.0, zorder=4)
axb.axvline(within.mean(), color=ACCENT_RED, linestyle="--", linewidth=1.0, zorder=4)
axb.annotate(f"mean {within.mean():.3f}",
             xy=(within.mean(), axb.get_ylim()[1] * 0.55),
             xytext=(0.72, axb.get_ylim()[1] * 0.92),
             fontsize=6.5, color=ACCENT_RED,
             arrowprops=dict(arrowstyle="-", color=ACCENT_RED, lw=0.7))
axb.set_xlim(-1.0, 1.0)
axb.set_xlabel("Within-chain MD Spearman $\\rho$ (RMSF vs $H_\\chi$)", fontsize=7)
axb.set_ylabel("Density", fontsize=7)
axb.text(-0.17, 1.08, "b", transform=axb.transAxes, fontsize=9, fontweight="bold", va="top")
axb.text(0.02, 0.97,
         f"n = {len(within):,} chains (both MD labels available)",
         transform=axb.transAxes, fontsize=6.5, color="#555555", va="top")

save_cns_figure(fig, os.path.join(OUT, "fig_flexppi_metric_orthogonality"))
print("saved", os.path.join(OUT, "fig_flexppi_metric_orthogonality"))
