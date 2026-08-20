# RMSF input-context ablation, same visual language as the H_chi panel (a):
# grouped bars -- chain role x input context (full complex / partner-removed
# pMHC-TCRab dimer / isolated chain), y = per-chain Spearman vs MD RMSF,
# error bars = complex-bootstrap 95% CI. No value labels, no grey notes.

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
from matplotlib.lines import Line2D

OUT = r"E:\Project\FlexPPI\output\flexppi_performance_figures"
os.makedirs(OUT, exist_ok=True)

PARTNER_CSV = r"E:\Project\FlexPPI\output\itrace_partner_ablation\per_chain_rmsf_metrics.csv"
ZERO_SHOT_CSV = r"E:\Project\FlexPPI\output\itrace_zero_shot_context\metrics\per_chain_metrics.csv"

pa = pd.read_csv(PARTNER_CSV)
zs = pd.read_csv(ZERO_SHOT_CSV)

ROLES = [
    ("peptide", "Peptide"),
    ("tcr_beta", "TCR beta"),
    ("tcr_alpha", "TCR alpha"),
    ("mhc", "MHC"),
    ("b2m", "beta2m"),
]

rng = np.random.default_rng(2026)

def boot_ci(values):
    boot = np.array([rng.choice(values, size=len(values), replace=True).mean()
                     for _ in range(2000)])
    return np.percentile(boot, 2.5), np.percentile(boot, 97.5)

groups = []
for role, lab in ROLES:
    fp = pa[pa.chain_role == role].groupby("pdb_id").full_spearman.first()
    pr = pa[pa.chain_role == role].groupby("pdb_id").partner_removed_spearman.first()
    iso = zs[zs.chain_role == role].groupby("pdb_id").isolated_spearman.first()
    full = (fp.mean(), *boot_ci(fp.values), len(fp))
    prm = (pr.mean(), *boot_ci(pr.values), len(pr))
    isl = (iso.mean(), *boot_ci(iso.values), len(iso))
    groups.append((lab, full, prm, isl))

# paired Wilcoxon (for the statistics report / caption)
for role, lab in ROLES:
    fp = pa[pa.chain_role == role].groupby("pdb_id").full_spearman.first()
    pr = pa[pa.chain_role == role].groupby("pdb_id").partner_removed_spearman.first()
    iso = zs[zs.chain_role == role].groupby("pdb_id").isolated_spearman.first()
    idx = fp.index.intersection(pr.index)
    d1 = (fp.loc[idx] - pr.loc[idx]).dropna().values
    w1 = stats.wilcoxon(d1, alternative="greater")
    idx2 = fp.index.intersection(iso.index)
    d2 = (fp.loc[idx2] - iso.loc[idx2]).dropna().values
    w2 = stats.wilcoxon(d2, alternative="greater")
    print(f"{role:8s} full-pr_mean_d={d1.mean():+.4f} p={w1.pvalue:.3e} | "
          f"full-iso_mean_d={d2.mean():+.4f} p={w2.pvalue:.3e}")

# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------
mm = 1 / 25.4
fig = plt.figure(figsize=(183 * mm, 78 * mm))

axa = fig.add_axes([0.07, 0.24, 0.86, 0.62])
x = np.arange(len(ROLES))
width = 0.20
offs = [-0.30, 0.0, 0.30]
colors = [ACCENT_RED, "#2166AC", GREY]
cond_names = ["full complex", "partner-removed (pMHC / TCRab dimer)",
              "isolated chain"]

for j, (off, col) in enumerate(zip(offs, colors)):
    vals = [g[j + 1][0] for g in groups]
    los = [g[j + 1][0] - g[j + 1][1] for g in groups]
    his = [g[j + 1][2] - g[j + 1][0] for g in groups]
    axa.bar(x + off, vals, width=width, color=col,
            yerr=[los, his], error_kw=dict(ecolor="#555555", lw=0.7,
                                           capsize=2.2),
            zorder=3)

axa.set_xticks(x)
axa.set_xticklabels([lab for _, lab in ROLES], fontsize=7)
axa.set_xlim(-0.55, len(ROLES) - 0.45)
axa.set_ylim(0.0, 0.95)
axa.set_yticks(np.arange(0.0, 1.0, 0.2))
axa.set_ylabel("Spearman $\\rho$ vs MD RMSF", fontsize=7)
axa.set_title("Partner context most affects peptide RMSF prediction",
              fontsize=8, loc="left", pad=4)
axa.legend(handles=[
    Line2D([0], [0], marker="s", color="none", markerfacecolor=ACCENT_RED,
           ms=5.0, label=cond_names[0]),
    Line2D([0], [0], marker="s", color="none", markerfacecolor="#2166AC",
           ms=5.0, label=cond_names[1]),
    Line2D([0], [0], marker="s", color="none", markerfacecolor=GREY,
           ms=5.0, label=cond_names[2]),
], loc="upper center", bbox_to_anchor=(0.5, 1.18), ncol=3, fontsize=6,
   frameon=False, handlelength=1.0, columnspacing=1.0)

save_cns_figure(fig, os.path.join(OUT, "fig_flexppi_rmsf_context_ablation"))
print("saved", os.path.join(OUT, "fig_flexppi_rmsf_context_ablation"))
