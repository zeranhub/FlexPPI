# TCR loop-region input-context ablation for RMSF and H_chi.
# For TCR alpha/beta chains, residues are split into CDR1/2/3 and framework
# (alpha/beta specific); per-chain Spearman vs MD is computed under three
# input contexts (full complex / partner-removed TCRab dimer / isolated
# chain) and summarized as equal-complex means with complex bootstrap 95% CI.
# Same grouped-bar visual language as the other context figures.

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

KEYS = ["pdb_id", "chain", "residue_number", "insertion_code"]
ALL_COND = r"E:\Project\FlexPPI\output\itrace_residue_biology\raw\all_condition_predictions.csv"
HCHI_CSV = r"E:\Project\FlexPPI\output\itrace_hchi_benchmark\raw\hchi_consensus_per_residue.csv"
CC_CSV = r"E:\Project\FlexPPI\output\itrace_zero_shot_context\raw\common_context_predictions.csv"

df = pd.read_csv(ALL_COND, usecols=KEYS + ["chain_role", "condition", "rmsf", "chi_entropy"],
                 low_memory=False)
h = pd.read_csv(HCHI_CSV, usecols=KEYS + ["region", "target_chi_entropy"], low_memory=False)
cc = pd.read_csv(CC_CSV, usecols=KEYS + ["target_chainfit_rmsf_angstrom"], low_memory=False)
df = df.merge(h, on=KEYS, how="left").merge(cc, on=KEYS, how="left")

iso_map = (df[df.condition.str.startswith("isolated")]
           .groupby(["pdb_id", "chain"]).condition.first().reset_index())
iso_map["condition"] = "isolated"

tcr = df[df.chain_role.isin(["tcr_alpha", "tcr_beta"])].copy()
full = tcr[tcr.condition == "full"].copy()
full["condition"] = "full"
pr = tcr[tcr.condition == "tcr"].copy()
pr["condition"] = "partner_removed"
isol = tcr[tcr.condition.str.startswith("isolated")].copy()
isol = isol.drop(columns=["condition"]).merge(iso_map, on=["pdb_id", "chain"])
master = pd.concat([full, pr, isol], ignore_index=True)

REGIONS = ["CDR1α", "CDR2α", "CDR3α", "TCR-FRα",
           "CDR1β", "CDR2β", "CDR3β", "TCR-FRβ"]
CONDITIONS = ["full", "partner_removed", "isolated"]
rng = np.random.default_rng(2026)

def region_stats(region, condition, descriptor):
    sub = master[(master.region == region) & (master.condition == condition)]
    target_col = "target_chainfit_rmsf_angstrom" if descriptor == "RMSF" else "target_chi_entropy"
    pred_col = "rmsf" if descriptor == "RMSF" else "chi_entropy"
    rows = []
    for (pid, ch), g in sub.groupby(["pdb_id", "chain"]):
        gg = g[[pred_col, target_col]].dropna()
        if len(gg) >= 4:
            sp, _ = stats.spearmanr(gg[pred_col], gg[target_col])
            rows.append((pid, ch, sp))
    if not rows:
        return None
    per = pd.DataFrame(rows, columns=["pdb_id", "chain", "sp"]).groupby("pdb_id").sp.first()
    boot = np.array([rng.choice(per, size=len(per), replace=True).mean()
                     for _ in range(2000)])
    return (per.mean(), np.percentile(boot, 2.5), np.percentile(boot, 97.5), len(per))

results = {}
for desc in ["RMSF", "Hchi"]:
    for region in REGIONS:
        for cond in CONDITIONS:
            s = region_stats(region, cond, desc)
            results[(desc, region, cond)] = s
            if s is None:
                print(f"{desc} {region} {cond}: NO DATA")

for desc in ["RMSF", "Hchi"]:
    print(f"\n===== {desc} =====")
    print("region       full       partner_removed  isolated")
    for region in REGIONS:
        cells = []
        for cond in CONDITIONS:
            s = results[(desc, region, cond)]
            cells.append(f"{s[0]:.3f} [{s[1]:.3f},{s[2]:.3f}]" if s else " n/a")
        print(f"{region:10s} {cells[0]:22s} {cells[1]:22s} {cells[2]}")

# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
mm = 1 / 25.4
colors = [ACCENT_RED, "#2166AC", GREY]
cond_names = ["full complex", "partner-removed (TCRab dimer)", "isolated chain"]
GROUP_LABELS = [(0, 3, "TCR$\\alpha$"), (4, 7, "TCR$\\beta$")]

def draw_grouped(desc, ylim, ylabel, title, fname):
    fig = plt.figure(figsize=(183 * mm, 80 * mm))
    ax = fig.add_axes([0.06, 0.24, 0.90, 0.62])
    x = np.arange(len(REGIONS))
    width = 0.20
    offs = [-0.30, 0.0, 0.30]
    for j, (off, col) in enumerate(zip(offs, colors)):
        vals, los, his = [], [], []
        for region in REGIONS:
            s = results[(desc, region, CONDITIONS[j])]
            vals.append(s[0]); los.append(s[0] - s[1]); his.append(s[2] - s[0])
        ax.bar(x + off, vals, width=width, color=col, yerr=[los, his],
               error_kw=dict(ecolor="#555555", lw=0.7, capsize=2.2), zorder=3)
    # group separators and labels
    ax.axvline(3.55, color="#CCCCCC", linewidth=0.7, zorder=1)
    for fx, lab in [(0.26, "TCR$\\alpha$"), (0.74, "TCR$\\beta$")]:
        ax.text(fx, -0.16, lab, ha="center", va="top", fontsize=7, color=BLACK,
                transform=ax.transAxes)
    ax.set_xticks(x)
    ax.set_xticklabels([r.replace("TCR-FR", "FR") for r in REGIONS], fontsize=6.5)
    ax.set_xlim(-0.75, len(REGIONS) - 0.25)
    ax.set_ylim(*ylim)
    ax.set_ylabel(ylabel, fontsize=7)
    ax.set_title(title, fontsize=8, loc="left", pad=4)
    ax.legend(handles=[
        Line2D([0], [0], marker="s", color="none", markerfacecolor=c, ms=5.0,
               label=n) for c, n in zip(colors, cond_names)
    ], loc="upper center", bbox_to_anchor=(0.5, 1.18), ncol=3, fontsize=6,
       frameon=False, handlelength=1.0, columnspacing=1.0)
    save_cns_figure(fig, os.path.join(OUT, fname))
    print("saved", os.path.join(OUT, fname))

# y-limits depend on observed values (choose after the summary print)
draw_grouped("RMSF", (0.70, 0.97),
             "Spearman $\\rho$ vs MD RMSF",
             "TCR loop-region RMSF prediction under different input contexts",
             "fig_flexppi_tcr_region_rmsf_context")
draw_grouped("Hchi", (0.70, 0.97),
             "Spearman $\\rho$ vs MD $H_\\chi$",
             "TCR loop-region $H_\\chi$ prediction under different input contexts",
             "fig_flexppi_tcr_region_hchi_context")
