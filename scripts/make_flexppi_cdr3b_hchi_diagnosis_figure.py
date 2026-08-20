# Mirror of the CDR3b RMSF diagnosis figure, evaluated on H_chi:
# (a) paired per-complex Spearman (full vs TCRab dimer)
# (b) within-chain prediction std: full input preserves H_chi spread
# (c) position profile along CDR3b: full input tracks MD
# (d) amino-acid-type bias: full vs dimer

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
ALL = r"E:\Project\FlexPPI\output\itrace_residue_biology\raw\all_condition_predictions.csv"
HCHI = r"E:\Project\FlexPPI\output\itrace_hchi_benchmark\raw\hchi_consensus_per_residue.csv"

df = pd.read_csv(ALL, usecols=KEYS + ["chain_role", "condition", "amino_acid_3", "chi_entropy"],
                 low_memory=False)
h = pd.read_csv(HCHI, usecols=KEYS + ["region_group", "target_chi_entropy",
                                      "is_tcr_pmhc_interface_8A"], low_memory=False)
df = df.merge(h, on=KEYS, how="left")

iso_map = df[df.condition.str.startswith("isolated")].groupby(["pdb_id", "chain"]).condition.first().reset_index()
iso_map["condition"] = "isolated"
tcr = df[df.chain_role == "tcr_beta"].copy()
full = tcr[tcr.condition == "full"].copy(); full["condition"] = "full"
pr = tcr[tcr.condition == "tcr"].copy(); pr["condition"] = "partner_removed"
isol = tcr[tcr.condition.str.startswith("isolated")].copy().drop(columns=["condition"]).merge(iso_map, on=["pdb_id", "chain"])
m = pd.concat([full, pr, isol], ignore_index=True)
cdr3 = m[(m.region_group == "cdr3") & m.target_chi_entropy.notna() & m.chi_entropy.notna()].copy()
cdr3["bias"] = cdr3.chi_entropy - cdr3.target_chi_entropy

# --- per-complex-chain spearman (a) and pred std (b) ---
sp_rows, std_rows = [], []
for cond in ["full", "partner_removed", "isolated"]:
    sub = cdr3[cdr3.condition == cond]
    for (pid, ch), g in sub.groupby(["pdb_id", "chain"]):
        gg = g[["chi_entropy", "target_chi_entropy"]].dropna()
        if len(gg) >= 4:
            sp, _ = stats.spearmanr(gg.chi_entropy, gg.target_chi_entropy)
            sp_rows.append((pid, ch, cond, sp))
        std_rows.append((pid, ch, cond, g.chi_entropy.std(), g.target_chi_entropy.std()))
per = pd.DataFrame(sp_rows, columns=["pdb_id", "chain", "condition", "sp"])
piv = per.pivot_table(index=["pdb_id", "chain"], columns="condition", values="sp")
st = pd.DataFrame(std_rows, columns=["pdb_id", "chain", "condition", "pred_std", "md_std"])

d = piv["full"] - piv["partner_removed"]
frac_above = (d > 0).mean()
w = stats.wilcoxon(d.dropna(), alternative="greater")
md_std_mean = st.groupby(["pdb_id", "chain"]).md_std.first().mean()

# --- position profile (c) ---
cdr3 = cdr3.sort_values(["pdb_id", "chain", "residue_number"])
cdr3["pos_rank"] = cdr3.groupby(["pdb_id", "chain"]).cumcount()
cdr3["loop_len"] = cdr3.groupby(["pdb_id", "chain"])["pos_rank"].transform("max")
cdr3["rel_pos"] = cdr3.pos_rank / cdr3.loop_len
cdr3["pos_bin"] = pd.cut(cdr3.rel_pos, bins=np.linspace(0, 1, 6),
                         labels=[f"{i}" for i in range(5)], include_lowest=True)
pos_md = cdr3.groupby("pos_bin").target_chi_entropy.mean()
pos_pred = cdr3.groupby(["pos_bin", "condition"]).chi_entropy.mean().unstack()
pos_contact = cdr3.groupby("pos_bin").is_tcr_pmhc_interface_8A.mean()

# --- amino-acid-type bias (d), exclude PRO (no rotamer entropy) ---
aa = cdr3[cdr3.amino_acid_3 != "PRO"].groupby(["amino_acid_3", "condition"]).agg(
    n=("bias", "size"), bias=("bias", "mean")).reset_index()
aa_n = aa.pivot_table(index="amino_acid_3", columns="condition", values="n")
aa_b = aa.pivot_table(index="amino_acid_3", columns="condition", values="bias")
aa = aa_b[aa_n.sum(axis=1) >= 50].dropna(subset=["full", "partner_removed"])
aa = aa.sort_values("full")

# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------
mm = 1 / 25.4
fig = plt.figure(figsize=(183 * mm, 142 * mm))

# ---- (a) paired Spearman ----------------------------------------------------
axa = fig.add_axes([0.065, 0.555, 0.38, 0.38])
axa_top = fig.add_axes([0.065, 0.955, 0.38, 0.035], sharex=axa)
axa_right = fig.add_axes([0.455, 0.555, 0.045, 0.38], sharey=axa)
axa_top.axis("off"); axa_right.axis("off")
x = piv["partner_removed"].values; y = piv["full"].values
axa.scatter(x, y, s=7, color="#4393C3", alpha=0.65, edgecolors="none", zorder=3)
lo, hi = 0.4, 1.02
axa.plot([lo, hi], [lo, hi], color=GREY, linestyle="--", linewidth=0.6, alpha=0.7, zorder=2)
axa.set_xlim(lo, hi); axa.set_ylim(lo, hi)
axa.set_xlabel("TCRab dimer: per-chain Spearman", fontsize=6.5)
axa.set_ylabel("Full complex: per-chain Spearman", fontsize=6.5)
axa_top.hist(x, bins=np.linspace(lo, hi, 18), color="#9ECAE1", edgecolor="white", linewidth=0.2)
axa_right.hist(y, bins=np.linspace(lo, hi, 18), color="#9ECAE1", orientation="horizontal", edgecolor="white", linewidth=0.2)
axa_top.set_ylim(0, axa_top.get_ylim()[1] * 1.4)
axa_right.set_xlim(0, axa_right.get_xlim()[1] * 1.4)
axa.text(0.03, 0.97,
         f"{frac_above:.0%} of complexes above diagonal\n"
         f"paired Wilcoxon p = {w.pvalue:.1e}",
         transform=axa.transAxes, va="top", ha="left", fontsize=6.5, color=BLACK)
axa.text(-0.15, 1.10, "a", transform=axa.transAxes, fontsize=9, fontweight="bold", va="top")

# ---- (b) prediction std ------------------------------------------------------
axb = fig.add_axes([0.58, 0.555, 0.36, 0.38])
order = ["full", "partner_removed", "isolated"]
colors_b = [ACCENT_RED, "#2166AC", GREY]
labels_b = ["full complex", "TCRab dimer", "isolated chain"]
for cond, col, lab in zip(order, colors_b, labels_b):
    vals = np.sort(st[st.condition == cond].pred_std.values)
    ecdf = np.arange(1, len(vals) + 1) / len(vals)
    axb.plot(vals, ecdf, color=col, linewidth=1.4 if cond == "full" else 0.9,
             label=lab)
axb.axvline(md_std_mean, color=BLACK, linestyle="--", linewidth=1.0, zorder=2)
axb.text(md_std_mean + 0.02, 0.08, f"MD std = {md_std_mean:.2f}",
         fontsize=6.5, color=BLACK, rotation=90, va="bottom")
axb.set_xlim(0.0, 1.8)
axb.set_ylim(0.0, 1.05)
axb.set_xlabel("Within-chain prediction std (nats)", fontsize=6.5)
axb.set_ylabel("Fraction of chains", fontsize=6.5)
axb.legend(loc="lower right", fontsize=6, frameon=False, handlelength=1.5)
axb.text(-0.14, 1.10, "b", transform=axb.transAxes, fontsize=9, fontweight="bold", va="top")

# ---- (c) position profile ----------------------------------------------------
axc = fig.add_axes([0.065, 0.10, 0.38, 0.38])
px = np.arange(5)
axc.plot(px, pos_md.values, color=BLACK, linestyle="--", linewidth=1.2,
         marker="o", ms=3.5, label="MD truth")
for cond, col, lab in zip(order, colors_b, labels_b):
    axc.plot(px, pos_pred[cond].values, color=col, linewidth=1.2, marker="s",
             ms=3.5, label=lab)
apex = np.where(pos_contact.values > 0.6)[0]
if len(apex):
    axc.axvspan(apex[0] - 0.5, apex[-1] + 0.5, color="#F0F0F0", zorder=0)
axc.text(2, 0.55, "apex (pMHC contact)", fontsize=6, color="#666666", ha="center")
axc.set_xticks(px)
axc.set_xticklabels(["N-term", "", "", "", "C-term"], fontsize=6.5)
axc.set_xlim(-0.5, 4.5)
axc.set_ylim(0.4, 2.1)
axc.set_xlabel("CDR3$\\beta$ loop position (N $\\rightarrow$ C)", fontsize=6.5)
axc.set_ylabel("Mean $H_\\chi$ (nats)", fontsize=6.5)
axc.legend(loc="upper left", fontsize=6, frameon=False, handlelength=1.5)
axc.text(-0.15, 1.10, "c", transform=axc.transAxes, fontsize=9, fontweight="bold", va="top")

# ---- (d) amino-acid-type bias -------------------------------------------------
axd = fig.add_axes([0.58, 0.10, 0.36, 0.38])
yy = np.arange(len(aa))[::-1]
for yi, (a, row) in zip(yy, aa.iterrows()):
    axd.barh(yi, row["full"], height=0.34, color="#F4A582", zorder=3)
    axd.barh(yi - 0.19, row["partner_removed"], height=0.20,
             color="#2166AC", alpha=0.75, zorder=3)
axd.axvline(0.0, color="#999999", linewidth=0.7, zorder=2)
axd.set_yticks(yy)
axd.set_yticklabels(aa.index, fontsize=6)
axd.set_xlim(-0.3, 0.6)
axd.set_xlabel("Mean bias pred $-$ MD (nats)", fontsize=6.5)
axd.legend(handles=[
    Line2D([0], [0], marker="s", color="none", markerfacecolor="#F4A582", ms=4.5, label="full complex"),
    Line2D([0], [0], marker="s", color="none", markerfacecolor="#2166AC", ms=4.5, label="TCRab dimer"),
], loc="lower left", fontsize=5.5, frameon=False, handlelength=1.0)
axd.text(-0.16, 1.10, "d", transform=axd.transAxes, fontsize=9, fontweight="bold", va="top")

save_cns_figure(fig, os.path.join(OUT, "fig_flexppi_cdr3b_hchi_diagnosis"))
print("saved", os.path.join(OUT, "fig_flexppi_cdr3b_hchi_diagnosis"))
print("frac above:", frac_above, "wilcoxon p:", w.pvalue, "md_std:", md_std_mean)
