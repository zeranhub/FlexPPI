"""Paper figure 3 (partner context ablation) draft.

Panels
------
a) Absolute RMSF Spearman per chain role and input context (full / partner removed /
   isolated target chain).
b) Paired full-minus-partner-removed Spearman difference per chain role with the
   existing bootstrap interval.
c) Same paired difference for Hchi.
d) Absolute Hchi Spearman per chain role and input context.

All numbers come from the frozen ITrace fixed-coordinate context run; this is an
input-environment test, not an experimental apo/holo response test.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
CONTEXT = ROOT / "artifacts/evaluation_v2/itrace_context"
OUTPUT = ROOT / "output/flexppi_paper_figs_v3"
ROLES = ("peptide", "tcr_alpha", "tcr_beta", "mhc", "b2m")
ROLE_LABELS = {"peptide": "peptide", "tcr_alpha": "TCR alpha", "tcr_beta": "TCR beta",
               "mhc": "MHC", "b2m": "b2m"}
CONTEXTS = ("full", "partner_removed", "isolated")
CONTEXT_COLORS = {"full": "#D55E00", "partner_removed": "#0072B2", "isolated": "#999999"}
CONTEXT_LABELS = {"full": "full complex", "partner_removed": "partner removed",
                  "isolated": "isolated chain"}


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    summary = pd.read_csv(CONTEXT / "context_summary.csv")
    paired = pd.read_csv(CONTEXT / "paired_context_differences.csv")
    spearman = summary[summary.statistic == "spearman"]
    paired_spearman = paired[paired.statistic == "spearman"]

    fig, axes = plt.subplots(2, 2, figsize=(11, 7.2))
    axis_a, axis_b, axis_c, axis_d = axes.ravel()

    width = 0.25
    for metric, axis, title in (("rmsf", axis_a, "a  absolute RMSF ranking"),
                                ("chi_entropy", axis_d, "d  absolute Hchi ranking")):
        block = spearman[spearman.metric == metric]
        for index, context in enumerate(CONTEXTS):
            values, errors = [], []
            for role in ROLES:
                row = block[(block.chain_role == role) & (block.context == context)]
                values.append(float(row["mean"].iloc[0]))
                errors.append(float(row["mean"].iloc[0] - row["ci_low"].iloc[0]))
            axis.bar(np.arange(len(ROLES)) + (index - 1) * width, values, width=width,
                     color=CONTEXT_COLORS[context], alpha=0.85, label=CONTEXT_LABELS[context])
            axis.errorbar(np.arange(len(ROLES)) + (index - 1) * width, values,
                          yerr=errors, fmt="none", ecolor="black", elinewidth=0.7, capsize=2)
        axis.set_xticks(np.arange(len(ROLES)))
        axis.set_xticklabels([ROLE_LABELS[role] for role in ROLES], fontsize=8)
        axis.set_ylim(0, 1)
        axis.set_ylabel("Spearman (equal-complex macro)")
        axis.set_title(title, fontsize=9, loc="left")
        axis.legend(fontsize=7, loc="lower right")

    for metric, axis, title in (("rmsf", axis_b, "b  benefit of partner context (RMSF)"),
                                ("chi_entropy", axis_c, "c  benefit of partner context (Hchi)")):
        block = paired_spearman[(paired_spearman.metric == metric)
                                & (paired_spearman.comparison == "full_minus_partner_removed")]
        values, lower, upper = [], [], []
        for role in ROLES:
            row = block[block.chain_role == role]
            values.append(float(row["mean_difference"].iloc[0]))
            lower.append(float(row["mean_difference"].iloc[0] - row["ci_low"].iloc[0]))
            upper.append(float(row["ci_high"].iloc[0] - row["mean_difference"].iloc[0]))
        axis.axhline(0, color="grey", lw=0.8)
        axis.errorbar(np.arange(len(ROLES)), values, yerr=[lower, upper], fmt="o",
                      color=CONTEXT_COLORS["full"], ecolor="black", elinewidth=0.8,
                      capsize=3, markersize=5)
        axis.set_xticks(np.arange(len(ROLES)))
        axis.set_xticklabels([ROLE_LABELS[role] for role in ROLES], fontsize=8)
        axis.set_ylabel("full minus partner-removed Spearman")
        axis.set_title(title, fontsize=9, loc="left")

    fig.suptitle("Partner context ablation on 244 ITrace complexes (fixed coordinates, 8 inputs)",
                 fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    for suffix in ("png", "pdf"):
        fig.savefig(OUTPUT / "fig3_partner_context.{}".format(suffix), dpi=300)
    plt.close(fig)
    print("wrote", OUTPUT / "fig3_partner_context.png")


if __name__ == "__main__":
    main()
