"""Paper figure 5 (biological interpretation and capability boundaries) draft.

Panels
------
a) Joint MD distribution of RMSF and Hchi (out-of-fold labels) with medians.
b) Region-resolved phenotype composition: fractions of the four median-split
   (low/high RMSF) x (low/high Hchi) phenotypes across sequence regions.
c) KIX partner-response boundary for Hchi (interface versus non-interface).
d) KIX partner-response boundary for RMSF (same layout).

Panels c/d are honest boundary results: interface-local Hchi response is stable,
RMSF response is not.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OOF = ROOT / "artifacts/evaluation_v2/analysis/oof_per_residue.csv.gz"
CONTEXT = ROOT / "output/itrace_residue_biology/analysis/tcr_residue_context.csv"
KIX = ROOT / "artifacts/external/kix_md_20260914/frame_evaluation/region_frame_delta_summary.csv"
OUTPUT = ROOT / "output/flexppi_paper_figs_v3"
QUADRANT_COLORS = {"low RMSF / low Hchi": "#9ecae1", "low RMSF / high Hchi": "#fdd0a2",
                   "high RMSF / low Hchi": "#a1d99b", "high RMSF / high Hchi": "#fc9272"}


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 7.8))
    axis_a, axis_b, axis_c, axis_d = axes.ravel()

    oof = pd.read_csv(OOF, usecols=["dataset", "target_rmsf", "target_chi_entropy",
                                    "valid_rmsf", "valid_chi_entropy"])
    valid = oof[(oof.valid_rmsf) & (oof.valid_chi_entropy)]
    hexbin = axis_a.hexbin(valid["target_rmsf"], valid["target_chi_entropy"],
                           gridsize=45, bins="log", cmap="viridis", mincnt=1)
    axis_a.axvline(valid["target_rmsf"].median(), color="white", ls="--", lw=0.8)
    axis_a.axhline(valid["target_chi_entropy"].median(), color="white", ls="--", lw=0.8)
    from flexppi.evaluation import _spearman
    rho = _spearman(valid["target_rmsf"].to_numpy(float),
                    valid["target_chi_entropy"].to_numpy(float))
    axis_a.text(0.03, 0.95, "n={:,}\nSpearman={:.2f}".format(len(valid), rho),
                transform=axis_a.transAxes, fontsize=8, va="top")
    fig.colorbar(hexbin, ax=axis_a, label="log10 residue count")
    axis_a.set_xlabel("MD RMSF (A)")
    axis_a.set_ylabel("MD Hchi (nats)")
    axis_a.set_title("a  two targets are not redundant", fontsize=9, loc="left")

    context = pd.read_csv(CONTEXT, usecols=["chain", "full_rmsf", "full_chi_entropy",
                                            "chi_defined", "region_group"])
    context = context[context.chi_defined.fillna(False)]
    medians = context.groupby("chain")[["full_rmsf", "full_chi_entropy"]].transform("median")
    flag_rmsf = np.where(context.full_rmsf >= medians.full_rmsf, "high RMSF", "low RMSF")
    flag_hchi = np.where(context.full_chi_entropy >= medians.full_chi_entropy,
                         "high Hchi", "low Hchi")
    context = context.assign(
        quadrant=pd.Series(flag_rmsf, index=context.index).astype(object)
        + " / " + pd.Series(flag_hchi, index=context.index).astype(object))
    groups = context.groupby("region_group")
    keep = [name for name, group in groups if len(group) >= 200]
    keep = sorted(keep, key=lambda name: -len(groups.get_group(name)))
    fractions = []
    for name in keep:
        group = groups.get_group(name)
        counts = group["quadrant"].value_counts(normalize=True)
        fractions.append([counts.get(quadrant, 0.0) for quadrant in QUADRANT_COLORS])
    fractions = np.asarray(fractions)
    bottom = np.zeros(len(keep))
    for index, (quadrant, color) in enumerate(QUADRANT_COLORS.items()):
        axis_b.bar(np.arange(len(keep)), fractions[:, index], bottom=bottom, color=color,
                   label=quadrant, alpha=0.9)
        bottom += fractions[:, index]
    axis_b.set_xticks(np.arange(len(keep)))
    axis_b.set_xticklabels(["{}\n(n={})".format(name, len(groups.get_group(name)))
                            for name in keep], fontsize=7, rotation=20, ha="right")
    axis_b.set_ylabel("fraction of residues")
    axis_b.set_ylim(0, 1)
    axis_b.legend(fontsize=6.5, ncol=2, loc="lower right")
    axis_b.set_title("b  regional phenotype composition (median split within chain)", fontsize=9,
                     loc="left")

    kix = pd.read_csv(KIX)
    states = [("myb-apo", "cMyb"), ("mll-apo", "MLL"), ("both-apo", "both")]
    interface_regions = {"myb-apo": "myb_interface", "mll-apo": "mll_interface",
                         "both-apo": "both_interface"}
    for metric, axis, title in (("hchi", axis_c, "c  KIX Hchi response boundary"),
                                ("rmsf", axis_d, "d  KIX RMSF response boundary")):
        positions = np.arange(len(states))
        width = 0.35
        interface_values, interface_sd, distal_values, distal_sd = [], [], [], []
        for contrast, _ in states:
            region = interface_regions[contrast]
            row = kix[(kix.region == region) & (kix.contrast == contrast)
                      & (kix.metric == metric)]
            interface_values.append(float(row["spearman_mean"].iloc[0]))
            interface_sd.append(float(row["spearman_sd"].iloc[0]))
            row = kix[(kix.region == "non_interface") & (kix.contrast == contrast)
                      & (kix.metric == metric)]
            distal_values.append(float(row["spearman_mean"].iloc[0]))
            distal_sd.append(float(row["spearman_sd"].iloc[0]))
        axis.axhline(0, color="grey", lw=0.8)
        axis.bar(positions - width / 2, interface_values, width, yerr=interface_sd,
                 color="#D55E00", alpha=0.85, label="bound-partner interface", capsize=3)
        axis.bar(positions + width / 2, distal_values, width, yerr=distal_sd,
                 color="#8A8F98", alpha=0.85, label="non-interface", capsize=3)
        axis.set_xticks(positions)
        axis.set_xticklabels([label for _, label in states], fontsize=8)
        axis.set_ylabel("delta Spearman (frame mean +/- frame SD)")
        axis.set_title(title, fontsize=9, loc="left")
        axis.legend(fontsize=7, loc="lower right")

    fig.suptitle("Biological interpretation and boundaries: two targets, regional phenotypes, "
                 "and where response prediction works or fails", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    for suffix in ("png", "pdf"):
        fig.savefig(OUTPUT / "fig5_biology_boundaries.{}".format(suffix), dpi=300)
    plt.close(fig)
    print("wrote", OUTPUT / "fig5_biology_boundaries.png")


if __name__ == "__main__":
    main()
