"""Paper figure 4 (independent generalization and external validation) draft.

Panels
------
a) Chen 2026 test-chain similarity to training chains (max identity per query chain).
b) Chen 2026 per-chain RMSF Spearman and top-20% AUROC for WT and S32M.
c) ITrace external benchmark (pretrained checkpoint, 244 dense-label complexes):
   per-region Spearman for FlexPPI versus FastProtFlex.
d) Cross-dataset RMSF Spearman summary with available baselines.

This figure deliberately mixes one sparse-reference external set (Chen 2026) with
one dense-reference external set (ITrace benchmark) and a non-TCR PPI set
(Zenodo 6638504); captions must state those differences.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
CHEN = ROOT / "artifacts/chen2026_tcr156"
ITRACE = ROOT / "output/itrace_benchmark/metrics/summary_metrics.csv"
ZENODO = ROOT / "artifacts/zenodo_6638504/evaluation/summary_metrics.csv"
OUTPUT = ROOT / "output/flexppi_paper_figs_v3"
COLORS = {"chen_wt": "#0072B2", "chen_s32m": "#D55E00", "flexppi": "#D55E00",
          "baseline": "#8A8F98", "chance": "#999999"}


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    chen_report = json.loads((CHEN / "report.json").read_text())
    overlap = json.loads((CHEN / "overlap_check.json").read_text())
    itrace = pd.read_csv(ITRACE)
    zenodo = pd.read_csv(ZENODO)

    fig, axes = plt.subplots(2, 2, figsize=(11.5, 7.6))
    axis_a, axis_b, axis_c, axis_d = axes.ravel()

    chains = ["A", "C", "D", "E"]
    chain_labels = {"A": "HLA\n(chain A)", "C": "peptide\n(chain C)",
                    "D": "TCR alpha\n(chain D)", "E": "TCR beta\n(chain E)"}
    best = []
    best_hit = []
    for chain in chains:
        hits = overlap["queries"][chain]["top_hits"]
        full = [hit for hit in hits if hit["coverage_query"] >= 0.9]
        if full:
            top = max(full, key=lambda hit: hit["identity"])
            best.append(top["identity"])
            best_hit.append("{} ({:.0%} q)".format(top["target"], top["coverage_query"]))
        else:
            best.append(0.0)
            best_hit.append("no full-length hit")
    axis_a.bar(np.arange(len(chains)), best, color=COLORS["baseline"], alpha=0.8)
    for position, (value, note) in enumerate(zip(best, best_hit)):
        axis_a.text(position, min(value + 0.03, 1.02), "{:.2f}".format(value),
                    ha="center", fontsize=8)
        axis_a.text(position, 0.02, note, ha="center", fontsize=6.5, rotation=90, va="bottom")
    rank = np.arange(len(chains))
    axis_a.axhline(0.30, color=COLORS["chance"], ls="--", lw=0.9)
    axis_a.text(0.02, 0.32, "family threshold (30%)", transform=axis_a.transAxes, fontsize=7)
    axis_a.set_xticks(rank)
    axis_a.set_xticklabels([chain_labels[chain] for chain in chains], fontsize=8)
    axis_a.set_ylabel("max identity to any training chain")
    axis_a.set_ylim(0, 1.05)
    axis_a.set_title("a  Chen 2026 novelty: combination-level, not family-level", fontsize=9,
                     loc="left")

    conditions = [("WT", "chen_wt"), ("S32M", "chen_s32m")]
    width = 0.35
    for index, (condition, color_key) in enumerate(conditions):
        block = chen_report["conditions"][condition]
        values = [block["chains"][chain]["rmsf"]["spearman"] for chain in chains]
        axis_b.bar(np.arange(len(chains)) + (index - 0.5) * width, values, width=width,
                   color=COLORS[color_key], alpha=0.85, label=condition)
        aurocs = [block["chains"][chain]["rmsf"].get("top20_auroc") for chain in chains]
        for position, value, auroc in zip(np.arange(len(chains)), values, aurocs):
            if auroc:
                axis_b.text(position + (index - 0.5) * width, 0.02,
                            "AUROC {:.2f}".format(auroc), rotation=90, fontsize=6,
                            ha="center", va="bottom")
    axis_b.axhline(0, color="grey", lw=0.6)
    axis_b.set_xticks(np.arange(len(chains)))
    axis_b.set_xticklabels([chain_labels[chain] for chain in chains], fontsize=8)
    axis_b.set_ylabel("RMSF Spearman (5-fold mean)")
    axis_b.set_ylim(0, 1)
    axis_b.legend(fontsize=8)
    axis_b.set_title("b  Chen 2026 external clonotype (10 ns-frame reference)", fontsize=9,
                     loc="left")

    regions = itrace[itrace.label_definition == "chain_fitted_internal_rmsf"]
    keep = ["tcr_alpha", "tcr_beta", "cdr3", "cdr_all", "mhc", "b2m", "pmhc", "overall"]
    positions = np.arange(len(keep))
    for index, (method, color, label) in enumerate(
            (("flexppi", COLORS["flexppi"], "FlexPPI"), ("fastprotflex", COLORS["baseline"],
                                                         "FastProtFlex"))):
        values = []
        for subset in keep:
            row = regions[(regions.method == method) & (regions.subset == subset)
                          & (regions.metric == "spearman")]
            values.append(float(row["macro_mean"].iloc[0]) if len(row) else np.nan)
        axis_c.bar(positions + (index - 0.5) * 0.35, values, width=0.35, color=color,
                   alpha=0.85, label=label)
    axis_c.set_xticks(positions)
    axis_c.set_xticklabels(keep, rotation=30, ha="right", fontsize=8)
    axis_c.set_ylabel("Spearman (equal-complex macro)")
    axis_c.set_ylim(0, 1)
    axis_c.legend(fontsize=8)
    axis_c.set_title("c  ITrace external benchmark, 244 dense-label complexes\n"
                     "(checkpoint not trained on ITrace)", fontsize=9, loc="left")

    rows = []
    block = chen_report["conditions"]["WT"]["chains"]
    rows.append({"dataset": "Chen 2026 WT (sparse)", "spearman": np.mean(
        [block[chain]["rmsf"]["spearman"] for chain in chains]), "note": "4 chains, 1 system"})
    row = regions[(regions.method == "flexppi") & (regions.subset == "overall")
                  & (regions.metric == "spearman")]
    rows.append({"dataset": "ITrace external (dense)", "spearman": float(row["macro_mean"].iloc[0]),
                 "note": "244 complexes"})
    row = zenodo[(zenodo.subset == "all_8") & (zenodo.task == "apo_state")
                 & (zenodo.metric == "rmsf") & (zenodo.method == "flexppi")
                 & (zenodo.measure == "spearman")]
    rows.append({"dataset": "Zenodo 6638504 (dense)", "spearman": float(row["macro_mean"].iloc[0]),
                 "note": "8 non-TCR complexes"})
    row = zenodo[(zenodo.subset == "all_8") & (zenodo.task == "apo_state")
                 & (zenodo.metric == "rmsf") & (zenodo.method == "gnm")
                 & (zenodo.measure == "spearman")]
    rows.append({"dataset": "GNM baseline", "spearman": float(row["macro_mean"].iloc[0]),
                 "note": "physics baseline"})
    summary = pd.DataFrame(rows)
    axis_d.bar(np.arange(len(summary)), summary["spearman"], color=COLORS["flexppi"], alpha=0.8)
    axis_d.set_xticks(np.arange(len(summary)))
    axis_d.set_xticklabels(["{}\n{}".format(row["dataset"], row["note"])
                            for _, row in summary.iterrows()], fontsize=7.5)
    axis_d.set_ylabel("RMSF Spearman")
    axis_d.set_ylim(0, 1)
    axis_d.set_title("d  external summary: sparse spike, dense ITrace, non-TCR PPI", fontsize=9,
                     loc="left")

    fig.suptitle("FlexPPI independent generalization: what the external tests do and do not show",
                 fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    for suffix in ("png", "pdf"):
        fig.savefig(OUTPUT / "fig4_independent_generalization.{}".format(suffix), dpi=300)
    plt.close(fig)
    print("wrote", OUTPUT / "fig4_independent_generalization.png")


if __name__ == "__main__":
    main()
