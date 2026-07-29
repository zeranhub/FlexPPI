"""Render a compact summary figure for the paired model experiments."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


COLORS = {
    "fastprotflex": "#4C78A8",
    "flexppi_pretrain": "#F58518",
    "flexppi_finetune": "#54A24B",
}
LABELS = {
    "fastprotflex": "FastProtFlex",
    "flexppi_pretrain": "FlexPPI static",
    "flexppi_finetune": "FlexPPI response-FT",
}


def _load(path: str):
    with Path(path).open("r", encoding="utf-8") as handle:
        return json.load(handle)


def _bar_with_ci(axis, x, value, ci, color, label=None, width=0.22):
    lower = value - ci[0]
    upper = ci[1] - value
    axis.bar(x, value, width=width, color=color, label=label, zorder=2)
    axis.errorbar(
        x, value, yerr=np.asarray([[lower], [upper]]), fmt="none",
        ecolor="#222222", capsize=2.5, linewidth=0.9, zorder=3,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--atlas",
        default="artifacts/evaluation/fastprotflex_comparison/atlas_comparison_report.json",
    )
    parser.add_argument(
        "--dynarepo",
        default=(
            "artifacts/evaluation/fastprotflex_comparison/dynarepo/"
            "dynarepo_comparison_report.json"
        ),
    )
    parser.add_argument(
        "--tcr3d",
        default=(
            "artifacts/evaluation/fastprotflex_comparison/tcr3d/"
            "tcr3d_comparison_report.json"
        ),
    )
    parser.add_argument(
        "--output",
        default="artifacts/evaluation/fastprotflex_comparison/comparison_summary.png",
    )
    arguments = parser.parse_args()

    atlas, dynarepo, tcr3d = (
        _load(arguments.atlas), _load(arguments.dynarepo), _load(arguments.tcr3d)
    )
    figure, axes = plt.subplots(1, 3, figsize=(13.5, 4.25), dpi=180)

    methods = ("fastprotflex", "flexppi_pretrain", "flexppi_finetune")
    for dataset_index, (dataset_label, report) in enumerate(
        (("ATLAS", atlas), ("DynaRepo", dynarepo))
    ):
        for method_index, method in enumerate(methods):
            metrics = report["methods"][method]["spearman"]
            x = dataset_index + (method_index - 1) * 0.23
            _bar_with_ci(
                axes[0], x, metrics["macro_mean"], metrics["protein_bootstrap_95_ci"],
                COLORS[method], LABELS[method] if dataset_index == 0 else None,
            )
    axes[0].set_xticks([0, 1], ["ATLAS\n(same-domain)", "DynaRepo\n(adaptation)"])
    axes[0].set_ylabel("Protein-macro Spearman")
    axes[0].set_ylim(0, 0.9)
    axes[0].set_title("A  Shared RMSF ranking")
    axes[0].legend(frameon=False, fontsize=8, loc="lower left")

    proxy_specs = (
        ("B-factor\nSpearman", "bfactor_rank_proxy", "spearman"),
        ("CDR3\nAUROC", "cdr3_vs_framework_retrieval", "auroc"),
        ("CDR3\nAP", "cdr3_vs_framework_retrieval", "average_precision"),
    )
    proxy_methods = ("fastprotflex", "flexppi_finetune")
    for proxy_index, (label, task, metric) in enumerate(proxy_specs):
        for method_index, method in enumerate(proxy_methods):
            metrics = tcr3d[task]["methods"][method][metric]
            x = proxy_index + (method_index - 0.5) * 0.25
            _bar_with_ci(
                axes[1], x, metrics["macro_mean"],
                metrics.get("protein_bootstrap_95_ci") or metrics["chain_bootstrap_95_ci"],
                COLORS[method], LABELS[method] if proxy_index == 0 else None,
                width=0.24,
            )
    axes[1].set_xticks(range(3), [item[0] for item in proxy_specs])
    axes[1].set_ylim(0, 0.7)
    axes[1].set_ylabel("Chain-macro score")
    axes[1].set_title("B  TCR3D external proxies")
    axes[1].legend(frameon=False, fontsize=8, loc="lower left")

    medians = [atlas["runtime"][method]["median_seconds"] for method in methods]
    positions = np.arange(len(methods))
    axes[2].bar(
        positions, medians, color=[COLORS[method] for method in methods], width=0.65
    )
    axes[2].set_xticks(positions, ["Fast", "Flex\nstatic", "Flex\nresponse-FT"])
    axes[2].set_ylabel("Median seconds / protein")
    axes[2].set_title("C  ATLAS CPU end-to-end")
    for position, value in zip(positions, medians):
        axes[2].text(position, value + 0.015, "{:.3f}".format(value), ha="center", fontsize=8)
    axes[2].set_ylim(0, max(medians) * 1.22)

    for axis in axes:
        axis.spines[["top", "right"]].set_visible(False)
        axis.grid(axis="y", alpha=0.18, linewidth=0.7, zorder=0)
    figure.suptitle(
        "FlexPPI vs FastProtFlex: paired common-task evaluation",
        fontsize=13, fontweight="bold", y=1.02,
    )
    figure.tight_layout()
    output = Path(arguments.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, bbox_inches="tight")
    print(output.resolve())


if __name__ == "__main__":
    main()
