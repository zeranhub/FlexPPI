"""Plot single-chain and complex-context comparison results."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


COLORS = {
    "fastprotflex": "#4C78A8",
    "fastprotflex_chain": "#4C78A8",
    "fastprotflex_complex": "#72A0CF",
    "flexppi_pretrain": "#F58518",
    "flexppi_pretrain_complex": "#F58518",
    "flexppi_finetune": "#54A24B",
    "flexppi_finetune_complex": "#54A24B",
}


def _bar(axis, x, metric, color, label, width):
    mean = metric["macro_mean"]
    ci = metric["protein_bootstrap_95_ci"]
    axis.bar(x, mean, width=width, color=color, label=label, zorder=2)
    axis.errorbar(
        x,
        mean,
        yerr=np.asarray([[mean - ci[0]], [ci[1] - mean]]),
        fmt="none",
        ecolor="#222222",
        capsize=2.5,
        linewidth=0.9,
        zorder=3,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--report",
        default=(
            "artifacts/evaluation/fastprotflex_comparison/chain_context/"
            "chain_context_comparison_report.json"
        ),
    )
    parser.add_argument(
        "--output",
        default=(
            "artifacts/evaluation/fastprotflex_comparison/chain_context/"
            "chain_context_summary.png"
        ),
    )
    arguments = parser.parse_args()
    with Path(arguments.report).open("r", encoding="utf-8") as handle:
        report = json.load(handle)

    figure, axes = plt.subplots(1, 3, figsize=(14.2, 4.4), dpi=180)
    single_methods = ("fastprotflex", "flexppi_pretrain", "flexppi_finetune")
    single_labels = ("FastProtFlex", "FlexPPI static", "FlexPPI response-FT")
    for dataset_index, dataset_name in enumerate(
        ("atlas", "dynarepo_single_protein_chain")
    ):
        methods = report["single_chain"][dataset_name]["methods"]
        for method_index, (method, label) in enumerate(zip(single_methods, single_labels)):
            _bar(
                axes[0],
                dataset_index + (method_index - 1) * 0.23,
                methods[method]["spearman"],
                COLORS[method],
                label if dataset_index == 0 else None,
                0.22,
            )
    axes[0].set_xticks([0, 1], ["ATLAS\n191 systems", "DynaRepo single-chain\n31 systems"])
    axes[0].set_title("A  Single protein chain")
    axes[0].set_ylabel("System-macro Spearman")
    axes[0].legend(frameon=False, fontsize=8, loc="lower left")

    complex_methods = (
        "fastprotflex_chain",
        "fastprotflex_complex",
        "flexppi_pretrain_complex",
        "flexppi_finetune_complex",
    )
    complex_labels = (
        "Fast / chain",
        "Fast / full complex",
        "Flex static / full",
        "Flex response-FT / full",
    )
    overall = report["protein_complex"]["equal_complex_methods"]
    for index, (method, label) in enumerate(zip(complex_methods, complex_labels)):
        _bar(
            axes[1], index, overall[method]["spearman"], COLORS[method], label, 0.7
        )
    axes[1].set_xticks(
        range(4), ["Fast\nchain", "Fast\nfull", "Flex static\nfull", "Flex FT\nfull"]
    )
    axes[1].set_title("B  Protein complexes (54 systems)")
    axes[1].set_ylabel("Equal-complex Spearman")

    region_methods = (
        "fastprotflex_chain",
        "fastprotflex_complex",
        "flexppi_pretrain_complex",
        "flexppi_finetune_complex",
    )
    for region_index, region in enumerate(("interface", "non_interface")):
        methods = report["protein_complex"]["by_region"][region]["methods"]
        for method_index, method in enumerate(region_methods):
            _bar(
                axes[2],
                region_index + (method_index - 1.5) * 0.19,
                methods[method]["spearman"],
                COLORS[method],
                complex_labels[method_index] if region_index == 0 else None,
                0.18,
            )
    axes[2].set_xticks([0, 1], ["Interface\n2,486 residues", "Non-interface\n16,667 residues"])
    axes[2].set_title("C  Where does context help?")
    axes[2].set_ylabel("Equal-complex Spearman")
    axes[2].legend(frameon=False, fontsize=7.5, loc="lower right")

    for axis in axes:
        axis.set_ylim(0, 0.9)
        axis.spines[["top", "right"]].set_visible(False)
        axis.grid(axis="y", alpha=0.18, linewidth=0.7, zorder=0)
    figure.suptitle(
        "FlexPPI vs FastProtFlex: separating chain identity from complex context",
        fontsize=13,
        fontweight="bold",
        y=1.02,
    )
    figure.tight_layout()
    output = Path(arguments.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, bbox_inches="tight")
    print(output.resolve())


if __name__ == "__main__":
    main()
