"""Aggregate performance-first v2 five-fold test evaluations.

Inputs are the ``fold_N_test.json`` files produced by
``evaluate_performance_first_v2_epoch40.sh``.  Aggregation is macro over folds:
each fold contributes its residue-pooled and its equal-system metrics, and the
final number is the mean and standard deviation across folds.  Source-stratified
metrics are reported per source with the same fold-macro treatment.  The H-chi
amino-acid-corrected report uses training-fold priors by construction.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np

METRICS = ["rmsf", "chi_entropy"]
SCALAR_KEYS = ("mae", "rmse", "bias", "pearson", "spearman", "r2",
               "concordance_correlation", "nmae_by_target_iqr")


def _report_scalars(report: Dict[str, Any]) -> Dict[str, Optional[float]]:
    return {key: report.get(key) for key in SCALAR_KEYS}


def _fold_mean(values: List[Optional[float]]) -> Optional[Dict[str, Optional[float]]]:
    finite = [float(value) for value in values if value is not None]
    if not finite:
        return None
    return {"mean": float(np.mean(finite)), "std": float(np.std(finite)),
            "n_folds": len(finite)}


def _fold_mean_of_scalars(
    reports: List[Dict[str, Any]],
) -> Dict[str, Optional[Dict[str, Optional[float]]]]:
    output: Dict[str, Optional[Dict[str, Optional[float]]]] = {}
    for key in SCALAR_KEYS:
        output[key] = _fold_mean([report.get(key) for report in reports])
    return output


def _aggregate_metric(metric_reports: List[Dict[str, Any]]) -> Dict[str, Any]:
    aggregated: Dict[str, Any] = {
        "overall": {
            "fold_values": [_report_scalars(report["overall"]) for report in metric_reports],
            "fold_mean": _fold_mean_of_scalars(
                [report["overall"] for report in metric_reports]
            ),
            "residue_count_per_fold": [report["overall"].get("count", 0)
                                       for report in metric_reports],
        },
        "macro_by_system": _fold_mean_of_scalars(
            [report["overall"].get("macro_by_system", {}) for report in metric_reports]
        ),
    }
    sources: Dict[str, List[Dict[str, Any]]] = {}
    for report in metric_reports:
        for source, source_report in report.get("by_dataset", {}).items():
            sources.setdefault(source, []).append(source_report)
    aggregated["by_dataset"] = {
        source: {
            "residue_pooled": _fold_mean_of_scalars(reports),
            "macro_by_system": _fold_mean_of_scalars(
                [report.get("macro_by_system", {}) for report in reports]
            ),
            "system_count_per_fold": [
                report.get("macro_by_system", {}).get("correlation_system_count", 0)
                for report in reports
            ],
        }
        for source, reports in sorted(sources.items())
    }
    return aggregated


def _aggregate_amino_acid_correction(
    metric_reports: List[Dict[str, Any]],
) -> Dict[str, Any]:
    corrected = [
        report.get("amino_acid_corrected") for report in metric_reports
        if report.get("amino_acid_corrected")
    ]
    if not corrected:
        return {"available": False}
    output: Dict[str, Any] = {
        "available": True,
        "definition": corrected[0].get("definition"),
        "residual_overall": _fold_mean_of_scalars(
            [report["residual_overall"] for report in corrected]
        ),
        "residual_macro_by_system": _fold_mean_of_scalars(
            [report["residual_overall"].get("macro_by_system", {}) for report in corrected]
        ),
        "within_amino_acid_type": {},
    }
    amino_acids: Dict[str, List[Dict[str, Any]]] = {}
    for report in corrected:
        for amino_acid, type_report in report.get("within_amino_acid_type", {}).items():
            amino_acids.setdefault(amino_acid, []).append(type_report)
    for amino_acid, reports in sorted(amino_acids.items()):
        output["within_amino_acid_type"][amino_acid] = {
            "residue_pooled": _fold_mean_of_scalars(reports),
            "macro_by_system": _fold_mean_of_scalars(
                [report.get("macro_by_system", {}) for report in reports]
            ),
        }
    return output


def aggregate(evaluation_dir: Path, output_path: Path) -> Dict[str, Any]:
    fold_reports = []
    for fold in range(5):
        path = evaluation_dir / "fold_{}_test.json".format(fold)
        if not path.exists():
            raise FileNotFoundError("Missing fold evaluation: {}".format(path))
        fold_reports.append(json.loads(path.read_text(encoding="utf-8")))
    summary: Dict[str, Any] = {
        "version": "performance_first_v2_epoch40",
        "protocol": (
            "Each fold trains from scratch for a fixed 40-epoch budget; the "
            "epoch-40 checkpoint evaluates its own untouched test fold. No "
            "model selection used test results."
        ),
        "checkpoints": [
            report.get("checkpoint") for report in fold_reports
        ],
        "checkpoint_epochs": [
            report.get("checkpoint_epoch") for report in fold_reports
        ],
        "per_fold_split": [
            report.get("split", {}) for report in fold_reports
        ],
        "metrics": {},
    }
    for metric in METRICS:
        metric_reports = [report["metrics"][metric] for report in fold_reports]
        aggregated = _aggregate_metric(metric_reports)
        if metric == "chi_entropy":
            aggregated["amino_acid_corrected"] = _aggregate_amino_acid_correction(
                metric_reports
            )
        summary["metrics"][metric] = aggregated
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    return summary


def _render_markdown(summary: Dict[str, Any]) -> str:
    lines: List[str] = [
        "# FlexPPI performance-first v2 - epoch-40 test evaluation (5-fold macro)",
        "",
        "| metric | scope | mae | rmse | pearson | spearman |",
        "|---|---|---|---|---|---|",
    ]

    def row(metric: str, scope: str, stats: Dict[str, Any]) -> str:
        def cell(key: str) -> str:
            value = stats.get(key)
            if value is None:
                return "-"
            return "{:.3f} ({:.3f})".format(value["mean"], value["std"])
        return "| {} | {} | {} | {} | {} | {} |".format(
            metric, scope, cell("mae"), cell("rmse"), cell("pearson"),
            cell("spearman"),
        )

    for metric in METRICS:
        block = summary["metrics"][metric]
        lines.append(row(metric, "residue-pooled", block["overall"]["fold_mean"]))
        lines.append(row(metric, "macro-by-system", block["macro_by_system"]))
        for source, source_block in sorted(block["by_dataset"].items()):
            lines.append(
                row(metric, "source {} pooled".format(source),
                    source_block["residue_pooled"])
            )
            lines.append(
                row(metric, "source {} macro-system".format(source),
                    source_block["macro_by_system"])
            )
    amino = summary["metrics"]["chi_entropy"].get("amino_acid_corrected", {})
    if amino.get("available"):
        lines += [
            "",
            "## H-chi amino-acid-type correction (priors fitted on each training fold)",
            "",
        ]
        lines.append(row("chi_entropy", "residual pooled", amino["residual_overall"]))
        lines.append(row("chi_entropy", "residual macro-system",
                         amino["residual_macro_by_system"]))
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--evaluation-dir",
        default="artifacts/performance_first_v2/evaluation",
    )
    parser.add_argument(
        "--output",
        default="artifacts/performance_first_v2/evaluation/summary_five_fold.json",
    )
    parser.add_argument(
        "--markdown",
        default="artifacts/performance_first_v2/evaluation/summary_five_fold.md",
    )
    args = parser.parse_args()
    summary = aggregate(Path(args.evaluation_dir), Path(args.output))
    markdown = _render_markdown(summary)
    Path(args.markdown).write_text(markdown, encoding="utf-8")
    print("wrote {} and {}".format(args.output, args.markdown))
    print(markdown)


if __name__ == "__main__":
    main()
