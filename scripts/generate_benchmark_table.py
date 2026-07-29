# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (table) publication benchmark table → cross-type inherit → param inherit
# RULE: "native run" = load pre-rendered PNG via Image.open().ax.imshow().
#       "param inherit" = drawing function below that copies Class A/B/C values.
#       If a panel says "native run" and you write a drawing function, you broke the contract.

"""Generate the manuscript Table 1 for FlexPPI versus FastProtFlex.

The table contains trajectory-derived RMSF experiments for ATLAS single-chain
proteins, DynaRepo single-chain proteins, and DynaRepo protein complexes. All
displayed values are read from evaluation reports.
"""

from __future__ import annotations

# Academic Figure Skill Typography Baseline — COPY VERBATIM, place at TOP of script
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

# Academic Figure Skill Export Baseline — COPY VERBATIM
mpl.rcParams.update({
    "pdf.fonttype": 42,         # TrueType font embedding
    "svg.fonttype": "none",     # editable text in SVG
    "savefig.bbox": "tight",    # trim whitespace
    "savefig.dpi": 300,
})
# Keep the tightly cropped vector master within the 183-mm double-column limit.
mpl.rcParams["savefig.pad_inches"] = 0.0

def save_cns_figure(fig, filename):
    """Standard Academic Figure Skill export: vector PDF + 300dpi PNG preview."""
    fig.savefig(f"{filename}.pdf", bbox_inches="tight", dpi=300)
    fig.savefig(f"{filename}.png", bbox_inches="tight", dpi=300)


import argparse
import csv
import json
import math
import shutil
from pathlib import Path
from typing import Any, Dict, List, Mapping, Sequence, Tuple

import matplotlib.pyplot as plt


METRICS: Tuple[Tuple[str, str, bool], ...] = (
    ("spearman", "ρ", True),
    ("pearson", "r", True),
    ("z_mae", "z-MAE", False),
    ("z_rmse", "z-RMSE", False),
    ("top20_auroc", "AUROC", True),
    ("top20_average_precision", "AP", True),
    ("top20_precision_at_k", "P@K", True),
)

GROUP_HEADERS: Tuple[Tuple[str, int, int], ...] = (
    ("FLEXIBILITY RANKING ↑", 0, 2),
    ("STANDARDIZED ERROR ↓", 2, 4),
    ("TOP-20% FLEXIBLE RESIDUES ↑", 4, 7),
)

COHORTS = ("ATLAS single-chain", "DynaRepo single-chain", "DynaRepo complexes")


def _cohort_label(dataset: str, sample_count: int, sample_unit: str) -> Tuple[str, str]:
    if dataset == "ATLAS single-chain":
        return "ATLAS", f"single-chain, n={sample_count} {sample_unit}"
    if dataset == "DynaRepo single-chain":
        return "DynaRepo", f"single-chain, n={sample_count} {sample_unit}"
    if dataset == "DynaRepo complexes":
        return "DynaRepo", f"complexes, n={sample_count} {sample_unit}"
    raise ValueError(f"Unexpected cohort: {dataset}")


def _load(path: Path) -> Dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def _value(metric: Mapping[str, Any]) -> float:
    value = float(metric["macro_mean"])
    if not math.isfinite(value):
        raise ValueError(f"Non-finite metric value: {metric}")
    return value


def _row(
    dataset: str,
    sample_count: int,
    sample_unit: str,
    model: str,
    context: str,
    metrics: Mapping[str, Any],
    source_key: str,
) -> Dict[str, Any]:
    return {
        "dataset": dataset,
        "sample_count": int(sample_count),
        "sample_unit": sample_unit,
        "model": model,
        "context": context,
        "source_key": source_key,
        **{key: _value(metrics[key]) for key, _, _ in METRICS},
    }


def build_rows(
    atlas: Mapping[str, Any],
    context_report: Mapping[str, Any],
) -> List[Dict[str, Any]]:
    """Build exactly two method rows for each paper-facing cohort."""
    rows: List[Dict[str, Any]] = []

    atlas_count = int(atlas["dataset"]["paired_evaluable_proteins"])
    rows.extend(
        [
            _row(
                "ATLAS single-chain", atlas_count, "proteins", "FastProtFlex", "single-chain",
                atlas["methods"]["fastprotflex"], "atlas.methods.fastprotflex",
            ),
            _row(
                "ATLAS single-chain", atlas_count, "proteins", "FlexPPI", "single-chain",
                atlas["methods"]["flexppi_pretrain"], "atlas.methods.flexppi_pretrain",
            ),
        ]
    )

    dyna_single = context_report["single_chain"]["dynarepo_single_protein_chain"]
    dyna_single_count = int(dyna_single["system_count"])
    dyna_single_methods = dyna_single["methods"]
    rows.extend(
        [
            _row(
                "DynaRepo single-chain", dyna_single_count, "proteins",
                "FastProtFlex", "single-chain",
                dyna_single_methods["fastprotflex"],
                "context.single_chain.dynarepo_single_protein_chain.methods.fastprotflex",
            ),
            _row(
                "DynaRepo single-chain", dyna_single_count, "proteins",
                "FlexPPI", "single-chain",
                dyna_single_methods["flexppi_pretrain"],
                "context.single_chain.dynarepo_single_protein_chain.methods.flexppi_pretrain",
            ),
        ]
    )

    complex_block = context_report["protein_complex"]
    complex_count = int(complex_block["system_count"])
    complex_methods = complex_block["equal_complex_methods"]
    rows.extend(
        [
            _row(
                "DynaRepo complexes", complex_count, "complexes", "FastProtFlex", "chain-wise",
                complex_methods["fastprotflex_chain"],
                "context.protein_complex.equal_complex_methods.fastprotflex_chain",
            ),
            _row(
                "DynaRepo complexes", complex_count, "complexes", "FlexPPI", "complex-aware",
                complex_methods["flexppi_pretrain_complex"],
                "context.protein_complex.equal_complex_methods.flexppi_pretrain_complex",
            ),
        ]
    )

    _validate_rows(rows)
    return rows


def _validate_rows(rows: Sequence[Mapping[str, Any]]) -> None:
    datasets = list(COHORTS)
    if [row["dataset"] for row in rows] != [name for name in datasets for _ in range(2)]:
        raise ValueError("Table rows must contain the three retained RMSF cohorts in order.")
    if any("FT" in str(row["model"]) or "finetune" in str(row["model"]).lower() for row in rows):
        raise ValueError("FlexPPI-FT must not appear in the presentation table.")
    for dataset in datasets:
        subset = [row for row in rows if row["dataset"] == dataset]
        if {row["model"] for row in subset} != {"FastProtFlex", "FlexPPI"}:
            raise ValueError(f"{dataset} does not contain the two required methods.")
    for row in rows:
        for key, _, _ in METRICS:
            if not math.isfinite(float(row[key])):
                raise ValueError(f"Non-finite {key} in {row['dataset']} / {row['model']}")


def _best_flags(rows: Sequence[Mapping[str, Any]]) -> Dict[Tuple[int, str], bool]:
    flags: Dict[Tuple[int, str], bool] = {}
    for dataset in COHORTS:
        group = [row for row in rows if row["dataset"] == dataset]
        for key, _, higher_is_better in METRICS:
            values = [float(row[key]) for row in group]
            best = max(values) if higher_is_better else min(values)
            for row in group:
                flags[(id(row), key)] = math.isclose(float(row[key]), best, rel_tol=0.0, abs_tol=5e-7)
    return flags


def write_csv(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(
            ["dataset", "sample_count", "sample_unit", "method", "context", "source_key"]
            + [label for _, label, _ in METRICS]
        )
        for row in rows:
            writer.writerow(
                [
                    row["dataset"], row["sample_count"], row["sample_unit"],
                    row["model"], row["context"], row["source_key"],
                ]
                + [f"{float(row[key]):.6f}" for key, _, _ in METRICS]
            )


def _latex_value(value: float, best: bool) -> str:
    text = f"{value:.3f}"
    return r"\bfseries " + text if best else text


def write_latex(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    flags = _best_flags(rows)
    groups = [(name, [row for row in rows if row["dataset"] == name]) for name in COHORTS]
    lines = [
        "% Requires: \\usepackage{booktabs,multirow,siunitx}",
        r"\sisetup{detect-weight=true,table-number-alignment=center}",
        r"\begin{table*}[t]",
        r"\centering",
        r"\caption{Static RMSF prediction by FlexPPI and FastProtFlex on single-chain proteins and protein complexes.}",
        r"\label{tab:flexppi_fastprotflex_benchmark}",
        r"\begin{tabular}{@{}ll*{7}{S[table-format=1.3]}@{}}",
        r"\toprule",
        r"\multirow{2}{*}{Dataset} & \multirow{2}{*}{Method (input context)} & "
        r"\multicolumn{2}{c}{Flexibility ranking $\uparrow$} & "
        r"\multicolumn{2}{c}{Standardized error $\downarrow$} & "
        r"\multicolumn{3}{c}{Top-20\% flexible residues $\uparrow$} \\",
        r"\cmidrule(lr){3-4}\cmidrule(lr){5-6}\cmidrule(lr){7-9}",
        r" & & {$\rho$} & {$r$} & {z-MAE} & {z-RMSE} & {AUROC} & {AP} & {P@$K$} \\",
        r"\midrule",
    ]
    for group_index, (dataset, group) in enumerate(groups):
        for row_index, row in enumerate(group):
            title, detail = _cohort_label(
                str(dataset),
                int(row["sample_count"]),
                str(row["sample_unit"]),
            )
            dataset_cell = (
                rf"\multirow{{2}}{{*}}{{\shortstack[l]{{{title}\\{detail}}}}}"
                if row_index == 0 else ""
            )
            method_cell = f"{row['model']} ({row['context']})"
            values = [
                _latex_value(float(row[key]), flags[(id(row), key)])
                for key, _, _ in METRICS
            ]
            lines.append(f"{dataset_cell} & {method_cell} & " + " & ".join(values) + r" \\")
        if group_index < len(groups) - 1:
            lines.append(r"\midrule")
    lines.extend(
        [
            r"\bottomrule",
            r"\end{tabular}",
            r"\vspace{2pt}",
            r"\begin{minipage}{0.99\textwidth}\footnotesize",
            r"Values are target-macro means; bold indicates the best method within each dataset. "
            r"Errors are computed after within-target standardization, and $K=\lceil0.20L\rceil$. "
            r"All targets use trajectory-derived MD RMSF. Both methods receive a single chain in "
            r"the single-chain cohorts. For DynaRepo complexes, FastProtFlex is applied to each "
            r"protein chain independently, whereas FlexPPI receives the full multichain complex.",
            r"\end{minipage}",
            r"\end{table*}",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def draw_table(path_stem: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    """Render a compact 183-mm, booktabs-style publication table."""
    flags = _best_flags(rows)
    width_in = 183.0 / 25.4
    height_in = 92.0 / 25.4
    figure = plt.figure(figsize=(width_in, height_in), facecolor="white")
    axis = figure.add_axes([0.0, 0.025, 1.0, 0.95])
    axis.set_axis_off()

    dataset_width = 0.155
    method_width = 0.215
    numeric_width = (1.0 - dataset_width - method_width) / len(METRICS)
    numeric_start = dataset_width + method_width
    centers = [numeric_start + numeric_width * (index + 0.5) for index in range(len(METRICS))]

    axis.text(0.0, 0.982, "Table 1.", ha="left", va="top", fontsize=9.0, fontweight="bold", color=BLACK)
    axis.text(
        0.078, 0.982,
        "Static RMSF prediction on single-chain proteins and protein complexes.",
        ha="left", va="top", fontsize=8.6, color=BLACK,
    )

    top = 0.875
    group_y = 0.826
    cmid_y = 0.790
    subheader_y = 0.742
    header_bottom = 0.695
    axis.plot([0.0, 1.0], [top, top], color=BLACK, linewidth=0.85, solid_capstyle="butt")
    axis.text(0.0, subheader_y, "Dataset", ha="left", va="center", fontsize=7.7, fontweight="bold")
    axis.text(dataset_width, subheader_y, "Method", ha="left", va="center", fontsize=7.7, fontweight="bold")

    for title, start, end in GROUP_HEADERS:
        left = numeric_start + numeric_width * start
        right = numeric_start + numeric_width * end
        axis.text((left + right) / 2.0, group_y, title, ha="center", va="center", fontsize=6.6, fontweight="bold")
        axis.plot([left + 0.006, right - 0.006], [cmid_y, cmid_y], color="#555555", linewidth=0.50)
    for center, (_, label, _) in zip(centers, METRICS):
        axis.text(
            center, subheader_y, label, ha="center", va="center", fontsize=7.5,
            fontstyle="italic" if label in {"ρ", "r", "z-MAE", "z-RMSE"} else "normal",
        )
    axis.plot([0.0, 1.0], [header_bottom, header_bottom], color="#555555", linewidth=0.62)

    y_positions = [0.640, 0.568, 0.488, 0.416, 0.336, 0.264]
    groups = [(name, [row for row in rows if row["dataset"] == name]) for name in COHORTS]
    row_index = 0
    for group_index, (dataset, group) in enumerate(groups):
        group_center = (y_positions[row_index] + y_positions[row_index + 1]) / 2.0
        title, detail = _cohort_label(
            str(dataset),
            int(group[0]["sample_count"]),
            str(group[0]["sample_unit"]),
        )
        axis.text(0.0, group_center + 0.012, title, ha="left", va="center", fontsize=7.7, fontweight="bold")
        axis.text(
            0.0, group_center - 0.022,
            detail,
            ha="left", va="center", fontsize=5.7, color="#666666",
        )
        for row in group:
            y = y_positions[row_index]
            axis.text(dataset_width, y + 0.010, str(row["model"]), ha="left", va="center", fontsize=7.6)
            axis.text(
                dataset_width, y - 0.019, f"({row['context']})",
                ha="left", va="center", fontsize=5.8, color="#666666",
            )
            for center, (key, _, _) in zip(centers, METRICS):
                axis.text(
                    center, y, f"{float(row[key]):.3f}", ha="center", va="center",
                    fontsize=7.55, color=BLACK,
                    fontweight="bold" if flags[(id(row), key)] else "normal",
                )
            row_index += 1
        if group_index < len(groups) - 1:
            boundary = (y_positions[row_index - 1] + y_positions[row_index]) / 2.0
            axis.plot([0.0, 1.0], [boundary, boundary], color="#888888", linewidth=0.42)

    bottom = 0.222
    axis.plot([0.0, 1.0], [bottom, bottom], color=BLACK, linewidth=0.85, solid_capstyle="butt")
    axis.text(
        0.0, 0.177,
        "Values are target-macro means; bold indicates the best method within each dataset. "
        "z-errors use within-target standardized profiles; K=ceil(0.20L).",
        ha="left", va="top", fontsize=6.15, color=BLACK,
    )
    axis.text(
        0.0, 0.128,
        "All three cohorts use trajectory-derived MD RMSF.",
        ha="left", va="top", fontsize=6.15, color=BLACK,
    )
    axis.text(
        0.0, 0.086,
        "On complexes, FastProtFlex is chain-wise and FlexPPI receives the full multichain structure.",
        ha="left", va="top", fontsize=6.15, color=BLACK,
    )
    axis.set_xlim(0.0, 1.0)
    axis.set_ylim(0.055, 1.0)

    path_stem.parent.mkdir(parents=True, exist_ok=True)
    save_cns_figure(figure, str(path_stem))
    figure.savefig(path_stem.with_suffix(".svg"), bbox_inches="tight")
    # A 600-dpi line-art preview is also supplied for journals that request it.
    figure.savefig(path_stem.with_name(path_stem.name + "_600dpi").with_suffix(".png"), bbox_inches="tight", dpi=600)
    plt.close(figure)


def write_manifest(
    path: Path,
    rows: Sequence[Mapping[str, Any]],
    inputs: Mapping[str, Path],
) -> None:
    payload = {
        "table": "table1_static_benchmark",
        "cohorts": ["ATLAS single-chain", "DynaRepo single-chain", "DynaRepo complexes"],
        "rows": list(rows),
        "input_reports": {name: str(path.resolve()) for name, path in inputs.items()},
        "aggregation": "macro mean over proteins or complexes as stated per row",
        "top_fraction_definition": "top 20%; K = ceil(0.20 * target length)",
        "paper_scope": {
            "included": [
                "ATLAS single-chain RMSF",
                "DynaRepo single-chain RMSF",
                "DynaRepo complex RMSF with chainwise versus full-complex context",
            ],
        },
        "scope_warning": (
            "ATLAS is a same-domain diagnostic. DynaRepo is an asymmetric practical comparison "
            "because FlexPPI was trained/adapted on DynaRepo while FastProtFlex was trained on ATLAS."
        ),
    }
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--atlas",
        default="artifacts/evaluation/fastprotflex_comparison/atlas_comparison_report.json",
    )
    parser.add_argument(
        "--context",
        default="artifacts/evaluation/fastprotflex_comparison/chain_context/chain_context_comparison_report.json",
    )
    parser.add_argument("--output-dir", default="artifacts/paper/benchmark_figures")
    parser.add_argument("--pdf-output-dir", default="output/pdf/flexppi_fastprotflex_benchmark")
    arguments = parser.parse_args()

    input_paths = {
        "atlas": Path(arguments.atlas),
        "context": Path(arguments.context),
    }
    reports = {name: _load(path) for name, path in input_paths.items()}
    rows = build_rows(reports["atlas"], reports["context"])

    output = Path(arguments.output_dir)
    stem = output / "table1_static_benchmark"
    write_csv(stem.with_suffix(".csv"), rows)
    write_latex(stem.with_suffix(".tex"), rows)
    draw_table(stem, rows)
    write_manifest(
        output / "table1_static_benchmark_manifest.json",
        rows,
        input_paths,
    )

    pdf_output = Path(arguments.pdf_output_dir)
    pdf_output.mkdir(parents=True, exist_ok=True)
    shutil.copy2(stem.with_suffix(".pdf"), pdf_output / stem.with_suffix(".pdf").name)
    print(stem.resolve())


if __name__ == "__main__":
    main()
