# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# -*- coding: utf-8 -*-
# (a) RMSF Spearman -> LineTrend/plot_sweep.py + user reference -> param inherit
# (b) RMSF MAE -> LineTrend/plot_sweep.py + user reference -> param inherit
# (c) H-chi Spearman -> LineTrend/plot_sweep.py + user reference -> param inherit
# (d) H-chi MAE -> LineTrend/plot_sweep.py + user reference -> param inherit
# User-requested visual reference overrides the asset's colors and layout:
# white canvas, narrow boxed facet headers, muted circle/line series.
# No native run: the asset has embedded optimization-sweep data and TeX fonts.
# RULE: "native run" = load pre-rendered PNG via Image.open().ax.imshow().
#       "param inherit" = drawing function below that copies Class A/B/C values.
#       If a panel says "native run" and you write a drawing function, you broke the contract.

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

def save_cns_figure(fig, filename):
    """Standard Academic Figure Skill export: vector PDF + 300dpi PNG preview."""
    fig.savefig(f"{filename}.pdf", bbox_inches="tight", dpi=300)
    fig.savefig(f"{filename}.png", bbox_inches="tight", dpi=300)


import argparse
import csv
import hashlib
import json
import math
from collections import Counter
from pathlib import Path
from statistics import mean, stdev

mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.ticker import FormatStrFormatter

ROOT = Path(__file__).resolve().parents[1]
STEM = "fig_flexppi_v2_five_fold_reference"
# Colors are deliberately taken from the supplied muted red/navy/teal palette.
# Every series is FlexPPI; series identify evaluation subsets, not methods.
SERIES = {
    "overall": {"label": "Overall", "color": "#D44D4C", "ls": "-", "lw": 1.6},
    "atlas": {"label": "ATLAS", "color": "#4E617F", "ls": (0, (4, 2)), "lw": 1.15},
    "dynarepo": {"label": "DynaRepo", "color": "#76A99B", "ls": (0, (6, 1.6)), "lw": 1.15},
    "itrace": {"label": "ITrace", "color": "#AAA56B", "ls": (0, (1.3, 1.5)), "lw": 1.15},
}
PANELS = [
    ("rmsf", "spearman", "RMSF · Spearman", "Spearman ρ ↑", (0.70, 1.00), [0.70, 0.80, 0.90, 1.00]),
    ("rmsf", "mae", "RMSF · MAE", "MAE (Å) ↓", (0.0, 0.85), [0.0, 0.2, 0.4, 0.6, 0.8]),
    ("chi_entropy", "spearman", r"$H_\chi$ · Spearman", "Spearman ρ ↑", (0.90, 1.00), [0.90, 0.925, 0.95, 0.975, 1.00]),
    ("chi_entropy", "mae", r"$H_\chi$ · MAE", "MAE (nats) ↓", (0.0, 0.21), [0.0, 0.05, 0.10, 0.15, 0.20]),
]


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def read_jsonl(path):
    return [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]


def read_data(evaluation_dir, split_dir):
    rows, provenance, test_sets, train_sets = [], [], [], []
    for fold in range(5):
        path = evaluation_dir / f"fold_{fold}_test.json"
        report = json.loads(path.read_text(encoding="utf-8"))
        split = report["split"]
        assert report["checkpoint_epoch"] == 40
        assert split["role"] == "untouched_internal_test"
        assert split["validation_samples"] == 0
        assert split["samples_before_limit"] == split["samples_evaluated"] == split["test_samples"]
        train_path = split_dir / f"fold_{fold}_train.jsonl"
        test_path = split_dir / f"fold_{fold}_test.jsonl"
        train_records, test_records = read_jsonl(train_path), read_jsonl(test_path)
        train_ids = {record["sample_id"] for record in train_records}
        test_ids = {record["sample_id"] for record in test_records}
        assert len(train_ids) == len(train_records) == split["train_samples"]
        assert len(test_ids) == len(test_records) == split["test_samples"]
        assert not train_ids.intersection(test_ids)
        assert all(record["metadata"]["fold"] == fold for record in test_records)
        test_sets.append(test_ids)
        train_sets.append(train_ids)
        dataset_counts = Counter(record["dataset"] for record in test_records)
        for task in ("rmsf", "chi_entropy"):
            block = report["metrics"][task]
            expected = {"atlas", "dynarepo", "itrace"} if task == "rmsf" else {"dynarepo", "itrace"}
            assert set(block["by_dataset"]) == expected
            groups = {"overall": block["overall"], **block["by_dataset"]}
            assert sum(item["system_count"] for item in block["by_dataset"].values()) == block["overall"]["system_count"]
            assert sum(item["count"] for item in block["by_dataset"].values()) == block["overall"]["count"]
            for source, values in groups.items():
                n_systems = int(values["system_count"])
                assert values["macro_by_system"]["correlation_system_count"] == n_systems
                if source != "overall":
                    assert n_systems <= dataset_counts[source]
                for metric in ("spearman", "mae"):
                    value = float(values["macro_by_system"][metric])
                    assert math.isfinite(value)
                    assert (-1 <= value <= 1) if metric == "spearman" else value >= 0
                    rows.append({
                        "task": task, "metric": metric, "subset": source,
                        "fold_file_index": fold, "fold_display": fold + 1,
                        "value": value, "n_systems": n_systems, "n_residues": int(values["count"]),
                        "train_systems": split["train_samples"], "heldout_systems": split["test_samples"],
                        "checkpoint_epoch": 40, "source_file": str(path.resolve()),
                        "json_field": f"metrics.{task}." + ("overall" if source == "overall" else f"by_dataset.{source}") + f".macro_by_system.{metric}",
                    })
            for metric in ("spearman", "mae"):
                weighted = sum(item["macro_by_system"][metric] * item["system_count"] for item in block["by_dataset"].values()) / block["overall"]["system_count"]
                assert abs(weighted - block["overall"]["macro_by_system"][metric]) < 1e-10
        config_path = evaluation_dir.parent / f"fold_{fold}" / "resolved_config.json"
        config = json.loads(config_path.read_text(encoding="utf-8"))
        provenance.append({
            "fold_file_index": fold, "checkpoint_reported": report["checkpoint"],
            "seed_recorded": config["seed"], "epoch": 40,
            "report_sha256": sha256(path), "train_manifest_sha256": sha256(train_path),
            "test_manifest_sha256": sha256(test_path), "config_sha256": sha256(config_path),
        })
    for i in range(5):
        for j in range(i):
            assert not test_sets[i].intersection(test_sets[j])
    assert len(set.union(*test_sets)) == 3068
    for fold in range(5):
        assert train_sets[fold] == set.union(*test_sets) - test_sets[fold]
    assert len(rows) == 70  # RMSF: 4 subsets x 2 metrics x 5 folds; Hchi: 3 x 2 x 5.
    summaries = []
    for task, metric, *_ in PANELS:
        for source in SERIES:
            selected = [row for row in rows if (row["task"], row["metric"], row["subset"]) == (task, metric, source)]
            if not selected:
                assert task == "chi_entropy" and source == "atlas"
                continue
            assert [row["fold_file_index"] for row in selected] == list(range(5))
            values = [row["value"] for row in selected]
            summaries.append({"task": task, "metric": metric, "subset": source,
                              "mean": mean(values), "sample_sd_ddof1": stdev(values),
                              "n_folds": 5, "total_systems": sum(row["n_systems"] for row in selected),
                              "total_valid_residues": sum(row["n_residues"] for row in selected)})
    return rows, summaries, provenance


def draw_panel(ax, panel, letter, rows, summaries, compact=False):
    task, metric, title, ylabel, ylim, yticks = panel
    ax.add_patch(Rectangle((0, 1.025), 1, 0.092, transform=ax.transAxes,
                           facecolor="white", edgecolor=BLACK, linewidth=0.75, clip_on=False))
    ax.text(0.5, 1.071, title, transform=ax.transAxes, ha="center", va="center", fontsize=8)
    ax.text(-0.17 if not compact else -0.25, 1.14, letter, transform=ax.transAxes,
            va="center", ha="left", fontsize=9, fontweight="bold")
    # Source series first, overall last: overall remains identifiable at crossings.
    for source in ("atlas", "dynarepo", "itrace", "overall"):
        selected = [row for row in rows if (row["task"], row["metric"], row["subset"]) == (task, metric, source)]
        if not selected:
            continue
        style = SERIES[source]
        xs, ys = [row["fold_display"] for row in selected], [row["value"] for row in selected]
        assert all(ylim[0] <= value <= ylim[1] for value in ys), "Point outside visible axes"
        ax.plot(xs, ys, color=style["color"], linestyle=style["ls"], linewidth=style["lw"],
                alpha=0.75, solid_capstyle="round", zorder=3 if source == "overall" else 2)
        if source == "atlas":
            # An open ring remains visible when ATLAS nearly coincides with Overall.
            ax.scatter(xs, ys, facecolor="none", edgecolor=style["color"], s=34,
                       linewidth=0.8, zorder=5)
        else:
            ax.scatter(xs, ys, color=style["color"], s=21 if source == "overall" else 17,
                       edgecolor="white", linewidth=0.22, zorder=6 if source == "overall" else 4)
    ax.set(xlim=(0.7, 5.3), ylim=ylim, xticks=[1, 2, 3, 4, 5], yticks=yticks,
           xlabel="Held-out fold", ylabel=ylabel)
    ax.yaxis.set_major_formatter(FormatStrFormatter("%.3f" if task == "chi_entropy" and metric == "spearman" else "%.2f"))
    ax.tick_params(length=2.6, pad=3)
    stat = next(row for row in summaries if (row["task"], row["metric"], row["subset"]) == (task, metric, "overall"))
    text = f"Overall: {stat['mean']:.3f} ± {stat['sample_sd_ddof1']:.3f}"
    ax.text(0.04, 0.055, text, transform=ax.transAxes, color=SERIES["overall"]["color"], fontsize=6.8)


def export_figure(fig, stem):
    # Fixed canvas preserves exactly 183 mm, without tight-bbox changing final size.
    with mpl.rc_context({"savefig.bbox": None}):
        for extension in ("pdf", "svg", "png"):
            fig.savefig(str(stem) + "." + extension, dpi=300, facecolor="white")


def compose_figure(rows, summaries, output_dir, layout):
    grid = layout == "grid"
    fig, axes = plt.subplots(2 if grid else 1, 2 if grid else 4,
                             figsize=(183 / 25.4, (141 if grid else 78) / 25.4))
    if grid:
        fig.subplots_adjust(left=0.095, right=0.97, bottom=0.145, top=0.855, wspace=0.35, hspace=0.59)
    else:
        fig.subplots_adjust(left=0.079, right=0.99, bottom=0.265, top=0.72, wspace=0.66)
    fig.text(0.095 if grid else 0.063, 0.975, "FlexPPI", fontsize=9, fontweight="bold", va="top")
    handles = [Line2D([0], [0], color=style["color"], ls=style["ls"], lw=style["lw"],
                      marker="o", markersize=4,
                      markerfacecolor="white" if source == "atlas" else style["color"],
                      markeredgecolor=style["color"] if source == "atlas" else "white",
                      markeredgewidth=0.8 if source == "atlas" else 0.22,
                      label=style["label"]) for source, style in SERIES.items()]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.60, 0.982),
               ncol=4, columnspacing=1.5, handlelength=2.3, handletextpad=0.6, fontsize=7.5,
               title="Evaluation subset", title_fontsize=6.6)
    for index, (ax, panel) in enumerate(zip(axes.flat, PANELS)):
        draw_panel(ax, panel, "abcd"[index], rows, summaries, compact=not grid)
    fig.text(0.5, 0.069 if grid else 0.118,
             "Points: equal-system means · annotations: five-fold mean ± sample SD (n = 5)",
             ha="center", fontsize=6.8)
    fig.text(0.5, 0.040 if grid else 0.067,
             r"Epoch 40 · RMSF: 3,068 systems · $H_\chi$: 1,139 systems · no ATLAS $H_\chi$ labels",
             ha="center", fontsize=6.6, color="#555555")
    filename = output_dir / (STEM + ("_row" if not grid else ""))
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    for artist in fig.findobj(mpl.text.Text):
        if not artist.get_visible() or not artist.get_text():
            continue
        box = artist.get_window_extent(renderer)
        assert box.x0 >= 0 and box.y0 >= 0 and box.x1 <= fig.bbox.width and box.y1 <= fig.bbox.height, (
            f"Text outside canvas in {layout}: {artist.get_text()} {box}")
    export_figure(fig, filename)
    bounds = [list(ax.get_position().bounds) for ax in axes.flat]
    plt.close(fig)
    return {"layout": layout, "width_mm": 183, "height_mm": 141 if grid else 78,
            "panel_bounds_figure_fraction": bounds, "stem": str(filename)}


def write_report(output_dir, summaries, provenance):
    overall = [row for row in summaries if row["subset"] == "overall"]
    report = {
        "model_version": "performance_first_v2", "statistics": summaries, "provenance": provenance,
        "n_definition": "A plotted point is one held-out fold; within a fold, whole-system metrics are averaged with equal system weight.",
        "center": "Arithmetic mean of five fold scores; folds have equal weight.",
        "spread": "Sample standard deviation, ddof=1 (denominator 4); descriptive fold variation, not a confidence interval.",
        "metrics": {"spearman": "Rank correlation between prediction and MD target within the entire system, then equal-system mean.",
                    "mae": "Mean absolute prediction error over valid residues within the entire system, then equal-system mean."},
        "raw_hchi": True, "statistical_test": "None; descriptive plot.", "multiple_comparison_correction": "Not applicable.",
        "baseline_definition": "No competing method plotted. Colored curves are evaluation subsets of FlexPPI.",
        "missing_data": "ATLAS does not supply Hchi labels. No ATLAS Hchi series is fabricated. All 70 selected aggregate metric values are displayed.",
        "validation": "Five disjoint test folds verified against manifests; train/test IDs disjoint within each fold; epoch 40 fixed; no validation split.",
        "scope": "RMSF and raw Hchi, Spearman and MAE; other report fields remain available in original JSONs and are not part of these four panels.",
    }
    (output_dir / "statistics_and_provenance.json").write_text(json.dumps(report, indent=2, ensure_ascii=False, allow_nan=False), encoding="utf-8")
    table = "\n".join(f"| {row['task']} | {row['metric']} | {row['mean']:.6f} | {row['sample_sd_ddof1']:.6f} | {row['total_systems']:,} | {row['total_valid_residues']:,} |" for row in overall)
    rmsf_count = next(row for row in overall if row["task"] == "rmsf")
    hchi_count = next(row for row in overall if row["task"] == "chi_entropy")
    caption = f"""# FlexPPI 五折性能图（用户参考图风格）

使用上一轮核实的 **performance-first v2，固定 epoch 40** 测试结果。主文件为 `{STEM}.pdf/.svg/.png`（2×2）；`{STEM}_row` 是相同四面板的横排版本，参考用户附图底部的排版。SVG保留文字、点、线和标题框可编辑，PDF嵌入TrueType字体，PNG为300 dpi，均为183 mm宽。

## 图注（可用于论文）

Five-fold cross-validation performance of FlexPPI. (a,b) RMSF rank correlation and absolute error. (c,d) Raw side-chain conformational entropy (Hχ) rank correlation and absolute error. Each point is the equal-system mean of the corresponding metric in one held-out test fold. All curves describe FlexPPI, evaluated overall or within ATLAS, DynaRepo and ITrace; ATLAS supplies no Hχ labels. Every system is evaluated by its own held-out model at the pre-specified epoch 40. Annotations report the arithmetic mean ± sample standard deviation of the five fold scores (ddof = 1). Lines connect categorical folds as visual guides and do not represent training trajectories. Overall results weight systems equally within each fold, not data sources equally. There are {rmsf_count['total_systems']:,} RMSF systems and {hchi_count['total_systems']:,} Hχ systems, contributing {rmsf_count['total_valid_residues']:,} and {hchi_count['total_valid_residues']:,} valid residues, respectively. RMSF is expressed in Å and Hχ in nats. Hχ here is the raw target, before subtracting amino-acid-type effects. No significance test was performed.

## 指标与统计口径

Spearman 衡量同一体系中各残基柔性高低的排序是否一致，越高越好；MAE 是预测值与MD真值之间绝对差的平均，越低越好。先对每个完整体系计算指标，再使各体系等权平均，得到图上的一个点；最后将五个折的点等权平均，计算图中均值及标准差。对于复合物，这是整个体系内的指标，不是先分链再平均的旧横评口径。两个指标不能互换：排序正确仍可能存在绝对量值偏差。

| Target | Metric | Five-fold mean | Sample SD | Systems | Valid residues |
|---|---|---:|---:|---:|---:|
{table}

每次训练使用其余四折，共2,454–2,455个体系；测试613–614个体系。固定40个epoch，无验证集，测试结果不用于挑选checkpoint。这五个模型对应五个留出划分，训练集之间存在重叠，因此这里的SD仅描述折间波动，不能称为五次独立重复的置信区间。ITrace属于v2训练来源，其留出折结果是内部OOF，不是零样本外部测试。来源间MD观测时长和体系构成不同，不能仅根据这些曲线的高低推断某一来源在一般意义上更容易。

旧 `summary_five_fold.md` 使用ddof=0，本图从五个原始报告重算ddof=1。因此均值相同，SD比旧表大约大11.8%。源数据CSV包含每个点的绝对来源路径、JSON字段、折号、体系数和有效残基数；JSON包含报告及划分文件SHA-256。图中Fold 1–5对应文件fold_0–fold_4。

## 复现

```powershell
.\\.venv-md\\Scripts\\python.exe scripts/plot_flexppi_five_fold_reference.py
```

参数 `--evaluation-dir`、`--split-dir`、`--output-dir` 可指定输入输出，`--layout grid/row/both` 可选择排版。脚本仅使用Python标准库和Matplotlib，不重新训练或推理。`source_data.csv`给出70个点；`five_fold_mean_sd.csv`给出14组均值与样本SD。
"""
    (output_dir / "README_and_caption.md").write_text(caption, encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description="Render four FlexPPI v2 held-out-fold panels in the user's reference style.")
    parser.add_argument("--evaluation-dir", type=Path, default=ROOT / "artifacts/performance_first_v2/evaluation")
    parser.add_argument("--split-dir", type=Path, default=ROOT / "data/splits/state_rmsf_hchi_v2_performance_first")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "output/flexppi_v2_five_fold_reference")
    parser.add_argument("--layout", choices=("grid", "row", "both"), default="both")
    args = parser.parse_args()
    rows, summaries, provenance = read_data(args.evaluation_dir, args.split_dir)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.output_dir / "source_data.csv", rows)
    write_csv(args.output_dir / "five_fold_mean_sd.csv", summaries)
    write_report(args.output_dir, summaries, provenance)
    layouts = ("grid", "row") if args.layout == "both" else (args.layout,)
    rendering = [compose_figure(rows, summaries, args.output_dir, layout) for layout in layouts]
    (args.output_dir / "rendering.json").write_text(json.dumps(rendering, indent=2), encoding="utf-8")
    print(f"Saved {len(rendering)} four-panel layouts and all 70 plotted values to {args.output_dir}")


if __name__ == "__main__":
    main()
