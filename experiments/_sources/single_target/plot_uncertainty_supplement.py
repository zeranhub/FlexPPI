# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a,b) coverage curves -> LineTrend/plot_trend.py -> param inherit (scientific axis/CI data differ from cumulative dates)
# (c,d) source coverage -> cross-type inherit from LineTrend -> param inherit
# (e,f) source interval width -> cross-type inherit from LineTrend -> param inherit
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
from pathlib import Path
import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from completion_common import sha256, write_json

KINDS = ["model", "model_plus_label", "group_crossfit_scaled_model"]
COLORS = dict(zip(KINDS, [CATEGORICAL[0], CATEGORICAL[3], CATEGORICAL[4]]))
LABELS = dict(zip(KINDS, ["Model variance", "Model + label variance", "Post-hoc group calibration"]))
MARKERS = dict(zip(KINDS, ["o", "s", "^"]))
SOURCE_ORDER = ["all", "atlas", "dynarepo", "itrace"]
SOURCE_LABELS = ["All", "ATLAS", "DynaRepo", "ITrace"]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--analysis", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--variant", default="full")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    base = pd.read_csv(args.analysis / "uncertainty_coverage_source_fold_macro.csv")
    calibrated = pd.read_csv(args.analysis / "group_calibration_coverage_source_fold_macro.csv")
    calibrated = calibrated[calibrated.variance_kind == "group_crossfit_scaled_model"]
    points = pd.concat([base, calibrated], ignore_index=True)
    points = points[points.variant == args.variant].copy()
    assert set(points.variance_kind) == set(KINDS)
    assert set(points.metric) == {"rmsf", "chi_entropy"}
    assert set(points.fold) == set(range(5))
    assert points.observed_coverage.between(0, 1).all()
    assert (points.mean_interval_width > 0).all() and np.isfinite(points.mean_interval_width).all()
    assert not points.duplicated(["variance_kind", "metric", "dataset", "fold", "nominal_coverage"]).any()
    points.to_csv(args.output / "all_coverage_and_width_plot_points.csv", index=False)
    means = points.groupby(["variance_kind", "metric", "dataset", "nominal_coverage"])[["observed_coverage", "mean_interval_width"]].agg(["mean", "std"])
    means.columns = [key + "_" + suffix for key, suffix in means.columns]
    means = means.reset_index()
    means.to_csv(args.output / "all_coverage_and_width_plot_summaries.csv", index=False)
    table = means[means.nominal_coverage == 0.95].copy()
    table["metric_unit"] = table.metric.map({"rmsf": "Å", "chi_entropy": "nats"})
    table["method_zh"] = table.variance_kind.map({"model": "原始模型方差", "model_plus_label": "模型方差加标签变异", "group_crossfit_scaled_model": "事后分组校准"})
    table["coverage_95_percent_mean"] = 100 * table.observed_coverage_mean
    table["coverage_95_percent_sd"] = 100 * table.observed_coverage_std
    table["coverage_95_formatted"] = ["{:.2f} ± {:.2f}".format(row.coverage_95_percent_mean, row.coverage_95_percent_sd) for row in table.itertuples()]
    table["width_95_formatted"] = ["{:.3f} ± {:.3f}".format(row.mean_interval_width_mean, row.mean_interval_width_std) for row in table.itertuples()]
    table.to_csv(args.output / "table_uncertainty_95_three_line_source.csv", index=False, encoding="utf-8-sig")
    write_json(args.output / "table_uncertainty_95_three_line_source.json", dict(
      title="冻结五折测试集中的 95% 高斯区间覆盖率与平均宽度",
      columns=["指标", "来源", "方差处理", "覆盖率（%，均值 ± SD）", "平均宽度（均值 ± SD）"],
      rows=[["RMSF" if row.metric == "rmsf" else "Hχ", SOURCE_LABELS[SOURCE_ORDER.index(row.dataset)], row.method_zh,
             row.coverage_95_formatted, row.width_95_formatted + " " + row.metric_unit] for row in table.itertuples()],
      notes="各记录先按有效残基求平均，再对记录等权汇总；误差为五折之间的样本标准差，不是置信区间。事后校准只在同一折模型未参与训练的另一半关联分量上拟合；标签变异不是标签均值的标准误。"))
    # Six equal-weight panels answer complementary calibration questions.
    fig, axes = plt.subplots(3, 2, figsize=(183 / 25.4, 166 / 25.4))
    fig.subplots_adjust(left=0.10, right=0.98, bottom=0.07, top=0.89, hspace=0.58, wspace=0.36)
    legend_handles = []
    for column, metric in enumerate(["rmsf", "chi_entropy"]):
        ax = axes[0, column]
        for kind in KINDS:
            data = means[(means.metric == metric) & (means.dataset == "all") & (means.variance_kind == kind)].sort_values("nominal_coverage")
            handle = ax.errorbar(data.nominal_coverage, data.observed_coverage_mean,
                                yerr=data.observed_coverage_std, color=COLORS[kind], marker=MARKERS[kind],
                                markersize=3.4, linewidth=1.2, capsize=2, label=LABELS[kind])
            if column == 0:
                legend_handles.append(handle)
        ax.plot([0, 1], [0, 1], linestyle="--", color=GREY, linewidth=0.7, zorder=0)
        ax.set_xlim(0.45, 1.0)
        ax.set_ylim(0.15, 1.02)
        ax.set_xticks([0.5, 0.68, 0.8, 0.9, 0.95])
        ax.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])
        ax.set_xlabel("Nominal Gaussian coverage")
        ax.set_ylabel("Observed coverage")
        ax.set_title("RMSF" if metric == "rmsf" else "Hχ", fontweight="bold", pad=7)
        for row_index, value, y_label in [(1, "observed_coverage", "95% observed coverage"),
                                        (2, "mean_interval_width", "95% interval width (" + ("Å" if metric == "rmsf" else "nats") + ")")]:
            ax = axes[row_index, column]
            for method_index, kind in enumerate(KINDS):
                data = points[(points.metric == metric) & (points.nominal_coverage == 0.95) & (points.variance_kind == kind)]
                for source_index, source in enumerate(SOURCE_ORDER):
                    fold_values = data[data.dataset == source].sort_values("fold")[value].to_numpy()
                    assert len(fold_values) == 5
                    x = source_index + (method_index - 1) * 0.23
                    ax.scatter(x + np.linspace(-0.035, 0.035, 5), fold_values, s=10, color=COLORS[kind], alpha=0.4,
                               marker=MARKERS[kind], linewidths=0, zorder=2)
                    ax.errorbar([x], [fold_values.mean()], yerr=[fold_values.std(ddof=1)], color=COLORS[kind],
                                marker=MARKERS[kind], markersize=4, capsize=2, linewidth=1.0, zorder=3)
            ax.set_xlim(-0.55, 3.55)
            ax.set_xticks(range(4), SOURCE_LABELS)
            ax.set_ylabel(y_label)
            if row_index == 1:
                ax.axhline(0.95, linestyle="--", color=GREY, linewidth=0.7, zorder=0)
                ax.set_ylim(0, 1.04)
                ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
            else:
                ax.set_ylim(bottom=0)
    for index, ax in enumerate(axes.flat):
        ax.text(-0.16, 1.09, chr(ord("a") + index), transform=ax.transAxes, fontsize=10, fontweight="bold")
    fig.legend(legend_handles, [LABELS[kind] for kind in KINDS], loc="upper center", bbox_to_anchor=(0.54, 0.995),
               ncol=3, frameon=False, fontsize=7.3, columnspacing=1.3)
    filename = args.output / "figS_uncertainty_calibration_corrected_v1"
    save_cns_figure(fig, filename)
    fig.savefig(str(filename) + ".svg", bbox_inches="tight", dpi=300)
    plt.close(fig)
    write_json(args.output / "figure_statistics_audit.json", dict(status="rendered_pending_visual_qa", variant=args.variant,
      records=3063, source_records={"atlas": 1929, "dynarepo": 890, "itrace": 244}, folds=5, seed=7, epoch=40,
      displayed_values="All five-fold macro points for each selected target/source/method; no subsampling of records",
      spread="Sample SD across five folds; raw points on source panels are individual fold macros",
      test="Descriptive interval diagnostics; no significance testing or multiplicity adjustment",
      source_data_sha256={path.name: sha256(path) for path in [args.analysis / "uncertainty_coverage_source_fold_macro.csv", args.analysis / "group_calibration_coverage_source_fold_macro.csv"]},
      code_sha256=sha256(__file__), intervals="Gaussian intervals for reference label means; coverage improvement must be read with width",
      structural_boundary="Original folds are not complete unseen-family splits; calibration halves are disjoint in frozen sequence-associated components"))


if __name__ == "__main__":
    main()
