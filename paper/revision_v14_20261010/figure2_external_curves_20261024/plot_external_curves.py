# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a) PR -> cross-type inherit from AUROC/plot_AUROC1.py -> param inherit
# (b) ROC -> AUROC/plot_AUROC1.py -> param inherit
# The generic asset has five fixed series and integrates a pooled display curve.
# These seven-method record-mean curves require independently averaged AP/AUROC.
# Preserve its explicit axes and reference line; retain native vector data paths.
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

from pathlib import Path
import argparse
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.transforms import Bbox

OUT = Path(__file__).resolve().parent
# Match the blue emphasis of Figure 2 and assign a distinct line/marker to each
# external method. Values are computed before plotting; styling never ranks data.
METHODS = [
    ("Full", "FlexPPI", CATEGORICAL[0], "-", None),
    ("backflip-2.1 (assembly)", "BackFlip-2.1", CATEGORICAL_EXTENDED[9], (0, (4, 1.6)), "o"),
    ("backflip-1.0 (assembly)", "BackFlip-1.0", CATEGORICAL_EXTENDED[7], (0, (4, 1.4, 1, 1.4)), "s"),
    ("DynaProt-M", "DynaProt-M", CATEGORICAL[2], (0, (1, 1.5)), "^"),
    ("FastProtFlex", "FastProtFlex", CATEGORICAL[4], (0, (6, 1.6)), "D"),
    ("ANM", "ANM", CATEGORICAL[5], (0, (2, 1.4)), "v"),
    ("GNM", "GNM", GREY, (0, (4, 1.4, 1, 1.4, 1, 1.4)), "P"),
]


def load_data():
    curves = pd.read_csv(OUT / "mean_curves.csv", float_precision="round_trip")
    summary = pd.read_csv(OUT / "summary_metrics.csv", float_precision="round_trip")
    audit = json.loads((OUT / "calculation_audit.json").read_text(encoding="utf-8"))
    curves = curves[curves.scope == "all"]
    summary = summary[summary.scope == "all"]
    assert set(curves.method) == set(m[0] for m in METHODS)
    assert set(summary.method) == set(m[0] for m in METHODS)
    for method, _, _, _, _ in METHODS:
        for kind in ("PR", "ROC"):
            part = curves[(curves.method == method) & (curves.curve_type == kind)].sort_values("x")
            assert len(part) == 1001
            assert np.array_equal(part.x.to_numpy(), np.linspace(0, 1, 1001)) or np.allclose(part.x, np.linspace(0, 1, 1001), atol=1e-15)
            assert np.isfinite(part["mean"]).all()
            assert part["mean"].between(-1e-12, 1 + 1e-12).all()
            if kind == "ROC":
                assert np.diff(part["mean"]).min() >= -1e-12
                assert part["mean"].iloc[-1] == 1
        assert (summary.loc[summary.method == method, "n_records"] == 615).all()
    # Compute the actual equal-record prevalence, never assume it is exactly .10.
    per = pd.read_csv(OUT / "per_record_metrics.csv", float_precision="round_trip")
    sub = per[per.method == "Full"]
    rate_col = next((c for c in ("positive_rate", "prevalence", "positive_fraction") if c in sub), None)
    if rate_col is not None:
        prevalence = float(sub[rate_col].mean())
    else:
        prevalence = float(audit["all_selected_precision"])
    assert 0 < prevalence < 1
    return curves, summary, audit, prevalence


def metric_value(summary, method, metric):
    # Keep the scalar AP/AUROC calculated from the exact per-record thresholds.
    # Never replace it with the area under a plotted record-average grid curve.
    row = summary[(summary.method == method) & (summary.metric.str.upper() == metric.upper())]
    assert len(row) == 1, (method, metric)
    return float(row.iloc[0]["mean"])


def panel(ax, kind, curves, summary, prevalence, label, with_band=False):
    handles, labels = [], []
    metric = "AP" if kind == "PR" else "AUROC"
    for method, name, color, style, marker in reversed(METHODS):
        part = curves[(curves.method == method) & (curves.curve_type == kind)].sort_values("x")
        x, y = part.x.to_numpy(), part["mean"].to_numpy()
        if with_band and method == "Full" and part.ci_low.notna().all():
            ax.fill_between(x, part.ci_low, part.ci_high, color=color, alpha=.12, linewidth=0, zorder=2)
        ax.plot(x, y, color=color, linestyle=style,
                linewidth=1.65 if method == "Full" else 1.05,
                marker=marker, markevery=[180, 380, 580, 780] if marker else None,
                markersize=2.7, markerfacecolor="white", markeredgewidth=.65,
                zorder=10 if method == "Full" else 4, clip_on=True)
        if kind == "PR":
            # Explicit all-selected endpoint; do not overwrite first-hit recall=1.
            ax.plot([1, 1], [y[-1], prevalence], color=color, linewidth=.7,
                    linestyle=style, clip_on=False, zorder=3)
        else:
            # Preserve the empirical starting vertical edge at FPR=0.
            ax.plot([0, 0], [0, y[0]], color=color, linewidth=.7,
                    linestyle=style, clip_on=False, zorder=3)
    for method, name, color, style, marker in METHODS:
        handles.append(Line2D([0], [0], color=color, lw=1.65 if method == "Full" else 1.05,
                              ls=style, marker=marker, markersize=2.8,
                              markerfacecolor="white", markeredgewidth=.65))
        labels.append(f"{name}  {metric_value(summary, method, metric):.3f}")
    if kind == "PR":
        ax.axhline(prevalence, color="#B7B7B7", ls=(0, (2, 2)), lw=.7, zorder=1)
        ax.text(.97, prevalence + .023, f"Positive fraction {prevalence:.2%}",
                ha="right", va="bottom", fontsize=6.5, color="#777777")
        ax.set_xlabel("Recall", labelpad=4)
        ax.set_ylabel("Precision", labelpad=4)
        title = "Precision-recall"
    else:
        ax.plot([0, 1], [0, 1], color="#B7B7B7", ls=(0, (3, 2)), lw=.7, zorder=1)
        ax.text(.72, .64, "Chance", color="#888888", rotation=45,
                ha="center", va="center", fontsize=6.5)
        ax.set_xlabel("False-positive rate", labelpad=4)
        ax.set_ylabel("True-positive rate", labelpad=4)
        title = "ROC"
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.005)
    ticks = np.linspace(0, 1, 6)
    ax.set_xticks(ticks, [f"{t:.1f}" for t in ticks])
    ax.set_yticks(ticks, [f"{t:.1f}" for t in ticks])
    ax.tick_params(length=2.5, width=.6, pad=2.5)
    ax.grid(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color("#535353")
    ax.set_title(title, loc="left", pad=7, fontsize=8.5)
    if label:
        ax.text(-.22, 1.035, label, transform=ax.transAxes, fontsize=10,
                fontweight="bold", color=BLACK, va="bottom")
    legend = ax.legend(handles, labels, loc="upper left", bbox_to_anchor=(-.08, -.22),
                      ncol=2, fontsize=6.8, handlelength=2.25, handletextpad=.55,
                      columnspacing=1.35, borderaxespad=0, labelspacing=.60,
                      title=f"Method     {metric}", title_fontsize=7.0)
    for text, model in zip(legend.get_texts(), METHODS):
        if model[0] == "Full":
            text.set_fontweight("bold")
    return legend


def export_exact(fig, stem):
    # Exact physical dimensions are kept even with baseline savefig.bbox='tight'.
    box = Bbox.from_bounds(0, 0, *fig.get_size_inches())
    for ext in ("pdf", "svg", "png"):
        fig.savefig(OUT / f"{stem}.{ext}", bbox_inches=box, dpi=600 if ext == "png" else 300,
                    facecolor="white")
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    bound = fig.bbox
    outside = []
    for artist in fig.findobj(mpl.text.Text):
        if not artist.get_visible() or not artist.get_text().strip():
            continue
        bb = artist.get_window_extent(renderer)
        if bb.x0 < bound.x0 - 1 or bb.y0 < bound.y0 - 1 or bb.x1 > bound.x1 + 1 or bb.y1 > bound.y1 + 1:
            outside.append(artist.get_text())
    return {"stem": stem, "width_mm": round(fig.get_figwidth() * 25.4, 5),
            "height_mm": round(fig.get_figheight() * 25.4, 5), "outside_text": outside}


def compose_figure(curves, summary, prevalence, with_band=False):
    results = []
    fig = plt.figure(figsize=(183 / 25.4, 110 / 25.4))
    fig.text(.5, .980, "High-RMSF retrieval | DynaRepo complexes (n = 615)",
             ha="center", va="top", fontsize=9.0)
    for x, kind, label in [(.085, "PR", "a"), (.585, "ROC", "b")]:
        ax = fig.add_axes([x, .295, .350, .5823])
        panel(ax, kind, curves, summary, prevalence, label, with_band=with_band)
    results.append(export_exact(fig, "external_RMSF_PR_ROC"))
    plt.close(fig)
    for kind, stem, label in [("PR", "fig2e_external_PR", "e"), ("ROC", "figS_external_ROC", "")]:
        fig = plt.figure(figsize=(89 / 25.4, 110 / 25.4))
        fig.text(.56, .980, "DynaRepo complexes (n = 615)", ha="center", va="top",
                 fontsize=8.2, color="#333333")
        ax = fig.add_axes([.205, .295, .720, .5825])
        panel(ax, kind, curves, summary, prevalence, label, with_band=with_band)
        results.append(export_exact(fig, stem))
        plt.close(fig)
    (OUT / "plot_layout_audit.json").write_text(json.dumps(results, indent=2), encoding="utf-8")
    assert not any(r["outside_text"] for r in results), results
    return results


def main():
    parser = argparse.ArgumentParser(description="Plot record-mean external flexibility benchmark curves.")
    parser.add_argument("--full-ci-band", action="store_true", help="Optional Full-only pointwise bootstrap band; all methods' scalar intervals are always in source tables.")
    args = parser.parse_args()
    curves, summary, audit, prevalence = load_data()
    result = compose_figure(curves, summary, prevalence, args.full_ci_band)
    print(json.dumps({"figures": result, "positive_fraction": prevalence}, indent=2))


if __name__ == "__main__":
    main()


