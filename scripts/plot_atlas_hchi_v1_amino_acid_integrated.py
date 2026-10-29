# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (single) Grouped point intervals -> GroupedBarChart/plot_GroupedBarChartv1.py
#     + LineTrend/plot_sweep.py -> param inherit (cross-type)
# Both assets embed unrelated results; the requested single categorical overview
# is not their grouped bar / optimization trajectory structure. Retain grouped
# point/interval emphasis, external legend and clean spines, not tests or lines.
# User's compactness request overrides the old four-facet geometry.
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


import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "artifacts/evaluation_atlas_hchi_v1/analysis"
OUT = ROOT / "output/flexppi_atlas_hchi_v1_amino_acid_integrated"
STEM = "fig04_amino_acid_integrated"
SOURCES = ["atlas", "dynarepo", "itrace"]
LABELS = {"atlas": "ATLAS", "dynarepo": "DynaRepo", "itrace": "ITrace"}
COLORS = {"atlas": "#4E617F", "dynarepo": "#76A99B", "itrace": "#AAA56B"}
MARKERS = {"atlas": "o", "dynarepo": "s", "itrace": "^"}
AA_ORDER = list("CDEFHIKLMNPQRSTVWY")
CONTROL_ITEMS = [
    ("rmsf", "AA-only", "AA\nonly", 0.),
    ("rmsf", "FlexPPI", "Flex\nPPI", 1.6),
    ("chi_entropy", "AA-only", "AA\nonly", 3.5),
    ("chi_entropy", "FlexPPI", "Flex\nPPI", 5.1),
    ("chi_entropy", "AA-corrected", "AA\nadj.", 6.7),
]
AA_X = 9.0 + np.arange(18)


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def prepare_data(folder=DATA):
    summary = json.loads((folder / "analysis_summary.json").read_text())
    if summary["status"] != "complete" or summary["n_systems"] != 3068:
        raise ValueError("Only complete new-model OOF analysis is accepted")
    frames = {}
    for name in ["source_fold_metrics", "within_aa_metrics"]:
        path = folder / (name + ".csv")
        if sha256(path) != summary["source_files"][path.name]:
            raise ValueError("Source data fingerprint changed: " + name)
        frames[name] = pd.read_csv(path)
    macro = frames["source_fold_metrics"]
    aa = frames["within_aa_metrics"]
    selected_aa = aa.loc[(aa.metric == "chi_entropy") & (aa.fold >= 0) & aa.dataset.isin(SOURCES)]
    if sorted(selected_aa.amino_acid.unique()) != AA_ORDER:
        raise ValueError("Unexpected amino-acid categories")
    records, missing = [], []
    for source in SOURCES:
        for metric, predictor, display, x in CONTROL_ITEMS:
            group = macro.loc[(macro.dataset == source) & (macro.fold >= 0)
                              & (macro.scope == "macro_by_system")
                              & (macro.metric == metric) & (macro.predictor == predictor)].sort_values("fold")
            if group.fold.tolist() != list(range(5)):
                raise ValueError("Missing or duplicate macro fold")
            for _, row in group.iterrows():
                records.append(dict(dataset=source, block="system_macro", metric=metric,
                                    condition=predictor, x=x, fold=int(row.fold),
                                    spearman=float(row.spearman), count=int(row["count"]),
                                    n_systems=int(row.n_systems), source_table="source_fold_metrics.csv"))
        for letter, x in zip(AA_ORDER, AA_X):
            group = selected_aa.loc[(selected_aa.dataset == source)
                                    & (selected_aa.amino_acid == letter)].sort_values("fold")
            if group.empty:
                missing.append(dict(dataset=source, amino_acid=letter, x=float(x),
                                    reason="No matching valid-residue stratum in the verified source table; not a zero correlation."))
                continue
            if group.fold.tolist() != list(range(5)):
                raise ValueError("Unexpected partially missing amino-acid folds")
            for _, row in group.iterrows():
                records.append(dict(dataset=source, block="within_aa_pooled", metric="chi_entropy",
                                    condition=letter, x=float(x), fold=int(row.fold),
                                    spearman=float(row.spearman), count=int(row["count"]),
                                    n_systems=None, source_table="within_aa_metrics.csv"))
    data = pd.DataFrame(records)
    if not np.isfinite(data.spearman).all() or not data.spearman.between(-1, 1).all():
        raise ValueError("Invalid correlation values")
    if data.duplicated(["dataset", "block", "metric", "condition", "fold"]).any():
        raise ValueError("Repeated scientific observation")
    if len(data) != 330 or len(missing) != 3:
        raise ValueError("Unexpected evidence inventory; review rather than silently remove rows")
    if {(r["dataset"], r["amino_acid"]) for r in missing} != {
            ("dynarepo", "H"), ("dynarepo", "I"), ("itrace", "I")}:
        raise ValueError("Missing-stratum inventory changed")
    groups = data.groupby(["dataset", "block", "metric", "condition", "x"], sort=False).spearman
    stats = groups.agg(mean="mean", sd="std", folds="count").reset_index()
    if len(stats) != 66 or not (stats.folds == 5).all():
        raise ValueError("Every displayed interval must summarize five folds")
    return data, stats, missing, summary


def compose_figure(data, stats, missing):
    """One categorical axis with explicitly separated statistical estimands."""
    mpl.rcParams.update({"xtick.labelsize": 8, "ytick.labelsize": 8})
    fig, ax = plt.subplots(figsize=(183 / 25.4, 95 / 25.4))
    fig.subplots_adjust(left=.09, right=.982, bottom=.235, top=.775)
    ax.set_xlim(-.8, 26.8)
    ax.set_ylim(-.075, 1.035)
    ax.set_ylabel("Spearman ρ", labelpad=6)
    ax.set_yticks([0, .25, .5, .75, 1.0], ["0", "0.25", "0.50", "0.75", "1.00"])
    positions = [item[3] for item in CONTROL_ITEMS] + list(AA_X)
    labels = [item[2] for item in CONTROL_ITEMS] + AA_ORDER
    ax.set_xticks(positions, labels)
    ax.tick_params(axis="x", length=0, pad=5)
    ax.tick_params(axis="y", length=3, pad=3)
    ax.axhline(0, c="#CCCCCC", lw=.55, ls=(0, (3, 3)), zorder=0)
    # Shared numeric units do not imply the same averaging scheme.
    ax.axvline(7.95, color="#C7CDD1", lw=.65, ls=(0, (3, 3)), zorder=0)
    for source_index, source in enumerate(SOURCES):
        offset = (source_index - 1) * .24
        part = data.loc[data.dataset == source]
        fold_jitter = (part.fold.to_numpy() - 2) * .027
        ax.scatter(part.x + offset + fold_jitter, part.spearman,
                   s=4, alpha=.42, color=COLORS[source], marker=MARKERS[source],
                   linewidths=0, zorder=2)
        group = stats.loc[stats.dataset == source]
        ax.errorbar(group.x + offset, group["mean"], yerr=group.sd,
                    color=COLORS[source], marker=MARKERS[source], markersize=3.1,
                    ls="none", elinewidth=.85, markeredgewidth=.35,
                    capsize=1.5, capthick=.65, zorder=3)
    # Brackets identify blocks without restoring four independent panels.
    trans = ax.get_xaxis_transform()
    for lo, hi, title, y in [(-.55, 7.3, "System-macro correlations", 1.18),
                              (8.55, 26.45, "Hχ within each amino acid (pooled residues)", 1.18),
                              (-.55, 2.15, "RMSF", 1.055),
                              (2.95, 7.3, "Hχ", 1.055)]:
        ax.plot([lo, hi], [y - .025, y - .025], transform=trans,
                color="#777777", lw=.55, clip_on=False)
        ax.text((lo + hi) / 2, y, title, transform=trans,
                ha="center", va="bottom", fontsize=8)
    handles = [Line2D([], [], color=COLORS[s], marker=MARKERS[s], lw=0,
                      ms=4, label=LABELS[s]) for s in SOURCES]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(.53, .994),
               ncol=3, handletextpad=.4, columnspacing=2, frameon=False)
    # Undefined combinations are outside the numeric plot, never assigned y=0.
    ax.text(8.1, -.22, "No valid labels:", transform=trans, fontsize=8,
            va="center", ha="right", color="#666666")
    for entry in missing:
        j = SOURCES.index(entry["dataset"])
        ax.plot(entry["x"] + (j - 1) * .24, -.22, transform=trans,
                marker=MARKERS[entry["dataset"]], ms=3, mfc="white",
                mec=COLORS[entry["dataset"]], mew=.8, ls="none", clip_on=False)
    # A short direct label explains the real near-zero group, not a data deletion.
    p_x = AA_X[AA_ORDER.index("P")]
    ax.annotate("P: near-constant labels", xy=(p_x, .008), xytext=(p_x + 1.1, .14),
                fontsize=8, color="#555555", ha="left", va="bottom",
                arrowprops={"arrowstyle": "-", "lw": .55, "color": "#AAAAAA"})
    fig.text(.09, .055, "Five folds: mean ± SD; faint marks: individual folds. Open marks below axis: no valid labels.", fontsize=8)
    fig.text(.09, .011, "AA = amino-acid type. AA-adj. = same model; prediction and target minus training-fold AA means.", fontsize=8)
    return fig


def export_and_audit(fig, stats, data, missing, summary):
    import fitz
    from PIL import Image
    OUT.mkdir(parents=True, exist_ok=True)
    stats.to_csv(OUT / "source_data_mean_sd.csv", index=False)
    data.to_csv(OUT / "source_data_all_folds.csv", index=False)
    fig.canvas.draw()
    # Check explicit text artists before export, including text wholly off-page.
    renderer = fig.canvas.get_renderer()
    checked_text = [t for t in fig.findobj(mpl.text.Text) if t.get_visible() and t.get_text()]
    outside = []
    for artist in checked_text:
        box = artist.get_window_extent(renderer)
        if box.x0 < -.5 or box.y0 < -.5 or box.x1 > fig.bbox.x1 + .5 or box.y1 > fig.bbox.y1 + .5:
            outside.append(artist.get_text())
    if outside:
        raise ValueError("Text outside fixed canvas: " + repr(outside))
    for ext in ["pdf", "svg", "png"]:
        with mpl.rc_context({"savefig.bbox": None}):
            fig.savefig(OUT / f"{STEM}.{ext}", bbox_inches=None, dpi=600,
                        facecolor="white", metadata={"Creator": "FlexPPI reproducible plotting"}
                        if ext == "pdf" else None)
    plt.close(fig)
    doc = fitz.open(OUT / f"{STEM}.pdf")
    page = doc[0]
    spans = [span for block in page.get_text("dict")["blocks"] if "lines" in block
             for line in block["lines"] for span in line["spans"] if span["text"].strip()]
    fonts = [dict(name=f[3], embedded_bytes=len(doc.extract_font(f[0])[3])) for f in page.get_fonts()]
    svg = (OUT / f"{STEM}.svg").read_text(encoding="utf-8")
    with Image.open(OUT / f"{STEM}.png") as im:
        png_size, png_dpi = im.size, im.info["dpi"]
    report = dict(model_run=summary["model_run"], status="complete",
                  layout="One categorical axis, separated system-macro and within-AA pooled blocks",
                  width_mm=page.rect.width * 25.4 / 72, height_mm=page.rect.height * 25.4 / 72,
                  old_height_mm=149, area_reduction_fraction=1 - 95 / 149,
                  unique_fold_observations=len(data), mean_sd_groups=len(stats),
                  duplicate_Hchi_raw_observations_removed=15, unique_observations_removed=0,
                  missing_strata=missing, mean="equal-fold arithmetic mean",
                  spread="sample SD, ddof=1, five folds", statistical_test="none",
                  selected_epoch=40, inference="one held-out model per system; no new inference/retraining",
                  baseline="AA target means estimated on four training folds only",
                  source_analysis_sha256=sha256(DATA / "analysis_summary.json"),
                  source_tables={n: sha256(DATA / n) for n in ["source_fold_metrics.csv", "within_aa_metrics.csv"]},
                  plot_script_sha256=sha256(Path(__file__)),
                  minimum_pdf_font_pt=min(s["size"] for s in spans), fonts=fonts,
                  svg_raster_images=svg.count("<image"), svg_editable_text="<text" in svg,
                  png_size=png_size, png_dpi=png_dpi, outside_canvas_text=outside)
    report["mechanical_qa_passed"] = (len(doc) == 1 and report["minimum_pdf_font_pt"] >= 7.9
        and all(f["embedded_bytes"] > 0 for f in fonts) and report["svg_raster_images"] == 0
        and report["svg_editable_text"] and abs(report["width_mm"] - 183) < .05
        and abs(report["height_mm"] - 95) < .05 and min(png_dpi) >= 599)
    doc.close()
    (OUT / "statistics_and_QA.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    if not report["mechanical_qa_passed"]:
        raise ValueError("Export QA failed")
    print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    data, stats, missing, summary = prepare_data()
    figure = compose_figure(data, stats, missing)
    export_and_audit(figure, stats, data, missing, summary)
