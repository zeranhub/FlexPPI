# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a) RMSF curves -> LineTrend/plot_posttraining.py -> param inherit
# (b) molecular structures -> no molecular asset; load authentic PyMOL renders
# The LineTrend asset has embedded benchmark arrays, categorical steps and
# opacity ramps, not dense residue profiles. Inherit its clean spines and
# unframed legend; use the user's reference colors and unsmoothed real data.
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
import json
import math
from pathlib import Path

mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator
import numpy as np
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
MM = 1 / 25.4
MD_COLOR = "#39348B"
PRED_COLOR = "#18B878"
STRUCTURE_CMAP = LinearSegmentedColormap.from_list(
    "rmsf_structure", ["#2D4EA1", "#F7F7F7", "#FF3026"]
)
DATASET_NAMES = {"atlas": "ATLAS", "dynarepo": "DynaRepo", "itrace": "ITrace"}


def absolute_path(value):
    path = Path(value)
    return path if path.is_absolute() else ROOT / path


def read_cases(path):
    obj = json.loads(path.read_text(encoding="utf-8-sig"))
    cases = obj if isinstance(obj, list) else obj["cases"]
    if not cases:
        raise ValueError(f"No cases supplied: {path}")
    return cases


def metric_column(array, index, n, name):
    array = np.asarray(array)
    if array.ndim != 2 or array.shape[0] != n or index >= array.shape[1]:
        raise ValueError(f"Invalid {name} shape: {array.shape}")
    return array[:, index]


def load_profile(case):
    """Keep every residue of the requested chain; missing labels remain gaps."""
    path = absolute_path(case["cache_path"])
    with np.load(path, allow_pickle=False) as data:
        chains = data["chain"].astype(str)
        metric_names = data["metric_names"].astype(str).tolist()
        if "rmsf" not in metric_names:
            raise ValueError(f"RMSF is missing from {path}")
        index = metric_names.index("rmsf")
        n = len(chains)
        selected = np.flatnonzero(chains == str(case["chain"]))
        if not len(selected):
            raise ValueError(f"Chain {case['chain']} not found in {path}")
        pred = metric_column(data["mean"], index, n, "mean")[selected].astype(float)
        target = metric_column(data["target"], index, n, "target")[selected].astype(float)
        valid = metric_column(data["valid"], index, n, "valid")[selected].astype(bool)
        variance = metric_column(data["label_variance"], index, n, "label_variance")[selected].astype(float)
        author_id = data["author_residue_id"].astype(str)[selected]
        residue_index = data["residue_index"][selected]
    if valid.sum() < 2:
        raise ValueError(f"Fewer than two valid residues for {case['case_id']}")
    for name, values in [("MD", target), ("prediction", pred), ("label variance", variance)]:
        if not np.isfinite(values[valid]).all():
            raise ValueError(f"Non-finite {name} at valid residues: {case['case_id']}")
        if (values[valid] < -1e-8).any():
            raise ValueError(f"Negative {name}: {case['case_id']}")
    # A predicted value at an invalid MD label cannot enter a comparison metric.
    if np.std(target[valid]) == 0 or np.std(pred[valid]) == 0:
        pearson = float("nan")
    else:
        pearson = float(np.corrcoef(target[valid], pred[valid])[0, 1])
    mae = float(np.abs(pred[valid] - target[valid]).mean())
    for key, computed in [("pearson", pearson), ("mae", mae)]:
        if key in case and case[key] is not None:
            expected = float(case[key])
            if not np.isclose(expected, computed, atol=5e-5, equal_nan=True):
                raise ValueError(f"{key} mismatch for {case['case_id']}: {expected} vs {computed}")
    return dict(case=case, x=np.arange(1, len(selected) + 1), pred=pred, target=target,
                valid=valid, variance=variance, author_id=author_id,
                residue_index=residue_index, cache_row=selected, pearson=pearson, mae=mae,
                n_valid=int(valid.sum()), n_invalid=int((~valid).sum()), path=str(path))


def display_name(case):
    explicit = case.get("display_name")
    if explicit:
        return str(explicit)
    sample = str(case["sample_id"])
    # Remove only literal dataset prefixes, keeping identifiers traceable.
    for prefix in (str(case.get("dataset", "")) + "_", str(case.get("dataset", "")) + ":"):
        if sample.lower().startswith(prefix.lower()) and len(prefix) > 1:
            sample = sample[len(prefix):]
            break
    chain = str(case['chain'])
    return sample if chain == '_' else f"{sample} · chain {chain}"


def dataset_name(case):
    name = str(case.get("dataset", ""))
    return DATASET_NAMES.get(name.lower(), name)


def format_r(value):
    return f"{value:.2f}" if np.isfinite(value) else "undefined"


def curve_panel(ax, profile, *, label_sd=False, compact=False, ylabel=True):
    valid = profile["valid"]
    x = profile["x"]
    target = np.where(valid, profile["target"], np.nan)
    pred = np.where(valid, profile["pred"], np.nan)
    sd = np.sqrt(np.maximum(profile["variance"], 0))
    if label_sd:
        lower = np.where(valid, np.maximum(target - sd, 0), np.nan)
        upper = np.where(valid, target + sd, np.nan)
        ax.fill_between(x, lower, upper, color=MD_COLOR, alpha=.16, linewidth=0)
    ax.plot(x, target, color=MD_COLOR, lw=.70, zorder=3, label="MD")
    ax.plot(x, pred, color=PRED_COLOR, lw=.80, zorder=4, label="FlexPPI-v3")
    max_y = float(max(np.nanmax(target + sd) if label_sd else np.nanmax(target), np.nanmax(pred)))
    top = max(.5, max_y * 1.08)
    ax.set_ylim(0, top)
    ax.set_xlim(.5, len(x) + .5)
    ax.xaxis.set_major_locator(MaxNLocator(nbins=4 if compact else 5, integer=True, min_n_ticks=3))
    ax.yaxis.set_major_locator(MaxNLocator(nbins=4, min_n_ticks=3))
    ax.set_xlabel("Residue (chain position)", fontsize=6.5 if compact else 7.5, labelpad=2)
    if ylabel:
        ax.set_ylabel("Local RMSF [Å]", fontsize=7 if compact else 8, labelpad=3)
    ax.tick_params(axis="both", labelsize=6 if compact else 7, pad=2, length=2.5)
    ax.spines["left"].set_color("#3A3A43")
    ax.spines["bottom"].set_color("#3A3A43")
    return dict(case_id=profile["case"]["case_id"], n_valid=profile["n_valid"],
                n_invalid=profile["n_invalid"], pearson=profile["pearson"], mae=profile["mae"],
                y_max=top, label_sd=label_sd, clipped_curve_points=0,
                md_range=[float(np.nanmin(target)), float(np.nanmax(target))],
                prediction_range=[float(np.nanmin(pred)), float(np.nanmax(pred))])


def structure_panel(ax, profile, structure_dir, role):
    path = structure_dir / f"{profile['case']['case_id']}_{role}.png"
    if not path.is_file():
        raise FileNotFoundError(f"Required authentic structure render is missing: {path}")
    with Image.open(path) as source:
        rgb = source.convert("RGBA")
        ax.imshow(rgb)
    ax.axis("off")
    return str(path)


def shared_legend(fig, *, y, label_sd=False, fontsize=8):
    handles = [Line2D([], [], color=PRED_COLOR, lw=.9, label="FlexPPI-v3"),
               Line2D([], [], color=MD_COLOR, lw=.8, label="MD")]
    if label_sd:
        handles.append(Patch(facecolor=MD_COLOR, alpha=.16, linewidth=0, label="MD label ± SD"))
    fig.legend(handles=handles, ncol=len(handles), loc="upper center",
               bbox_to_anchor=(.54, y), frameon=False, fontsize=fontsize,
               handlelength=2.5, columnspacing=2.2)


def colorbar(fig, bbox, vmax):
    ax = fig.add_axes(bbox)
    bar = fig.colorbar(mpl.cm.ScalarMappable(norm=Normalize(0, vmax), cmap=STRUCTURE_CMAP),
                       cax=ax, orientation="horizontal", ticks=[0, vmax / 2, vmax])
    bar.outline.set_visible(False)
    bar.ax.tick_params(length=0, pad=2, labelsize=7)
    bar.set_label("Local RMSF [Å]", fontsize=8, labelpad=1)


def export_exact_size(fig, output_dir, stem):
    # A fixed canvas preserves the actual 183-mm print width. The skill baseline's
    # tight bounding box would otherwise change it as title lengths change.
    paths = []
    with mpl.rc_context({"savefig.bbox": None}):
        for extension in ("pdf", "svg", "png"):
            path = output_dir / f"{stem}.{extension}"
            fig.savefig(path, dpi=300, bbox_inches=None, facecolor="white")
            paths.append(str(path))
    plt.close(fig)
    return paths


def main_figure(profiles, structure_dir, output_dir, vmax):
    if len(profiles) != 3:
        raise ValueError(f"The main figure needs exactly 3 representative cases, got {len(profiles)}")
    fig = plt.figure(figsize=(183 * MM, 180 * MM), facecolor="white")
    fig.text(.52, .987, "Representative held-out predictions", ha="center", va="top", fontsize=9)
    shared_legend(fig, y=.966)
    lefts = [.095, .402, .709]
    audit = []
    for profile, left in zip(profiles, lefts):
        ax = fig.add_axes([left, .701, .245, .180])
        item = curve_panel(ax, profile)
        ax.set_title(display_name(profile["case"]), fontsize=8.2, pad=20)
        ax.text(.5, 1.045,
                f"{dataset_name(profile['case'])}   r = {format_r(profile['pearson'])}   MAE = {profile['mae']:.2f} Å",
                transform=ax.transAxes, ha="center", va="bottom", fontsize=6.6)
        for bottom, role in [(.382, "md"), (.085, "prediction")]:
            sx = fig.add_axes([left - .013, bottom, .286, .276])
            item[role + "_render"] = structure_panel(sx, profile, structure_dir, role)
        audit.append(item)
    for letter, name, y in [("A", "Residue profiles", .922), ("B", "MD", .658), ("C", "FlexPPI-v3", .361)]:
        fig.text(.021, y, letter, fontsize=12, fontweight="bold", va="top")
        if name != "Residue profiles":
            fig.text(.095, y, name, fontsize=8.5, va="top")
    colorbar(fig, [.37, .058, .30, .014], vmax)
    paths = export_exact_size(fig, output_dir, "fig_v3_representative_profiles_structures")
    return {"paths": paths, "cases": audit, "dimensions_mm": [183, 180]}


def showcase_figure(profiles, structure_dir, output_dir, vmax, label_sd):
    if len(profiles) != 8:
        raise ValueError(f"The showcase gallery needs exactly 8 cases, got {len(profiles)}")
    fig = plt.figure(figsize=(183 * MM, 174 * MM), facecolor="white")
    fig.text(.52, .99, "Selected high-agreement examples", ha="center", va="top", fontsize=9)
    shared_legend(fig, y=.976, label_sd=label_sd, fontsize=7.5)
    lefts = [.078, .319, .560, .801]
    audit = []
    for i, profile in enumerate(profiles):
        left = lefts[i % 4]
        top = i < 4
        sy = .711 if top else .063
        cy = .535 if top else .267
        ty = .917 if top else .459
        sx = fig.add_axes([left - .019, sy, .220, .185])
        render = structure_panel(sx, profile, structure_dir, "prediction")
        ax = fig.add_axes([left, cy, .183, .135])
        item = curve_panel(ax, profile, compact=True, label_sd=label_sd, ylabel=i % 4 == 0)
        fig.text(left + .0915, ty, display_name(profile["case"]), ha="center", fontsize=7.3)
        fig.text(left + .0915, ty - .016, dataset_name(profile["case"]), ha="center", fontsize=6.2)
        ax.text(.5, 1.035, f"r = {format_r(profile['pearson'])}   MAE = {profile['mae']:.2f} Å",
                transform=ax.transAxes, fontsize=6.3, ha="center", va="bottom")
        item["prediction_render"] = render
        audit.append(item)
    colorbar(fig, [.385, .041, .270, .010], vmax)
    paths = export_exact_size(fig, output_dir, "fig_v3_high_agreement_gallery")
    return {"paths": paths, "cases": audit, "dimensions_mm": [183, 174]}


def all_profiles_pdf(cases, output_dir, label_sd):
    """Full coverage in stable input order: no selection by model agreement."""
    path = output_dir / "all_test_chain_profiles.pdf"
    rows = []
    with PdfPages(path) as pdf:
        for page in range(math.ceil(len(cases) / 12)):
            fig, axes = plt.subplots(4, 3, figsize=(183 * MM, 235 * MM))
            fig.subplots_adjust(left=.09, right=.98, top=.92, bottom=.075, hspace=.88, wspace=.41)
            fig.text(.52, .985, "All held-out chain profiles", fontsize=10, ha="center", va="top")
            shared_legend(fig, y=.97, label_sd=label_sd, fontsize=7.5)
            batch = cases[page * 12:(page + 1) * 12]
            for ax, case in zip(axes.flat, batch):
                profile = load_profile(case)
                rows.append(curve_panel(ax, profile, compact=True, label_sd=label_sd))
                ax.set_title(display_name(case), fontsize=7.2, pad=16)
                ax.text(.5, 1.035,
                        f"{dataset_name(case)}   r = {format_r(profile['pearson'])}   MAE = {profile['mae']:.2f} Å",
                        ha="center", transform=ax.transAxes, fontsize=5.7, va="bottom")
            for ax in list(axes.flat)[len(batch):]:
                ax.axis("off")
            fig.text(.52, .025, f"Page {page + 1} / {math.ceil(len(cases) / 12)} · all valid residues, no smoothing",
                     fontsize=7, ha="center")
            with mpl.rc_context({"savefig.bbox": None}):
                pdf.savefig(fig, dpi=300, bbox_inches=None)
            plt.close(fig)
            if (page + 1) % 20 == 0:
                print(f"Full profile PDF: page {page + 1}", flush=True)
    return {"path": str(path), "pages": math.ceil(len(cases) / 12), "n_chains": len(cases), "cases": rows}


def source_csv(profiles, output_dir):
    path = output_dir / "figure_source_residues.csv"
    columns = ["case_id", "selection_group", "sample_id", "dataset", "fold", "chain",
               "chain_position", "cache_row", "residue_index", "author_residue_id", "valid_rmsf",
               "md_rmsf_A", "prediction_rmsf_A", "md_label_variance_A2", "md_label_sd_A", "cache_path"]
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for profile in profiles:
            case = profile["case"]
            for j, x in enumerate(profile["x"]):
                row = {k: case.get(k, "") for k in columns[:6]}
                row.update(chain_position=int(x), cache_row=int(profile["cache_row"][j]),
                           residue_index=int(profile["residue_index"][j]),
                           author_residue_id=profile["author_id"][j], valid_rmsf=bool(profile["valid"][j]),
                           md_rmsf_A=float(profile["target"][j]), prediction_rmsf_A=float(profile["pred"][j]),
                           md_label_variance_A2=float(profile["variance"][j]),
                           md_label_sd_A=float(np.sqrt(max(0, profile["variance"][j]))), cache_path=profile["path"])
                writer.writerow(row)
    return str(path)


def overall_figure(metrics_path, output_dir):
    """Empirical cumulative distributions retain every system with a valid metric."""
    import pandas as pd
    data = pd.read_csv(metrics_path)
    data = data.loc[data.metric == 'rmsf']
    if len(data) != 3068 or data.sample_id.duplicated().any():
        raise ValueError('Overview requires all 3,068 unique RMSF test systems')
    fig, axes = plt.subplots(1, 2, figsize=(183*MM, 84*MM))
    fig.subplots_adjust(left=.10, right=.985, bottom=.21, top=.79, wspace=.34)
    fig.text(.53, .985, 'Complete held-out RMSF test set · 3,068 protein samples',
             ha='center', va='top', fontsize=9)
    styles = [('atlas', '#2166AC', '-'), ('dynarepo', '#C78124', '--'),
              ('itrace', '#762A83', '-.')]
    handles, audit = [], []
    for dataset, color, linestyle in styles:
        group = data.loc[data.dataset == dataset]
        handles.append(Line2D([], [], color=color, ls=linestyle,
                       label=f'{DATASET_NAMES[dataset]} (n={len(group):,})'))
        for ax, metric in zip(axes, ['pearson', 'mae']):
            values = group[metric].to_numpy(float)
            if not np.isfinite(values).all():
                raise ValueError('Undefined test statistic in complete overview')
            values = np.sort(values)
            ax.step(values, np.arange(1, len(values)+1)/len(values), where='post',
                    color=color, linestyle=linestyle, linewidth=1)
            audit.append(dict(dataset=dataset, statistic=metric, count=len(values),
                              median=float(np.median(values)), minimum=float(values[0]),
                              maximum=float(values[-1])))
    r_min = min(0.0, float(np.floor(data.pearson.min()*10)/10))
    axes[0].set(xlim=(r_min, 1), xlabel='Per-sample Pearson r')
    axes[1].set(xscale='log', xlabel='Per-sample MAE [Å] · log scale')
    for ax, letter in zip(axes, ['A', 'B']):
        ax.set(ylim=(0, 1.02), ylabel='Fraction of samples ≤ value')
        ax.text(-.17, 1.04, letter, transform=ax.transAxes, fontsize=10)
    fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(.53, .93),
               ncol=3, fontsize=7.5)
    fig.text(.53, .035, 'All samples contribute once; curves describe the full distribution of test performance.',
             fontsize=6.6, ha='center')
    return dict(paths=export_exact_size(fig, output_dir, 'fig_v3_all_test_performance'), statistics=audit)


def write_report(output_dir, cases_path, audit, profiles, vmax, label_sd):
    stats = []
    for profile in profiles:
        valid = profile["valid"]
        stats.append({"case": profile["case"], "n_valid_rmsf_residues": profile["n_valid"],
                      "n_invalid_rmsf_residues_shown_as_gaps": profile["n_invalid"],
                      "pearson": profile["pearson"], "mae_A": profile["mae"],
                      "md_residues_above_color_max": int((profile["target"][valid] > vmax).sum()),
                      "prediction_residues_above_color_max": int((profile["pred"][valid] > vmax).sum())})
    audit["statistics"] = stats
    audit["source_cases"] = str(cases_path)
    audit["structure_color_limits_A"] = [0, vmax]
    audit["visual_review_status"] = "Rendered files require visual inspection; automated checks are not a substitute."
    (output_dir / "figure_statistics_and_audit.json").write_text(json.dumps(audit, indent=2, ensure_ascii=False), encoding="utf-8")
    total = sum(p["n_valid"] + p["n_invalid"] for p in profiles)
    invalid = sum(p["n_invalid"] for p in profiles)
    text = f"""# FlexPPI-v3 test figure captions and reproducibility notes

The main figure pairs the residue-wise MD RMSF labels and FlexPPI-v3 held-out predictions with static structures colored by each quantity. Three representative cases were selected by the analysis script using its recorded per-dataset chain-correlation rule. The gallery contains eight explicitly selected high-agreement examples, and must not be interpreted as the average performance of the benchmark. Overall conclusions must use the complete test-set metrics. The source selection file is `{cases_path}`; all recorded selection fields are included in the statistics JSON.

**Main figure caption.** (A) Original, unsmoothed per-residue RMSF values from MD labels (purple) and FlexPPI-v3 predictions (green). Pearson r and mean absolute error (MAE, Å) are computed over paired valid RMSF residues in each displayed chain. (B) Input structure colored by MD RMSF. (C) The same input structure, camera and geometry colored by predicted RMSF. These are static color mappings, not predicted conformations or an MD ensemble. The same structure color scale, 0–{vmax:g} Å, is used throughout. Values above {vmax:g} Å use the maximum color; the exact counts are documented in the statistics JSON. Curves are never clipped to that color scale.

**Gallery caption.** Eight selected high-agreement held-out examples, each showing a predicted-RMSF structure color map and the paired residue curves. All source cases and their selection ranks are recorded. {'The purple band denotes MD label mean ± the square root of the cached MD label variance, with the lower boundary truncated at zero. It is label standard deviation, not a confidence interval or predictive uncertainty.' if label_sd else 'No uncertainty band is drawn.'} The colorbar applies to structure colors only; curve axes adapt to show all values.

**Complete-test overview caption.** Empirical cumulative distributions of (A) per-sample Pearson correlation and (B) per-sample MAE, with one observation from every held-out protein sample. ATLAS includes 1,929 samples, DynaRepo 895 and ITrace 244. Higher correlations shift a distribution to the right in A; lower errors shift it to the left in B. The MAE axis is logarithmic and includes the entire observed range. No samples are excluded from either panel, and every source is labelled with its sample count.

**Statistics.** Here n is the number of residues with valid paired MD labels and predictions in a chain, not the number of independent proteins or independent replicates. Pearson r compares the two residue profiles; MAE is the arithmetic mean of absolute predicted-minus-MD deviations. No hypothesis test or p-value is used, and no multiple-comparison correction applies. The model prediction is the cached mean output; the MD target and optional label variance are read directly from the evaluation cache. Fold identifiers are preserved for every case. The precise train/validation/test audit and checkpoint provenance are reported by the evaluation analysis; these plots do not create a fresh split or refit the model.

**Source data.** The selected cases contain {total:,} cache residues, including {invalid:,} invalid RMSF label positions. Every selected-chain residue is retained in `figure_source_residues.csv`; invalid positions are visible as line gaps and do not enter paired metrics. The horizontal axis is the 1-based position within the cached chain. The original author residue identifier, cache residue index, cache row, chain and sample ID are retained alongside every plotted value. No residue downsampling, curve smoothing or normalization is applied. Cases outside the illustration panels remain part of the overall evaluation; when `--all-cases` is supplied, `all_test_chain_profiles.pdf` shows every supplied chain in input order.

**Design and QA.** The reference illustrations determine the purple/green comparison palette and curve-plus-cartoon arrangement. The academic LineTrend asset was inspected but contains categorical synthetic benchmark steps; its clean spines and unframed legend were inherited rather than adapting the benchmark data. Typography, palette and export baselines are present in the script. Figure width is exactly 183 mm; previews are 300 dpi, and PDF/SVG preserve vector lines and editable text while embedding authentic PyMOL raster renders. Structure color stops match the renderer exactly: #2D4EA1, #F7F7F7, #FF3026. Gray structure residues indicate missing mapped values. Automated checks reject invalid shapes, nonfinite values at valid residues, negative RMSF/variance, missing structure assets and inconsistencies with the analysis metrics. Bounds include all line data. Visual inspection of final PNGs must confirm text separation, image crops and colorbar placement; this document does not claim that inspection has occurred automatically.
"""
    (output_dir / "figure_captions_and_QA.md").write_text(text, encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description="Create real held-out RMSF curve/structure figures.")
    parser.add_argument("--cases", type=Path, default=ROOT / "artifacts/evaluation_v3/analysis/selected_cases.json")
    parser.add_argument("--structure-dir", type=Path, default=ROOT / "artifacts/evaluation_v3/structure_renders")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "output/flexppi_v3_test_figures")
    parser.add_argument("--all-cases", type=Path, help="Optional JSON of every chain for a complete profile PDF.")
    parser.add_argument("--system-metrics", type=Path,
                        default=ROOT / "artifacts/evaluation_v3/analysis/per_system_metrics.csv")
    parser.add_argument("--label-sd", action="store_true", help="Show genuine cached MD label SD in gallery/PDF.")
    parser.add_argument("--vmax", type=float, default=5.0, help="Must equal the structure renderer color maximum.")
    args = parser.parse_args()
    if not args.vmax > 0:
        parser.error("--vmax must be positive")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    cases = read_cases(args.cases)
    profiles = [load_profile(case) for case in cases]
    representative = [p for p in profiles if p["case"]["selection_group"] == "representative"]
    showcase = [p for p in profiles if p["case"]["selection_group"] == "showcase"]
    audit = {"representative": main_figure(representative, args.structure_dir, args.output_dir, args.vmax),
             "showcase": showcase_figure(showcase, args.structure_dir, args.output_dir, args.vmax, args.label_sd)}
    audit["source_csv"] = source_csv(profiles, args.output_dir)
    audit['overall'] = overall_figure(args.system_metrics, args.output_dir)
    if args.all_cases:
        audit["full_profiles"] = all_profiles_pdf(read_cases(args.all_cases), args.output_dir, args.label_sd)
    write_report(args.output_dir, args.cases, audit, profiles, args.vmax, args.label_sd)
    print(f"Generated RMSF figures in {args.output_dir}", flush=True)


if __name__ == "__main__":
    main()
