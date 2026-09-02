# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a) apo absolute RMSF profile -> assets/figures/LineTrend/plot_posttraining.py (cross-type inherit) -> param inherit
# (b) holo absolute RMSF profile -> assets/figures/LineTrend/plot_posttraining.py (cross-type inherit) -> param inherit
# (c) delta RMSF profile -> assets/figures/LineTrend/plot_posttraining.py (cross-type inherit) -> param inherit
# (d) apo absolute H_chi profile -> assets/figures/LineTrend/plot_posttraining.py (cross-type inherit) -> param inherit
# (e) holo absolute H_chi profile -> assets/figures/LineTrend/plot_posttraining.py (cross-type inherit) -> param inherit
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
import json
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
from scipy.stats import spearmanr

# Keep the tight-bounding-box export within the 183 mm double-column contract.
mpl.rcParams["savefig.pad_inches"] = 0.02


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = ROOT / "artifacts" / "apo_holo_validation" / "multi_system_8" / "residue_predictions.csv"
DEFAULT_OVERVIEW = ROOT / "artifacts" / "apo_holo_validation" / "multi_system_8" / "complex_overview.csv"
DEFAULT_OUTPUT = ROOT / "artifacts" / "apo_holo_validation" / "multi_system_8" / "per_complex_profiles"

PDB_ORDER = ("1AK4", "1AY7", "1BRS", "1EMV", "1GCQ", "1PVH", "2OOB", "3SGB")
SHORT_NAMES = {
    "1AK4": "CypA–HIV-1 capsid NTD",
    "1AY7": "RNase Sa–barstar",
    "1BRS": "barnase–barstar",
    "1EMV": "colicin E9 DNase–Im9",
    "1GCQ": "Vav SH3–Grb2 SH3",
    "1PVH": "LIF–gp130",
    "2OOB": "Cbl-b UBA–ubiquitin",
    "3SGB": "SGPB–ovomucoid inhibitor",
}

SIDE_COLORS = {"A": CATEGORICAL[0], "B": CATEGORICAL[3]}
TRUE_STYLE = {"linestyle": "-", "linewidth": 1.25, "alpha": 0.98}
PRED_STYLE = {"linestyle": "--", "linewidth": 1.05, "alpha": 0.92}
BOUNDARY_COLOR = "#BDBDBD"
ANNOTATION_COLOR = "#555555"

USED_COLUMNS = (
    "pdb_id", "side", "pair_id", "overlap_class", "apo_residue_index",
    "holo_residue_index", "author_residue_id", "amino_acid",
    "true_apo_rmsf", "true_holo_rmsf", "true_delta_rmsf",
    "flexppi_apo_rmsf", "flexppi_holo_full_rmsf",
    "flexppi_delta_response_rmsf", "true_apo_chi_entropy",
    "true_holo_chi_entropy", "flexppi_apo_chi_entropy",
    "flexppi_holo_full_chi_entropy",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Draw per-complex apo/holo RMSF, delta RMSF, and H_chi profiles."
    )
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--overview", type=Path, default=DEFAULT_OVERVIEW)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def finite_pair(true_values: Iterable[float], pred_values: Iterable[float]) -> Tuple[np.ndarray, np.ndarray]:
    true_array = np.asarray(true_values, dtype=float)
    pred_array = np.asarray(pred_values, dtype=float)
    valid = np.isfinite(true_array) & np.isfinite(pred_array)
    return true_array[valid], pred_array[valid]


def agreement_stats(true_values: Iterable[float], pred_values: Iterable[float]) -> Dict[str, float]:
    true_array, pred_array = finite_pair(true_values, pred_values)
    n = int(true_array.size)
    if n < 2 or np.nanstd(true_array) == 0 or np.nanstd(pred_array) == 0:
        rho = float("nan")
    else:
        rho = float(spearmanr(true_array, pred_array).statistic)
    mae = float(np.mean(np.abs(true_array - pred_array))) if n else float("nan")
    return {"n": n, "spearman_rho": rho, "mae": mae}


def format_metric_annotation(stats: Mapping[str, float]) -> str:
    rho = stats["spearman_rho"]
    rho_text = "NA" if not np.isfinite(rho) else f"{rho:.2f}"
    return f"ρ = {rho_text}   MAE = {stats['mae']:.2f}   n = {int(stats['n'])}"


def add_concat_index(data: pd.DataFrame) -> pd.DataFrame:
    pieces: List[pd.DataFrame] = []
    offset = 0
    for side in ("A", "B"):
        side_data = data.loc[data["side"].eq(side)].sort_values("apo_residue_index").copy()
        side_data["matched_residue_concat_index"] = np.arange(len(side_data), dtype=int) + offset
        pieces.append(side_data)
        offset += len(side_data)
    return pd.concat(pieces, ignore_index=True)


def plot_profile_pair(
    ax: plt.Axes,
    data: pd.DataFrame,
    true_column: str,
    pred_column: str,
    title: str,
    y_label: str,
    y_limits: Tuple[float, float],
    mask_prediction_by_truth: bool = False,
    zero_line: bool = False,
) -> Dict[str, float]:
    if zero_line:
        ax.axhline(0.0, color=BLACK, linewidth=0.65, alpha=0.75, zorder=0)

    for side in ("A", "B"):
        side_data = data.loc[data["side"].eq(side)].sort_values("matched_residue_concat_index")
        x = side_data["matched_residue_concat_index"].to_numpy(dtype=float)
        truth = side_data[true_column].to_numpy(dtype=float)
        prediction = side_data[pred_column].to_numpy(dtype=float)
        if mask_prediction_by_truth:
            prediction = np.where(np.isfinite(truth), prediction, np.nan)
        ax.plot(x, truth, color=SIDE_COLORS[side], **TRUE_STYLE)
        ax.plot(x, prediction, color=SIDE_COLORS[side], **PRED_STYLE)

    stats = agreement_stats(data[true_column], data[pred_column])
    ax.set_title(title, loc="left", pad=4, fontweight="bold")
    ax.set_ylabel(y_label)
    ax.set_ylim(*y_limits)
    ax.set_xlim(-1, len(data))
    ax.text(
        0.99, 0.97, format_metric_annotation(stats), transform=ax.transAxes,
        ha="right", va="top", fontsize=6, color=ANNOTATION_COLOR,
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.78, "pad": 1.0},
    )
    return stats


def add_chain_context(ax: plt.Axes, n_a: int, n_total: int) -> None:
    boundary = n_a - 0.5
    ax.axvline(boundary, color=BOUNDARY_COLOR, linestyle=":", linewidth=0.75, zorder=0)
    a_center = (n_a - 1) / 2 if n_a else 0
    b_center = n_a + (n_total - n_a - 1) / 2 if n_total > n_a else n_a
    label_box = {"facecolor": "white", "edgecolor": "none", "alpha": 0.72, "pad": 0.7}
    ax.text(a_center, 0.87, "chain A", transform=ax.get_xaxis_transform(),
            ha="center", va="top", fontsize=6, color=SIDE_COLORS["A"], bbox=label_box)
    ax.text(b_center, 0.87, "chain B", transform=ax.get_xaxis_transform(),
            ha="center", va="top", fontsize=6, color=SIDE_COLORS["B"], bbox=label_box)


def shared_absolute_limits(data: pd.DataFrame, columns: Sequence[str]) -> Tuple[float, float]:
    values = data.loc[:, columns].to_numpy(dtype=float)
    finite = values[np.isfinite(values)]
    upper = float(np.max(finite)) if finite.size else 1.0
    return 0.0, max(0.2, upper * 1.08)


def symmetric_delta_limits(data: pd.DataFrame, columns: Sequence[str]) -> Tuple[float, float]:
    values = data.loc[:, columns].to_numpy(dtype=float)
    finite = values[np.isfinite(values)]
    magnitude = float(np.max(np.abs(finite))) if finite.size else 1.0
    magnitude = max(0.2, magnitude * 1.08)
    return -magnitude, magnitude


def panel_statistics(data: pd.DataFrame) -> Dict[str, Dict[str, float]]:
    return {
        "apo_rmsf": agreement_stats(data["true_apo_rmsf"], data["flexppi_apo_rmsf"]),
        "holo_rmsf": agreement_stats(data["true_holo_rmsf"], data["flexppi_holo_full_rmsf"]),
        "delta_rmsf": agreement_stats(data["true_delta_rmsf"], data["flexppi_delta_response_rmsf"]),
        "apo_hchi": agreement_stats(data["true_apo_chi_entropy"], data["flexppi_apo_chi_entropy"]),
        "holo_hchi": agreement_stats(data["true_holo_chi_entropy"], data["flexppi_holo_full_chi_entropy"]),
    }


def draw_complex_figure(
    data: pd.DataFrame,
    pdb_id: str,
    overview_row: Mapping[str, object],
) -> Tuple[plt.Figure, Dict[str, object]]:
    data = add_concat_index(data)
    n_a = int(data["side"].eq("A").sum())
    n_b = int(data["side"].eq("B").sum())
    n_total = n_a + n_b

    rmsf_limits = shared_absolute_limits(
        data,
        ("true_apo_rmsf", "flexppi_apo_rmsf", "true_holo_rmsf", "flexppi_holo_full_rmsf"),
    )
    hchi_limits = shared_absolute_limits(
        data,
        ("true_apo_chi_entropy", "flexppi_apo_chi_entropy",
         "true_holo_chi_entropy", "flexppi_holo_full_chi_entropy"),
    )
    delta_limits = symmetric_delta_limits(
        data, ("true_delta_rmsf", "flexppi_delta_response_rmsf")
    )

    mm_to_inch = 1 / 25.4
    fig = plt.figure(figsize=(183 * mm_to_inch, 160 * mm_to_inch))
    grid = fig.add_gridspec(
        nrows=3, ncols=2, height_ratios=(1.0, 1.18, 1.0),
        left=0.085, right=0.985, bottom=0.105, top=0.865,
        wspace=0.22, hspace=0.55,
    )
    axes = [
        fig.add_subplot(grid[0, 0]),
        fig.add_subplot(grid[0, 1]),
        fig.add_subplot(grid[1, :]),
        fig.add_subplot(grid[2, 0]),
        fig.add_subplot(grid[2, 1]),
    ]

    stats = panel_statistics(data)
    plot_profile_pair(
        axes[0], data, "true_apo_rmsf", "flexppi_apo_rmsf",
        "apo A/B — absolute RMSF", "RMSF (Å)", rmsf_limits,
    )
    plot_profile_pair(
        axes[1], data, "true_holo_rmsf", "flexppi_holo_full_rmsf",
        "holo AB — absolute RMSF", "RMSF (Å)", rmsf_limits,
    )
    plot_profile_pair(
        axes[2], data, "true_delta_rmsf", "flexppi_delta_response_rmsf",
        "ΔRMSF — holo AB minus apo A/B", "ΔRMSF (Å)", delta_limits,
        zero_line=True,
    )
    plot_profile_pair(
        axes[3], data, "true_apo_chi_entropy", "flexppi_apo_chi_entropy",
        "apo A/B — absolute $H_\\chi$", "$H_\\chi$ (nats)", hchi_limits,
        mask_prediction_by_truth=True,
    )
    plot_profile_pair(
        axes[4], data, "true_holo_chi_entropy", "flexppi_holo_full_chi_entropy",
        "holo AB — absolute $H_\\chi$", "$H_\\chi$ (nats)", hchi_limits,
        mask_prediction_by_truth=True,
    )

    for label, ax in zip("abcde", axes):
        add_chain_context(ax, n_a=n_a, n_total=n_total)
        ax.text(-0.11, 1.08, label, transform=ax.transAxes, fontsize=9,
                fontweight="bold", va="top", ha="left")
        ax.set_xlabel("Matched residue (chain A then chain B)")
        ax.tick_params(length=2.5)

    handles = [
        Line2D([0], [0], color=SIDE_COLORS["A"], linewidth=2.0, label="chain A"),
        Line2D([0], [0], color=SIDE_COLORS["B"], linewidth=2.0, label="chain B"),
        Line2D([0], [0], color=BLACK, linestyle="-", linewidth=1.25, label="MD truth"),
        Line2D([0], [0], color=BLACK, linestyle="--", linewidth=1.05, label="FlexPPI"),
    ]
    trajectory_ns = int(overview_row["trajectory_ns_per_state"])
    frames = 12_501 if trajectory_ns == 500 else 7_501
    fig.suptitle(
        f"{pdb_id}  {SHORT_NAMES[pdb_id]} — apo/holo residue profiles",
        x=0.5, y=0.976, fontsize=11, fontweight="bold",
    )
    fig.text(
        0.5, 0.943,
        f"A: {n_a} matched residues   B: {n_b} matched residues   "
        f"MD: {trajectory_ns} ns/state, {frames:,} frames/state",
        ha="center", va="center", fontsize=7, color=ANNOTATION_COLOR,
    )
    fig.legend(
        handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.925),
        ncol=4, handlelength=2.8, columnspacing=1.7,
    )
    fig.text(
        0.5, 0.025,
        "Solid lines: MD; dashed lines: FlexPPI. "
        "ΔRMSF uses the direct response head. Residues without a defined side-chain rotamer entropy are omitted.",
        ha="center", va="bottom", fontsize=6, color=ANNOTATION_COLOR,
    )

    output_stats: Dict[str, object] = {
        "pdb_id": pdb_id,
        "proteins": str(overview_row["proteins"]),
        "trajectory_ns_per_state": trajectory_ns,
        "frames_per_state": frames,
        "n_chain_A": n_a,
        "n_chain_B": n_b,
        "n_matched_residues": n_total,
        "n_hchi_defined": int(data["true_apo_chi_entropy"].notna().sum()),
        "panels": stats,
        "definitions": {
            "true_delta_rmsf": "MD holo RMSF minus MD apo RMSF for matched residues",
            "predicted_delta_rmsf": "FlexPPI direct binding-response head",
            "hchi": "natural-log Shannon entropy over side-chain rotamer states (nats)",
        },
    }
    return fig, output_stats


def flatten_stats(statistics: Sequence[Mapping[str, object]]) -> pd.DataFrame:
    rows: List[Dict[str, object]] = []
    for complex_stats in statistics:
        for panel_name, values in complex_stats["panels"].items():
            rows.append({
                "pdb_id": complex_stats["pdb_id"],
                "panel": panel_name,
                "n": values["n"],
                "spearman_rho": values["spearman_rho"],
                "mae": values["mae"],
            })
    return pd.DataFrame(rows)


def write_readme(output_dir: Path, statistics: Sequence[Mapping[str, object]]) -> None:
    lines = [
        "# Zenodo 6638504 per-complex apo/holo profiles",
        "",
        "Each PDB entry has a PDF and SVG vector master plus a 300-dpi PNG preview. "
        "The multipage PDF contains all eight vector figures.",
        "",
        "The x axis concatenates the matched residues of chain A followed by chain B. "
        "Absolute apo panels correspond to separately simulated apo A and apo B; holo panels "
        "correspond to the A and B chains within the independently simulated holo AB complex. "
        "Delta RMSF is holo minus apo for MD, while the FlexPPI curve is the model's direct "
        "binding-response head, matching the original validation plot convention.",
        "",
        "H_chi is the natural-log Shannon entropy of side-chain rotamer states. Predictions are "
        "masked wherever the MD H_chi target is undefined, so missing targets are not silently "
        "replaced by model output.",
        "",
        "## Data and statistics",
        "",
        "- Source data: `../residue_predictions.csv`",
        "- Figure-specific source data: `per_complex_profile_source_data.csv`",
        "- Panel metrics: `per_complex_profile_statistics.csv` and `.json`",
        "- Metrics shown in panels: per-residue Spearman rho and mean absolute error (MAE).",
        "- No hypothesis test or confidence interval is claimed; each residue is a sequence "
        "position, not an independent biological replicate.",
        "- MD sampling: one trajectory per state. The 7,501/12,501 saved frames are correlated "
        "time samples, not independent replicates.",
        "",
        "## Files",
        "",
    ]
    for item in statistics:
        pdb_id = item["pdb_id"]
        lines.append(f"- `{pdb_id}_apo_holo_profiles.[pdf|svg|png]` — {item['proteins']}")
    lines.extend([
        "- `zenodo_6638504_all_complex_profiles.pdf` — eight-page vector collection",
        "- `zenodo_6638504_profile_contact_sheet.png` — raster overview for quick inspection",
        "",
    ])
    (output_dir / "README.md").write_text("\n".join(lines), encoding="utf-8")


def make_contact_sheet(png_paths: Sequence[Path], output_path: Path) -> None:
    from PIL import Image, ImageOps, ImageDraw

    previews = [Image.open(path).convert("RGB") for path in png_paths]
    target_width = 1050
    scaled = []
    for image in previews:
        target_height = int(round(image.height * target_width / image.width))
        scaled.append(image.resize((target_width, target_height), Image.Resampling.LANCZOS))
    gutter = 36
    rows = 4
    cols = 2
    cell_height = max(image.height for image in scaled)
    sheet = Image.new(
        "RGB",
        (cols * target_width + (cols + 1) * gutter, rows * cell_height + (rows + 1) * gutter),
        "white",
    )
    draw = ImageDraw.Draw(sheet)
    for index, image in enumerate(scaled):
        row, col = divmod(index, cols)
        x = gutter + col * (target_width + gutter)
        y = gutter + row * (cell_height + gutter)
        framed = ImageOps.expand(image, border=1, fill="#DDDDDD")
        sheet.paste(framed, (x, y))
        draw.text((x + 8, y + 8), PDB_ORDER[index], fill="#222222")
    sheet.save(output_path, dpi=(300, 300))


def main() -> None:
    args = parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    residue = pd.read_csv(args.input, usecols=list(USED_COLUMNS))
    overview = pd.read_csv(args.overview).set_index("pdb_id")
    missing = [pdb_id for pdb_id in PDB_ORDER if pdb_id not in set(residue["pdb_id"])]
    if missing:
        raise ValueError(f"Missing expected complexes: {missing}")
    if not set(residue["side"]).issubset({"A", "B"}):
        raise ValueError("Only chain sides A and B are supported.")

    source_pieces: List[pd.DataFrame] = []
    statistics: List[Dict[str, object]] = []
    png_paths: List[Path] = []
    multipage_path = args.output / "zenodo_6638504_all_complex_profiles.pdf"

    with PdfPages(multipage_path) as multipage:
        for pdb_id in PDB_ORDER:
            complex_data = residue.loc[residue["pdb_id"].eq(pdb_id)].copy()
            source_pieces.append(add_concat_index(complex_data))
            figure, complex_stats = draw_complex_figure(
                complex_data, pdb_id=pdb_id, overview_row=overview.loc[pdb_id]
            )
            base = args.output / f"{pdb_id}_apo_holo_profiles"
            save_cns_figure(figure, base)
            figure.savefig(f"{base}.svg", bbox_inches="tight", dpi=300)
            multipage.savefig(figure, bbox_inches="tight")
            plt.close(figure)
            statistics.append(complex_stats)
            png_paths.append(Path(f"{base}.png"))

    source_data = pd.concat(source_pieces, ignore_index=True)
    source_data.to_csv(args.output / "per_complex_profile_source_data.csv", index=False)
    flatten_stats(statistics).to_csv(
        args.output / "per_complex_profile_statistics.csv", index=False
    )
    (args.output / "per_complex_profile_statistics.json").write_text(
        json.dumps(statistics, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    write_readme(args.output, statistics)
    make_contact_sheet(
        png_paths, args.output / "zenodo_6638504_profile_contact_sheet.png"
    )
    print(f"Wrote {len(statistics)} per-complex profile figures to {args.output}")


if __name__ == "__main__":
    main()
