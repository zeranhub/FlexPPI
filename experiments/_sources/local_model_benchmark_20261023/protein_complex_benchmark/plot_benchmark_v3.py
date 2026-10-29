# Academic Figure Skill Typography Baseline 鈥?COPY VERBATIM, place at TOP of script
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


# Academic Figure Skill Export Baseline 鈥?COPY VERBATIM
mpl.rcParams.update({
    "pdf.fonttype": 42,         # TrueType font embedding
    "svg.fonttype": "none",     # editable text in SVG
    "savefig.bbox": "tight",    # trim whitespace
    "savefig.dpi": 300,
})

"""Same two rows as v2, but every structure-consuming model reads the whole
assembly. BackFlip is shown in its official whole-complex mode, which is the
only mode that receives the same input information as FlexPPI; the chain-split
run stays in figures_v2 as an ablation. Scores are the frozen per-record values.
"""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.transforms import Bbox

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "figures_v3"
MEMBERSHIP = ROOT / "classification" / "membership_frozen.csv"
PER_RECORD = ROOT / "evaluation" / "per_record_metrics.csv"
SUMMARY = ROOT / "evaluation" / "group_summary.csv"

METHODS = ["Full", "backflip-2.1 (assembly)", "backflip-1.0 (assembly)", "DynaProt-M",
           "FastProtFlex", "ANM", "GNM", "Static", "AA"]
ABLATION = ["Full", "backflip-2.1 (chains)", "backflip-1.0 (chains)"]
DISPLAY = {"Full": "FlexPPI",
           "backflip-2.1 (assembly)": "BackFlip-2.1 \u00b7 assembly",
           "backflip-1.0 (assembly)": "BackFlip-1.0 \u00b7 assembly",
           "DynaProt-M": "DynaProt-M",
           "FastProtFlex": "FastProtFlex",
           "ANM": "ANM", "GNM": "GNM",
           "Static": "Static features", "AA": "Residue-type mean"}
BLUE = "#31688E"; OCHRE = "#B78038"; TEAL = "#348980"; INK = "#253846"
RELATIVE = {"FastProtFlex", "ANM", "GNM"}
METRICS = ["spearman", "top10_ap", "z_mae"]
SPECS = [("spearman", "Spearman \u03c1 \u2191", (-0.1, 1.0), [0, .5, 1]),
         ("top10_ap", "High-RMSF average precision \u2191", (0, 1), [0, .5, 1]),
         ("z_mae", "Standardized MAE \u2193", (0, 1.1), [0, .5, 1.0])]
SEED = 20261024
BOOTSTRAP = 2000
ROWS = [("multichain_protein", "Multichain proteins"), ("protein_complex", "Protein complexes")]


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as fh:
        for block in iter(lambda: fh.read(1024 ** 2), b""):
            h.update(block)
    return h.hexdigest()


def save_exact_canvas(fig, stem):
    box = Bbox.from_bounds(0, 0, *fig.get_size_inches())
    for suffix in ("pdf", "png", "svg"):
        fig.savefig(stem.with_suffix("." + suffix), bbox_inches=box, dpi=300)


def visual(name):
    if name == "Full":
        return BLUE, "s", 3.6
    if name in ("ANM", "GNM", "Static", "AA"):
        return TEAL, "^", 3.1
    return OCHRE, "o", 3.0


def family_bootstrap(frame, methods):
    records = sorted(frame[frame.method == "Full"].sample_id)
    fam = frame[frame.method == "Full"].set_index("sample_id").loc[records, "family"]
    names = sorted(fam.unique())
    pos = {f: i for i, f in enumerate(names)}
    codes = np.array([pos[x] for x in fam], dtype=np.int32)
    sizes = np.bincount(codes, minlength=len(names))
    rng = np.random.default_rng(SEED)
    if len(names) > 1:
        weights = rng.multinomial(len(names), np.full(len(names), 1 / len(names)), size=BOOTSTRAP)
        denom = weights @ sizes
    else:
        weights = denom = None
    out = {}
    for method in methods:
        block = frame[frame.method == method].set_index("sample_id").loc[records]
        for metric in METRICS:
            values = block[metric].to_numpy(float)
            if not np.isfinite(values).any():
                out[(method, metric)] = (np.nan, np.nan, np.nan, 0)
                continue
            assert np.isfinite(values).all(), f"non-finite {metric} for {method}"
            lo = hi = np.nan
            if weights is not None:
                sums = np.bincount(codes, weights=values, minlength=len(names))
                draws = (weights @ sums) / denom
                lo, hi = np.quantile(draws, [.025, .975])
            out[(method, metric)] = (float(values.mean()), float(lo), float(hi), len(values))
    return out, len(names)


def load():
    mem = pd.read_csv(MEMBERSHIP, encoding="utf-8-sig")
    per = pd.read_csv(PER_RECORD)
    summ = pd.read_csv(SUMMARY)
    multi_ids = mem.loc[(mem.benchmark_group == "protein") & (mem.actual_input_n_chains >= 2), "sample_id"]
    assert len(multi_ids) == 10, f"expected 10 multichain protein records, got {len(multi_ids)}"
    per_multi = per[per.sample_id.isin(set(multi_ids))]
    assert per_multi.groupby("method").size().eq(10).all()

    points, counts, fams = [], {}, {}
    boot, nfam = family_bootstrap(per_multi, METHODS)
    counts["multichain_protein"] = 10
    fams["multichain_protein"] = nfam
    for method in METHODS:
        for metric in METRICS:
            m, lo, hi, n = boot[(method, metric)]
            points.append(dict(row="multichain_protein", method=method, metric=metric,
                               mean=m, ci_low=lo, ci_high=hi, n_records=n, n_families=nfam))

    sub = summ[(summ.group == "protein_complex") & (summ.method.isin(METHODS)) & (summ.metric.isin(METRICS))]
    counts["protein_complex"] = int(sub.n_records.max())
    fams["protein_complex"] = int(sub.n_families.max())
    assert counts["protein_complex"] == 860 and fams["protein_complex"] == 261
    for row in sub.itertuples():
        points.append(dict(row="protein_complex", method=row.method, metric=row.metric,
                           mean=row.mean, ci_low=row.ci_low, ci_high=row.ci_high,
                           n_records=row.n_records, n_families=row.n_families))
    pts = pd.DataFrame(points)
    assert len(pts) == 54, f"expected 54 cells (standardized error for every method), got {len(pts)}"
    assert pts.loc[pts.metric == "z_mae", "mean"].notna().all(), "standardized MAE must be defined for every method"

    # ablation numbers kept beside the main output, not drawn
    OUT.mkdir(exist_ok=True)
    abl = per_multi[per_multi.method.isin(ABLATION)].pivot_table(
        index="sample_id", columns="method", values=METRICS)
    abl.to_csv(OUT / "chain_split_ablation_n10.csv")
    return dict(points=pts, counts=counts, families=fams,
                hashes={str(p): digest(p) for p in (MEMBERSHIP, PER_RECORD, SUMMARY)})


def draw(data):
    pts = data["points"]
    fig, axes = plt.subplots(2, 3, figsize=(183 / 25.4, 138 / 25.4), sharey=True)
    fig.subplots_adjust(left=.245, right=.985, bottom=.215, top=.90, wspace=.20, hspace=.55)
    for ri, (row, title) in enumerate(ROWS):
        position = axes[ri, 0].get_position()
        fig.text(.012, position.y1 + .040, title, fontsize=7, fontweight="bold", color=INK)
        fig.text(.985, position.y1 + .040,
                 f'n = {data["counts"][row]:,} records \u00b7 {data["families"][row]:,} families',
                 fontsize=6.5, ha="right", color=INK)
        for ci, (metric, label, limits, ticks) in enumerate(SPECS):
            ax = axes[ri, ci]
            for yi, method in enumerate(METHODS):
                if metric == "mae" and method in RELATIVE:
                    ax.text(.97, yi, "relative score", transform=ax.get_yaxis_transform(),
                            ha="right", va="center", fontsize=6.3, color="#6C797F")
                    continue
                r = pts[(pts.row == row) & (pts.method == method) & (pts.metric == metric)].iloc[0]
                color, marker, size = visual(method)
                if np.isfinite(r.ci_low) and np.isfinite(r.ci_high):
                    ax.plot([r.ci_low, r.ci_high], [yi, yi], color=color, lw=.9,
                            solid_capstyle="round", zorder=2)
                ax.plot(r["mean"], yi, marker=marker, ms=size, color=color, mec="white",
                        mew=.3, ls="", zorder=3)
            ax.axhspan(-.45, .45, color=BLUE, alpha=.055, lw=0, zorder=0)
            ax.set_ylim(len(METHODS) - .45, -.65)
            ax.set_yticks(range(len(METHODS)), [DISPLAY[m] for m in METHODS], fontsize=6.5)
            ax.tick_params(axis="y", length=0, pad=5, labelleft=ci == 0)
            ax.tick_params(axis="x", labelsize=6.3, length=2.5, pad=2)
            ax.spines["left"].set_visible(False)
            ax.spines["bottom"].set_color("#80909A")
            ax.set_xlim(*limits); ax.set_xticks(ticks)
            ax.set_xlabel(label, fontsize=6.7, labelpad=4)
            ax.text(-.06, 1.0, chr(97 + ri * 3 + ci), transform=ax.transAxes,
                    fontsize=8.5, fontweight="bold", color=INK)
    handles = [Line2D([], [], color=BLUE, marker="o", lw=0, ms=3.2, label="FlexPPI"),
               Line2D([], [], color=OCHRE, marker="o", lw=0, ms=3.2, label="Published learned models"),
               Line2D([], [], color=TEAL, marker="o", lw=0, ms=3.2, label="Graph / simple controls")]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(.5, .012), ncol=3,
               fontsize=6.3, handletextpad=.45, columnspacing=1.6)
    fig.text(.012, .112, "All structure-consuming models read the whole assembly, including BackFlip",
             fontsize=6.3, color=INK)
    fig.text(.012, .088, "Multichain row excludes single-chain inputs; complexes unchanged",
             fontsize=6.3, color=INK)
    fig.text(.012, .064, "Points: mean across records; intervals: 95% family-cluster bootstrap.",
             fontsize=6.3, color=INK)
    OUT.mkdir(exist_ok=True)
    stem = OUT / "benchmark_multichain_v4"
    save_exact_canvas(fig, stem)
    plt.close(fig)
    pts.to_csv(OUT / "benchmark_multichain_v4_source.csv", index=False)
    return stem


def main():
    data = load()
    stem = draw(data)
    for suffix in ("pdf", "png", "svg"):
        assert stem.with_suffix("." + suffix).is_file()
    assert all(digest(Path(p)) == h for p, h in data["hashes"].items())
    audit = dict(generated_utc=datetime.now(timezone.utc).isoformat(),
                 script_sha256=digest(__file__), input_sha256=data["hashes"],
                 methods=METHODS, rows=dict(ROWS), row_counts=data["counts"],
                 row_families=data["families"], cells=54,
                 error_metric="standardized MAE (z-MAE): per-record standardised absolute error, defined for every method; no Angstrom MAE is claimed for relative-output models",
                 input_protocol="whole assembly for every structure-consuming method",
                 ablation="BackFlip chain-split per-record values in chain_split_ablation_n10.csv and figures_v2/",
                 protein_row_definition="benchmark_group=protein with actual_input_n_chains>=2",
                 complex_row_source="frozen evaluation/group_summary.csv, unchanged",
                 statistics=f"record-weighted mean, 95% family-cluster percentile, seed={SEED}, draws={BOOTSTRAP}",
                 no_rescoring=True, no_score_based_selection=True,
                 artifact_sha256={f"benchmark_multichain_v4.{s}": digest(stem.with_suffix("." + s))
                                  for s in ("pdf", "png", "svg")})
    (OUT / "plot_audit_v4.json").write_text(json.dumps(audit, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(dict(counts=data["counts"], families=data["families"]), ensure_ascii=False))


if __name__ == "__main__":
    main()

