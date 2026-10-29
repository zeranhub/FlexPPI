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
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
    "savefig.bbox": "tight",
    "savefig.dpi": 300,
})

"""Paired difference against FlexPPI, whole-assembly protocol.

Every comparator reads the full multi-chain coordinates. A row whose interval
covers zero is a tie by construction; a row entirely on the FlexPPI side is a
lead. This is the companion to benchmark_multichain_v3.
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
OUT = ROOT / "figures_v5"
MEMBERSHIP = ROOT / "classification" / "membership_frozen.csv"
PER_RECORD = ROOT / "evaluation" / "per_record_metrics.csv"
PAIRED = ROOT / "evaluation" / "paired_differences.csv"

COMPARATORS = ["backflip-2.1 (assembly)", "backflip-1.0 (assembly)", "DynaProt-M",
               "FastProtFlex", "ANM", "GNM"]
DISPLAY = {"backflip-2.1 (assembly)": "BackFlip-2.1",
           "backflip-1.0 (assembly)": "BackFlip-1.0",
           "DynaProt-M": "DynaProt-M", "FastProtFlex": "FastProtFlex",
           "ANM": "ANM", "GNM": "GNM", "Static": "Static features", "AA": "Residue-type mean"}
RELATIVE = {"FastProtFlex", "ANM", "GNM"}
METRICS = [("spearman", "Spearman \u03c1: FlexPPI \u2212 competitor", True),
           ("top10_ap", "High-RMSF AP: FlexPPI \u2212 competitor", True),
           ("z_mae", "Standardized MAE: FlexPPI \u2212 competitor", False)]
BLUE = "#31688E"; BROWN = "#A6611A"; NEUTRAL = "#7C8B94"; INK = "#253846"
SEED = 20261024
BOOTSTRAP = 2000
ROWS = [("multichain_protein", "Multichain proteins"), ("protein_complex", "Protein complexes")]


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as fh:
        for block in iter(lambda: fh.read(1024 ** 2), b""):
            h.update(block)
    return h.hexdigest()


def boot_row(frame, comparators):
    records = sorted(frame[frame.method == "Full"].sample_id)
    fam = frame[frame.method == "Full"].set_index("sample_id").loc[records, "family"]
    names = sorted(fam.unique())
    pos = {f: i for i, f in enumerate(names)}
    codes = np.array([pos[x] for x in fam], dtype=np.int32)
    sizes = np.bincount(codes, minlength=len(names))
    rng = np.random.default_rng(SEED)
    weights = rng.multinomial(len(names), np.full(len(names), 1 / len(names)), size=BOOTSTRAP)
    denom = weights @ sizes
    out = []
    for comp in comparators:
        for metric, _, _ in METRICS:
            a = frame[frame.method == "Full"].set_index("sample_id").loc[records, metric].to_numpy(float)
            b = frame[frame.method == comp].set_index("sample_id").loc[records, metric].to_numpy(float)
            if not (np.isfinite(a).all() and np.isfinite(b).all()):
                out.append(dict(group="multichain_protein", comparator=comp, metric=metric,
                                delta=np.nan, lo=np.nan, hi=np.nan, n=len(records), n_families=len(names)))
                continue
            d = a - b
            sums = np.bincount(codes, weights=d, minlength=len(names))
            draws = (weights @ sums) / denom
            out.append(dict(group="multichain_protein", comparator=comp, metric=metric,
                            delta=float(d.mean()), lo=float(np.quantile(draws, .025)),
                            hi=float(np.quantile(draws, .975)), n=len(records), n_families=len(names)))
    return out


def load():
    mem = pd.read_csv(MEMBERSHIP, encoding="utf-8-sig")
    per = pd.read_csv(PER_RECORD)
    paired = pd.read_csv(PAIRED)
    ids = list(mem.loc[(mem.benchmark_group == "protein") & (mem.actual_input_n_chains >= 2), "sample_id"])
    assert len(ids) == 10
    frame = per[per.sample_id.isin(set(ids))]

    rows = boot_row(frame, COMPARATORS)
    sel = paired[(paired.group == "protein_complex") & (paired.metric.isin([m for m, _, _ in METRICS]))]
    for comp in COMPARATORS:
        for metric, _, _ in METRICS:
            hit = sel[(sel.comparator == comp) & (sel.metric == metric)]
            if hit.empty:
                rows.append(dict(group="protein_complex", comparator=comp, metric=metric,
                                 delta=np.nan, lo=np.nan, hi=np.nan, n=860, n_families=261))
            else:
                r = hit.iloc[0]
                rows.append(dict(group="protein_complex", comparator=comp, metric=metric,
                                 delta=float(r.full_minus_comparator), lo=float(r.ci_low),
                                 hi=float(r.ci_high), n=int(r.n_records), n_families=int(r.n_families)))
    df = pd.DataFrame(rows)
    df["verdict"] = np.where(df.lo.isna(), "n/a",
                    np.where((df.lo > 0) | (df.hi < 0), "separated", "tie"))
    return df


def colour(metric, lo, hi):
    if not np.isfinite(lo):
        return NEUTRAL
    if lo <= 0 <= hi:
        return NEUTRAL
    favours_full = (hi < 0) if metric in ("mae", "z_mae") else (lo > 0)
    return BLUE if favours_full else BROWN


def draw(df):
    fig, axes = plt.subplots(2, 3, figsize=(183 / 25.4, 124 / 25.4), sharey=True)
    fig.subplots_adjust(left=.215, right=.985, bottom=.205, top=.885, wspace=.20, hspace=.52)
    limits = {}
    for metric, _, _ in METRICS:
        lo = float(np.nanmin(df.loc[df.metric == metric, "lo"]))
        hi = float(np.nanmax(df.loc[df.metric == metric, "hi"]))
        pad = 0.08 * (hi - lo)
        limits[metric] = (lo - pad, hi + pad)
    for ri, (grp, title) in enumerate(ROWS):
        for ci, (metric, label, higher_better) in enumerate(METRICS):
            ax = axes[ri, ci]
            ax.axvline(0, color="#5C6B73", lw=.8, zorder=1)
            for yi, comp in enumerate(COMPARATORS):
                r = df[(df.group == grp) & (df.comparator == comp) & (df.metric == metric)]
                if r.empty or not np.isfinite(r.iloc[0].delta):
                    ax.text(.5, yi, "\u2014", transform=ax.get_yaxis_transform(),
                            ha="center", va="center", fontsize=7, color="#A8B4BA")
                    continue
                r = r.iloc[0]
                c = colour(metric, r.lo, r.hi)
                ax.plot([r.lo, r.hi], [yi, yi], color=c, lw=1.15, solid_capstyle="round", zorder=2)
                ax.plot(r.delta, yi, "o", ms=3.3, color=c, mec="white", mew=.3, zorder=3)
            ax.set_ylim(len(COMPARATORS) - .4, -.6)
            ax.set_yticks(range(len(COMPARATORS)), [DISPLAY[c] for c in COMPARATORS], fontsize=6.4)
            ax.tick_params(axis="y", length=0, pad=5, labelleft=ci == 0)
            ax.tick_params(axis="x", labelsize=6.2, length=2.5, pad=2)
            ax.spines["left"].set_visible(False)
            ax.spines["bottom"].set_color("#80909A")
            ax.set_xlim(*limits[metric])
            ax.set_xlabel(label, fontsize=6.4, labelpad=4)
            ax.text(-.06, 1.0, chr(97 + ri * 3 + ci), transform=ax.transAxes,
                    fontsize=8.5, fontweight="bold", color=INK)
        axes[ri, 0].figure.text(.012, axes[ri, 0].get_position().y1 + .045, title,
                                fontsize=7, fontweight="bold", color=INK)
        n = int(df[df.group == grp].n.max()); k = int(df[df.group == grp].n_families.max())
        axes[ri, 0].figure.text(.985, axes[ri, 0].get_position().y1 + .045,
                                f"n = {n:,} records \u00b7 {k:,} families",
                                fontsize=6.5, ha="right", color=INK)
    handles = [Line2D([], [], color=BLUE, marker="o", lw=1.1, ms=3.3, label="FlexPPI better (interval excludes 0)"),
               Line2D([], [], color=NEUTRAL, marker="o", lw=1.1, ms=3.3, label="tie (interval covers 0)"),
               Line2D([], [], color=BROWN, marker="o", lw=1.1, ms=3.3, label="competitor better")]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(.5, .030), ncol=3,
               fontsize=6.3, handletextpad=.5, columnspacing=1.8)
    fig.text(.012, .118, "Whole-assembly input for every method; paired record-wise differences, 95% family-cluster bootstrap",
             fontsize=6.2, color=INK)
    fig.text(.012, .094, "Chain-split BackFlip rows are reported separately and are not shown here",
             fontsize=6.2, color=INK)
    OUT.mkdir(exist_ok=True)
    stem = OUT / "benchmark_paired_vs_flexppi_v5"
    box = Bbox.from_bounds(0, 0, *fig.get_size_inches())
    for suffix in ("pdf", "png", "svg"):
        fig.savefig(stem.with_suffix("." + suffix), bbox_inches=box, dpi=300)
    plt.close(fig)
    df.to_csv(OUT / "benchmark_paired_vs_flexppi_v5_source.csv", index=False)
    return stem


def main():
    df = load()
    stem = draw(df)
    (OUT / "paired_plot_audit_v5.json").write_text(json.dumps(dict(
        generated_utc=datetime.now(timezone.utc).isoformat(),
        script_sha256=digest(__file__),
        inputs={str(p): digest(p) for p in (MEMBERSHIP, PER_RECORD, PAIRED)},
        comparators=COMPARATORS, metrics=[m for m, _, _ in METRICS],
        rows=dict(ROWS), statistics=f"paired family-cluster bootstrap, seed={SEED}, draws={BOOTSTRAP}",
        complex_row_source="frozen evaluation/paired_differences.csv",
        multichain_row_source="recomputed from per_record_metrics.csv on the frozen 10 records",
        artifacts={f"benchmark_paired_vs_flexppi_v5.{s}": digest(stem.with_suffix("." + s)) for s in ("pdf", "png", "svg")}
    ), indent=2, ensure_ascii=False), encoding="utf-8")
    print(df.verdict.value_counts().to_string())


if __name__ == "__main__":
    main()





