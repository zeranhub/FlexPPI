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

"""Two biological strata, three metrics, seven methods.

Every structure-consuming method receives the whole assembly; the chain-split
BackFlip run is reported separately. Values are the frozen per-record scores.
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
SUMMARY = ROOT / "evaluation" / "group_summary.csv"

# (column key, display name, family, marker, face)
METHODS = [
    ("Full", "FlexPPI", "own", "s"),
    ("backflip-2.1 (assembly)", "BackFlip-2.1", "learned", "o"),
    ("backflip-1.0 (assembly)", "BackFlip-1.0", "learned", "o"),
    ("DynaProt-M", "DynaProt-M", "learned", "D"),
    ("FastProtFlex", "FastProtFlex", "learned", "^"),
    ("ANM", "ANM", "geometry", "v"),
    ("GNM", "GNM", "geometry", "v"),
]
FAMILY_COLOUR = {"own": "#1F4E79", "learned": "#C0842B", "geometry": "#3F8F83"}
METRICS = [("spearman", "Spearman \u03c1 \u2191", (-0.05, 1.0), [0, 0.5, 1.0]),
           ("top10_ap", "High-RMSF AP \u2191", (0, 1), [0, 0.5, 1.0]),
           ("z_mae", "Standardized MAE \u2193", (0, 0.72), [0, 0.2, 0.4, 0.6])]
ROWS = [("multichain_protein", "Multichain proteins"), ("protein_complex", "Protein complexes")]
SEED = 20261024
BOOTSTRAP = 2000
INK = "#22313A"
HAIR = "#C9D2D6"


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as fh:
        for block in iter(lambda: fh.read(1024 ** 2), b""):
            h.update(block)
    return h.hexdigest()


def family_bootstrap(frame, methods):
    records = sorted(frame[frame.method == "Full"].sample_id)
    fam = frame[frame.method == "Full"].set_index("sample_id").loc[records, "family"]
    names = sorted(fam.unique())
    pos = {f: i for i, f in enumerate(names)}
    codes = np.array([pos[x] for x in fam], dtype=np.int32)
    sizes = np.bincount(codes, minlength=len(names))
    rng = np.random.default_rng(SEED)
    weights = rng.multinomial(len(names), np.full(len(names), 1 / len(names)), size=BOOTSTRAP)
    denom = weights @ sizes
    out = {}
    for method in methods:
        block = frame[frame.method == method].set_index("sample_id").loc[records]
        for metric, _, _, _ in METRICS:
            values = block[metric].to_numpy(float)
            assert np.isfinite(values).all(), f"non-finite {metric} for {method}"
            sums = np.bincount(codes, weights=values, minlength=len(names))
            draws = (weights @ sums) / denom
            out[(method, metric)] = (float(values.mean()),
                                     float(np.quantile(draws, .025)),
                                     float(np.quantile(draws, .975)))
    return out, len(names)


def load():
    mem = pd.read_csv(MEMBERSHIP, encoding="utf-8-sig")
    per = pd.read_csv(PER_RECORD)
    summ = pd.read_csv(SUMMARY)
    ids = list(mem.loc[(mem.benchmark_group == "protein") & (mem.actual_input_n_chains >= 2), "sample_id"])
    assert len(ids) == 10
    frame = per[per.sample_id.isin(set(ids))]
    keys = [m[0] for m in METHODS]
    assert frame.groupby("method").size().eq(10).all()

    boot, nfam = family_bootstrap(frame, keys)
    points, counts, fams = [], {"multichain_protein": 10, "protein_complex": 860}, \
        {"multichain_protein": nfam, "protein_complex": 261}
    for method in keys:
        for metric, _, _, _ in METRICS:
            m, lo, hi = boot[(method, metric)]
            points.append(dict(row="multichain_protein", method=method, metric=metric,
                               mean=m, ci_low=lo, ci_high=hi, n_records=10, n_families=nfam))
    sub = summ[(summ.group == "protein_complex") & (summ.method.isin(keys))
               & (summ.metric.isin([m for m, _, _, _ in METRICS]))]
    assert int(sub.n_records.max()) == 860 and int(sub.n_families.max()) == 261
    for r in sub.itertuples():
        points.append(dict(row="protein_complex", method=r.method, metric=r.metric,
                           mean=r.mean, ci_low=r.ci_low, ci_high=r.ci_high,
                           n_records=860, n_families=261))
    pts = pd.DataFrame(points)
    assert len(pts) == 42, f"expected 42 cells, got {len(pts)}"
    assert pts["mean"].notna().all()
    return dict(points=pts, counts=counts, families=fams,
                hashes={str(p): digest(p) for p in (MEMBERSHIP, PER_RECORD, SUMMARY)})


def draw(data):
    pts = data["points"]
    fig, axes = plt.subplots(2, 3, figsize=(183 / 25.4, 104 / 25.4), sharey=True)
    fig.subplots_adjust(left=.205, right=.988, bottom=.175, top=.845, wspace=.14, hspace=.50)
    n_method = len(METHODS)

    for ri, (row, title) in enumerate(ROWS):
        top = axes[ri, 0].get_position().y1
        fig.text(.008, top + .052, f"{title}  (n = {data['counts'][row]:,})", fontsize=7.4, fontweight="bold", color=INK)
        for ci, (metric, label, limits, ticks) in enumerate(METRICS):
            ax = axes[ri, ci]
            ax.set_facecolor("white")
            # subtle band behind the proposed model
            ax.axhspan(-0.5, 0.5, color="#1F4E79", alpha=.045, lw=0, zorder=0)
            for yi, (key, name, family, marker) in enumerate(METHODS):
                r = pts[(pts.row == row) & (pts.method == key) & (pts.metric == metric)].iloc[0]
                colour = FAMILY_COLOUR[family]
                ax.errorbar(r["mean"], yi, xerr=[[r["mean"] - r.ci_low], [r.ci_high - r["mean"]]],
                            fmt="none", ecolor=colour, elinewidth=.9, capsize=1.6,
                            capthick=.9, zorder=2)
                ax.plot(r["mean"], yi, marker=marker, ms=4.2 if family == "own" else 3.6,
                        color=colour, mec="white", mew=.5, ls="", zorder=3)
            ax.set_ylim(n_method - .45, -.55)
            ax.set_xlim(*limits)
            ax.set_xticks(ticks)
            ax.set_yticks(range(n_method), [m[1] for m in METHODS], fontsize=6.7)
            for tick, m in zip(ax.get_yticklabels(), METHODS):
                if m[2] == "own":
                    tick.set_fontweight("bold")
                    tick.set_color(FAMILY_COLOUR["own"])
            ax.tick_params(axis="y", length=0, pad=6, labelleft=(ci == 0))
            ax.tick_params(axis="x", labelsize=6.2, length=2.4, pad=2, colors="#3C4A52")
            ax.spines["left"].set_visible(False)
            ax.spines["bottom"].set_color("#93A1A8")
            ax.spines["bottom"].set_linewidth(.6)
            if ci == 0:
                for sep in (0.5, 4.5):
                    ax.axhline(sep, color=HAIR, lw=.6, ls="-", zorder=1)
            if ri == 0:
                ax.text(.5, 1.16, label, transform=ax.transAxes, ha="center", va="bottom",
                        fontsize=6.9, color=INK)
            else:
                ax.set_xlabel(label, fontsize=6.6, labelpad=3.5, color=INK)
            ax.text(-.03, 1.06, chr(97 + ri * 3 + ci), transform=ax.transAxes,
                    fontsize=8.2, fontweight="bold", color=INK, va="bottom")
    handles = [Line2D([], [], color=FAMILY_COLOUR["own"], marker="s", lw=0, ms=3.8, label="FlexPPI"),
               Line2D([], [], color=FAMILY_COLOUR["learned"], marker="o", lw=0, ms=3.4,
                      label="Published learned models"),
               Line2D([], [], color=FAMILY_COLOUR["geometry"], marker="v", lw=0, ms=3.4,
                      label="Elastic-network baselines")]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(.5, .005), ncol=3,
               fontsize=6.3, handletextpad=.55, columnspacing=2.0)
    fig.text(.008, .068, "Every method reads the full multi-chain coordinates; points are record means, bars are 95% family-cluster bootstrap intervals",
             fontsize=6.1, color="#4A5A62")
    OUT.mkdir(exist_ok=True)
    stem = OUT / "benchmark_v5"
    box = Bbox.from_bounds(0, 0, *fig.get_size_inches())
    for suffix in ("pdf", "png", "svg"):
        fig.savefig(stem.with_suffix("." + suffix), bbox_inches=box, dpi=300)
    plt.close(fig)
    pts.to_csv(OUT / "benchmark_v5_source.csv", index=False)
    return stem


def main():
    data = load()
    stem = draw(data)
    assert all(digest(Path(p)) == h for p, h in data["hashes"].items())
    audit = dict(generated_utc=datetime.now(timezone.utc).isoformat(),
                 script_sha256=digest(__file__), input_sha256=data["hashes"],
                 methods=[m[1] for m in METHODS], rows=dict(ROWS),
                 row_counts=data["counts"], row_families=data["families"], cells=42,
                 input_protocol="whole assembly for every structure-consuming method",
                 protein_row_definition="benchmark_group=protein with actual_input_n_chains>=2",
                 complex_row_source="frozen evaluation/group_summary.csv",
                 statistics=f"record-weighted mean, 95% family-cluster percentile, seed={SEED}, draws={BOOTSTRAP}",
                 no_rescoring=True, no_score_based_selection=True,
                 artifacts={f"benchmark_v5.{s}": digest(stem.with_suffix("." + s)) for s in ("pdf", "png", "svg")})
    (OUT / "plot_audit_v5.json").write_text(json.dumps(audit, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(dict(counts=data["counts"], families=data["families"]), ensure_ascii=False))


if __name__ == "__main__":
    main()


