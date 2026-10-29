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
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
    "savefig.bbox": "tight",
    "savefig.dpi": 300,
})

"""Same bar layout as v6, but the comparison is split by input composition:
single-chain proteins versus multichain inputs. This keeps the two groups
within a factor of 2.5 and avoids a 10-versus-860 panel.
"""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.transforms import Bbox

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "figures_v7"
MEMBERSHIP = ROOT / "classification" / "membership_frozen.csv"
PER_RECORD = ROOT / "evaluation" / "per_record_metrics.csv"

METHODS = [
    ("Full", "FlexPPI", "own"),
    ("backflip-2.1 (assembly)", "BackFlip-2.1", "learned"),
    ("backflip-1.0 (assembly)", "BackFlip-1.0", "learned"),
    ("DynaProt-M", "DynaProt-M", "learned"),
    ("FastProtFlex", "FastProtFlex", "learned"),
    ("ANM", "ANM", "geometry"),
    ("GNM", "GNM", "geometry"),
]
METRICS = [("spearman", "Spearman \u03c1 \u2191", (0, 1.0), [0, 0.5, 1.0]),
           ("top10_ap", "High-RMSF AP \u2191", (0, 1.0), [0, 0.5, 1.0]),
           ("z_mae", "Standardized MAE \u2193", (0, 0.72), [0, 0.2, 0.4, 0.6])]
ROWS = [("single_chain", "Single-chain proteins", "#A7C3D9", "#7D97A8"),
        ("multi_chain", "Multichain inputs", "#1F4E79", "#1B2C36")]
SEED = 20261024
BOOTSTRAP = 2000
INK = "#22313A"


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
            out[(method, metric)] = (float(values.mean()), float(np.quantile(draws, .025)),
                                     float(np.quantile(draws, .975)))
    return out, len(names)


def load():
    mem = pd.read_csv(MEMBERSHIP, encoding="utf-8-sig")
    per = pd.read_csv(PER_RECORD)
    m = mem[["sample_id", "actual_input_n_chains"]]
    df = per.merge(m, on="sample_id", how="left")
    df = df[df.method.isin([x[0] for x in METHODS])]
    df["row"] = np.where(df.actual_input_n_chains >= 2, "multi_chain", "single_chain")
    assert df.groupby("row").sample_id.nunique().to_dict() == {"single_chain": 2193, "multi_chain": 870}, \
        df.groupby("row").sample_id.nunique().to_dict()

    points, meta = [], {}
    for row, _, _, _ in ROWS:
        frame = df[df.row == row]
        assert frame.groupby("method").size().eq(len(frame) // 7).all()
        boot, nfam = family_bootstrap(frame, [x[0] for x in METHODS])
        nrec = frame.sample_id.nunique()
        meta[row] = dict(records=int(nrec), families=int(nfam))
        for method in [x[0] for x in METHODS]:
            for metric, _, _, _ in METRICS:
                mean, lo, hi = boot[(method, metric)]
                points.append(dict(row=row, method=method, metric=metric, mean=mean,
                                   ci_low=lo, ci_high=hi, n_records=nrec, n_families=nfam))
    pts = pd.DataFrame(points)
    assert len(pts) == 42 and pts["mean"].notna().all()
    sub = df[(df.row == "multi_chain")].groupby("group").sample_id.nunique().to_dict()
    return dict(points=pts, meta=meta, multi_split=sub,
                hashes={str(p): digest(p) for p in (MEMBERSHIP, PER_RECORD)})


def draw(data):
    pts = data["points"]
    fig, axes = plt.subplots(1, 3, figsize=(183 / 25.4, 74 / 25.4))
    fig.subplots_adjust(left=.068, right=.992, bottom=.275, top=.865, wspace=.26)
    xs = np.arange(len(METHODS))
    width = .34
    offsets = {"single_chain": -width / 2, "multi_chain": width / 2}
    for ci, (metric, label, limits, ticks) in enumerate(METRICS):
        ax = axes[ci]
        for row, _, colour, edge in ROWS:
            vals, los, his = [], [], []
            for key, name, family in METHODS:
                r = pts[(pts.row == row) & (pts.method == key) & (pts.metric == metric)].iloc[0]
                vals.append(r["mean"]); los.append(r["mean"] - r.ci_low); his.append(r.ci_high - r["mean"])
            ax.bar(xs + offsets[row], vals, width=width, color=colour, edgecolor=edge,
                   linewidth=.45, zorder=2, yerr=[los, his],
                   error_kw=dict(ecolor="#31424C", elinewidth=.6, capsize=1.4, capthick=.6, zorder=3))
        ax.set_ylim(*limits)
        ax.set_yticks(ticks)
        ax.set_xticks(xs)
        ax.set_xticklabels([m[1] for m in METHODS], rotation=32, ha="right",
                           rotation_mode="anchor", fontsize=6.4)
        for tick, m in zip(ax.get_xticklabels(), METHODS):
            if m[2] == "own":
                tick.set_fontweight("bold")
        ax.tick_params(axis="y", labelsize=6.2, length=2.2, pad=2)
        ax.tick_params(axis="x", length=2.2, pad=2)
        ax.set_ylabel(label, fontsize=6.8, labelpad=3)
        ax.spines["left"].set_color("#93A1A8")
        ax.grid(axis="y", color="#E3E8EA", lw=.5, zorder=0)
        ax.set_axisbelow(True)
        ax.text(-.10, 1.03, chr(97 + ci), transform=ax.transAxes, fontsize=8.2,
                fontweight="bold", color=INK, va="bottom")
    handles = [Patch(facecolor=c, edgecolor=e, lw=.45,
                     label=f"{name} (n = {data['meta'][row]['records']:,})")
               for row, name, c, e in ROWS]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(.53, .995), ncol=2,
               fontsize=6.6, handlelength=1.2, handletextpad=.5, columnspacing=2.4)
    fig.text(.008, .028, "Bars are record means; whiskers are 95% family-cluster bootstrap intervals. "
                         "Multichain inputs = 10 homooligomers + 860 heteromeric complexes.",
             fontsize=6.0, color="#4A5A62")
    OUT.mkdir(exist_ok=True)
    stem = OUT / "benchmark_v7"
    box = Bbox.from_bounds(0, 0, *fig.get_size_inches())
    for suffix in ("pdf", "png", "svg"):
        fig.savefig(stem.with_suffix("." + suffix), bbox_inches=box, dpi=300)
    plt.close(fig)
    pts.to_csv(OUT / "benchmark_v7_source.csv", index=False)
    return stem


def main():
    data = load()
    stem = draw(data)
    assert all(digest(Path(p)) == h for p, h in data["hashes"].items())
    (OUT / "plot_audit_v7.json").write_text(json.dumps(dict(
        generated_utc=datetime.now(timezone.utc).isoformat(),
        script_sha256=digest(__file__), input_sha256=data["hashes"],
        rows={r: data["meta"][r] for r, _, _, _ in ROWS},
        multichain_split=data["multi_split"],
        methods=[m[1] for m in METHODS], metrics=[m[0] for m in METRICS], cells=42,
        input_protocol="whole assembly for every structure-consuming method",
        statistics=f"record-weighted mean, 95% family-cluster percentile, seed={SEED}, draws={BOOTSTRAP}",
        artifacts={f"benchmark_v7.{s}": digest(stem.with_suffix("." + s)) for s in ("pdf", "png", "svg")}
    ), indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(dict(rows=data["meta"], split=data["multi_split"]), ensure_ascii=False))


if __name__ == "__main__":
    main()
