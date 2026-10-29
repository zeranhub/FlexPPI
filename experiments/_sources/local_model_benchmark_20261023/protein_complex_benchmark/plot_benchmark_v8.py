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

"""Stratification by the number of protein chains in the simulated input.

Four strata (1, 2, 3, >=4 chains), three metrics, seven methods, one row of
three panels. Bars are record means with family-cluster bootstrap intervals.
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
OUT = ROOT / "figures_v8"
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
STRATA = [("1", "1 chain", "#DCE7F0"),
          ("2", "2 chains", "#9BB9D2"),
          ("3", "3 chains", "#4E7FA6"),
          ("4+", "\u22654 chains", "#1B3A57")]
SEED = 20261024
BOOTSTRAP = 2000
INK = "#22313A"


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as fh:
        for block in iter(lambda: fh.read(1024 ** 2), b""):
            h.update(block)
    return h.hexdigest()


def strata_key(n):
    n = int(n)
    return "1" if n <= 1 else "2" if n == 2 else "3" if n == 3 else "4+"


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
    df = per.merge(mem[["sample_id", "actual_input_n_chains"]], on="sample_id", how="left")
    df = df[df.method.isin([m[0] for m in METHODS])].copy()
    df["stratum"] = df.actual_input_n_chains.map(strata_key)
    sizes = df.groupby("stratum").sample_id.nunique().to_dict()
    assert sizes == {"1": 2193, "2": 422, "3": 131, "4+": 317}, sizes

    points, meta = [], {}
    for key, _, _ in STRATA:
        frame = df[df.stratum == key]
        boot, nfam = family_bootstrap(frame, [m[0] for m in METHODS])
        nrec = int(frame.sample_id.nunique())
        meta[key] = dict(records=nrec, families=int(nfam))
        for method in [m[0] for m in METHODS]:
            for metric, _, _, _ in METRICS:
                mean, lo, hi = boot[(method, metric)]
                points.append(dict(stratum=key, method=method, metric=metric, mean=mean,
                                   ci_low=lo, ci_high=hi, n_records=nrec, n_families=nfam))
    pts = pd.DataFrame(points)
    assert len(pts) == 84 and pts["mean"].notna().all()
    return dict(points=pts, meta=meta,
                hashes={str(p): digest(p) for p in (MEMBERSHIP, PER_RECORD)})


def draw(data):
    pts = data["points"]
    fig, axes = plt.subplots(1, 3, figsize=(183 / 25.4, 84 / 25.4))
    fig.subplots_adjust(left=.068, right=.992, bottom=.265, top=.855, wspace=.26)
    xs = np.arange(len(METHODS))
    n_str = len(STRATA)
    width = .185
    offs = {k: (i - (n_str - 1) / 2) * width for i, (k, _, _) in enumerate(STRATA)}
    for ci, (metric, label, limits, ticks) in enumerate(METRICS):
        ax = axes[ci]
        for key, _, colour in STRATA:
            vals, los, his = [], [], []
            for mkey, name, fam in METHODS:
                r = pts[(pts.stratum == key) & (pts.method == mkey) & (pts.metric == metric)].iloc[0]
                vals.append(r["mean"]); los.append(r["mean"] - r.ci_low); his.append(r.ci_high - r["mean"])
            ax.bar(xs + offs[key], vals, width=width, color=colour, edgecolor="#1B2C36",
                   linewidth=.30, zorder=2, yerr=[los, his],
                   error_kw=dict(ecolor="#31424C", elinewidth=.5, capsize=1.0, capthick=.5, zorder=3))
        ax.set_ylim(*limits)
        ax.set_yticks(ticks)
        ax.set_xticks(xs)
        ax.set_xticklabels([m[1] for m in METHODS], rotation=32, ha="right",
                           rotation_mode="anchor", fontsize=6.2)
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
    handles = [Patch(facecolor=c, edgecolor="#1B2C36", lw=.35,
                     label=f"{name} (n = {data['meta'][k]['records']:,})")
               for k, name, c in STRATA]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(.5, .995), ncol=4,
               fontsize=6.3, handlelength=1.1, handletextpad=.45, columnspacing=1.6)
    fig.text(.008, .026, "Bars are record means over the residues with a valid RMSF label; whiskers are 95% family-cluster bootstrap intervals. "
                         "Every method reads the full input structure.",
             fontsize=5.9, color="#4A5A62")
    OUT.mkdir(exist_ok=True)
    stem = OUT / "benchmark_v8"
    box = Bbox.from_bounds(0, 0, *fig.get_size_inches())
    for suffix in ("pdf", "png", "svg"):
        fig.savefig(stem.with_suffix("." + suffix), bbox_inches=box, dpi=300)
    plt.close(fig)
    pts.to_csv(OUT / "benchmark_v8_source.csv", index=False)
    return stem


def main():
    data = load()
    stem = draw(data)
    assert all(digest(Path(p)) == h for p, h in data["hashes"].items())
    (OUT / "plot_audit_v8.json").write_text(json.dumps(dict(
        generated_utc=datetime.now(timezone.utc).isoformat(),
        script_sha256=digest(__file__), input_sha256=data["hashes"],
        strata={k: data["meta"][k] for k, _, _ in STRATA},
        methods=[m[1] for m in METHODS], metrics=[m[0] for m in METRICS], cells=84,
        stratum_definition="number of protein chains in the simulated input; 1/2/3/>=4",
        input_protocol="whole structure for every structure-consuming method",
        statistics=f"record-weighted mean, 95% family-cluster percentile, seed={SEED}, draws={BOOTSTRAP}",
        artifacts={f"benchmark_v8.{s}": digest(stem.with_suffix("." + s)) for s in ("pdf", "png", "svg")}
    ), indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(data["meta"], ensure_ascii=False))


if __name__ == "__main__":
    main()

