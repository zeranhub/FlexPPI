# Academic Figure Skill Typography Baseline — COPY VERBATIM, place at TOP of script
import matplotlib as mpl
mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans"],
    "font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 8,
    "figure.titlesize": 9, "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 0.6, "xtick.direction": "out", "ytick.direction": "out",
    "xtick.major.width": 0.6, "ytick.major.width": 0.6, "legend.frameon": False,
})
mpl.rcParams.update({"pdf.fonttype": 42, "svg.fonttype": "none",
                     "savefig.bbox": "tight", "savefig.dpi": 300})

"""Figure 2 panel d: DynaRepo complexes grouped by inter-chain interface count.

Palette follows the blue/orange/teal of panels a-c; bars use the blue ramp so
the panel reads as part of the same figure.
"""
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.transforms import Bbox

ROOT = Path(r"E:\Project\FlexPPI\experiments\_sources\local_model_benchmark_20261023\protein_complex_benchmark")
OUT = Path(__file__).resolve().parent
MEMBERSHIP = ROOT / "classification" / "membership_frozen.csv"
PER_RECORD = ROOT / "evaluation" / "per_record_metrics.csv"
PPI = ROOT / "classification" / "ppi_interfaces_20261024.csv"

METHODS = [
    ("Full", "FlexPPI", True),
    ("backflip-2.1 (assembly)", "BackFlip-2.1", False),
    ("backflip-1.0 (assembly)", "BackFlip-1.0", False),
    ("DynaProt-M", "DynaProt-M", False),
    ("FastProtFlex", "FastProtFlex", False),
    ("ANM", "ANM", False),
    ("GNM", "GNM", False),
]
METRICS = [("spearman", "Spearman \u03c1 \u2191", (0, 1.0), [0, 0.5, 1.0]),
           ("top10_ap", "High-RMSF AP \u2191", (0, 1.0), [0, 0.5, 1.0]),
           ("z_mae", "Standardized MAE \u2193", (0, 0.72), [0, 0.2, 0.4, 0.6])]
STRATA = [("i1", "1 interface", "#B9CEE3"),
          ("i23", "2\u20133 interfaces", "#5E8FBF"),
          ("i4", "\u22654 interfaces", "#2060A8")]
SEED = 20261024; BOOTSTRAP = 2000; INK = "#1F2C33"


def family_bootstrap(frame, methods):
    records = sorted(frame[frame.method == "Full"].sample_id)
    fam = frame[frame.method == "Full"].set_index("sample_id").loc[records, "family"]
    names = sorted(fam.unique()); pos = {f: i for i, f in enumerate(names)}
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
            assert np.isfinite(values).all()
            sums = np.bincount(codes, weights=values, minlength=len(names))
            draws = (weights @ sums) / denom
            out[(method, metric)] = (float(values.mean()), float(np.quantile(draws, .025)),
                                     float(np.quantile(draws, .975)))
    return out, len(names)


def load():
    mem = pd.read_csv(MEMBERSHIP, encoding="utf-8-sig")
    per = pd.read_csv(PER_RECORD)
    ppi = pd.read_csv(PPI)
    df = per.merge(mem[["sample_id", "actual_input_n_chains", "benchmark_group"]], on="sample_id", how="left")
    df = df[(df.method.isin([m[0] for m in METHODS]))
            & (df.dataset == "dynarepo")
            & (df.benchmark_group == "protein_complex")
            & (df.actual_input_n_chains >= 2)].merge(ppi[["sample_id", "n_interfaces"]], on="sample_id")
    df = df[df.n_interfaces >= 1].copy()
    n = df.n_interfaces.astype(int)
    df["stratum"] = np.select([n == 1, (n >= 2) & (n <= 3), n >= 4], ["i1", "i23", "i4"], default="")
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
    assert len(pts) == 63 and pts["mean"].notna().all()
    assert sum(v["records"] for v in meta.values()) == 615
    return pts, meta


def draw(pts, meta):
    fig, axes = plt.subplots(1, 3, figsize=(150 / 25.4, 52 / 25.4))
    fig.subplots_adjust(left=.085, right=.990, bottom=.330, top=.795, wspace=.32)
    xs = np.arange(len(METHODS)); width = .20
    offs = {k: (i - 1) * width for i, (k, _, _) in enumerate(STRATA)}
    for ci, (metric, label, limits, ticks) in enumerate(METRICS):
        ax = axes[ci]
        for key, _, colour in STRATA:
            vals, los, his = [], [], []
            for mkey, name, _ in METHODS:
                r = pts[(pts.stratum == key) & (pts.method == mkey) & (pts.metric == metric)].iloc[0]
                vals.append(r["mean"]); los.append(r["mean"] - r.ci_low); his.append(r.ci_high - r["mean"])
            ax.bar(xs + offs[key], vals, width=width, color=colour, edgecolor="#2B4A63",
                   linewidth=.25, zorder=2, yerr=[los, his],
                   error_kw=dict(ecolor="#3C5A6E", elinewidth=.45, capsize=.9, capthick=.45, zorder=3))
        ax.set_ylim(*limits); ax.set_yticks(ticks)
        ax.set_xticks(xs)
        ax.set_xticklabels([m[1] for m in METHODS], rotation=30, ha="right",
                           rotation_mode="anchor", fontsize=6.0)
        for tick, m in zip(ax.get_xticklabels(), METHODS):
            if m[2]:
                tick.set_fontweight("bold")
        ax.tick_params(axis="y", labelsize=6.2, length=2.0, pad=1.5)
        ax.tick_params(axis="x", length=2.0, pad=1.5)
        ax.set_ylabel(label, fontsize=6.6, labelpad=2.5)
        ax.spines["left"].set_color("#8E9BA3")
        ax.grid(axis="y", color="#E6EBEE", lw=.45, zorder=0)
        ax.set_axisbelow(True)
        if ci == 0:
            ax.text(-.135, 1.04, "d", transform=ax.transAxes,
                    fontsize=9.5, fontweight="bold", color=INK, va="bottom")
    handles = [Patch(facecolor=c, edgecolor="#2B4A63", lw=.25,
                     label=f"{name} (n = {meta[k]['records']}, {meta[k]['families']} fam)")
               for k, name, c in STRATA]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(.5, 1.0), ncol=3,
               fontsize=6.1, handlelength=.85, handletextpad=.4, columnspacing=2.4)
    OUT.mkdir(exist_ok=True)
    stem = OUT / "fig2_panel_d"
    box = Bbox.from_bounds(0, 0, *fig.get_size_inches())
    for suffix in ("pdf", "png", "svg"):
        fig.savefig(stem.with_suffix("." + suffix), bbox_inches=box, dpi=450)
    plt.close(fig)
    pts.to_csv(OUT / "fig2_panel_d_source.csv", index=False)
    return stem


if __name__ == "__main__":
    pts, meta = load()
    stem = draw(pts, meta)
    print("panel:", stem, meta)
