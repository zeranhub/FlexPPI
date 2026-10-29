# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a,b) BarAblation/plot_bar_ablation.py -> cross-type inherit -> param inherit
# (c,d) Forest from BarAblation -> cross-type inherit -> param inherit
# (supplement) heatmap/plot_composition.py -> cross-type inherit -> param inherit
# The assets contain aggregate example scores/counts, not paired protein errors.
# Inherit direct labels, reference lines, clean spines and restrained emphasis;
# use native vector artists to retain all observations and valid paired intervals.
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
import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
METHODS = ["flexppi", "nn_geometry", "nn_intrachain", "nn_sequence"]
ALL_METHODS = ["aa", "sequence", "geometry", "sequence_intra", "sequence_complex"] + METHODS
LABELS = {"flexppi": "Full FlexPPI", "nn_geometry": "Identity masked",
          "nn_intrachain": "Intrachain only", "nn_sequence": "Sequence only",
          "aa": "AA-type prior", "sequence": "Sequence probe", "geometry": "Geometry probe",
          "sequence_intra": "Sequence + local geometry", "sequence_complex": "Sequence + geometry + contacts"}
COLORS = {"flexppi": "#0F4D92", "nn_geometry": "#F1A340", "nn_intrachain": "#762A83", "nn_sequence": "#777777"}
MARKERS = {"flexppi": "o", "nn_geometry": "s", "nn_intrachain": "D", "nn_sequence": "^"}
METRICS = ["rmsf", "chi_entropy"]
METRIC_NAMES = {"rmsf": "Backbone flexibility · RMSF", "chi_entropy": r"Side-chain entropy · $H_\chi$"}
UNITS = {"rmsf": "Å", "chi_entropy": "nats"}
SEED = 20260922
REPEATS = 5000

def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

def bootstrap_means(values, seed=SEED):
    """Resample whole protein rows; all methods share each resampling draw."""
    values = np.asarray(values, float)
    assert values.ndim == 2 and np.isfinite(values).all()
    rng = np.random.default_rng(seed)
    out = np.empty((REPEATS, values.shape[1]))
    for start in range(0, REPEATS, 100):
        size = min(100, REPEATS - start)
        ix = rng.integers(0, len(values), size=(size, len(values)))
        out[start:start + size] = values[ix].mean(axis=1)
    return out

def statistics(table):
    performance, contrasts = [], []
    for metric in METRICS:
        block = table[table.metric.eq(metric)]
        values = block[ALL_METHODS].to_numpy(float)
        means = values.mean(axis=0)
        boot = bootstrap_means(values)
        intervals = np.quantile(boot, [.025, .975], axis=0)
        full_idx = ALL_METHODS.index("flexppi")
        for j, method in enumerate(ALL_METHODS):
            performance.append(dict(metric=metric, method=method, n_proteins=len(block),
                                    mean_mae=means[j], ci_low=intervals[0, j], ci_high=intervals[1, j], unit=UNITS[metric]))
            if method in METHODS[1:]:
                delta = values[:, j] - values[:, full_idx]
                lo, hi = np.quantile(boot[:, j] - boot[:, full_idx], [.025, .975])
                contrasts.append(dict(metric=metric, method=method, n_proteins=len(block),
                                      mean_delta=delta.mean(), ci_low=lo, ci_high=hi,
                                      median_delta=np.median(delta), full_better_fraction=(delta > 0).mean(),
                                      ties=int((delta == 0).sum()),
                                      relative_error_increase_pct=100 * delta.mean() / means[full_idx],
                                      unit=UNITS[metric]))
    return pd.DataFrame(performance), pd.DataFrame(contrasts)

def panel_title(ax, letter, title, subtitle):
    ax.text(-.04, 1.16, letter, transform=ax.transAxes, weight="bold", fontsize=10,
            ha="right", va="bottom", color=BLACK)
    ax.text(.01, 1.17, title, transform=ax.transAxes, fontsize=8.5, weight="bold", color=BLACK)
    ax.text(.01, 1.07, subtitle, transform=ax.transAxes, fontsize=6.5, color="#666666")

def export(fig, name, out):
    """Fixed canvas preserves the declared 183 mm export width."""
    for suffix in ("pdf", "svg", "png"):
        with mpl.rc_context({"savefig.bbox": None}):
            fig.savefig(out / (name + "." + suffix), dpi=300, facecolor="white", bbox_inches=None)
    plt.close(fig)

def compose_figure(performance, contrasts, out, archetype="quantitative_grid"):
    assert archetype == "quantitative_grid"
    fig = plt.figure(figsize=(183 / 25.4, 140 / 25.4))
    fig.text(.045, .969, "Information ablations reveal distinct error profiles", fontsize=11, weight="bold", color=BLACK)
    fig.text(.045, .94, "Matched five-fold OOF predictions · 3,068 records · one training seed", fontsize=7, color="#666666")
    positions = [[.195, .59, .30, .235], [.665, .59, .30, .235],
                 [.195, .21, .30, .23], [.665, .21, .30, .23]]
    for col, metric in enumerate(METRICS):
        data = performance[performance.metric.eq(metric)].set_index("method")
        ax = fig.add_axes(positions[col])
        n = int(data.loc["flexppi", "n_proteins"])
        panel_title(ax, "ab"[col], METRIC_NAMES[metric], f"{n:,} protein identities · lower error is better")
        xmax = .83 if metric == "rmsf" else .78
        for y, method in enumerate(METHODS):
            row = data.loc[method]
            ax.barh(y, row.mean_mae, height=.51, color=COLORS[method], alpha=.88, edgecolor="none", zorder=2)
            ax.errorbar(row.mean_mae, y, xerr=[[row.mean_mae-row.ci_low], [row.ci_high-row.mean_mae]],
                        fmt="none", ecolor=BLACK, elinewidth=.75, capsize=2, capthick=.7, zorder=3)
            ax.text(row.ci_high + xmax*.027, y, f"{row.mean_mae:.3f}", va="center", fontsize=7, color=BLACK)
        ax.set(yticks=range(4), yticklabels=[LABELS[m] for m in METHODS], xlim=(0, xmax), ylim=(3.65, -.6),
               xlabel=f"Mean absolute error ({UNITS[metric]})")
        ax.set_xticks([0, .2, .4, .6, .8] if col == 0 else [0, .2, .4, .6])
        ax.tick_params(axis="y", length=0, pad=6)
        ax.spines["left"].set_visible(False)
        ax.xaxis.label.set_size(7)

        ax = fig.add_axes(positions[col+2])
        panel_title(ax, "cd"[col], "Cost of each information restriction", "ΔMAE = ablated − full; positive favours full")
        diffs = contrasts[contrasts.metric.eq(metric)].set_index("method")
        xmin, xmax = (-.035, .29) if col == 0 else (-.06, .82)
        ax.axvline(0, color="#999999", linestyle="--", linewidth=.7, zorder=0)
        for y, method in enumerate(METHODS[1:]):
            row = diffs.loc[method]
            ax.errorbar(row.mean_delta, y, xerr=[[row.mean_delta-row.ci_low], [row.ci_high-row.mean_delta]],
                        fmt=MARKERS[method], color=COLORS[method], ms=4.5, elinewidth=1.3,
                        capsize=2.4, capthick=.8, zorder=3)
            label = f"{row.mean_delta:+.3f}"
            ax.text(row.ci_high+(xmax-xmin)*.045, y, label, fontsize=7, va="center", color=BLACK)
            ax.text(.99, y, f"{row.full_better_fraction:.1%}",
                    transform=ax.get_yaxis_transform(), fontsize=6.5, color="#666666", ha="right", va="center")
        ax.text(.99, 1.005, "Full wins", transform=ax.transAxes, fontsize=6.5, ha="right", color="#666666")
        ax.set(yticks=range(3), yticklabels=[LABELS[m] for m in METHODS[1:]],
               ylim=(2.65, -.6), xlim=(xmin, xmax), xlabel=f"Paired error difference ({UNITS[metric]})")
        ax.set_xticks([0, .1, .2] if col == 0 else [0, .2, .4, .6, .8])
        ax.tick_params(axis="y", length=0, pad=6)
        ax.spines["left"].set_visible(False)
        ax.xaxis.label.set_size(7)
    fig.text(.045, .105, "Bars / points: protein-equal means. Whiskers: 95% protein-bootstrap intervals (5,000 resamples).", fontsize=6.4, color="#555555")
    fig.text(.045, .081, "Full wins: fraction of protein identities with lower error in Full FlexPPI. No observations are trimmed.", fontsize=6.4, color="#555555")
    fig.text(.045, .057, "Sequence-only capacity = 61.3% of full; identity masking retains structural clues to residue type.", fontsize=6.4, color="#555555")
    fig.text(.045, .033, "Record-level OOF: repeated sequences can span folds. Intervals describe test-set variation, not seed stability.", fontsize=6.4, color="#555555")
    export(fig, "fig_neural_ablation_main", out)

def figure_distributions(table, contrasts, out):
    fig = plt.figure(figsize=(183/25.4, 93/25.4))
    fig.text(.05, .955, "The full distribution of paired protein errors", fontsize=11, weight="bold", color=BLACK)
    fig.text(.05, .912, "All protein pairs retained · one training seed · symmetric logarithmic scale", fontsize=7, color="#666666")
    for col, metric in enumerate(METRICS):
        ax = fig.add_axes([.105 + col*.48, .245, .365, .49])
        block = table[table.metric.eq(metric)]
        panel_title(ax, "ab"[col], METRIC_NAMES[metric], f"{len(block):,} paired protein identities")
        ax.axvline(0, color="#999999", linewidth=.7, linestyle="--")
        for method in METHODS[1:]:
            delta = np.sort((block[method]-block.flexppi).to_numpy())
            ax.step(delta, np.arange(1,len(delta)+1)/len(delta), where="post", color=COLORS[method], lw=1.2,
                    label=LABELS[method])
        ax.set_xscale("symlog", linthresh=.01)
        ax.set_ylim(0, 1)
        ax.set_yticks([0, .25, .5, .75, 1])
        ax.set_xlabel(f"ΔMAE = ablated − full ({UNITS[metric]})", fontsize=7)
        ax.set_ylabel("Cumulative fraction of proteins", fontsize=7)
        ax.set_xticks([-10,-1,-.1,0,.1,1,10] if metric=="rmsf" else [-.1,0,.1,1])
        ax.set_xticklabels(["−10","−1","−0.1","0","0.1","1","10"] if metric=="rmsf" else ["−0.1","0","0.1","1"])
        delta_all = (block[METHODS[1:]].to_numpy()-block.flexppi.to_numpy()[:,None]).ravel()
        ax.set_xlim(delta_all.min()*1.2, delta_all.max()*1.2)
    handles = [Line2D([0],[0],color=COLORS[m],lw=1.4,label=LABELS[m]) for m in METHODS[1:]]
    fig.legend(handles=handles,loc="lower center",bbox_to_anchor=(.52,.094),ncol=3,fontsize=7)
    fig.text(.05,.042,"Right of zero: full has lower MAE. Linear within ±0.01; logarithmic beyond that range. All protein pairs retained.",fontsize=6.5,color="#555555")
    export(fig,"fig_neural_ablation_distributions",out)

def figure_all_controls(performance, out):
    fig=plt.figure(figsize=(183/25.4,112/25.4))
    fig.text(.055,.96,"Neural ablations and feature probes in context",fontsize=11,weight="bold",color=BLACK)
    fig.text(.055,.921,"Different model families are separated; probe-to-neural gaps are descriptive comparisons",fontsize=7,color="#666666")
    order=METHODS+["aa","sequence","geometry","sequence_intra","sequence_complex"]
    yvalues=np.array([0,1,2,3,5,6,7,8,9])
    for col,metric in enumerate(METRICS):
        ax=fig.add_axes([.28+col*.405,.23,.29,.55])
        block=performance[performance.metric.eq(metric)].set_index("method")
        panel_title(ax,"ab"[col],METRIC_NAMES[metric],f"{int(block.iloc[0].n_proteins):,} protein identities")
        for y,method in zip(yvalues,order):
            row=block.loc[method]
            color=COLORS.get(method,"#999999")
            ax.errorbar(row.mean_mae,y,xerr=[[row.mean_mae-row.ci_low],[row.ci_high-row.mean_mae]],fmt=MARKERS.get(method,"o"),
                        ms=4,elinewidth=1.1,capsize=2,color=color)
            ax.text(row.ci_high+.017,y,f"{row.mean_mae:.3f}",fontsize=6.5,va="center",color=BLACK)
        ax.axhline(4,color="#dddddd",lw=.7)
        ax.set(yticks=yvalues,yticklabels=[LABELS[m] for m in order] if col==0 else [""]*len(order),
               ylim=(9.7,-.7),xlim=(.40,.98) if col==0 else (.11,.84),xlabel=f"Protein-equal MAE ({UNITS[metric]})")
        ax.tick_params(axis="y",length=0,pad=6)
        ax.spines["left"].set_visible(False)
        ax.xaxis.label.set_size(7)
    fig.text(.055,.108,"Top: full neural model and three retrained neural restrictions. Bottom: AA prior and fixed-budget boosting probes.",fontsize=6.5,color="#555555")
    fig.text(.055,.080,"Points: protein-equal mean MAE; whiskers: 95% protein-bootstrap intervals. Smaller error is better.",fontsize=6.5,color="#555555")
    fig.text(.055,.052,"Model capacity and representation differ across families; these gaps do not identify an isolated causal contribution.",fontsize=6.5,color="#555555")
    fig.text(.055,.024,"One training seed; intervals describe test-set variation conditional on fitted models. Sequences may span folds.",fontsize=6.5,color="#555555")
    export(fig,"fig_ablation_all_controls",out)

def figure_heterogeneity(audit_dir, out):
    """Retain mixed-source identities; fold panels use complete residue tables."""
    source_path = audit_dir / "subgroups/nn_subgroup_summary.csv"
    fold_path = audit_dir / "recomputed_fold_mae.csv"
    grouped = pd.read_csv(source_path)
    folds = pd.read_csv(fold_path)
    groups = ["atlas", "dynarepo", "itrace", "atlas|dynarepo", "dynarepo|itrace"]
    group_names = ["ATLAS", "DynaRepo", "iTrace", "ATLAS + DynaRepo", "DynaRepo + iTrace"]
    fig = plt.figure(figsize=(183/25.4, 155/25.4))
    fig.text(.045, .973, "Source heterogeneity and consistency across folds", fontsize=11, weight="bold", color=BLACK)
    fig.text(.045, .948, "ΔMAE = ablated − full; positive values favour the full model", fontsize=7, color="#666666")
    cmap = LinearSegmentedColormap.from_list("paired_delta", DIVERGING)
    fold_markers = ["o", "s", "^", "D", "v"]
    exported = []
    for col, metric in enumerate(METRICS):
        ax = fig.add_axes([.195+col*.47, .635, .30, .215])
        block = grouped[grouped.group_type.eq("source") & grouped.metric.eq(metric)]
        values = block.pivot(index="group", columns="model", values="mean_delta").loc[groups, METHODS[1:]]
        counts = block.groupby("group").n.first().loc[groups]
        limit = float(np.abs(values.to_numpy()).max())
        norm = TwoSlopeNorm(vmin=-limit, vcenter=0, vmax=limit)
        im = ax.imshow(values, cmap=cmap, norm=norm, aspect="auto")
        panel_title(ax, "ab"[col], METRIC_NAMES[metric], "Protein-equal mean difference by data source")
        for y in range(5):
            for x in range(3):
                value = float(values.iloc[y, x])
                ax.text(x, y, f"{value:+.4f}", ha="center", va="center", fontsize=6.5,
                        color="white" if abs(value) > limit*.63 else BLACK)
        ax.set_xticks(range(3), ["Identity\nmasked", "Intrachain\nonly", "Sequence\nonly"], fontsize=6.5)
        ax.set_yticks(range(5), [f"{name} (n={int(n):,})" for name, n in zip(group_names, counts)], fontsize=6.2)
        ax.set_xticks(np.arange(-.5, 3, 1), minor=True)
        ax.set_yticks(np.arange(-.5, 5, 1), minor=True)
        ax.grid(which="minor", color="white", linewidth=1)
        ax.tick_params(which="both", length=0)
        for spine in ax.spines.values():
            spine.set_visible(False)
        cax = fig.add_axes([.195+col*.47, .548, .30, .010])
        cb = fig.colorbar(im, cax=cax, orientation="horizontal", ticks=[-limit,0,limit])
        cb.ax.tick_params(labelsize=6, length=2, width=.5)
        cb.outline.set_visible(False)
        cb.set_label(f"Paired ΔMAE ({UNITS[metric]})", fontsize=6.5, labelpad=2)

        ax = fig.add_axes([.195+col*.47, .205, .30, .235])
        panel_title(ax, "cd"[col], "Five-fold directional consistency", "Residue-pooled MAE; one fitted model per fold")
        wide = folds[folds.metric.eq(metric)].pivot(index="fold", columns="method", values="residue_pooled_mae")
        ax.axvline(0, color="#999999", linestyle="--", linewidth=.7)
        for y, method in enumerate(METHODS[1:]):
            delta = wide[method] - wide.flexppi
            for fold in range(5):
                ax.scatter(delta.loc[fold], y+(fold-2)*.075, marker=fold_markers[fold],
                           s=16, color=COLORS[method], linewidths=.4, edgecolors="white", zorder=3)
                exported.append(dict(metric=metric, method=method, fold=fold,
                                     mean_full=wide.loc[fold,"flexppi"], mean_ablated=wide.loc[fold,method],
                                     delta=float(delta.loc[fold]), aggregation="residue_pooled"))
        all_deltas = wide[METHODS[1:]].to_numpy()-wide.flexppi.to_numpy()[:,None]
        low, high = all_deltas.min(), all_deltas.max()
        padding = (high-low)*.14
        ax.set(xlim=(min(low-padding, -.01), high+padding), ylim=(2.6,-.6),
               yticks=range(3), yticklabels=[LABELS[m] for m in METHODS[1:]], xlabel=f"Paired ΔMAE ({UNITS[metric]})")
        ax.xaxis.label.set_size(7)
        ax.tick_params(axis="y", length=0, pad=6)
        ax.spines["left"].set_visible(False)
    handles = [Line2D([0],[0], marker=marker, linestyle="none", color="#555555", ms=4,
                      label=f"Fold {i}") for i, marker in enumerate(fold_markers)]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(.55,.119), ncol=5, fontsize=6.5,
               handletextpad=.4, columnspacing=1.2)
    fig.text(.045,.096,"All five source memberships are retained, including mixed sources. Source panels use equal protein weights.",fontsize=6.4,color="#555555")
    fig.text(.045,.074,"Fold panels retain all valid residues, including protein identities spanning folds; they use a different weighting.",fontsize=6.4,color="#555555")
    fig.text(.045,.052,"The five folds share training data and use one seed. The points show directional consistency, not five independent repeats.",fontsize=6.4,color="#555555")
    fig.text(.045,.030,"Source composition and chain context are associated; subgroup patterns are descriptive comparisons.",fontsize=6.4,color="#555555")
    pd.DataFrame(exported).to_csv(out/"source_data_fold_deltas.csv",index=False)
    block = grouped[grouped.group_type.eq("source")]
    block.to_csv(out/"source_data_source_groups.csv",index=False)
    export(fig,"fig_neural_ablation_heterogeneity",out)

def figure_dependence_sensitivity(contrasts, audit_dir, out):
    sensitivity = pd.read_csv(audit_dir/"cluster_bootstrap_sensitivity.csv")
    fig = plt.figure(figsize=(183/25.4, 105/25.4))
    fig.text(.045,.963,"Dependence-aware resampling widens uncertainty",fontsize=11,weight="bold",color=BLACK)
    fig.text(.045,.925,"Identical protein-equal point estimates; paired 95% intervals under two resampling units",fontsize=7,color="#666666")
    rows=[]
    for col,metric in enumerate(METRICS):
        ax=fig.add_axes([.205+col*.465,.295,.285,.425])
        panel_title(ax,"ab"[col],METRIC_NAMES[metric],"Protein rows versus connected biological groups")
        ordinary=contrasts[contrasts.metric.eq(metric)].set_index("method")
        grouped=sensitivity[sensitivity.metric.eq(metric)].set_index("ablation")
        ax.axvline(0,color="#999999",linestyle="--",linewidth=.7)
        for y,method in enumerate(METHODS[1:]):
            row=ordinary.loc[method]
            cluster=grouped.loc[method]
            assert abs(row.mean_delta-cluster.mean_delta_mae)<1e-12
            ax.errorbar(row.mean_delta,y+.10,xerr=[[row.mean_delta-row.ci_low],[row.ci_high-row.mean_delta]],
                        fmt=MARKERS[method],color=COLORS[method],mfc="white",ms=4,elinewidth=.8,capsize=2,alpha=.75)
            ax.errorbar(row.mean_delta,y-.10,xerr=[[row.mean_delta-cluster.cluster_ci_low],[cluster.cluster_ci_high-row.mean_delta]],
                        fmt=MARKERS[method],color=COLORS[method],ms=4.5,elinewidth=1.5,capsize=2)
            rows.append(dict(metric=metric,method=method,mean_delta=row.mean_delta,
                             protein_ci_low=row.ci_low,protein_ci_high=row.ci_high,
                             cluster_ci_low=cluster.cluster_ci_low,cluster_ci_high=cluster.cluster_ci_high,
                             n_proteins=int(row.n_proteins),n_connected_units=int(cluster.n_connected_units),
                             largest_unit_proteins=int(cluster.largest_unit_proteins)))
        ax.set(yticks=range(3),yticklabels=[LABELS[m] for m in METHODS[1:]],ylim=(2.6,-.6),
               xlim=(-.045,.225) if col==0 else (-.065,.58),xlabel=f"ΔMAE = ablated − full ({UNITS[metric]})")
        ax.tick_params(axis="y",length=0,pad=6)
        ax.spines["left"].set_visible(False)
        ax.xaxis.label.set_size(7)
    handles=[Line2D([0],[0],marker="o",color="#666666",mfc="white",lw=.8,ms=4,label="Resample protein rows"),
             Line2D([0],[0],marker="o",color="#666666",lw=1.5,ms=4,label="Resample connected groups")]
    fig.legend(handles=handles,loc="lower center",bbox_to_anchor=(.54,.164),ncol=2,fontsize=7)
    fig.text(.045,.131,"5,000 paired resamples. Each sampled group contributes all its proteins; the final mean keeps equal protein weights.",fontsize=6.4,color="#555555")
    fig.text(.045,.101,"2,035 connected groups; the largest contains 625 proteins (16.2%). Group dependence increases several interval widths.",fontsize=6.4,color="#555555")
    fig.text(.045,.071,"Groups merge existing biological split clusters linked by shared exact sequences; full homology clustering was not rerun.",fontsize=6.4,color="#555555")
    fig.text(.045,.041,"This sensitivity analysis cannot remove cross-fold sequence overlap or measure variability from new training seeds.",fontsize=6.4,color="#555555")
    pd.DataFrame(rows).to_csv(out/"source_data_dependence_intervals.csv",index=False)
    export(fig,"fig_neural_ablation_dependence",out)

def main():
    parser=argparse.ArgumentParser(description="Reproducible, audited neural ablation figures.")
    parser.add_argument("--source",type=Path,default=ROOT/"output/information_ablation_figures/source_data_protein_error.csv")
    parser.add_argument("--output",type=Path,default=ROOT/"output/information_ablation_figures/verified")
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    table=pd.read_csv(args.source)
    assert not table.duplicated(["protein_id","metric"]).any()
    assert table.groupby("metric").size().to_dict()=={"chi_entropy":3867,"rmsf":3870}
    assert np.isfinite(table[ALL_METHODS].to_numpy()).all()
    assert (table[ALL_METHODS].to_numpy()>=0).all()
    performance,contrasts=statistics(table)
    performance.to_csv(args.output/"performance_protein_equal.csv",index=False)
    contrasts.to_csv(args.output/"paired_ablation_effects.csv",index=False)
    table.to_csv(args.output/"source_data_all_proteins.csv",index=False)
    compose_figure(performance,contrasts,args.output)
    figure_distributions(table,contrasts,args.output)
    figure_all_controls(performance,args.output)
    audit_dir=ROOT/"output/information_ablation_analysis/verified_audit"
    figure_heterogeneity(audit_dir,args.output)
    figure_dependence_sensitivity(contrasts,audit_dir,args.output)
    report={"status":"rendered_pending_visual_QA","source":str(args.source.resolve()),"source_sha256":sha256(args.source),
            "script_sha256":sha256(__file__),"n_definition":"exact protein sequence identity, including short peptides",
            "aggregation":"residue absolute error -> chain mean -> mean over chain copies within record -> mean over records of one protein -> equal protein mean",
            "n_proteins":{"rmsf":3870,"chi_entropy":3867},"dropped_protein_rows":0,"downsampling":False,
            "interval":"percentile 95% bootstrap over exact-sequence protein rows; paired methods share each draw",
            "resampling_seed":SEED,"bootstrap_repeats":REPEATS,"tests":"None: effect estimation, no significance stars",
            "multiplicity":"Not applied; intervals are marginal descriptive intervals, not simultaneous familywise inference",
            "training_seeds":[7],"n_folds":5,"checkpoint":"fixed epoch 40 (metadata claims; ablation checkpoints unavailable locally)",
            "limitations":["Exact sequences can span folds; this is record OOF rather than protein-disjoint validation",
                           "Protein bootstrap does not establish independence among homologs, chains in one complex, or training runs",
                           "Sequence-only model has 836676 parameters vs 1365812 full/identity-masked/intrachain",
                           "Identity masking retains geometric clues to amino-acid type",
                           "Frozen valid-label masks are preserved; no unverified legacy label-coverage claim is inherited"],
            "figure_width_mm":183,"png_dpi":300,"fonts":"Arial with TrueType embedding; SVG editable text",
            "ecdf_scale":"symlog, linear inside +/-0.01; all tails retained"}
    report["diagnostic_source_files"]={str(p.resolve()):sha256(p) for p in [
        audit_dir/"subgroups/nn_subgroup_summary.csv", audit_dir/"recomputed_fold_mae.csv",
        audit_dir/"cluster_bootstrap_sensitivity.csv"]}
    report["heterogeneity"]={"source_aggregation":"protein_equal; all five source membership categories retained",
                             "fold_aggregation":"residue_pooled; complete five held-out folds, no arbitrary protein-fold assignment",
                             "folds_independent_repeats":False}
    (args.output/"statistics_and_reproducibility.json").write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding="utf-8")
    print(contrasts.drop(columns="unit").to_string(index=False))
    print("Saved figures to",args.output)

if __name__=="__main__":
    main()
