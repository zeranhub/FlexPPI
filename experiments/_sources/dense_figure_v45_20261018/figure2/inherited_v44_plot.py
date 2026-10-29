# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a,b) record-weighted density -> KernelDensity/plot_CorrelationDensity1.R -> param inherit
# (c,d) precision-recall curves -> cross-type from AUROC/plot_AUROC2.py -> param inherit
# (e,f) categorical paired gains -> verified v43 plot_model_capability.py -> param inherit
# RULE: no native-run panels. Production density uses unweighted R KDE; PR is not ROC.
from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parent
EXP = ROOT.parent.parent
sys.path.insert(0, str(EXP / 'completed_results_20261012' / 'plot_runtime'))
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
import shutil
import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LogNorm
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter, NullLocator

# Source canvas 183 mm, manuscript image 167.64 mm; 9 pt -> 8.24 pt.
mpl.rcParams.update({'font.size': 9, 'axes.titlesize': 9.5, 'axes.labelsize': 9,
                     'xtick.labelsize': 9, 'ytick.labelsize': 9, 'legend.fontsize': 9,
                     'axes.unicode_minus': False})
ANALYSIS = ROOT.parent / 'analysis'
OLD = EXP / 'figure_merge_v43_20261018' / 'figure2'
ARCHIVE = ROOT / 'source_archive'
ARCHIVE.mkdir(exist_ok=True)
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
paired_path = OLD / 'source_type_paired_gains_all80_rows.csv'
shutil.copy2(paired_path, ARCHIVE / paired_path.name)
assert sha(paired_path) == sha(ARCHIVE / paired_path.name)
paired = pd.read_csv(paired_path, float_precision='round_trip')
summary = json.loads((ANALYSIS / 'summary.json').read_text(encoding='utf8'))
pr = pd.read_csv(ANALYSIS / 'macro_pr_curve.csv', float_precision='round_trip')
endpoints = pd.read_csv(ANALYSIS / 'all_selected_endpoint.csv', float_precision='round_trip')
trend = pd.read_csv(ANALYSIS / 'reference_binned_trend.csv', float_precision='round_trip')
TASKS = ['rmsf', 'chi_entropy']
TN = {'rmsf': 'RMSF', 'chi_entropy': 'Hχ'}
UN = {'rmsf': 'Å', 'chi_entropy': 'nats'}
METHODS = ['Full', 'Intra', 'Static', 'AA']
MC = dict(zip(METHODS, [CATEGORICAL[0], CATEGORICAL[3], CATEGORICAL[2], CATEGORICAL[5]]))
LS = dict(zip(METHODS, ['-', '--', '-.', ':']))
CMAP = LinearSegmentedColormap.from_list('record_weighted_blue', ['#E4EDF5', '#6BAED6', '#08306B'])
CMAP.set_bad('white')
AA = list('AGCSTVNDHILFPWYQEMRK')
points = {t: np.load(ANALYSIS / f'prediction_reference_{t}.npz', allow_pickle=False) for t in TASKS}
nonzero = np.concatenate([points[t]['histogram_mass'][points[t]['histogram_mass'] > 0] for t in TASKS])
norm = LogNorm(vmin=nonzero.min(), vmax=nonzero.max())
fig = plt.figure(figsize=(183 / 25.4, 190 / 25.4))
xpos = [.107, .601]
width = .365
rects = [[x, y, width, h] for y, h in [(.703, .199), (.395, .202), (.075, .164)] for x in xpos]
axs = [fig.add_axes(r) for r in rects]

def label(ax, letter, title):
    # The first row reserves a distinct subtitle line for the numerical metrics.
    y = 1.15 if letter in 'ab' else 1.08
    ax.text(0, y, title, transform=ax.transAxes, fontsize=9.5, ha='left', va='bottom')
    ax.text(-.18, y, letter, transform=ax.transAxes, fontsize=11, fontweight='bold', ha='left', va='bottom')

def fmt(ax, axis):
    obj = getattr(ax, axis + 'axis')
    obj.set_major_formatter(FuncFormatter(lambda v, p: format(v, 'g')))
    obj.set_minor_locator(NullLocator())

panel_info = {}
type_points = []
for j, target in enumerate(TASKS):
    z = points[target]
    meta = summary['metrics'][target]['Full']
    reference, prediction = z['reference'], z['prediction']
    assert len(reference) == len(prediction) == meta['n_positions']
    assert np.isfinite(reference).all() and np.isfinite(prediction).all()
    assert np.isclose(z['histogram_mass'].sum(), 1.0)
    ax = axs[j]
    mesh = ax.pcolormesh(z['x_edges'], z['y_edges'], np.ma.masked_equal(z['histogram_mass'].T, 0),
                         cmap=CMAP, norm=norm, shading='flat', rasterized=False, linewidth=0)
    lo = min(float(reference.min()), float(prediction.min()), float(z['x_edges'][0]), float(z['y_edges'][0]))
    hi = max(float(reference.max()), float(prediction.max()), float(z['x_edges'][-1]), float(z['y_edges'][-1]))
    ax.plot([lo, hi], [lo, hi], color='#888888', linestyle='--', linewidth=.8, zorder=3)
    tr = trend[trend.target.eq(target)]
    ax.fill_between(tr.x_mean, tr.ci_low, tr.ci_high, color=BLACK, alpha=.12, linewidth=0, zorder=4)
    ax.plot(tr.x_mean, tr.prediction_mean, color=BLACK, linewidth=1.2, zorder=5)
    ax.set(xlim=(lo, hi), ylim=(lo, hi), xlabel=f'MD reference ({UN[target]})', ylabel=f'Full prediction ({UN[target]})')
    ax.set_aspect('equal', adjustable='box', anchor='W')
    if target == 'rmsf':
        thresh = float(z['density_linthresh']) if 'density_linthresh' in z else 1.0
        ax.set_xscale('symlog', linthresh=thresh)
        ax.set_yscale('symlog', linthresh=thresh)
        ticks = [v for v in [0, 1, 10, 30, 100, 1000] if lo <= v <= hi]
        ax.set_xticks(ticks); ax.set_yticks(ticks); fmt(ax, 'x'); fmt(ax, 'y')
    else:
        ax.set_xticks([0, 1, 2, 3, 4]); ax.set_yticks([0, 1, 2, 3, 4])
    label(ax, 'ab'[j], TN[target] + ': numerical agreement')
    ax.text(0, 1.025, f'MAE = {meta["macro_mae"]:.3f} {UN[target]}   ρ = {meta["mean_spearman"]:.3f}',
            transform=ax.transAxes, fontsize=9, ha='left', va='bottom', color='#444444')
    panel_info['ab'[j]] = {'positions': len(reference), 'records': meta['n_records'],
                           'density_sum': float(z['histogram_mass'].sum()),
                           'range_reference': [float(reference.min()), float(reference.max())],
                           'range_prediction': [float(prediction.min()), float(prediction.max())],
                           'xy_limits': [lo, hi], 'nonzero_density_cells': int((z['histogram_mass'] > 0).sum())}

    ax = axs[2 + j]
    legend = []
    for method in METHODS:
        curve = pr[pr.target.eq(target) & pr.method.eq(method)]
        m = summary['metrics'][target][method]
        assert len(curve) == 101
        # Preserve the 101-point first-reached-recall grid, then append the original
        # all-selected endpoint at the same recall=1 (all score thresholds exhausted).
        endpoint = endpoints[endpoints.target.eq(target) & endpoints.method.eq(method)].iloc[0]
        assert endpoint.recall == 1 and np.isclose(endpoint.precision, m['prevalence_mean'])
        xx = np.r_[curve.recall.to_numpy(), endpoint.recall]
        yy = np.r_[curve.precision.to_numpy(), endpoint.precision]
        ax.step(xx, yy, where='post', color=MC[method], ls=LS[method], lw=1.5)
        # AP is averaged over original per-record ranking metrics, not the area under this mean curve.
        legend.append(Line2D([], [], color=MC[method], ls=LS[method], lw=1.5,
                              label=f'{method} AP {m["mean_ap"]:.3f}'))
    base = float(pr[pr.target.eq(target)].no_skill_precision.iloc[0])
    ax.axhline(base, color='#AAAAAA', ls=(0, (2, 2)), lw=.7, zorder=0)
    ax.set(xlim=(-.005, 1.005), ylim=(0, 1.02), xlabel='Recall', ylabel='Precision')
    ax.set_xticks([0, .5, 1]); ax.set_yticks([0, .5, 1])
    ax.grid(axis='both', color='#E0E0E0', linewidth=.3, alpha=.5, zorder=0)
    label(ax, 'cd'[j], TN[target] + ': high-dynamics retrieval')
    ax.legend(handles=legend, loc='upper center', bbox_to_anchor=(.51, -.25), ncol=2,
              handlelength=1.5, columnspacing=.8, handletextpad=.4, borderaxespad=0, labelspacing=.25)
    panel_info['cd'[j]] = {'methods': METHODS, 'recall_grid_vertices_per_method': 101, 'all_selected_endpoint_appended': True,
                           'all_selected_endpoint_recall': 1.0, 'all_selected_endpoint_precision': base, 'no_skill_precision': base,
                           'n_records': int(curve.n_records.iloc[0]), 'all_records_same_eligibility': True}

    ax = axs[4 + j]
    for k, method in enumerate(['intrachain', 'static_sequence_complex']):
        d = paired[paired.metric.eq(target) & paired.comparator.eq(method)].set_index('amino_acid').reindex(AA)
        valid = d.relative_MAE_reduction_percent.notna().to_numpy()
        x = np.arange(20) + (k - .5) * .38
        y = d.relative_MAE_reduction_percent.to_numpy()
        key = 'Intra' if k == 0 else 'Static'
        ax.scatter(x[valid], y[valid], s=17, marker='o' if k == 0 else '^', color=MC[key], edgecolors='none', zorder=3)
        for aa, xx, yy in zip(np.array(AA)[valid], x[valid], y[valid]):
            type_points.append({'panel': 'ef'[j], 'metric': target, 'method': method, 'amino_acid': aa, 'x': xx, 'y': yy})
    ax.axhline(0, color=GREY, lw=.7, ls=':')
    ax.set_xlim(-.65, 19.65); ax.set_xticks(range(20), AA)
    ax.tick_params(axis='x', length=2.5, pad=4)
    ax.set_xlabel('Amino acid'); ax.set_ylabel('Full MAE reduction (%)')
    if target == 'rmsf':
        ax.set_ylim(-1, 21); ax.set_yticks([0, 5, 10, 15, 20])
    else:
        ax.set_yscale('symlog', linthresh=10); ax.set_ylim(-5.5, 130)
        ax.set_yticks([-5, 0, 5, 10, 30, 100]); fmt(ax, 'y')
        ax.text(.025, .9, 'A/G: no χ', transform=ax.transAxes, ha='left', va='top', fontsize=9, color='#777777')
        for tick in ax.get_xticklabels()[:2]: tick.set_color('#999999')
    label(ax, 'ef'[j], TN[target] + ': within-type gain')

# Shared density key and source line key outside the data regions.
cbax = fig.add_axes([.727, .979, .22, .006])
cb = fig.colorbar(mesh, cax=cbax, orientation='horizontal')
cb.solids.set_rasterized(False)
cb.ax.xaxis.set_ticks_position('bottom')
cb.ax.tick_params(axis='x', pad=0, length=2)
cb.locator = mpl.ticker.LogLocator(base=10, numticks=3); cb.update_ticks()
cb.ax.xaxis.set_minor_locator(NullLocator()); cb.outline.set_visible(False)
fig.text(.708, .982, 'Cell probability', fontsize=9, ha='right', va='center')
fig.legend(handles=[Line2D([], [], color='#888888', ls='--', lw=.8, label='Ideal: prediction = reference'),
                    Line2D([], [], color=BLACK, lw=1.2, label='Reference-bin trend (95% CI)')],
           loc='center', bbox_to_anchor=(.548, .643), ncol=2, handlelength=2, columnspacing=1, handletextpad=.5)
fig.legend(handles=[Line2D([], [], marker='o', ls='none', ms=4, color=MC['Intra'], label='vs Intra'),
                    Line2D([], [], marker='^', ls='none', ms=4, color=MC['Static'], label='vs Static')],
           loc='center', bbox_to_anchor=(.55, .287), ncol=2, columnspacing=1.75, handletextpad=.35)

type_points = pd.DataFrame(type_points)
assert type_points.panel.value_counts().to_dict() == {'e': 40, 'f': 36}
assert (type_points.loc[type_points.panel.eq('f'), 'y'] < 0).sum() == 2
type_points.to_csv(ROOT / 'source_type_visible_points.csv', index=False)
fig.canvas.draw()
renderer = fig.canvas.get_renderer()
ignore = set()
for ax in fig.axes:
    for axis in [ax.xaxis, ax.yaxis]:
        lo, hi = sorted(axis.get_view_interval())
        for tick in axis.get_major_ticks() + axis.get_minor_ticks():
            if tick.get_loc() < lo - 1e-10 or tick.get_loc() > hi + 1e-10:
                ignore.update([tick.label1, tick.label2])
texts = [t for t in fig.findobj(mpl.text.Text) if t.get_visible() and t.get_text() and t not in ignore]
outside = []
for text in texts:
    bb = text.get_window_extent(renderer)
    if bb.x0 < -.5 or bb.y0 < -.5 or bb.x1 > fig.bbox.width + .5 or bb.y1 > fig.bbox.height + .5:
        outside.append(text.get_text())
assert not outside, outside
assert min(t.get_fontsize() for t in texts) >= 9
assert not any('\u4e00' <= c <= '\u9fff' for t in texts for c in t.get_text())
for ext in ['png', 'pdf', 'svg']:
    fig.savefig(ROOT / f'fig2_predictive_capability.{ext}', bbox_inches=fig.bbox_inches, dpi=450)
svg = ET.parse(ROOT / 'fig2_predictive_capability.svg')
assert not svg.findall('.//{http://www.w3.org/2000/svg}image')
visible_text = [t.get_text() for t in texts]
(ROOT / 'visible_text.json').write_text(json.dumps(visible_text, ensure_ascii=False, indent=2), encoding='utf8')
audit = {'status': 'numerical_and_canvas_passed_visual_pending', 'width_mm': 183, 'height_mm': 190,
         'dpi': 450, 'minimum_font_pt': min(t.get_fontsize() for t in texts),
         'minimum_font_at_manuscript_width_pt': min(t.get_fontsize() for t in texts) * 167.64 / 183,
         'panel_info': panel_info, 'type_points': len(type_points), 'type_negative_points': 2,
         'source_type_csv_sha256': sha(paired_path), 'type_csv_archive_sha256': sha(ARCHIVE / paired_path.name),
         'svg_raster_image_count': 0, 'canvas_text_outside': outside, 'english_only_visible_text': True,
         'no_downsampling': True, 'no_axis_tail_clipping': True,
         'source_summary_sha256': sha(ANALYSIS / 'summary.json'), 'source_pr_sha256': sha(ANALYSIS / 'macro_pr_curve.csv'),
         'source_endpoint_sha256': sha(ANALYSIS / 'all_selected_endpoint.csv'),
         'source_trend_sha256': sha(ANALYSIS / 'reference_binned_trend.csv')}
(ROOT / 'data_audit.json').write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding='utf8')
print(json.dumps(audit, ensure_ascii=False, indent=2))
