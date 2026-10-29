# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a,b) record-weighted full density -> KernelDensity/plot_CorrelationDensity1.R -> param inherit
# (c,d,e) macro precision-recall -> cross-type AUROC/plot_AUROC2.py -> param inherit
# (f) annotated source/task AP matrix -> heatmap/plot_composition.py -> param inherit
# RULE: no native-run panels; unweighted R KDE, ROC and raw-count heatmap have incompatible statistical semantics.
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
from matplotlib.colors import LinearSegmentedColormap, LogNorm, Normalize
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter, NullLocator

mpl.rcParams.update({'font.size': 9, 'axes.titlesize': 10, 'axes.labelsize': 9,
                     'xtick.labelsize': 9, 'ytick.labelsize': 9, 'legend.fontsize': 9,
                     'axes.unicode_minus': False})
ANALYSIS = ROOT.parent / 'analysis'
V44 = EXP / 'predictive_figure_v44_20261018' / 'analysis'
ARCHIVE = ROOT / 'source_archive'
ARCHIVE.mkdir(exist_ok=True)
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
for name in ['summary.json', 'macro_pr_curve.csv', 'all_selected_endpoint.csv',
             'prediction_reference_rmsf.npz', 'prediction_reference_chi_entropy.npz']:
    shutil.copy2(V44 / name, ARCHIVE / name)
    assert sha(V44 / name) == sha(ARCHIVE / name)
for name in ['four_method_binned_trend.csv', 'within_type_macro_pr.csv',
             'within_type_all_selected_endpoint.csv', 'source_ap_matrix.csv']:
    shutil.copy2(ANALYSIS / name, ARCHIVE / name)
    assert sha(ANALYSIS / name) == sha(ARCHIVE / name)
summary = json.loads((ARCHIVE / 'summary.json').read_text(encoding='utf8'))
pr = pd.read_csv(ARCHIVE / 'macro_pr_curve.csv', float_precision='round_trip')
endpoints = pd.read_csv(ARCHIVE / 'all_selected_endpoint.csv', float_precision='round_trip')
within_pr = pd.read_csv(ARCHIVE / 'within_type_macro_pr.csv', float_precision='round_trip')
within_end = pd.read_csv(ARCHIVE / 'within_type_all_selected_endpoint.csv', float_precision='round_trip')
trend = pd.read_csv(ARCHIVE / 'four_method_binned_trend.csv', float_precision='round_trip')
ap_source = pd.read_csv(ARCHIVE / 'source_ap_matrix.csv', float_precision='round_trip')
METHODS = ['Full', 'Intra', 'Static', 'AA']
TASKS = ['rmsf', 'chi_entropy', 'within_type_hchi']
MC = dict(zip(METHODS, [CATEGORICAL[0], CATEGORICAL[3], CATEGORICAL[2], CATEGORICAL[5]]))
LS = dict(zip(METHODS, ['-', '--', '-.', ':']))
points = {t: np.load(ARCHIVE / f'prediction_reference_{t}.npz', allow_pickle=False) for t in TASKS[:2]}
nonzero = np.concatenate([points[t]['histogram_mass'][points[t]['histogram_mass'] > 0] for t in TASKS[:2]])
density_cmap = LinearSegmentedColormap.from_list('neutral_density', ['#EEEEEE', '#C6C6C6', '#5A5A5A'])
density_cmap.set_bad('white')
density_norm = LogNorm(vmin=nonzero.min(), vmax=nonzero.max())
ap_cmap = LinearSegmentedColormap.from_list('mean_AP', SEQUENTIAL)
W, H = 183., 190.
fig = plt.figure(figsize=(W / 25.4, H / 25.4))
def axes_mm(x, y, w, h):
    return fig.add_axes([x / W, y / H, w / W, h / H])
def text_mm(x, y, text, **kwargs):
    return fig.text(x / W, y / H, text, **kwargs)
def heading(x, y, letter, title):
    text_mm(x - 10, y, letter, fontsize=11, fontweight='bold', ha='left', va='bottom')
    text_mm(x, y, title, fontsize=10, ha='left', va='bottom')
def fmt(ax, axis):
    obj = getattr(ax, axis + 'axis')
    obj.set_major_formatter(FuncFormatter(lambda v, p: format(v, 'g')))
    obj.set_minor_locator(NullLocator())
panel_info = {}
visible_sources = []
density_axes = []
for j, target in enumerate(TASKS[:2]):
    x = [15, 107][j]
    ax = axes_mm(x, 112, 65, 65)
    density_axes.append(ax)
    z = points[target]
    meta = summary['metrics'][target]['Full']
    ref, pred = z['reference'], z['prediction']
    assert len(ref) == len(pred) == meta['n_positions']
    assert np.isfinite(ref).all() and np.isfinite(pred).all()
    assert np.isclose(z['histogram_mass'].sum(), 1.0)
    mesh = ax.pcolormesh(z['x_edges'], z['y_edges'], np.ma.masked_equal(z['histogram_mass'].T, 0),
                         cmap=density_cmap, norm=density_norm, shading='flat', rasterized=False,
                         linewidth=0, zorder=1)
    lo = min(float(ref.min()), float(pred.min()), float(z['x_edges'][0]), float(z['y_edges'][0]))
    hi = max(float(ref.max()), float(pred.max()), float(z['x_edges'][-1]), float(z['y_edges'][-1]))
    ax.plot([lo, hi], [lo, hi], color='#333333', linestyle=(0, (3, 3)), linewidth=.75, zorder=2)
    for method in METHODS[::-1]:
        tr = trend[trend.target.eq(target) & trend.method.eq(method)].sort_values('bin_index')
        assert len(tr) == 20
        populated = tr.n_positions.gt(0)
        assert np.isfinite(tr.loc[populated, 'prediction_mean']).all()
        assert tr.loc[~populated, 'prediction_mean'].isna().all()
        if method == 'Full':
            assert np.isfinite(tr.loc[populated, 'ci_low']).all() and np.isfinite(tr.loc[populated, 'ci_high']).all()
            ax.fill_between(tr.x_mean, tr.ci_low, tr.ci_high, color=MC[method], alpha=.15, linewidth=0, zorder=3)
        ax.plot(tr.x_mean, tr.prediction_mean, color=MC[method], ls=LS[method], linewidth=1.45,
                zorder=4 if method != 'Full' else 5)
    un = 'Å' if j == 0 else 'nats'
    ax.set(xlim=(lo, hi), ylim=(lo, hi), xlabel=f'MD reference ({un})', ylabel=f'Prediction ({un})')
    ax.set_aspect('equal', adjustable='box')
    if j == 0:
        ax.set_xscale('symlog', linthresh=1.0)
        ax.set_yscale('symlog', linthresh=1.0)
        ax.set_xticks([0, 1, 10, 30]); ax.set_yticks([0, 1, 10, 30]); fmt(ax, 'x'); fmt(ax, 'y')
    else:
        ax.set_xticks([0, 1, 2, 3, 4]); ax.set_yticks([0, 1, 2, 3, 4])
    heading(x, 179.3, 'ab'[j], ['RMSF: agreement', 'Hχ: agreement'][j])
    ax.text(.035, .97, f'Full MAE {meta["macro_mae"]:.3f} {un}\nρ {meta["mean_spearman"]:.3f} · n = {len(ref):,}',
            transform=ax.transAxes, ha='left', va='top', fontsize=9, linespacing=1.4,
            bbox=dict(facecolor='white', edgecolor='none', alpha=.9, pad=1.5))
    panel_info['ab'[j]] = {'positions': len(ref), 'records': meta['n_records'],
                           'density_sum': float(z['histogram_mass'].sum()),
                           'nonzero_density_cells': int((z['histogram_mass'] > 0).sum()),
                           'reference_range': [float(ref.min()), float(ref.max())],
                           'prediction_range': [float(pred.min()), float(pred.max())],
                           'xy_limits': [lo, hi], 'trend_methods': METHODS, 'trend_vertices': 80}

# A single semantic method key is reused in density trends, all PR panels and the source matrix.
fig.legend(handles=[Line2D([], [], color=MC[m], ls=LS[m], lw=1.5, label=m) for m in METHODS],
           loc='center', bbox_to_anchor=(70 / W, 187 / H), ncol=4, columnspacing=1.1,
           handlelength=2, handletextpad=.45, borderaxespad=0)
text_mm(132, 187, 'Mass', fontsize=9, ha='right', va='center', color='#555555')
cbax = axes_mm(136, 186, 36, 1.6)
cb = fig.colorbar(mesh, cax=cbax, orientation='horizontal')
cb.solids.set_rasterized(False)
cb.outline.set_visible(False)
cb.ax.xaxis.set_ticks_position('bottom')
cb.set_ticks([1e-5, 1e-3, 1e-2]); cb.ax.xaxis.set_minor_locator(NullLocator())
cb.ax.tick_params(labelsize=9, pad=.5, length=1.5)

pr_axes = []
for j, target in enumerate(TASKS):
    x = [15, 72, 129][j]
    ax = axes_mm(x, 54, 46, 38)
    pr_axes.append(ax)
    frame = within_pr if j == 2 else pr[pr.target.eq(target)]
    ends = within_end if j == 2 else endpoints[endpoints.target.eq(target)]
    mm = summary['within_type'] if j == 2 else summary['metrics'][target]
    for method in METHODS:
        cur = frame[frame.method.eq(method)].sort_values('recall')
        end = ends[ends.method.eq(method)].iloc[0]
        assert len(cur) == 101 and end.recall == 1
        expected = mm[method]['mean_record_type_prevalence'] if j == 2 else mm[method]['prevalence_mean']
        assert np.isclose(end.precision, expected)
        xx, yy = np.r_[cur.recall.to_numpy(), end.recall], np.r_[cur.precision.to_numpy(), end.precision]
        ax.step(xx, yy, where='post', color=MC[method], ls=LS[method], lw=1.4,
                zorder=6 - METHODS.index(method))
        visible_sources.extend({'panel': 'cde'[j], 'task': target, 'method': method,
                                'recall': float(r), 'precision': float(p), 'is_final_endpoint': i == 101}
                               for i, (r, p) in enumerate(zip(xx, yy)))
    base = float(frame.no_skill_precision.iloc[0])
    assert np.isclose(base, expected)
    ax.axhline(base, color='#AAAAAA', ls=(0, (2, 2)), lw=.7, zorder=0)
    ax.set(xlim=(0, 1), ylim=(0, 1.02), xlabel='Recall', ylabel='Precision' if j == 0 else '')
    ax.set_xticks([0, .5, 1]); ax.set_yticks([0, .5, 1])
    if j > 0:
        ax.set_yticklabels([])
    ax.tick_params(pad=2.5)
    heading(x, 98.4, 'cde'[j], ['RMSF', 'Hχ', 'Hχ within type'][j])
    text_mm(x, 94.2, ['MD ≥ 90th percentile', 'MD ≥ 90th percentile', 'MD ≥ within-type P90'][j],
            ha='left', va='bottom', fontsize=9, color='#555555')
    # Values use the shared Full/Intra/Static/AA ordering. These are means of exact APs,
    # never areas under these macro curves.
    text_mm(x - 2, 41.7, 'AP', fontsize=9, ha='left', va='center', color='#444444')
    for k, method in enumerate(METHODS):
        text_mm(x + 8 + k * 11.7, 41.7, f'{mm[method]["mean_ap"]:.3f}', color=MC[method],
                fontsize=9, ha='center', va='center')
    panel_info['cde'[j]] = {'task': target, 'methods': METHODS, 'n_records': 3063,
                           'vertices_per_method': 102, 'original_endpoint_appended': True,
                           'all_selected_precision': base,
                           'n_eligible_record_types': 42108 if j == 2 else None,
                           'mean_AP': {m: mm[m]['mean_ap'] for m in METHODS}}

# A compact 4 x 9 matrix answers source consistency, with exactly one metric and
# a common 0..1 scale. This uses source strata of existing out-of-fold records.
axh = axes_mm(24, 3.5, 147, 20)
SOURCES = ['atlas', 'dynarepo', 'itrace']
SOURCE_NAMES = ['ATLAS', 'DynaRepo', 'ITrace']
matrix = np.empty((4, 9))
for i, method in enumerate(METHODS):
    for j, task in enumerate(TASKS):
        for k, source in enumerate(SOURCES):
            cell = ap_source[ap_source.task.eq(task) & ap_source.dataset.eq(source) & ap_source.method.eq(method)]
            assert len(cell) == 1
            matrix[i, 3 * j + k] = float(cell.mean_ap.iloc[0])
assert np.isfinite(matrix).all() and ((matrix >= 0) & (matrix <= 1)).all()
im = axh.pcolormesh(np.arange(10), np.arange(5), matrix, cmap=ap_cmap, norm=Normalize(0, 1),
                    linewidth=.8, edgecolor='white', rasterized=False)
axh.set(xlim=(0, 9), ylim=(4, 0))
axh.set_xticks([]); axh.set_yticks(np.arange(4) + .5, METHODS)
axh.tick_params(axis='y', length=0, pad=5)
for t, method in zip(axh.get_yticklabels(), METHODS):
    t.set_color(MC[method]); t.set_fontweight('bold' if method == 'Full' else 'normal')
for spine in axh.spines.values(): spine.set_visible(False)
for i in range(4):
    for j in range(9):
        axh.text(j + .5, i + .5, f'{matrix[i, j]:.3f}', ha='center', va='center', fontsize=9,
                 color='white' if matrix[i, j] >= .63 else '#222222')
for xx in [3, 6]:
    axh.axvline(xx, color='white', lw=2)
text_mm(5, 33.3, 'f', fontsize=11, fontweight='bold', ha='left', va='bottom')
text_mm(15, 33.3, 'Mean AP', fontsize=9, ha='left', va='bottom')
for j, task_title in enumerate(['RMSF', 'Hχ', 'Hχ within type']):
    text_mm(24 + (j * 3 + 1.5) * 147 / 9, 33.3, task_title,
            fontsize=9.5, ha='center', va='bottom')
    for k, name in enumerate(SOURCE_NAMES):
        source = SOURCES[k]
        nvals = ap_source[ap_source.task.eq(TASKS[j]) & ap_source.dataset.eq(source)].n_records.unique()
        assert len(nvals) == 1
        text_mm(24 + (j * 3 + k + .5) * 147 / 9, 24.9,
                f'{name}\nn = {int(nvals[0]):,}', fontsize=9, ha='center', va='bottom', linespacing=1.2)
cb2ax = axes_mm(174, 3.5, 1.9, 20)
cb2 = fig.colorbar(im, cax=cb2ax)
cb2.solids.set_rasterized(False); cb2.outline.set_visible(False)
cb2.set_ticks([0, .5, 1]); cb2.ax.tick_params(length=1.5, pad=2, labelsize=9)
panel_info['f'] = {'cells': 36, 'methods': METHODS, 'tasks': TASKS, 'sources': SOURCES,
                    'n_records_by_source': dict(zip(SOURCES, [1929, 890, 244])),
                    'metric': 'mean exact per-record AP', 'color_limits': [0, 1],
                    'not_source_held_out': True}
pd.DataFrame(visible_sources).to_csv(ROOT / 'source_visible_pr_vertices.csv', index=False)

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
for t in texts:
    bb = t.get_window_extent(renderer)
    if bb.x0 < -.5 or bb.y0 < -.5 or bb.x1 > fig.bbox.width + .5 or bb.y1 > fig.bbox.height + .5:
        outside.append(t.get_text())
assert not outside, outside
assert min(t.get_fontsize() for t in texts) >= 9
assert not any('\u4e00' <= c <= '\u9fff' for t in texts for c in t.get_text())
for ext in ['png', 'pdf', 'svg']:
    fig.savefig(ROOT / f'fig2_dense_capability.{ext}', bbox_inches=fig.bbox_inches, dpi=450)
svg = ET.parse(ROOT / 'fig2_dense_capability.svg')
assert not svg.findall('.//{http://www.w3.org/2000/svg}image')
(ROOT / 'visible_text.json').write_text(json.dumps([t.get_text() for t in texts], ensure_ascii=False, indent=2), encoding='utf8')
audit = {'status': 'numerical_and_canvas_passed_visual_pending', 'width_mm': W, 'height_mm': H,
         'dpi': 450, 'minimum_font_pt': min(t.get_fontsize() for t in texts),
         'minimum_font_at_manuscript_width_pt': min(t.get_fontsize() for t in texts) * 167.64 / W,
         'panel_info': panel_info, 'svg_raster_image_count': 0, 'canvas_text_outside': outside,
         'english_only_visible_text': True, 'no_downsampling': True, 'no_axis_tail_clipping': True,
         'removed_original_type_dot_panels': True,
         'source_hashes': {p.name: sha(p) for p in ARCHIVE.iterdir() if p.is_file()},
         'output_hashes': {ext: sha(ROOT / f'fig2_dense_capability.{ext}') for ext in ['png', 'pdf', 'svg']}}
(ROOT / 'data_audit.json').write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding='utf8')
print(json.dumps(audit, ensure_ascii=False, indent=2))
