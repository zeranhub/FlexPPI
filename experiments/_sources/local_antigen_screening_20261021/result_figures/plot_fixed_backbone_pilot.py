# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a,b) Held-position-out prediction scatter -> PairedBoxScatter/plot_PairedBoxScatter1.R -> param inherit
# (c) Seven-method MAE dots -> BarComparison/plot_comparison_GeneRegulatory.py -> param inherit
# (d) Signed annotated matrix -> heatmap/plot_comparison.py -> param inherit
# RULE: no native-run panels. Assets contain unrelated synthetic/gene-regulatory
#       data and inference comparisons, so their whole panels are not executable
#       with this candidate-level heldout schema. Inherit only clean spines,
#       individual-point connectors and annotated-cell contrast; no synthetic data.

import os
_dll_handle = os.add_dll_directory(r'D:\Anaconda\envs\pvb_cpu\Library\bin')

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

from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.text import Text
from PIL import Image

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
RESULTS = ROOT / 'evaluation_review/results_fixed_protocol'
ANALYSIS = 'within_system_leave_position_out'
MAIN_LABEL = 'reported_hybrid'
BASE = 'sequence_plus_static'
JOINT = 'sequence_plus_static_plus_RH'
METHODS = ['zero_change', 'sequence_only', 'static_only', BASE,
           'sequence_plus_static_plus_R', 'sequence_plus_static_plus_H', JOINT]
METHOD_LABELS = ['Zero change', 'Sequence', 'Static', 'Sequence + static',
                 'Baseline + R', 'Baseline + H', 'Baseline + R + H']
SYSTEMS = ['nyeso', 'tax']
SYSTEM_COLOR = {'nyeso': '#CF8C27', 'tax': '#762A83'}
SYSTEM_NAME = {'nyeso': '1G4', 'tax': 'A6'}
TEMPLATES = ['2BNQ', '2BNR', '1AO7', '1QSE', '1QSF']
LABELS = ['reported_hybrid', 'all_constrained', 'all_fitted']
mpl.rcParams.update({key: 9 for key in [
    'font.size', 'axes.titlesize', 'axes.labelsize', 'xtick.labelsize',
    'ytick.labelsize', 'legend.fontsize', 'figure.titlesize']})
mpl.rcParams['savefig.bbox'] = None  # Preserve exact requested 183 mm canvas.


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')


def main():
    source = HERE / 'source_data'
    source.mkdir(exist_ok=True)
    originals = {name: RESULTS / name for name in [
        'candidate_predictions.csv', 'candidate_level_metrics.csv', 'dynamic_increment_comparisons.csv']}
    archived = {}
    for name, path in originals.items():
        target = source / name
        shutil.copy2(path, target)
        assert sha(path) == sha(target)
        archived[name] = {'original': str(path), 'archive': str(target), 'SHA256': sha(path)}
    predictions = pd.read_csv(originals['candidate_predictions.csv'], float_precision='round_trip')
    metrics = pd.read_csv(originals['candidate_level_metrics.csv'], float_precision='round_trip')
    increment = pd.read_csv(originals['dynamic_increment_comparisons.csv'], float_precision='round_trip')
    primary_predictions = predictions[(predictions.analysis == ANALYSIS) &
        (predictions.label_variant == MAIN_LABEL) & predictions.is_primary_template].copy()
    primary_metrics = metrics[(metrics.analysis == ANALYSIS) &
        (metrics.label_variant == MAIN_LABEL) & metrics.is_primary_template].copy()
    sensitivity = increment[increment.analysis == ANALYSIS].copy()
    scatter = primary_predictions[primary_predictions.model.isin([BASE, JOINT])].copy()
    assert len(primary_predictions) == 26 * 7 and len(primary_metrics) == 14
    assert len(scatter) == 52 and len(sensitivity) == 15
    assert set(primary_metrics.model) == set(METHODS)
    assert set(sensitivity.template_pdb) == set(TEMPLATES)
    assert set(sensitivity.label_variant) == set(LABELS)
    assert (primary_predictions.analysis_set == 'primary_single_substitution').all()
    assert primary_predictions.groupby(['system', 'model']).size().eq(13).all()
    assert primary_predictions.groupby(['system', 'model']).candidate_id.nunique().eq(13).all()
    assert set(primary_predictions.system) == set(SYSTEMS)
    assert (primary_metrics.n_candidates == 13).all() and (sensitivity.n_candidates == 13).all()
    assert not primary_metrics.formal_significance_test_performed.any()
    leakage = []
    for row in primary_predictions.itertuples():
        train = json.loads(row.training_candidate_ids)
        assert row.candidate_id not in train
        assert json.loads(row.mutation_positions) == [int(row.heldout_mutation_position)]
        leakage.append({'candidate_id': row.candidate_id, 'model': row.model,
            'heldout_mutation_position': int(row.heldout_mutation_position), 'candidate_not_in_training': True})
    for (system, model), group in primary_predictions.groupby(['system', 'model']):
        published = primary_metrics[(primary_metrics.system == system) & (primary_metrics.model == model)].iloc[0]
        mae = np.abs(group.predicted_delta_pKD.to_numpy() - group.true_delta_pKD.to_numpy()).mean()
        assert abs(float(mae) - published.MAE_delta_pKD) < 1e-12
        assert np.isfinite(group[['predicted_delta_pKD', 'true_delta_pKD']].to_numpy()).all()
    assert np.allclose(sensitivity.joint_increment_vs_baseline,
                       sensitivity.baseline_MAE - sensitivity.joint_MAE, rtol=0, atol=1e-14)
    scatter.to_csv(source / 'panels_ab_all52_predictions.csv', index=False, float_format='%.17g')
    primary_metrics.to_csv(source / 'panel_c_all14_method_MAE.csv', index=False, float_format='%.17g')
    sensitivity.to_csv(source / 'panel_d_all15_sensitivity_cells.csv', index=False, float_format='%.17g')
    pairs = []
    fig = plt.figure(figsize=(183 / 25.4, 196 / 25.4), facecolor='white')
    axes = []
    ticks = [-2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5]
    limits = [-2.7, .55]
    for k, system in enumerate(SYSTEMS):
        left = .105 if k == 0 else .605
        ax = fig.add_axes([left, .610, .355, .331])
        axes.append(ax)
        ax.set_aspect('equal', adjustable='box')
        data = scatter[scatter.system == system]
        wide = data.pivot(index=['candidate_id', 'sequence', 'true_delta_pKD'],
                          columns='model', values='predicted_delta_pKD').reset_index()
        assert len(wide) == 13
        color = SYSTEM_COLOR[system]
        for row in wide.itertuples(index=False):
            observed = row.true_delta_pKD
            baseline = getattr(row, BASE)
            joint = getattr(row, JOINT)
            ax.plot([observed, observed], [baseline, joint], color='#B5B5B5', lw=.65, zorder=1)
            pairs.append({'system': system, 'candidate_id': row.candidate_id, 'sequence': row.sequence,
                'observed_delta_pKD': observed, 'baseline_prediction': baseline, 'RH_prediction': joint})
        ax.scatter(wide.true_delta_pKD, wide[BASE], s=26, facecolors='white',
                   edgecolors=color, linewidths=1.05, marker='o', zorder=3)
        ax.scatter(wide.true_delta_pKD, wide[JOINT], s=20, color=color,
                   marker='s', linewidths=.7, zorder=4)
        ax.plot(limits, limits, color='#A3A3A3', ls=(0, (3, 3)), lw=.65, zorder=0)
        ax.axhline(0, color='#DDDDDD', lw=.45, zorder=0)
        ax.axvline(0, color='#DDDDDD', lw=.45, zorder=0)
        ax.set(xlim=limits, ylim=limits, xticks=ticks[::2], yticks=ticks[::2],
               xlabel='Measured ΔpKD', ylabel='Held-position-out ΔpKD')
        ax.tick_params(length=2.7, pad=3)
        label = 'a' if system == 'nyeso' else 'b'
        title = '1G4 · NY-ESO-1' if system == 'nyeso' else 'A6 · Tax'
        fig.text(left-.077, .979, label, weight='bold', fontsize=11, va='top')
        fig.text(left, .979, title, color=color, fontsize=10, va='top')
        template = '2BNQ' if system == 'nyeso' else '1AO7'
        ax.text(.015, .975, template + ' · n = 13', transform=ax.transAxes, va='top', fontsize=9)
    legend_handles = [
        Line2D([], [], marker='o', mfc='white', mec='#555555', color='none', ms=5,
               label='Sequence + static (baseline)'),
        Line2D([], [], marker='s', mfc='#555555', mec='#555555', color='none', ms=4.5,
               label='Baseline + R + H')]
    fig.legend(handles=legend_handles, loc='center', bbox_to_anchor=(.53, .540),
        ncol=2, columnspacing=2.0, handletextpad=.45, borderaxespad=0)
    fig.text(.53, .513, 'R = RMSF; H = Hχ     Dashed line: prediction = measurement',
             ha='center', color='#555555', fontsize=9)
    pd.DataFrame(pairs).to_csv(source / 'panels_ab_all26_candidate_pairs.csv', index=False, float_format='%.17g')
    # Complete prespecified method comparison: no selection of the best model.
    fig.text(.028, .475, 'c', weight='bold', fontsize=11, va='top')
    fig.text(.105, .475, 'All seven prespecified methods', fontsize=9.5, va='top')
    cax = fig.add_axes([.258, .180, .222, .250])
    axes.append(cax)
    for j, system in enumerate(SYSTEMS):
        data = primary_metrics[primary_metrics.system == system].set_index('model').loc[METHODS]
        y = np.arange(len(METHODS)) + (-.12 if j == 0 else .12)
        cax.scatter(data.MAE_delta_pKD, y, s=24, color=SYSTEM_COLOR[system],
                    marker='o' if j == 0 else 's', linewidths=.6, zorder=3)
    cax.set(yticks=np.arange(len(METHODS)), yticklabels=METHOD_LABELS,
            ylim=(6.55, -.55), xlim=(.30, 1.06), xticks=[.4, .6, .8, 1.0],
            xlabel='MAE (ΔpKD)\nLower is better')
    cax.tick_params(axis='y', length=0, pad=5)
    cax.tick_params(axis='x', length=2.7)
    cax.axhline(3.5, color='#DDDDDD', lw=.6)
    for yindex in [3, 6]:
        cax.axhspan(yindex-.43, yindex+.43, color='#EFEFEF', lw=0, zorder=0)
    system_handles = [Line2D([], [], color=SYSTEM_COLOR[s], marker='o' if j == 0 else 's',
        ls='', ms=4.5, label=SYSTEM_NAME[s] + ' (n = 13)') for j, s in enumerate(SYSTEMS)]
    fig.legend(handles=system_handles, loc='center', bbox_to_anchor=(.286, .090),
        ncol=2, columnspacing=1.0, handletextpad=.35, borderaxespad=0)
    # Every fit x template combination; sign is explicitly baseline minus joint.
    fig.text(.535, .475, 'd', weight='bold', fontsize=11, va='top')
    fig.text(.585, .475, 'R + H gain across all label fits', fontsize=9.5, va='top')
    dax = fig.add_axes([.655, .195, .305, .210])
    axes.append(dax)
    matrix = sensitivity.pivot(index='template_pdb', columns='label_variant',
        values='joint_increment_vs_baseline').loc[TEMPLATES, LABELS].to_numpy()
    assert matrix.shape == (5, 3) and np.isfinite(matrix).all()
    vmax = .08
    assert np.abs(matrix).max() <= vmax
    cmap = LinearSegmentedColormap.from_list('MAE_gain', ['#B2182B', '#F7F7F7', '#2166AC'])
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0., vmax=vmax)
    mesh = dax.pcolormesh(np.arange(4), np.arange(6), matrix, cmap=cmap, norm=norm,
                          edgecolors='white', linewidth=.8, rasterized=False)
    dax.set(xlim=(0, 3), ylim=(5, 0), xticks=np.arange(3)+.5,
            xticklabels=['Reported\nhybrid', 'All\nconstrained', 'All\nBmax-fitted'],
            yticks=np.arange(5)+.5,
            yticklabels=['1G4  2BNQ*', '1G4  2BNR', 'A6  1AO7*', 'A6  1QSE', 'A6  1QSF'])
    dax.xaxis.set_ticks_position('top')
    dax.tick_params(axis='both', length=0, pad=5)
    for spine in dax.spines.values():
        spine.set_visible(False)
    dax.axhline(2, color='white', lw=2)
    for (i, j), value in np.ndenumerate(matrix):
        text_color = 'white' if abs(value) > .045 else '#222222'
        dax.text(j+.5, i+.5, f'{value:+.4f}', ha='center', va='center',
                 fontsize=9, color=text_color)
    color_ax = fig.add_axes([.655, .125, .305, .013])
    colorbar = fig.colorbar(mesh, cax=color_ax, orientation='horizontal', ticks=[-.08, 0., .08])
    colorbar.solids.set_rasterized(False)
    colorbar.ax.set_xticklabels(['−0.08', '0', '+0.08'])
    colorbar.ax.tick_params(length=2, pad=2, labelsize=9)
    colorbar.outline.set_visible(False)
    fig.text(.8075, .084, 'MAE(baseline) − MAE(+R+H)', ha='center', fontsize=9)
    fig.text(.8075, .057, 'Negative: worse · Positive: better', ha='center', color='#555555', fontsize=9)
    fig.text(.655, .168, '* Primary template', fontsize=9, color='#555555')
    fig.text(.028, .018, 'Historical 2431 model · Fixed-backbone retrospective pilot',
             fontsize=9, color='#555555')
    # Canvas/font/coordinate audit; checks visible text only and leaves data intact.
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    width, height = fig.canvas.get_width_height()
    visible_text = []
    for item in fig.findobj(Text):
        if not item.get_visible() or not item.get_text().strip():
            continue
        bbox = item.get_window_extent(renderer)
        if bbox.width == 0 or bbox.height == 0:
            continue
        visible_text.append({'text': item.get_text(), 'font_pt': item.get_fontsize(),
            'bounds_fraction': [bbox.x0/width, bbox.y0/height, bbox.x1/width, bbox.y1/height]})
        assert item.get_fontsize() >= 9, item.get_text()
        assert bbox.x0 >= -.5 and bbox.y0 >= -.5 and bbox.x1 <= width+.5 and bbox.y1 <= height+.5, item.get_text()
        assert not any('\u4e00' <= char <= '\u9fff' for char in item.get_text())
    for ax in axes[:2]:
        assert ax.get_xlim() == tuple(limits) and ax.get_ylim() == tuple(limits)
    for system in SYSTEMS:
        values = scatter[scatter.system == system][['true_delta_pKD', 'predicted_delta_pKD']].to_numpy()
        assert values.min() > limits[0] and values.max() < limits[1]
    stem = HERE / 'fixed_backbone_antigen_pilot'
    fig.savefig(str(stem)+'.pdf', dpi=450, bbox_inches=None)
    fig.savefig(str(stem)+'.svg', dpi=450, bbox_inches=None)
    fig.savefig(str(stem)+'.png', dpi=450, bbox_inches=None)
    with Image.open(str(stem)+'.png') as preview:
        png_size = list(preview.size)
    primary_gains = sensitivity[(sensitivity.label_variant == MAIN_LABEL) &
        sensitivity.template_pdb.isin(['2BNQ', '1AO7'])]
    result = {'status': 'data_and_export_checks_passed_pending_personal_visual_QA',
        'dimensions_mm': [183, 196], 'raster_dpi': 450, 'PNG_pixels': png_size,
        'minimum_visible_font_pt': min(row['font_pt'] for row in visible_text),
        'all_text_inside_canvas': True, 'all_figure_labels_English': True,
        'archetype': 'quantitative_grid',
        'source_files': archived, 'source_row_counts': {'candidate_predictions': len(predictions),
            'candidate_level_metrics': len(metrics), 'dynamic_increment_comparisons': len(increment)},
        'plotted_scope': {'panels_a_b': 'All 26 eligible primary single substitutions, two held-position-out predictions per candidate, total52 markers.',
            'panel_c': 'All7 prespecified methods in both primary systems,14 MAE values.',
            'panel_d': 'All3 label fits x5 predeclared templates;15 MAE increments.'},
        'plot_filters': {'analysis': ANALYSIS, 'primary_label': MAIN_LABEL,
            'main_templates': ['2BNQ', '1AO7'], 'scatter_methods': [BASE, JOINT]},
        'all_requested_candidate_points_methods_and_sensitivity_cells_retained': True,
        'fitted_training_predictions_used': False,
        'heldout_candidate_not_in_training_verified_rows': len(leakage),
        'candidate_MAE_recomputed_matches_metric_source': True,
        'heatmap_increment_recomputed_matches_source': True,
        'identical_axis_limits_a_b': limits, 'scatter_reference': 'identity y=x; not fitted regression or a mean trend',
        'heatmap_zero_symmetric_limits': [-vmax, vmax],
        'heatmap_sign': 'positive baseline-MAE minus joint-MAE improves; negative worsens',
        'heatmap_positive_cells': int((matrix > 0).sum()), 'heatmap_negative_cells': int((matrix < 0).sum()),
        'primary_label_primary_template_gains': primary_gains[['system', 'template_pdb', 'baseline_MAE', 'joint_MAE', 'joint_increment_vs_baseline']].to_dict(orient='records'),
        'zero_change_ranking_displayed': False,
        'statistical_tests': None, 'confidence_intervals': None,
        'all_original_source_hashes_unchanged_after_plotting': all(sha(originals[name]) == row['SHA256'] for name, row in archived.items()),
        'visible_text': visible_text,
        'export_SHA256': {ext: sha(str(stem)+'.'+ext) for ext in ['png', 'pdf', 'svg']}}
    write_json(HERE / 'data_audit.json', result)
    write_json(HERE / 'figure_contract.json', {
        'question': 'Do RMSF and Hchi reliably add held-position-out antigen-affinity ranking information beyond sequence/static features in this fixed-backbone pilot?',
        'conclusion': 'The joint readouts provide no stable gain: primary MAE increases in both systems and only2 of15 label/template combinations improve.',
        'archetype': 'quantitative_grid', 'dimensions_mm': [183,196], 'font_min_pt': 9,
        'panels': {'a': 'All candidate-level paired predictions in1G4.', 'b': 'All candidate-level paired predictions inA6.',
            'c': 'Absolute MAE for every prespecified method; one-readout ablations are retained.',
            'd': 'Sign and size of joint increment under all predeclared fits/templates.'},
        'scope': 'Independent figure, no manuscript insertion, historical2431 not manuscript2430 model.',
        'inherited_assets': {'PairedBoxScatter1.R': 'neutral within-candidate connectors,individual points,clean spines; no boxes or t-test because not applicable',
            'BarComparison/plot_comparison_GeneRegulatory.py': 'complete method ordering,clean spines; dots replace bars to compare two MAE values without error bars',
            'heatmap/plot_comparison.py': 'cell annotations and luminance-adaptive text; a single zero-symmetric scale replaces per-column normalization'}})
    write_json(HERE / 'statistics_and_reproducibility.json', {
        'independent_observation': 'A biological peptide candidate, not a residue, template, conformer or fit variant.',
        'primary_n': {'1G4':13, 'A6':13, 'total':26}, 'independent_TCR_HLA_contexts':2,
        'HLA': 'Both HLA-A*02:01; related systems, not a broad-generalization sample.',
        'base_model': 'Historical strict_full2431 ATLAS/DynaRepo no-ITrace model,seed7,epoch40; no fine-tuning. A different model from current manuscript2430.',
        'model_sha256': 'b3eed060919df15f075329bc892f959d9355232be5feb561417edd49ea82b39f',
        'split': 'Within each system, leave one single-mutation position out; six position folds. All substitutions at the held position are tested together.',
        'calibration': 'Prespecified ridge,lambda1,no intercept; source training-fold RMS scaling only. No hyperparameter selection or refitting for figure.',
        'label': 'Reported-hybrid SPR pKD candidate-minus-parent; all-constrained and all-Bmax-fitted are sensitivity fits,not independent samples.',
        'center_statistic': 'MAE=mean absolute candidate-level heldout delta_pKD error across13 candidates per system.',
        'spread_or_interval': 'None plotted; finite retrospective panels without confidence intervals.',
        'statistical_test': 'None', 'multiple_comparison_adjustment': 'Not applicable',
        'baseline': 'Sequence plus static BSA/contact features; R/H are unchanged fixed-environment mean predicted RMSF/Hchi changes.',
        'zero_change': 'Shown only as its MAE; no Spearman or lexical tie-resolved top3 ranking is plotted.',
        'scatter_source_columns': ['true_delta_pKD', 'predicted_delta_pKD', 'candidate_id', 'system', 'template_pdb', 'model'],
        'method_source_column': 'candidate_level_metrics.MAE_delta_pKD',
        'heatmap_source_column': 'dynamic_increment_comparisons.joint_increment_vs_baseline',
        'complete_sources': archived,
        'unplotted_source_scope': 'All original2037 prediction rows,210 metric rows and30 increment rows are archived. This figure addresses the predeclared primary single-mutation within-system question; secondary cross-system transport and double-substitution evaluations remain in the full source files,not silently discarded.',
        'limits': 'Retrospective known panels with fixed backbone and deterministic sidechain completion; no prospective design,activation,affinity-free-energy or causal Hchi-to-RMSF claim.'})
    print(json.dumps({'status': 'rendered', 'figure': str(stem)+'.png',
        'sources': len(predictions), 'scatter_candidates': 26, 'sensitivity_cells':15,
        'heatmap_positive': int((matrix>0).sum())}), flush=True)


if __name__ == '__main__':
    main()
