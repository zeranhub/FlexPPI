# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a) Method comparison dot plot -> BarComparison/plot_comparison_Trajectory.py -> param inherit
# (b) Method comparison dot plot -> BarComparison/plot_comparison_Trajectory.py -> param inherit
# (c) Time-block heatmap -> heatmap/plot_comparison.py -> param inherit
# RULE: "native run" = load pre-rendered PNG via Image.open().ax.imshow().
#       "param inherit" = drawing function below that copies Class A/B/C values.
#       If a panel says "native run" and you write a drawing function, you broke the contract.
# The assets were read and their companion PNGs inspected. Their multi-dataset
# bar grid and normalized model-benchmark heatmap have different data structures
# from one matched-residue cohort and a physical time-block matrix. Cross-type
# parameter inheritance preserves direct labels, minimal axes, method ordering,
# and signed cells; journal typography/line widths supersede screen-size values.

import os
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'
os.environ['OPENBLAS_NUM_THREADS'] = '1'

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

mpl.use('Agg')
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.text import Text
from matplotlib import font_manager

OUT = Path(__file__).resolve().parent
SOURCE = OUT.parent/'all_qualifying_pairs_with_physical_verification.csv'
METHODS = ['full', 'intra', 'aa_sasa_contacts', 'static47_sasa']
LABELS = ['Full', 'Intra', 'AA + SASA + contacts', 'Static47 + SASA']
COLORS = [CATEGORICAL[0], CATEGORICAL[3], GREY, CATEGORICAL[5]]
MARKERS = ['o', '^', 's', 'D']


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def data_section():
    source = pd.read_csv(SOURCE)
    p = source[(source.sample_id == 'A00C9') & (source.chain == 'A')].copy()
    assert len(p) == 211
    assert len(set(p.low_residue_index) | set(p.high_residue_index)) == 122
    assert p.true_gap.gt(1e-7).all()
    assert p.sasa_abs_gap_A2.le(10).all() and p.rmsf_abs_gap_A.le(.1).all()
    assert p.minimum_end_distance.ge(3).all()
    summaries = []
    for method, label in zip(METHODS, LABELS):
        delta = p[f'{method}_gap'].to_numpy()
        accuracy = np.where(delta > 1e-7, 1., np.where(delta < -1e-7, 0., .5)).mean()*100
        summaries.append(dict(method=method, label=label, n_pairs=len(p),
                              correct_order_percent=float(accuracy),
                              mean_absolute_gap_error=float(p[f'{method}_gap_abs_error'].mean())))
    metrics = pd.DataFrame(summaries)
    high = p[p.true_gap >= .4].sort_values(
        ['amino_acid', 'low_residue_index', 'high_residue_index'], kind='stable').copy()
    assert len(high) == 46
    matrix = high[[f'hchi_block_{i}_gap' for i in range(1, 6)]].to_numpy()
    assert np.isfinite(matrix).all()
    assert ((matrix > 1e-7).sum(axis=1) >= 4).all()
    assert int(((matrix > 1e-7).sum(axis=1) == 5).sum()) == 42
    assert matrix.min() >= -2 and matrix.max() <= 2
    high.insert(0, 'display_row', np.arange(1, len(high)+1))
    p.to_csv(OUT/'source_data_all_211_pairs.csv', index=False)
    high.to_csv(OUT/'source_data_heatmap_all_46_pairs.csv', index=False)
    metrics.to_csv(OUT/'source_data_method_metrics.csv', index=False)
    return metrics, high, matrix


def panel_dot(ax, metrics, metric, limits, ticks, xlabel, fmt):
    for y, (value, color, marker) in enumerate(zip(metrics[metric], COLORS, MARKERS)):
        ax.scatter(value, y, s=30, color=color, marker=marker,
                   edgecolors='white', linewidths=.35, zorder=3)
        ax.annotate(fmt.format(value), (value, y), xytext=(5, 0),
                    textcoords='offset points', ha='left', va='center',
                    fontsize=8, color=BLACK)
    ax.set_yticks(range(4), LABELS)
    ax.tick_params(axis='y', length=0, pad=7, labelsize=8)
    ax.tick_params(axis='x', length=3, labelsize=8)
    ax.get_yticklabels()[0].set_fontweight('bold')
    ax.set_ylim(3.5, -.5)
    ax.set_xlim(*limits)
    ax.set_xticks(ticks)
    ax.set_xlabel(xlabel, labelpad=7, fontsize=8)
    ax.spines['left'].set_visible(False)


def compose_figure(metrics, high, matrix):
    # asymmetric_mixed: the complete 46x5 matrix needs more vertical space
    # than the two aggregate-method comparisons. Figure size is final size.
    fig = plt.figure(figsize=(7.2, 8.9), facecolor='white')
    fig.text(.065, .967, 'DraG: a whole-chain comparison of matched sidechains',
             fontsize=9, fontweight='bold', color=BLACK)
    fig.text(.065, .943, '3O5T, chain A  |  211 matched pairs  |  122 distinct residues',
             fontsize=8, color=BLACK)
    axa = fig.add_axes([.265, .648, .245, .205])
    axb = fig.add_axes([.265, .358, .245, .205])
    axc = fig.add_axes([.685, .170, .255, .710])

    fig.text(.065, .886, 'a', fontsize=9, fontweight='bold')
    fig.text(.105, .886, 'Correct ordering', fontsize=8, fontweight='bold')
    panel_dot(axa, metrics, 'correct_order_percent', (46, 70), [50, 60, 70],
              'Correct order (%)  ↑', '{:.1f}')
    axa.axvline(50, ls=(0, (3, 3)), lw=.5, color=GREY, zorder=0)

    fig.text(.065, .594, 'b', fontsize=9, fontweight='bold')
    fig.text(.105, .594, 'Accuracy of contrast magnitude', fontsize=8, fontweight='bold')
    panel_dot(axb, metrics, 'mean_absolute_gap_error', (.235, .286), [.24, .26, .28],
              'Mean |predicted gap − MD gap| (nats)  ↓', '{:.3f}')
    axb.set_xticklabels(['0.24', '0.26', '0.28'])

    fig.text(.580, .916, 'c', fontsize=9, fontweight='bold')
    fig.text(.620, .916, 'MD contrast across time', fontsize=8, fontweight='bold')
    fig.text(.620, .897, 'All 46 pairs with mean gap ≥0.4', fontsize=8)
    # All cells use one physical color scale; no per-row normalization.
    signed = LinearSegmentedColormap.from_list('orange_white_blue',
                                               [CATEGORICAL[3], '#FFFFFF', CATEGORICAL[0]])
    image = axc.imshow(matrix, aspect='auto', interpolation='nearest',
                       cmap=signed, norm=TwoSlopeNorm(vmin=-2, vcenter=0, vmax=2))
    labels = [f'{lo} / {hi}' for lo, hi in zip(high.low_residue, high.high_residue)]
    axc.set_yticks(np.arange(len(labels)), labels, fontsize=8)
    axc.tick_params(axis='y', length=0, pad=4)
    axc.set_xticks(np.arange(5), ['0–100', '100–200', '200–300', '300–400', '400–500'],
                  rotation=45, ha='right', rotation_mode='anchor', fontsize=8)
    axc.tick_params(axis='x', length=0, pad=4)
    axc.set_xlabel('Consecutive time blocks (ns)', fontsize=8, labelpad=5)
    for spine in axc.spines.values():
        spine.set_visible(False)
    aa = high.amino_acid.to_numpy()
    for boundary in np.flatnonzero(aa[1:] != aa[:-1]):
        axc.axhline(boundary+.5, color='white', lw=.65)
    # Negative cells are additionally marked so sign is visible in greyscale.
    rev = np.argwhere(matrix < -1e-7)
    axc.scatter(rev[:, 1], rev[:, 0], s=6, c=BLACK, marker='o', linewidths=0)
    cbax = fig.add_axes([.685, .101, .255, .012])
    colorbar = fig.colorbar(image, cax=cbax, orientation='horizontal', ticks=[-2, 0, 2])
    colorbar.outline.set_visible(False)
    colorbar.ax.tick_params(labelsize=8, length=2, width=.6)
    colorbar.set_label('MD Hχ gap (nats)', fontsize=8, labelpad=2)

    fig.text(.065, .284, 'Matched within one chain and amino acid', fontsize=8, fontweight='bold')
    fig.text(.065, .260, '|ΔSASA| ≤10 Å²; |ΔRMSF| ≤0.10 Å', fontsize=8)
    fig.text(.065, .238, 'Four or more sidechains of each type;', fontsize=8)
    fig.text(.065, .217, 'first and last three residues excluded.', fontsize=8)
    fig.text(.065, .178, 'Time-block consistency', fontsize=8, fontweight='bold')
    fig.text(.065, .154, '42/46 pairs keep the order in all five blocks;', fontsize=8)
    fig.text(.065, .133, '46/46 keep the order in at least four blocks.', fontsize=8)
    fig.text(.065, .088, 'Rows: low / high mean Hχ; fixed AA and residue order.', fontsize=8)
    fig.text(.065, .067, 'Orange cells and dots mark a reversed MD order.', fontsize=8)
    fig.text(.065, .036, 'Pairs share residues. Time blocks are not independent repeats.', fontsize=8,
             color=CATEGORICAL[5])
    return fig


def main():
    metrics, high, matrix = data_section()
    fig = compose_figure(metrics, high, matrix)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    limits = fig.bbox
    outside = []
    for txt in fig.findobj(Text):
        if not txt.get_visible() or not txt.get_text().strip():
            continue
        box=txt.get_window_extent(renderer)
        if box.x0 < limits.x0-1 or box.y0 < limits.y0-1 or box.x1 > limits.x1+1 or box.y1 > limits.y1+1:
            outside.append(txt.get_text())
    assert not outside, outside
    stem=OUT/'drag_whole_chain_matched_sidechain_comparison'
    # bbox_inches=None alone falls back to the global "tight" baseline;
    # override that rcParam in a context to preserve exact physical dimensions.
    with mpl.rc_context({'savefig.bbox': None}):
        for suffix in ['pdf', 'svg', 'png']:
            fig.savefig(stem.with_suffix('.'+suffix), dpi=300, bbox_inches=None, pad_inches=0,
                        facecolor='white', metadata={'Creator':'FlexPPI numerical case audit'})
    from PIL import Image
    with Image.open(stem.with_suffix('.png')) as preview:
        actual_preview_pixels=list(preview.size)
    assert actual_preview_pixels == [2160,2670], actual_preview_pixels
    plt.close(fig)
    audit=dict(status='passed', sample='A00C9', pdb='3O5T', chain='A', archetype='asymmetric_mixed',
               figure_size_inches=[7.2,8.9], preview_pixels=actual_preview_pixels, font_size_pt=8,
               resolved_font=font_manager.findfont(mpl.font_manager.FontProperties(family='Arial')),
               n_pairs=211,n_unique_residues=122,n_large_gap_pairs=46,n_time_blocks=5,
               negative_cells=int((matrix < -1e-7).sum()), all5_positive_pairs=int(((matrix>1e-7).sum(1)==5).sum()),
               minimum_positive_blocks_per_pair=int((matrix>1e-7).sum(1).min()),
               heatmap_range=[float(matrix.min()),float(matrix.max())],color_scale=[-2,2],
               heatmap_unique_values=len(np.unique(matrix)),heatmap_mean_row_sd=float(np.std(matrix,axis=1).mean()),
               generic_heatmap_sd_threshold='Not applied: this is a time-consistency audit, not expression clustering; low time variation is meaningful and was not changed.',
               text_outside_canvas=outside, data_downsampling=False, model_performance_selection=False,
               intervals='None: shared-residue pairs and time blocks are not independent replicates.',
               source_sha256=sha(SOURCE),script_sha256=sha(Path(__file__)),
               metrics=metrics.to_dict('records'))
    (OUT/'figure_data_and_layout_audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
    print(json.dumps(audit,indent=2))


if __name__ == '__main__':
    main()
