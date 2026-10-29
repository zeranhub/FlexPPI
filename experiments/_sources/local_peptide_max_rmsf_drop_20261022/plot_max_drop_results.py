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

from pathlib import Path
from datetime import datetime, timezone
import csv
import hashlib
import json
import shutil
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.font_manager import findfont, FontProperties

ROOT = Path(__file__).resolve().parent
DATA = ROOT / 'functional_screening'
OUT = ROOT / 'result_figures'
OUT.mkdir(exist_ok=True)
PANELS = [f'MSKCC_TCR{i}' for i in range(1, 8)]
LABELS = ['TCR1 · CMV*', 'TCR2 · CMV', 'TCR3 · CMV', 'TCR4 · gp100',
          'TCR5 · gp100', 'TCR6 · gp100', 'TCR7 · GRL']
METHODS = ['max_rmsf_drop_A', 'raw_BLOSUM62', 'HERMES', 'ProteinMPNN020', 'ProteinMPNN002']
METHOD_LABELS = ['Max drop', 'BLOSUM62', 'HERMES', 'MPNN\n020', 'MPNN\n002']
BLUE, ORANGE = CATEGORICAL[0], CATEGORICAL[3]

def read(path):
    with path.open(encoding='utf-8-sig', newline='') as handle:
        return list(csv.DictReader(handle))

def write_csv(path, rows):
    with path.open('w', encoding='utf-8-sig', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def export(fig, stem):
    fig.canvas.draw()
    save_cns_figure(fig, OUT / stem)
    fig.savefig(OUT / (stem + '.svg'), bbox_inches='tight')
    plt.close(fig)

def overview(rows):
    lookup = {(r['panel_id'], r['score']): r for r in rows}
    assert len(rows) == 42
    assert set(lookup) == {(p, m) for p in PANELS for m in METHODS + ['second_largest_drop_A']}
    source_rows = []
    gains = np.zeros((7, 5))
    for i, panel in enumerate(PANELS):
        for j, method in enumerate(METHODS):
            r = lookup[panel, method]
            assert int(r['budget']) == 30 and int(r['n_candidates']) == 171
            gains[i, j] = float(r['expected_hits']) - float(r['random_expected_hits'])
            source_rows.append(dict(panel_id=panel, score=method, budget=30,
                                    expected_hits=float(r['expected_hits']),
                                    pool_random_expected_hits=float(r['random_expected_hits']),
                                    pool_random_hits_p10=float(r['random_hits_p10']),
                                    pool_random_hits_p90=float(r['random_hits_p90']),
                                    position_matched_random_expected_hits=float(r['mutation_position_matched_random_hits']),
                                    additional_hits_vs_pool_random=float(gains[i, j])))
    write_csv(OUT / 'max_drop_overview_source.csv', source_rows)
    fig = plt.figure(figsize=(183 / 25.4, 108 / 25.4))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.03, 1.12], left=.14, right=.975,
                          bottom=.23, top=.88, wspace=.27)
    ax, heat = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])
    fig.text(.025, .96, 'a', weight='bold', fontsize=10)
    fig.text(.54, .96, 'b', weight='bold', fontsize=10)
    ax.set_title('Max-drop rule: yield from 30 candidates', loc='left', pad=9)
    for i, panel in enumerate(PANELS):
        r = lookup[panel, METHODS[0]]
        rand, match, observed = map(float, [r['random_expected_hits'],
                                          r['mutation_position_matched_random_hits'], r['expected_hits']])
        lo, hi = float(r['random_hits_p10']), float(r['random_hits_p90'])
        ax.plot([lo, hi], [i, i], color='#CECECE', linewidth=3, solid_capstyle='round', zorder=1)
        ax.plot(rand, i, marker='|', markersize=8, color='#666666', markeredgewidth=1, zorder=2)
        ax.plot(match, i + .13, marker='D', color=ORANGE, markersize=4.3,
                markerfacecolor='white', markeredgewidth=1, zorder=3)
        ax.plot(observed, i - .13, marker='o', color=BLUE, markersize=5, zorder=4)
    ax.set_xlim(0, 31)
    ax.set_xticks([0, 10, 20, 30])
    ax.set_ylim(6.55, -.55)
    ax.set_yticks(range(7), LABELS)
    ax.set_xlabel('Functional hits among selected candidates')
    ax.tick_params(axis='y', length=0, pad=5)
    handles = [Line2D([], [], color=BLUE, marker='o', linestyle='none', markersize=4.7, label='Max-drop rule'),
               Line2D([], [], color=ORANGE, marker='D', markerfacecolor='white', linestyle='none', markersize=4.3,
                      label='Position-matched random'),
               Line2D([], [], color='#666666', marker='|', linestyle='-', linewidth=2,
                      markeredgewidth=1, markersize=7, label='Pool random: mean / 10–90%')]
    ax.legend(handles=handles, loc='upper left', bbox_to_anchor=(-.02, -.17), fontsize=7,
              handletextpad=.5, labelspacing=.45, borderaxespad=0)
    heat.set_title('Zero-label comparators: gain over pool random', loc='left', pad=9)
    limit = float(np.ceil(np.abs(gains).max()))
    cmap = LinearSegmentedColormap.from_list('additional_hits', DIVERGING[::-1])
    im = heat.imshow(gains, cmap=cmap, vmin=-limit, vmax=limit, aspect='auto', interpolation='nearest')
    heat.add_patch(Rectangle((-.5, -.5), 1, 7, fill=False, edgecolor=BLUE, linewidth=1))
    heat.set_yticks(range(7), [])
    heat.set_xticks(range(5), METHOD_LABELS)
    heat.tick_params(length=0, pad=5, labelsize=7)
    heat.spines[['left', 'bottom']].set_visible(False)
    for i in range(7):
        for j in range(5):
            value = gains[i, j]
            heat.text(j, i, f'{value:+.1f}', ha='center', va='center', fontsize=7,
                      color='white' if abs(value) > .62 * limit else BLACK)
    cax = heat.inset_axes([.10, -.26, .80, .055])
    cb = fig.colorbar(im, cax=cax, orientation='horizontal')
    cb.set_ticks([-limit, 0, limit])
    cb.set_label('Additional functional hits per 30 selections', fontsize=7, labelpad=3)
    cb.outline.set_visible(False)
    cb.ax.tick_params(length=2, labelsize=6.5, pad=2)
    export(fig, 'max_drop_overview')

def scatter(candidates, ranking, metrics):
    assert len(candidates) == len(ranking) == 1197
    ranked = {r['candidate_id']: r for r in ranking}
    lookup = {r['panel_id']: r for r in metrics if r['score'] == METHODS[0]}
    source = []
    all_y = np.asarray([float(r['functional_response']) for r in candidates])
    assert np.isfinite(all_y).all()
    ylo, yhi = min(0., all_y.min()), max(1., all_y.max())
    yrange = yhi - ylo
    ylim = (ylo - .045 * yrange, yhi + .23 * yrange)
    fig, axes = plt.subplots(2, 4, figsize=(183 / 25.4, 110 / 25.4), sharey=True)
    fig.subplots_adjust(left=.085, right=.985, bottom=.16, top=.87, hspace=.41, wspace=.24)
    for i, (panel, ax) in enumerate(zip(PANELS, axes.flat)):
        rr = [r for r in candidates if r['panel_id'] == panel]
        assert len(rr) == 171
        x = np.asarray([float(r['max_rmsf_drop_A']) for r in rr])
        y = np.asarray([float(r['functional_response']) for r in rr])
        selected = np.asarray([float(ranked[r['candidate_id']]['top30_expected_membership']) for r in rr])
        assert abs(selected.sum() - 30) < 1e-9
        assert set(selected).issubset({0., 1.}), 'Boundary ties need fractional-membership visual encoding.'
        assert np.isfinite(x).all()
        xmin, xmax = min(0., x.min()), max(0., x.max())
        xspan = xmax - xmin
        ax.axvline(0, color='#999999', linewidth=.55, linestyle=(0, (3, 2)), zorder=0)
        ax.axhline(.3, color='#999999', linewidth=.55, linestyle=(0, (3, 2)), zorder=0)
        ax.scatter(x[selected == 0], y[selected == 0], s=8, facecolors='#BBBBBB',
                   edgecolors='white', linewidths=.2, alpha=.72, zorder=2)
        ax.scatter(x[selected == 1], y[selected == 1], s=12, facecolors=BLUE,
                   edgecolors='white', linewidths=.25, alpha=.9, zorder=3)
        ax.set_xlim(xmin - .06 * xspan, xmax + .06 * xspan)
        ax.set_ylim(ylim)
        ax.set_title(LABELS[i], loc='left', pad=6)
        ax.text(-.18, 1.10, 'abcdefg'[i], transform=ax.transAxes, weight='bold', fontsize=9)
        metric = lookup[panel]
        ax.text(.04, .96, f"Spearman ρ = {float(metric['spearman_response']):+.2f}\nAUROC = {float(metric['AUROC']):.2f}",
                transform=ax.transAxes, ha='left', va='top', fontsize=6.5, linespacing=1.4)
        ax.tick_params(length=2.5, pad=2.5)
        ax.xaxis.set_major_locator(mpl.ticker.MaxNLocator(3, min_n_ticks=3))
        ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(.5))
        for r, membership in zip(rr, selected):
            source.append(dict(panel_id=panel, candidate_id=r['candidate_id'],
                               mutation=r['mutation'], max_rmsf_drop_A=float(r['max_rmsf_drop_A']),
                               functional_response=float(r['functional_response']), is_hit=int(r['is_hit']),
                               top30_expected_membership=float(membership),
                               spearman_response=float(metric['spearman_response']), AUROC=float(metric['AUROC'])))
    legend = axes.flat[-1]
    legend.axis('off')
    legend.set_title('How to read this test', loc='left', pad=6)
    handles = [Line2D([], [], marker='o', linestyle='none', markersize=4,
                      color=BLUE, label='Selected top 30'),
               Line2D([], [], marker='o', linestyle='none', markersize=3.5,
                      color='#BBBBBB', label='Other 141 candidates'),
               Line2D([], [], color='#999999', linewidth=.7, linestyle=(0, (3, 2)),
                      label='Zero drop / hit threshold')]
    legend.legend(handles=handles, loc='upper left', bbox_to_anchor=(-.04, .98),
                  fontsize=6.8, handletextpad=.5, labelspacing=.7)
    legend.text(0., .46, 'Drop = parent − mutant RMSF\nMaximum over all 9 peptide sites\nHigher score → earlier selection\nHit: functional response ≥ 0.3\nAll 171 mutants shown per context',
                transform=legend.transAxes, fontsize=6.5, va='top', linespacing=1.6)
    fig.text(.52, .063, 'Maximum predicted peptide RMSF decrease (Å)', ha='center', fontsize=8)
    fig.text(.012, .52, 'Highest-dose normalized 4-1BB response', ha='center', va='center', rotation=90, fontsize=8)
    fig.text(.085, .967, 'Single-mutant screening · No functional-label calibration', fontsize=9)
    fig.text(.085, .015, '* TCR1 construct sequence unresolved. GRL: GRLKALCQR neoantigen.', fontsize=6.7)
    write_csv(OUT / 'max_drop_all_contexts_source.csv', source)
    export(fig, 'max_drop_all_contexts')

def main():
    files = ['primary_top30_metrics.csv', 'all_candidate_scores.csv', 'max_drop_rankings_all_candidates.csv']
    metrics, candidates, ranking = [read(DATA / name) for name in files]
    overview(metrics)
    scatter(candidates, ranking, metrics)
    # Keep exact input tables beside the figures for direct reproduction.
    for name in files:
        shutil.copy2(DATA / name, OUT / ('input_' + name))
    captions = {
        'max_drop_overview': 'Panel a: Expected functional hits among 30 candidates selected without functional-label calibration. Pool-random bar spans the exact hypergeometric 10th to 90th percentiles; grey tick is its mean. Position-matched random expectation holds fixed the max-drop rule\'s selected mutation-position counts, randomizing amino-acid choices within each position. Panel b: Additional hits versus pool-random expectation for five predeclared zero-label rankings; all methods select 30. Exact score ties have exchangeable expected membership, so comparator counts can be fractional. Each context includes all 171 nonparent single mutants. A functional hit is a highest-dose normalized CD8 4-1BB response ≥ 0.3; this endpoint is not binding affinity. TCR1*: construct-sequence discrepancy remains unresolved. CMV: NLVPMVATV; gp100: IMDQVPFSV; GRL: GRLKALCQR.',
        'max_drop_all_contexts': 'Each point is one of all 171 nonparent single mutants in a fixed TCR context. The score is the largest signed predicted RMSF decrease, parent minus mutant, among the nine peptide residues on the same parent backbone. Blue highlights the predeclared top 30; grey shows the remaining 141. Dashed lines mark zero decrease and functional-hit threshold 0.3. Spearman correlation and AUROC use all 171 candidates; no response-dependent score direction or site exclusion was chosen. Every measured response is shown without clipping. This is an exploratory retrospective functional-response test, not direct affinity validation. TCR1*: construct-sequence discrepancy remains unresolved.'
    }
    (OUT / 'figure_captions.txt').write_text('\n\n'.join(f'{k}\n{v}' for k, v in captions.items()), encoding='utf-8')
    manifest = {
        'created_utc': datetime.now(timezone.utc).isoformat(),
        'script': str(Path(__file__).resolve()),
        'script_sha256': sha(Path(__file__)),
        'question': 'Does the largest predicted decrease at any peptide RMSF site enrich functional candidates?',
        'source': [{'path': str(DATA / name), 'sha256': sha(DATA / name)} for name in files],
        'figures': [{'stem': stem, 'caption': caption,
                     'files': [{'path': str(OUT / (stem + ext)), 'sha256': sha(OUT / (stem + ext))}
                               for ext in ['.pdf', '.svg', '.png']],
                     'source_csv': str(OUT / (stem + '_source.csv'))} for stem, caption in captions.items()],
        'font_path': findfont(FontProperties(family='Arial')),
        'style': {'width_mm': 183, 'raster_dpi': 300, 'pdf_fonttype': 42, 'svg_fonttype': 'none'},
        'data_integrity': {'candidate_rows': len(candidates), 'contexts': 7, 'candidates_per_context': 171,
                           'scatter_rows_shown': 1197, 'scatter_rows_dropped': 0, 'functional_responses_clipped': False,
                           'overview_primary_methods': METHODS, 'diagnostic_not_promoted': 'second_largest_drop_A'},
        'visual_qa': 'Pending final PNG inspection.'
    }
    (OUT / 'figure_manifest.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding='utf-8')

if __name__ == '__main__':
    main()
