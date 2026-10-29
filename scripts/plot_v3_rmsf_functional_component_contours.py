# -*- coding: utf-8 -*-
# Academic Figure Skill Typography Baseline -- COPY VERBATIM, place at TOP of script
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
# Academic Figure Skill Export Baseline -- COPY VERBATIM
mpl.rcParams.update({
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
    "savefig.bbox": "tight",
    "savefig.dpi": 300,
})

def save_cns_figure(fig, filename):
    """Standard Academic Figure Skill export: vector PDF + 300dpi PNG preview."""
    fig.savefig(f"{filename}.pdf", bbox_inches="tight", dpi=300)
    fig.savefig(f"{filename}.png", bbox_inches="tight", dpi=300)

# Copy-first adaptation of plot_v3_rmsf_chain_count_contours.py.
# Academic asset KernelDensity/plot_CorrelationDensity1.R and its preview were
# inspected in the preceding figure workflow. Reuse that verified 2-D density
# design and the exact histogram/Gaussian estimator; no R dependency is needed.
# This is a PARTIAL-ANNOTATION PREVIEW, not a complete functional classification.
import json
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib import patheffects
from matplotlib.lines import Line2D
from matplotlib.path import Path as PlotPath
from scipy.ndimage import gaussian_filter

import plot_v3_overall_capability as base
import plot_v3_overall_capability_contours as previous

ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = ROOT / 'artifacts/evaluation_v3/analysis'
PRIOR = ROOT / 'output/flexppi_v3_overall_capability'
OUT = ROOT / 'output/flexppi_v3_overall_capability_functional_component_contours'
MANIFEST = ROOT / 'artifacts/evaluation_v3/oof/prediction_manifest.jsonl'
ROLES = ROOT / 'output/itrace_benchmark/raw/itrace_chain_roles.json'
GROUPS = ['1', '3plus', 'unassigned']
LABELS = {'1': '1 component*', '3plus': '\u22653 components\u2020', 'unassigned': 'Unassigned'}
COLORS = {'1': '#B2182B', '3plus': '#762A83', 'unassigned': '#777777'}
STYLES = {'1': '--', '3plus': '-.', 'unassigned': ':'}
EXPECTED = {'1': (2193, 496093), '3plus': (244, 199668), 'unassigned': (631, 263586)}


def load_groups():
    chains = pd.read_csv(ANALYSIS / 'per_chain_metrics.csv', dtype={'sample_id': str, 'chain': str})
    chains = chains.loc[chains.metric.eq('rmsf')].copy()
    if chains.duplicated(['sample_id', 'chain']).any() or not chains.valid_fraction.eq(1).all():
        raise ValueError('Expected unique chains and full RMSF coverage')
    chain_lookup = chains.groupby('sample_id').chain.apply(list).to_dict()
    members = chains.groupby(['sample_id', 'dataset'], sort=False).agg(
        chain_count=('chain', 'nunique'), n_residues=('n_residues', 'sum')).reset_index()
    records = [json.loads(line)['record'] for line in MANIFEST.read_text(encoding='utf-8').splitlines() if line]
    records = {r['sample_id']: r for r in records}
    roles = json.loads(ROLES.read_text(encoding='utf-8'))
    annotations = []
    for row in members.itertuples(index=False):
        record = records[row.sample_id]
        if record['dataset'] != row.dataset:
            raise ValueError('Dataset mismatch')
        pdb_id = record['metadata']['pdb_id'].lower()
        all_chains = chain_lookup[row.sample_id]
        annotation = dict(sample_id=row.sample_id, pdb_id=pdb_id,
                          functional_component_count=None, functional_group='unassigned',
                          annotation_status='unassigned_multichain',
                          rule='No validated functional chain-to-component mapping',
                          component_chain_groups_json='', chain_roles_json='', evidence='')
        if row.dataset == 'itrace':
            chain_roles = roles[pdb_id]
            if set(chain_roles) != set(all_chains) or sorted(chain_roles.values()) != sorted(
                    ['mhc', 'b2m', 'peptide', 'tcr_alpha', 'tcr_beta']):
                raise ValueError('Incomplete or ambiguous ITrace chain roles')
            components = {
                'TCR': sorted(c for c, role in chain_roles.items() if role in ['tcr_alpha', 'tcr_beta']),
                'MHC': sorted(c for c, role in chain_roles.items() if role in ['mhc', 'b2m']),
                'peptide': sorted(c for c, role in chain_roles.items() if role == 'peptide'),
            }
            annotation.update(functional_component_count=3, functional_group='3plus',
                annotation_status='role_supported_user_defined_grouping',
                rule='TCR alpha+beta = 1; MHC heavy chain+b2m = 1; peptide = 1',
                component_chain_groups_json=json.dumps(components, sort_keys=True),
                chain_roles_json=json.dumps(chain_roles, sort_keys=True), evidence=str(ROLES))
        elif row.chain_count == 1:
            annotation.update(functional_component_count=1, functional_group='1',
                annotation_status='operational_single_chain_not_functionally_audited',
                rule='One input polypeptide counted as one component; fusion proteins not audited',
                component_chain_groups_json=json.dumps({'input_polypeptide': all_chains}),
                evidence=str(ANALYSIS / 'per_chain_metrics.csv'))
        annotations.append(annotation)
    members = members.merge(pd.DataFrame(annotations), on='sample_id', validate='one_to_one')
    members['functional_component_count'] = members.functional_component_count.astype('Int64')
    data = pd.read_csv(ANALYSIS / 'all_residue_predictions.csv.gz',
        usecols=['sample_id', 'dataset', 'valid_rmsf', 'md_rmsf_angstrom', 'pred_rmsf_angstrom'],
        dtype={'sample_id': str})
    data = data.merge(members[['sample_id', 'dataset', 'functional_group']],
                      on=['sample_id', 'dataset'], how='left', validate='many_to_one')
    if len(members) != 3068 or len(records) != 3068 or len(data) != 959347 or data.functional_group.isna().any():
        raise ValueError('Incomplete assignments or residue coverage')
    values = data[['md_rmsf_angstrom', 'pred_rmsf_angstrom']].to_numpy(float)
    if not data.valid_rmsf.all() or not np.isfinite(values).all() or values.min() <= 0:
        raise ValueError('Invalid log-plot observations')
    observed = data.groupby('sample_id').size().sort_index()
    expected = members.set_index('sample_id').n_residues.sort_index()
    if not observed.index.equals(expected.index) or not np.array_equal(observed.to_numpy(), expected.to_numpy()):
        raise ValueError('Per-sample residue counts do not reconcile')
    for group, (n_samples, n_residues) in EXPECTED.items():
        if int(members.functional_group.eq(group).sum()) != n_samples or int(data.functional_group.eq(group).sum()) != n_residues:
            raise ValueError('Unexpected functional-group coverage')
    return data, members


def estimate_regions(data):
    plot_range = (np.log10(.2), np.log10(36))
    edges = np.linspace(plot_range[0] - .25, plot_range[1] + .25, 257)
    centers = (edges[:-1] + edges[1:]) / 2
    sigma = .035 / (edges[1] - edges[0])
    regions = {}
    for group in GROUPS:
        subset = data.loc[data.functional_group.eq(group)]
        points = np.log10(subset[['md_rmsf_angstrom', 'pred_rmsf_angstrom']].to_numpy(float))
        counts, _, _ = np.histogram2d(points[:, 0], points[:, 1], bins=(edges, edges))
        if int(counts.sum()) != len(subset):
            raise ValueError('Lost residues in density estimate')
        probability = gaussian_filter(counts, sigma=sigma, mode='constant')
        probability /= probability.sum()
        ordered = np.sort(probability.ravel())[::-1]
        threshold = ordered[np.searchsorted(np.cumsum(ordered), .95)]
        enclosed = probability >= threshold
        if enclosed[0].any() or enclosed[-1].any() or enclosed[:, 0].any() or enclosed[:, -1].any():
            raise ValueError('Clipped density region')
        fig, ax = base.plt.subplots()
        contour = ax.contour(centers, centers, probability.T, levels=[threshold])
        paths = [p.copy() for p in contour.allsegs[0] if len(p) >= 4]
        base.plt.close(fig)
        inside = np.zeros(len(points), dtype=bool)
        for p in paths:
            if not np.allclose(p[0], p[-1]) or p.min() <= plot_range[0] or p.max() >= plot_range[1]:
                raise ValueError('Open or visibly clipped contour')
            inside ^= PlotPath(p).contains_points(points)
        if not .93 <= inside.mean() <= .98:
            raise ValueError('Unexpected contour coverage')
        regions[group] = dict(paths=paths, counts=counts, probability=probability,
            n_samples=int(subset.sample_id.nunique()), n_residues=len(subset),
            estimated_probability_enclosed=float(probability[enclosed].sum()),
            observed_residue_fraction_enclosed=float(inside.mean()), observed_residues_enclosed=int(inside.sum()),
            threshold_probability_per_cell=float(threshold), contour_components=len(paths),
            color=COLORS[group], linestyle=STYLES[group],
            pooled_metrics=base.stats(subset.pred_rmsf_angstrom, subset.md_rmsf_angstrom))
    return regions, edges


def add_regions(ax, regions, enlarged=False):
    width = 1.25 if enlarged else 1.05
    handles = {}
    for group in ['unassigned', '1', '3plus']:
        halo = [patheffects.Stroke(linewidth=width + 1.1, foreground='white'), patheffects.Normal()]
        for p in regions[group]['paths']:
            line, = ax.plot(p[:, 0], p[:, 1], color=COLORS[group], lw=width,
                            ls=STYLES[group], zorder=4 if group == 'unassigned' else 5)
            line.set_path_effects(halo)
        label = LABELS[group] if enlarged else {'1': '1*', '3plus': '\u22653\u2020', 'unassigned': 'Unassigned'}[group]
        handles[group] = Line2D([], [], color=COLORS[group], lw=width, ls=STYLES[group],
                               label=label, path_effects=halo)
    # Pending is a text-only legend entry: no data or outline is fabricated.
    pending = Line2D([], [], lw=0, label='2: pending')
    legend = ax.legend(handles=[handles['1'], pending, handles['3plus'], handles['unassigned']],
        title='95% density regions' if enlarged else '95% regions', loc='upper left',
        bbox_to_anchor=(.012, .83), borderaxespad=0, handlelength=2.7,
        fontsize=6.5 if enlarged else 5.2, title_fontsize=6.5 if enlarged else 5.2,
        labelspacing=.32, frameon=True, facecolor='white', edgecolor='none', framealpha=.92)
    legend.get_texts()[1].set_color('#777777')
    legend.set_zorder(10)


def check_legend_clearance(fig, ax, regions):
    fig.canvas.draw()
    box = ax.get_legend().get_window_extent(fig.canvas.get_renderer())
    for region in regions.values():
        for vertices in region['paths']:
            screen_path = PlotPath(ax.transData.transform(vertices))
            if screen_path.intersects_bbox(box, filled=False):
                raise ValueError('Legend obscures a density contour')
    return True


def plot(data, members, regions):
    fig, ax = base.plt.subplots(figsize=(110 * base.MM, 120 * base.MM))
    fig.subplots_adjust(left=.16, right=.89, bottom=.25, top=.86)
    base.heading(ax, 'A', 'Absolute-value agreement', 'Functional components (partial) \u00b7 both axes log-scaled')
    audit = base.density(fig, ax, data.pred_rmsf_angstrom, data.md_rmsf_angstrom, 'RMSF [\u00c5]', True)
    add_regions(ax, regions, enlarged=True)
    audit['legend_does_not_obscure_contours'] = check_legend_clearance(fig, ax, regions)
    fig.text(.52, .105,
        'Samples: 1 component*, 2,193 \u00b7 3 components, 244 \u00b7 Unassigned, 631\n'
        '* Single-chain inputs count as one; fusion proteins have not been audited.\n'
        '\u2020 All annotated \u22653 entries have 3 components: TCR + MHC + peptide.\n'
        'Blue background: all 959,347 held-out residues; grey outline: unassigned.',
        fontsize=5.8, ha='center', va='center', linespacing=1.65)
    panel_paths = base.export(fig, OUT, 'fig_A_rmsf_functional_component_contours_enlarged')
    systems = pd.read_csv(PRIOR / 'all_protein_rmsf_metrics_and_baseline.csv')
    folds = pd.read_csv(PRIOR / 'paired_five_fold_comparison.csv')
    overlay_hook, export_hook = previous.add_regions, previous.export

    def export_full(fig, out, name):
        ax = fig.axes[0]
        ax.texts[0].set_text('Functional components (partial) \u00b7 log-scaled axes')
        ax.text(.5, -.235, '* 1 = single-chain input; \u2020 all \u22653 entries = TCR + MHC + peptide.\n'
                'Unassigned: 631 samples. Two-component annotations pending.',
                transform=ax.transAxes, ha='center', va='top', fontsize=5.1, linespacing=1.4)
        check_legend_clearance(fig, ax, regions)
        return base.export(fig, out, 'fig_v3_overall_rmsf_capability_functional_component_contours')

    try:
        previous.add_regions, previous.export = add_regions, export_full
        full = previous.rmsf_figure(data, systems, folds, regions, OUT)
    finally:
        previous.add_regions, previous.export = overlay_hook, export_hook
    return dict(density=audit, panel_a_paths=panel_paths, full_figure=full)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    protected = [PRIOR / 'fig_v3_overall_rmsf_capability.png',
        ROOT / 'output/flexppi_v3_overall_capability_chain_count_contours/fig_A_rmsf_chain_count_contours_enlarged.png',
        ROOT / 'output/flexppi_v3_overall_capability_dataset_contours/fig_A_rmsf_dataset_contours_enlarged.png']
    before = {str(p): previous.sha256(p) for p in protected}
    data, members = load_groups()
    regions, edges = estimate_regions(data)
    members.to_csv(OUT / 'sample_functional_component_assignments.csv', index=False)
    crossing = members.groupby(['functional_group', 'dataset'], sort=False).agg(
        n_samples=('sample_id', 'size'), n_residues=('n_residues', 'sum')).reset_index()
    crossing.to_csv(OUT / 'functional_groups_by_dataset.csv', index=False)
    audit = plot(data, members, regions)
    summaries = {g: {k: v for k, v in r.items() if k not in ['paths', 'counts', 'probability']} for g, r in regions.items()}
    pd.DataFrame([dict(functional_group=g, **{k: v for k, v in r.items() if k not in ['pooled_metrics', 'n_residues']},
                       **r['pooled_metrics']) for g, r in summaries.items()]).to_csv(OUT / 'functional_group_statistics.csv', index=False)
    coordinates = []
    for group, region in regions.items():
        for component, path in enumerate(region['paths']):
            for vertex, (x, y) in enumerate(path):
                coordinates.append(dict(functional_group=group, component=component, vertex=vertex,
                    md_log10=x, pred_log10=y, md_rmsf_angstrom=10 ** x, pred_rmsf_angstrom=10 ** y))
    pd.DataFrame(coordinates).to_csv(OUT / 'functional_component_contour_coordinates.csv', index=False)
    np.savez_compressed(OUT / 'functional_component_density_grids.npz', edges_log10=edges,
                       **{g + '_' + k: r[k] for g, r in regions.items() for k in ['counts', 'probability']})
    if before != {str(p): previous.sha256(p) for p in protected}:
        raise ValueError('A preceding figure was changed')
    from PIL import Image
    original = np.asarray(Image.open(PRIOR / 'fig_v3_overall_rmsf_capability.png'))
    modified = np.asarray(Image.open(OUT / 'fig_v3_overall_rmsf_capability_functional_component_contours.png'))
    y, x = np.where(np.any(original != modified, axis=2))
    # The new A caption ends above y=1050; the unchanged C header starts below.
    if x.max() >= 1100 or y.max() >= 1050:
        raise ValueError('Changes extend beyond panel A and its caption')
    audit.update(groups=summaries, coverage=dict(n_samples=len(members), n_residues=len(data),
        operational_single_chain_samples=2193, role_supported_samples=244, unassigned_samples=631),
        method=dict(grouping='Functional components under a stated user convention; partial annotation',
            target_probability=.95, grid_bins=256, gaussian_sigma_log10=.035,
            coordinate_space='log10(MD RMSF), log10(predicted RMSF)',
            normalization='independent per plotted group, including the unassigned collection',
            unit_weight='one weight per residue; each residue inherits its sample-level group',
            missing_two_component_group='Annotation pending, NOT evidence that the group is absent',
            caveat='Single-chain inputs operationally count as one; fusion proteins not audited. All currently annotated >=3 samples are ITrace TCR-pMHC with exactly three components. No causal complexity claim.'),
        originals_unchanged=True, original_sha256=before, only_panel_a_changed=True,
        changed_pixel_bbox_xyxy=[int(x.min()), int(y.min()), int(x.max()), int(y.max())],
        input_sha256={str(p): previous.sha256(p) for p in [MANIFEST, ROLES, ANALYSIS / 'per_chain_metrics.csv',
            ANALYSIS / 'all_residue_predictions.csv.gz', Path(__file__)]})
    (OUT / 'functional_component_statistics_and_QA.json').write_text(json.dumps(audit, indent=2), encoding='utf-8')
    (OUT / 'README_functional_component_contours.md').write_text('''# Functional-component contour preview

This copy changes A only. B/C/D retain their original source-based grouping and are pixel-identical to the original. The original overall, dataset-contour, and chain-count figures are preserved.

## Classification and its limits

This is a partial-annotation preview, not an exhaustive biological classification. A single input polypeptide is operationally counted as one functional component (2,193 samples, 496,093 residues); fusion proteins were not audited and no absent biological assembly partners were inferred. The 244 ITrace inputs (199,668 residues) have complete chain-role maps checked against the evaluated chain IDs. Following the user's convention, TCR alpha+beta form one component, MHC heavy chain+b2m form one, and the antigen peptide is one, hence exactly three per input. These are the currently annotated members of the >=3 group.

The remaining 631 multichain DynaRepo inputs (263,586 residues) remain unassigned; their grey dotted outline is a distribution of unassigned inputs, not a fourth complexity level. No confirmed two-component group has yet been assembled from uniform annotations. 'Pending' does not mean two-component complexes are absent. PPI contact labels, source names, and sequence clusters were not converted into functional-component counts.

Every one of the 959,347 valid held-out residues from 3,068 sample IDs remains in the blue background. No observations were sampled, filtered, or assigned from model errors. Each entire input, including its non-interface residues, inherits a single sample-level group.

## Contours

Each group is independently normalized. All its residues enter a 256 x 256 histogram in log10(MD RMSF), log10(predicted RMSF), then a Gaussian filter with sigma 0.035 log10 units. The outline is the 95% highest-density region of that smoothed estimate, not a confidence interval, protein-level coverage, or an extreme-value envelope. Settings match the previous contour versions. Actual original-residue coverage, polygons, grids, metrics, and source composition are exported alongside the figures.

## Interpretation and QA

The three-component group is currently entirely ITrace, while the unassigned group is entirely DynaRepo. Differences cannot be attributed to functional-component count. The red one-component group comprises ATLAS 1,929 and DynaRepo 264 inputs. Exact group membership, chain-to-component mapping, evidence, and annotation status are retained in sample_functional_component_assignments.csv.

Numerical checks enforce all sample/residue totals, positive finite plotting coordinates, full chain-role coverage, closed unclipped contours, and plausible empirical contour coverage. Original-file hashes and a pixel comparison verify that preceding figures and panels B/C/D are unchanged. PNG is 300 dpi; PDF/SVG retain vector text/lines and rasterize only the density layer. The enlarged 110 mm preview follows the user's existing comparison format; the four-panel master is 183 mm wide. Internal legends occupy the otherwise empty upper-left region and are visually inspected.
''', encoding='utf-8')
    print(json.dumps(dict(output=str(OUT), groups=summaries, only_A_changed=True), indent=2))


if __name__ == '__main__':
    main()
