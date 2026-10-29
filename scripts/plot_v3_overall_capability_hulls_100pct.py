# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (A) All-residue density -> MarginalDensity/plot_VariableCorrelation.py -> param inherit
# (B,C) Full protein distributions -> Violin/plot_Violin.py -> param inherit
# (D) Paired fold baseline comparison -> cross-type inherit -> param inherit
# Assets inspected with companion PNGs. The genomic filtering/log regression and
# categorical significance tests do not answer held-out prediction accuracy.
# Inherit clean axes, bounded violins and explicit count labels; no data filters.

# Academic Figure Skill Typography Baseline \u2014 COPY VERBATIM, place at TOP of script
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
# Academic Figure Skill Export Baseline \u2014 COPY VERBATIM
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

# Dataset-contour copy of the completed RMSF evaluation figure.
# Academic Figure Skill asset: KernelDensity/plot_CorrelationDensity1.R + PNG.
# Parameter inheritance: clean axes and distinct source outlines, implemented with a
# convex hull so all 959,347 observations are enclosed without sampling.
# Typography, palette, and export settings are inherited from the original script.
import argparse
import hashlib
import json
import shutil
from pathlib import Path

import plot_v3_overall_capability as original
from plot_v3_overall_capability import (
    ROOT, MM, SOURCES, NAMES, COLORS, BLACK, plt, np, pd,
    Line2D, heading, density, distribution, export,
)
from matplotlib import patheffects
from matplotlib.path import Path as PlotPath
from scipy.spatial import ConvexHull

STYLES = {'atlas': '-', 'dynarepo': '--', 'itrace': '-.'}


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def estimate_regions(data):
    """Enclose all observed points separately by source, in the displayed log space."""
    plot_range = (np.log10(.2), np.log10(36))
    tolerance = 1e-12
    result = {}
    for source in SOURCES:
        group = data.loc[data.dataset.eq(source)]
        values = group[['md_rmsf_angstrom', 'pred_rmsf_angstrom']].to_numpy(float)
        points = np.log10(values)
        hull = ConvexHull(points)
        vertices = points[hull.vertices]
        closed = np.vstack([vertices, vertices[0]])
        if not np.isfinite(closed).all() or closed.min() <= plot_range[0] or closed.max() >= plot_range[1]:
            raise ValueError('Non-finite or clipped hull')
        # A point lies inside/on this convex polygon if it is on the inner side
        # of every edge. Use original points, including the boundary vertices.
        maximum_violation = np.full(len(points), -np.inf)
        for nx, ny, offset in hull.equations:
            signed_distance = points[:, 0]*nx + points[:, 1]*ny + offset
            maximum_violation = np.maximum(maximum_violation, signed_distance)
        outside = maximum_violation > tolerance
        if outside.any():
            raise ValueError(f'{source}: hull does not enclose all observations')
        result[source] = dict(
            paths=[closed], n_residues=len(group), n_samples=int(group.sample_id.nunique()),
            observed_residue_fraction_enclosed=1.0, observed_residues_enclosed=len(group),
            observations_outside=0, hull_vertices=len(vertices), contour_components=1,
            maximum_signed_edge_distance=float(maximum_violation.max()),
            md_range_angstrom=[float(values[:, 0].min()), float(values[:, 0].max())],
            pred_range_angstrom=[float(values[:, 1].min()), float(values[:, 1].max())],
            color=COLORS[source], linestyle=STYLES[source],
        )
    return result, dict(
        target_observed_point_coverage=1.0,
        method='convex hull of every observed point; scipy.spatial.ConvexHull',
        coordinate_space='log10(MD RMSF), log10(predicted RMSF)',
        plotted_range_log10=list(plot_range), boundary_included=True,
        containment_tolerance_log10=tolerance, smoothing=False, points_discarded=0,
        interpretation='Full range of observed valid held-out residues, not a probability contour or confidence region; interior may include empty space',
    )


def add_regions(ax, regions, enlarged=False):
    handles=[]
    linewidth=1.25 if enlarged else 1.05
    for source in SOURCES:
        halo=[patheffects.Stroke(linewidth=linewidth+1.1, foreground='white'), patheffects.Normal()]
        for vertices in regions[source]['paths']:
            line, = ax.plot(vertices[:, 0], vertices[:, 1], color=COLORS[source],
                            lw=linewidth, ls=STYLES[source], zorder=5)
            line.set_path_effects(halo)
        handle=Line2D([], [], color=COLORS[source], lw=linewidth, ls=STYLES[source],
                      label=NAMES[source], path_effects=halo)
        handles.append(handle)
    legend=ax.legend(handles=handles, title='100% observed-point convex hulls', loc='upper left',
                     bbox_to_anchor=(.015, .995), borderaxespad=0, handlelength=1.7, ncol=3, columnspacing=.7, handletextpad=.4,
                     fontsize=7 if enlarged else 5.6, title_fontsize=7 if enlarged else 5.6,
                     labelspacing=.35, frameon=True, facecolor='white', edgecolor='none', framealpha=.9)
    legend.set_zorder(10)
    for label in ax.texts:
        if label.get_text().startswith('r = '):
            label.set_position((.98, .025))
            label.set_ha('right')
            label.set_va('bottom')


def rmsf_figure(data, systems, folds, regions, out):
    # Layout and panels B-D are copied verbatim from the original figure function.
    fig=plt.figure(figsize=(183*MM,163*MM))
    gs=fig.add_gridspec(2,2,left=.095,right=.98,bottom=.095,top=.83,hspace=.66,wspace=.40)
    axes=[fig.add_subplot(gs[i,j]) for i in range(2) for j in range(2)]
    fig.text(.53,.984,'FlexPPI-v3 \xb7 overall RMSF prediction',ha='center',va='top',fontsize=11)
    fig.text(.53,.949,'3,068 held-out protein samples \xb7 959,347 residues \xb7 fixed epoch-40 models',
             ha='center',va='top',fontsize=7.8)
    audit={}
    heading(axes[0],'A','Absolute-value agreement','All valid residues \xb7 both axes log-scaled')
    audit['density']=density(fig,axes[0],data.pred_rmsf_angstrom,data.md_rmsf_angstrom,'RMSF [\xc5]',True)
    add_regions(axes[0], regions)
    heading(axes[1],'B','Profile agreement by protein','One dot per sample; white: median, black: IQR')
    groups=[systems.loc[systems.dataset==s,'pearson'].to_numpy() for s in SOURCES]
    audit['per_protein_r']=distribution(axes[1],groups,[NAMES[s] for s in SOURCES],
                                      [COLORS[s] for s in SOURCES],limits=(0,1.04))
    axes[1].set_ylabel('Per-protein Pearson r')
    heading(axes[2],'C','Absolute error by protein','All samples shown \xb7 logarithmic error axis')
    groups=[systems.loc[systems.dataset==s,'mae'].to_numpy() for s in SOURCES]
    audit['per_protein_mae']=distribution(axes[2],groups,[NAMES[s] for s in SOURCES],
                                        [COLORS[s] for s in SOURCES],logarithmic=True)
    axes[2].set_ylabel('Per-protein MAE [\xc5]')
    heading(axes[3],'D','Improvement over a simple baseline','Each connected pair: one held-out fold')
    for _,row in folds.iterrows():
        offset=(int(row.fold)-2)*.045
        axes[3].plot([offset,1+offset],[row.baseline_mae,row.mae],color='#A5A9AF',lw=.75,zorder=1)
        axes[3].scatter([offset],[row.baseline_mae],color='#7A7F87',s=24,zorder=3)
        axes[3].scatter([1+offset],[row.mae],color='#1B9474',s=24,zorder=3)
    means=folds[['baseline_mae','mae']].mean().to_numpy()
    for x,value in enumerate(means):
        axes[3].plot([x-.16,x+.16],[value,value],color=BLACK,lw=1.5,zorder=4)
        axes[3].text(x,value+.065,f'{value:.3f} \xc5',ha='center',fontsize=7.5)
    improvement=float(1-means[1]/means[0])
    axes[3].text(.98,.98,f'{improvement:.1%} lower MAE',ha='right',va='top',
                 transform=axes[3].transAxes,color='#15785E',fontsize=9)
    axes[3].set_xticks([0,1],['Training-fold\nmean predictor','FlexPPI-v3'],fontsize=7.5)
    axes[3].set(xlim=(-.45,1.45),ylim=(0,1.15),ylabel='Mean per-protein MAE [\xc5]')
    axes[3].tick_params(axis='x',length=0)
    audit['baseline']=dict(fold_mean_baseline_mae=float(means[0]),fold_mean_model_mae=float(means[1]),
                           relative_reduction=improvement,folds=folds.to_dict('records'))
    fig.text(.53,.025,'Every sample is predicted only by its held-out fold model. No illustrative-case selection.',
             fontsize=7,ha='center')
    audit['paths']=export(fig,out,'fig_v3_overall_rmsf_capability_dataset_hulls_100pct')
    return audit


def panel_a(data, regions, out):
    fig, ax=plt.subplots(figsize=(110*MM,112*MM))
    fig.subplots_adjust(left=.16, right=.89, bottom=.18, top=.82)
    heading(ax, 'A', 'Absolute-value agreement', '959,347 held-out residues \xb7 both axes log-scaled')
    density(fig,ax,data.pred_rmsf_angstrom,data.md_rmsf_angstrom,'RMSF [\xc5]',True)
    add_regions(ax,regions,enlarged=True)
    fig.text(.52,.055,'Colored outlines enclose 100% of each dataset\u2019s observed test residues.\nConvex hulls may include empty space and are sensitive to extreme values.',
             fontsize=6.2,ha='center',va='center',linespacing=1.5)
    return export(fig,out,'fig_A_rmsf_dataset_hulls_100pct_enlarged')


def main():
    from PIL import Image
    import fitz
    parser = argparse.ArgumentParser(description='100% observed-point hull version; preserves existing figures')
    parser.add_argument('--analysis', type=Path, default=ROOT/'artifacts/evaluation_v3/analysis')
    parser.add_argument('--original-output', type=Path, default=ROOT/'output/flexppi_v3_overall_capability')
    parser.add_argument('--output', type=Path, default=ROOT/'output/flexppi_v3_overall_capability_dataset_hulls_100pct')
    args = parser.parse_args()
    previous = ROOT/'output/flexppi_v3_overall_capability_dataset_contours'
    protected_dirs = [args.original_output.resolve(), previous.resolve()]
    if args.output.resolve() in protected_dirs:
        raise ValueError('Save the 100% version in a separate directory')
    protected = [p for folder in protected_dirs for p in folder.iterdir() if p.is_file()]
    before = {str(p): sha256(p) for p in protected}
    args.output.mkdir(parents=True, exist_ok=True)
    raw = args.analysis/'all_residue_predictions.csv.gz'
    data = pd.read_csv(raw, usecols=['sample_id', 'dataset', 'valid_rmsf', 'md_rmsf_angstrom', 'pred_rmsf_angstrom'], dtype={'sample_id': str})
    expected = {'atlas': 433207, 'dynarepo': 326472, 'itrace': 199668}
    if len(data) != 959347 or data.dataset.value_counts().to_dict() != expected or not data.valid_rmsf.all():
        raise ValueError('Unexpected valid RMSF residue coverage')
    values = data[['md_rmsf_angstrom', 'pred_rmsf_angstrom']].to_numpy()
    if not np.isfinite(values).all() or values.min() <= 0:
        raise ValueError('Invalid RMSF coordinates')
    systems = pd.read_csv(args.original_output/'all_protein_rmsf_metrics_and_baseline.csv')
    folds = pd.read_csv(args.original_output/'paired_five_fold_comparison.csv')
    if len(systems) != 3068 or len(folds) != 5:
        raise ValueError('Missing samples/folds')
    regions, method = estimate_regions(data)
    audit = rmsf_figure(data, systems, folds, regions, args.output)
    audit['enlarged_panel_paths'] = panel_a(data, regions, args.output)
    audit['method'] = method
    audit['sources'] = {s: {k: v for k, v in r.items() if k != 'paths'} for s, r in regions.items()}
    coords = []
    for source, region in regions.items():
        for index, (x, y) in enumerate(region['paths'][0]):
            coords.append(dict(dataset=source, vertex=index, md_log10=float(x), pred_log10=float(y),
                               md_rmsf_angstrom=float(10**x), pred_rmsf_angstrom=float(10**y)))
    pd.DataFrame(coords).to_csv(args.output/'dataset_hull_coordinates.csv', index=False)
    # Only panel A may change relative to the original full PNG.
    a = np.asarray(Image.open(args.original_output/'fig_v3_overall_rmsf_capability.png'))
    b = np.asarray(Image.open(args.output/'fig_v3_overall_rmsf_capability_dataset_hulls_100pct.png'))
    if a.shape != b.shape:
        raise ValueError('Changed full-figure dimensions')
    y, x = np.where(np.any(a != b, axis=2))
    if not (x.max() < 1100 and y.max() < 900):
        raise ValueError('Unexpected edits outside panel A')
    audit['pixel_verification'] = dict(panels_B_C_D_identical=True,
                                      changed_pixel_bbox_xyxy=[int(x.min()), int(y.min()), int(x.max()), int(y.max())])
    after = {str(p): sha256(p) for p in protected}
    if before != after:
        raise ValueError('Previously delivered files changed')
    audit['all_previously_delivered_files_unchanged'] = True
    audit['protected_sha256'] = before
    audit['input_sha256'] = {str(raw): sha256(raw), str(Path(original.__file__)): sha256(Path(original.__file__)), str(Path(__file__)): sha256(Path(__file__))}
    pdf_checks = []
    for path in args.output.glob('*.pdf'):
        with fitz.open(path) as document:
            page = document[0]
            assert len(document) == 1 and '100%' in page.get_text() and '95%' not in page.get_text()
            pdf_checks.append(dict(path=str(path), width_mm=page.rect.width*25.4/72,
                                   height_mm=page.rect.height*25.4/72, contains_100_percent_label=True))
    audit['pdf_verification'] = pdf_checks
    (args.output/'dataset_hull_statistics_and_QA.json').write_text(json.dumps(audit, indent=2, ensure_ascii=False), encoding='utf-8')
    notes = """# A \u56fe\uff1a100% \u5df2\u89c2\u6d4b\u6d4b\u8bd5\u6b8b\u57fa\u8303\u56f4

\u6b64\u7248\u672c\u4fdd\u7559\u84dd\u8272\u5e95\u56fe\u7684\u5168\u90e8 959,347 \u4e2a\u6709\u6548\u6d4b\u8bd5\u6b8b\u57fa\uff0c\u7528\u84dd\u8272\u5b9e\u7ebf\uff08ATLAS\uff09\u3001\u6a59\u8272\u865a\u7ebf\uff08DynaRepo\uff09\u548c\u7d2b\u8272\u70b9\u5212\u7ebf\uff08ITrace\uff09\u5206\u522b\u56f4\u4f4f\u5404\u6765\u6e90\u7684\u5168\u90e8\u70b9\u3002\u539f\u59cb\u7248\u672c\u548c 95% \u7248\u672c\u4fdd\u7559\uff1bB\u3001C\u3001D \u4e0e\u539f\u59cb\u7248\u672c\u9010\u50cf\u7d20\u4e00\u81f4\u3002

\u8fd9\u91cc\u7684 100% \u6307\u5f53\u524d\u4e94\u6298\u6d4b\u8bd5\u9884\u6d4b\u4e2d\u5b9e\u9645\u89c2\u6d4b\u5230\u7684\u5168\u90e8\u6709\u6548\u6b8b\u57fa\uff0c\u4e0d\u4ee3\u8868\u672a\u6765\u6837\u672c\u6216\u603b\u4f53\u5206\u5e03\u7684 100%\u3002\u5bf9\u6bcf\u4e2a\u6765\u6e90\uff0c\u5728\u539f\u56fe\u7684 log10 \u6a2a\u7eb5\u5750\u6807\u4e2d\u6784\u5efa\u51f8\u5305\uff1a\u53ef\u4ee5\u60f3\u8c61\u7528\u4e00\u6761\u6a61\u76ae\u7b4b\u5305\u4f4f\u6240\u6709\u70b9\uff0c\u6536\u7d27\u540e\u53ea\u63a5\u89e6\u6700\u5916\u56f4\u7684\u70b9\uff0c\u8fb9\u754c\u7531\u8fd9\u4e9b\u70b9\u8fde\u6210\u3002\u8fb9\u754c\u4e0a\u7684\u70b9\u4e5f\u7b97\u5728\u5185\uff0c\u6ca1\u6709\u5220\u9664\u6781\u7aef\u503c\uff0c\u4e5f\u6ca1\u6709\u4e3a\u4e86\u7f8e\u89c2\u5e73\u6ed1\u6216\u79fb\u52a8\u8fb9\u754c\u3002

\u51f8\u5305\u53ef\u80fd\u5305\u5165\u6ca1\u6709\u70b9\u7684\u7a7a\u767d\u533a\u57df\uff1b\u4e24\u4e2a\u51f8\u5305\u91cd\u53e0\u4e0d\u7b49\u4e8e\u91cd\u53e0\u533a\u57df\u5904\u5904\u90fd\u6709\u4e24\u5957\u6570\u636e\u3002\u8fd9\u79cd 100% \u8303\u56f4\u5bf9\u5c11\u6570\u6781\u7aef\u503c\u654f\u611f\uff0c\u9002\u5408\u5c55\u793a\u5b8c\u6574\u8de8\u5ea6\uff1b95% \u7248\u672c\u9002\u5408\u5c55\u793a\u4e3b\u8981\u96c6\u4e2d\u533a\u57df\u3002100% \u7248\u672c\u4e0d\u5e94\u79f0\u4f5c\u5bc6\u5ea6\u7b49\u503c\u7ebf\u6216\u7f6e\u4fe1\u533a\u95f4\u3002

\u5168\u90e8\u771f\u5b9e\u70b9\u5747\u901a\u8fc7\u51f8\u5305\u6bcf\u6761\u8fb9\u7684\u5185\u4fa7\u68c0\u67e5\uff0c\u5305\u542b\u8fb9\u754c\uff0c\u6d6e\u70b9\u5bb9\u5dee\u4e3a 1e-12\uff08log10 \u5750\u6807\u5355\u4f4d\uff09\u3002\u6bcf\u4e2a\u6765\u6e90\u7684\u8986\u76d6\u7387\u3001\u8fb9\u754c\u9876\u70b9\u548c\u539f\u6587\u4ef6\u6821\u9a8c\u503c\u89c1\u968f\u9644 CSV\u3001JSON\u3002
"""
    for source, region in regions.items():
        notes += f"\n{NAMES[source]}\uff1a{region['n_residues']:,} \u4e2a\u6b8b\u57fa\uff0c\u8986\u76d6 100%\uff0c\u8303\u56f4\u5916 0 \u4e2a\u3002\n"
    (args.output/'README_100pct_hulls.md').write_text(notes, encoding='utf-8')
    print(json.dumps({'output': str(args.output), 'sources': audit['sources'], 'previous_files_unchanged': True}, indent=2))


if __name__ == '__main__':
    main()
