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
# Parameter inheritance: clean axes and density outlines, implemented with a
# binned Gaussian estimate so all 959,347 observations contribute without sampling.
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
from scipy.ndimage import gaussian_filter

STYLES = {'atlas': '--', 'dynarepo': '--', 'itrace': '-.'}
CONTOUR_COLORS = dict(COLORS, atlas='#B2182B')


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def estimate_regions(data, mass=.95):
    """Estimate each source separately in log10(MD), log10(prediction) space."""
    # Use identical bandwidth and cells for every source. Extend the calculation
    # beyond the visible axes so the Gaussian tails cannot truncate the contours.
    plot_range = (np.log10(.2), np.log10(36))
    bandwidth_log10 = .035
    edges = np.linspace(plot_range[0]-.25, plot_range[1]+.25, 257)
    centers = (edges[:-1]+edges[1:])/2
    sigma_cells = bandwidth_log10/(edges[1]-edges[0])
    regions = {}
    for source in SOURCES:
        group = data.loc[data.dataset.eq(source)]
        points = np.log10(group[['md_rmsf_angstrom', 'pred_rmsf_angstrom']].to_numpy(float))
        counts, _, _ = np.histogram2d(points[:, 0], points[:, 1], bins=(edges, edges))
        if int(counts.sum()) != len(group):
            raise ValueError('Contour calculation lost observations')
        probability = gaussian_filter(counts, sigma=sigma_cells, mode='constant')
        probability /= probability.sum()
        ordered = np.sort(probability.ravel())[::-1]
        threshold = ordered[np.searchsorted(np.cumsum(ordered), mass)]
        enclosed = probability >= threshold
        if enclosed[0].any() or enclosed[-1].any() or enclosed[:, 0].any() or enclosed[:, -1].any():
            raise ValueError('Contour is truncated by calculation grid')
        # Extract the same interpolated polygons that will appear in the figure.
        temporary, axes = plt.subplots()
        contour = axes.contour(centers, centers, probability.T, levels=[threshold])
        paths = [segment.copy() for segment in contour.allsegs[0] if len(segment) >= 4]
        plt.close(temporary)
        inside = np.zeros(len(points), dtype=bool)
        for vertices in paths:
            if not np.allclose(vertices[0], vertices[-1]):
                raise ValueError('Contour is not closed')
            if vertices.min() <= plot_range[0] or vertices.max() >= plot_range[1]:
                raise ValueError('Contour would be clipped by visible axes')
            # Odd-even membership also handles any nested holes correctly.
            inside ^= PlotPath(vertices).contains_points(points)
        fraction = float(inside.mean())
        if not .93 <= fraction <= .98:
            raise ValueError(f'Unexpected empirical coverage for {source}: {fraction}')
        regions[source] = dict(
            paths=paths, counts=counts, probability=probability,
            n_residues=len(group), n_samples=int(group.sample_id.nunique()),
            contour_threshold_probability_per_cell=float(threshold),
            estimated_probability_enclosed=float(probability[enclosed].sum()),
            observed_residue_fraction_enclosed=fraction,
            observed_residues_enclosed=int(inside.sum()),
            contour_components=len(paths), color=CONTOUR_COLORS[source], linestyle=STYLES[source],
        )
    return regions, dict(
        target_probability=mass, coordinate_space='log10(MD RMSF), log10(predicted RMSF)',
        normalization='independent per dataset; each residue has equal weight within its dataset',
        estimator='256 x 256 count histogram followed by an isotropic Gaussian filter',
        gaussian_sigma_log10=bandwidth_log10, gaussian_sigma_cells=float(sigma_cells),
        grid_edges_log10=edges.tolist(), plotted_range_log10=list(plot_range),
        interpretation='95% highest-density region of the smoothed distribution, not a confidence interval or the full range',
    )


def add_regions(ax, regions, enlarged=False):
    handles=[]
    linewidth=1.25 if enlarged else 1.05
    for source in SOURCES:
        halo=[patheffects.Stroke(linewidth=linewidth+1.1, foreground='white'), patheffects.Normal()]
        for vertices in regions[source]['paths']:
            line, = ax.plot(vertices[:, 0], vertices[:, 1], color=CONTOUR_COLORS[source],
                            lw=linewidth, ls=STYLES[source], zorder=5)
            line.set_path_effects(halo)
        handle=Line2D([], [], color=CONTOUR_COLORS[source], lw=linewidth, ls=STYLES[source],
                      label=NAMES[source], path_effects=halo)
        handles.append(handle)
    legend=ax.legend(handles=handles, title='95% density regions', loc='upper left',
                     bbox_to_anchor=(.015, .79), borderaxespad=0, handlelength=2.9,
                     fontsize=7 if enlarged else 5.6, title_fontsize=7 if enlarged else 5.6,
                     labelspacing=.35, frameon=True, facecolor='white', edgecolor='none', framealpha=.9)
    legend.set_zorder(10)


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
    audit['paths']=export(fig,out,'fig_v3_overall_rmsf_capability_dataset_contours')
    return audit


def panel_a(data, regions, out):
    fig, ax=plt.subplots(figsize=(110*MM,112*MM))
    fig.subplots_adjust(left=.16, right=.89, bottom=.18, top=.82)
    heading(ax, 'A', 'Absolute-value agreement', '959,347 held-out residues \xb7 both axes log-scaled')
    density(fig,ax,data.pred_rmsf_angstrom,data.md_rmsf_angstrom,'RMSF [\xc5]',True)
    add_regions(ax,regions,enlarged=True)
    fig.text(.52,.055,'Colored outlines: each dataset\u2019s 95% density region.\nBlue background: all valid test residues, including those outside the outlines.',
             fontsize=6.2,ha='center',va='center',linespacing=1.5)
    return export(fig,out,'fig_A_rmsf_dataset_contours_enlarged')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--analysis',type=Path,default=ROOT/'artifacts/evaluation_v3/analysis')
    parser.add_argument('--original-output',type=Path,default=ROOT/'output/flexppi_v3_overall_capability')
    parser.add_argument('--output',type=Path,default=ROOT/'output/flexppi_v3_overall_capability_dataset_contours')
    args=parser.parse_args()
    if args.output.resolve() == args.original_output.resolve():
        raise ValueError('Use a new directory to preserve the original figure')
    args.output.mkdir(parents=True,exist_ok=True)
    original_files=[args.original_output/f'fig_v3_overall_rmsf_capability.{ext}' for ext in ['png','pdf','svg']]
    hashes_before={str(path):sha256(path) for path in original_files}
    for path in original_files:
        shutil.copy2(path,args.output/f'{path.stem}_original{path.suffix}')
    raw=args.analysis/'all_residue_predictions.csv.gz'
    data=pd.read_csv(raw,usecols=['sample_id','dataset','valid_rmsf','md_rmsf_angstrom','pred_rmsf_angstrom'],dtype={'sample_id':str})
    expected={'atlas':433207,'dynarepo':326472,'itrace':199668}
    if len(data)!=959347 or data.dataset.value_counts().to_dict()!=expected or not data.valid_rmsf.all():
        raise ValueError('Unexpected valid RMSF residue coverage')
    values=data[['md_rmsf_angstrom','pred_rmsf_angstrom']].to_numpy()
    if not np.isfinite(values).all() or values.min()<=0:
        raise ValueError('Invalid RMSF coordinates')
    systems=pd.read_csv(args.original_output/'all_protein_rmsf_metrics_and_baseline.csv')
    folds=pd.read_csv(args.original_output/'paired_five_fold_comparison.csv')
    if len(systems)!=3068 or len(folds)!=5:
        raise ValueError('Missing previously evaluated samples or folds')
    regions,method=estimate_regions(data)
    audit=rmsf_figure(data,systems,folds,regions,args.output)
    audit['enlarged_panel_paths']=panel_a(data,regions,args.output)
    audit['method']=method
    audit['sources']={source:{k:v for k,v in region.items() if k not in ['paths','counts','probability']}
                      for source,region in regions.items()}
    coordinates=[]
    for source,region in regions.items():
        for component,vertices in enumerate(region['paths']):
            for vertex,(x,y) in enumerate(vertices):
                coordinates.append(dict(dataset=source,component=component,vertex=vertex,
                                        md_log10=float(x),pred_log10=float(y),
                                        md_rmsf_angstrom=float(10**x),pred_rmsf_angstrom=float(10**y)))
    pd.DataFrame(coordinates).to_csv(args.output/'dataset_contour_coordinates.csv',index=False)
    np.savez_compressed(args.output/'dataset_contour_density_grids.npz',
                        edges_log10=np.array(method['grid_edges_log10']),
                        **{f'{s}_{k}':regions[s][k] for s in SOURCES for k in ['counts','probability']})
    hashes_after={str(path):sha256(path) for path in original_files}
    if hashes_before!=hashes_after:
        raise ValueError('An original figure was changed')
    audit['original_files_unchanged']=True
    audit['original_sha256']=hashes_before
    audit['input_sha256']={str(raw):sha256(raw),str(Path(original.__file__)):sha256(Path(original.__file__)),
                           str(Path(__file__)):sha256(Path(__file__))}
    (args.output/'dataset_contour_statistics_and_QA.json').write_text(json.dumps(audit,indent=2,ensure_ascii=False),encoding='utf-8')
    notes='''# A \u56fe\u6570\u636e\u96c6\u8f6e\u5ed3\u526f\u672c\n\n\u539f\u56fe\u4fdd\u7559\uff1b\u65b0\u56db\u8054\u56fe\u53ea\u5728 A \u56fe\u6dfb\u52a0\u5404\u6570\u636e\u96c6\u7684\u8f6e\u5ed3\u548c\u56fe\u4f8b\uff0cB\u3001C\u3001D \u590d\u7528\u539f\u56fe\u6570\u503c\u4e0e\u6392\u7248\u3002\u53e6\u6709 A \u56fe\u653e\u5927\u7248\u3002\u84dd\u8272\u5e95\u56fe\u4ecd\u4f7f\u7528\u4e94\u6298\u6d4b\u8bd5\u9884\u6d4b\u5408\u5e76\u540e\u7684\u5168\u90e8 959,347 \u4e2a\u6709\u6548\u6b8b\u57fa\uff0c\u6bcf\u4e2a\u6837\u672c\u53ea\u4f7f\u7528\u5176\u6240\u5c5e\u6d4b\u8bd5\u6298\u7684\u6a21\u578b\u9884\u6d4b\u3002\n\n\u84dd\u8272\u5b9e\u7ebf\u4ee3\u8868 ATLAS\uff0c\u6a59\u8272\u865a\u7ebf\u4ee3\u8868 DynaRepo\uff0c\u7d2b\u8272\u70b9\u5212\u7ebf\u4ee3\u8868 ITrace\uff1b\u4e0e B\u3001C \u56fe\u914d\u8272\u4e00\u81f4\uff0c\u767d\u8272\u886c\u8fb9\u5e2e\u52a9\u533a\u5206\u84dd\u8272\u5e95\u56fe\u4e0a\u7684\u8f6e\u5ed3\u3002\u4e09\u6761\u7ebf\u5708\u51fa\u5404\u6765\u6e90\u7ea6 95% \u7684\u4e3b\u8981\u5206\u5e03\u533a\u57df\uff0c\u91cd\u53e0\u90e8\u5206\u8868\u793a\u8fd9\u4e9b\u6570\u636e\u96c6\u5728\u6b64\u5904\u5747\u6709\u6b8b\u57fa\u3002\u7ebf\u5916\u7684\u6b8b\u57fa\u4ecd\u5728\u84dd\u8272\u5e95\u56fe\u4e2d\u3002\n\n\u8f6e\u5ed3\u662f\u5728 log10(MD RMSF)\u3001log10(\u9884\u6d4b RMSF) \u5750\u6807\u5185\u4f30\u8ba1\u7684\u6700\u9ad8\u5bc6\u5ea6\u533a\u57df\uff0c\u4e0d\u662f 95% \u7f6e\u4fe1\u533a\u95f4\uff0c\u4e5f\u4e0d\u662f\u5168\u90e8\u6781\u7aef\u70b9\u7684\u5916\u8fb9\u754c\u3002\u6bcf\u4e2a\u6765\u6e90\u5404\u81ea\u5f52\u4e00\u5316\uff0c\u4f7f\u7528\u76f8\u540c\u7684 256\xd7256 \u7f51\u683c\u548c 0.035 log10 \u5355\u4f4d\u7684\u9ad8\u65af\u5e73\u6ed1\u6807\u51c6\u5dee\uff1b\u6309\u5e73\u6ed1\u540e\u5bc6\u5ea6\u4ece\u9ad8\u5230\u4f4e\u7d2f\u52a0\u81f3 95%\uff0c\u518d\u63d0\u53d6\u5c01\u95ed\u7b49\u503c\u7ebf\u3002\u6bcf\u4e2a\u6b8b\u57fa\u7b49\u6743\uff0c\u56e0\u6b64\u957f\u86cb\u767d\u5bf9\u8f6e\u5ed3\u8d21\u732e\u66f4\u591a\u3002\n\n\u6570\u503c QA\uff1a\u5168\u90e8\u6b8b\u57fa\u8fdb\u5165\u5e95\u56fe\u4e0e\u5404\u81ea\u6765\u6e90\u7684\u5bc6\u5ea6\u4f30\u8ba1\uff1b\u68c0\u67e5\u5c01\u95ed\u8f6e\u5ed3\u53ca\u663e\u793a\u8303\u56f4\uff0c\u7edf\u8ba1\u771f\u5b9e\u6b8b\u57fa\u843d\u5165\u8f6e\u5ed3\u7684\u6bd4\u4f8b\uff0c\u5e76\u4ee5 SHA-256 \u6838\u5bf9\u539f\u59cb PNG\u3001PDF\u3001SVG \u672a\u53d1\u751f\u53d8\u5316\u3002\u9608\u503c\u3001\u771f\u5b9e\u8986\u76d6\u7387\u3001\u6570\u636e\u6821\u9a8c\u4fe1\u606f\u3001\u8f6e\u5ed3\u5750\u6807\u548c\u5bc6\u5ea6\u7f51\u683c\u5747\u968f\u56fe\u4fdd\u5b58\u3002\n\n\u5404\u6765\u6e90\u5b9e\u9645\u843d\u5165\u8f6e\u5ed3\u7684\u6b8b\u57fa\u6bd4\u4f8b\uff1a\n'''
    for source in SOURCES:
        region=regions[source]
        notes+=f"\n- {NAMES[source]}\uff1a{region['n_residues']:,} \u4e2a\u6b8b\u57fa\uff0c\u5b9e\u9645\u8986\u76d6 {region['observed_residue_fraction_enclosed']:.2%}\u3002\n"
    notes=notes.replace('\u84dd\u8272\u5b9e\u7ebf\u4ee3\u8868 ATLAS', '\u7ea2\u8272\u865a\u7ebf\u4ee3\u8868 ATLAS')
    notes=notes.replace('\u4e0e B\u3001C \u56fe\u914d\u8272\u4e00\u81f4\uff0c', '')
    (args.output/'README_dataset_contours.md').write_text(notes,encoding='utf-8')
    print(json.dumps({'output':str(args.output),'sources':audit['sources'],'original_unchanged':True},indent=2))


if __name__=='__main__':
    main()
