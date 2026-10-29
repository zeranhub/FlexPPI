"""A-panel copy stratified by 1, 2, or >=3 input protein/peptide chains."""
# Reuse the verified academic typography, blue density and 95% contour estimator
# settings. The original dataset-based figures and their source data are retained.
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

ROOT=Path(__file__).resolve().parents[1]
ANALYSIS=ROOT/'artifacts/evaluation_v3/analysis'
PRIOR=ROOT/'output/flexppi_v3_overall_capability'
OUT=ROOT/'output/flexppi_v3_overall_capability_chain_count_contours'
GROUPS=['1','2','3plus']
LABELS={'1':'1 chain','2':'2 chains','3plus':'\u22653 chains'}
COLORS={'1':'#B2182B','2':'#C78124','3plus':'#762A83'}
STYLES={'1':'--','2':'--','3plus':'-.'}
EXPECTED={'1':(2193,496093),'2':(424,149016),'3plus':(451,314238)}


def load_groups():
    chains=pd.read_csv(ANALYSIS/'per_chain_metrics.csv',dtype={'sample_id':str,'chain':str})
    chains=chains.loc[chains.metric.eq('rmsf')].copy()
    if chains.duplicated(['sample_id','chain']).any() or not chains.valid_fraction.eq(1).all():
        raise ValueError('Expected unique chains with complete RMSF coverage')
    members=chains.groupby(['sample_id','dataset'],sort=False).agg(
        chain_count=('chain','nunique'),n_residues=('n_residues','sum')).reset_index()
    members['chain_group']=members.chain_count.map(lambda n:'1' if n==1 else ('2' if n==2 else '3plus'))
    data=pd.read_csv(ANALYSIS/'all_residue_predictions.csv.gz',
                     usecols=['sample_id','dataset','valid_rmsf','md_rmsf_angstrom','pred_rmsf_angstrom'],
                     dtype={'sample_id':str})
    data=data.merge(members[['sample_id','dataset','chain_count','chain_group']],
                    on=['sample_id','dataset'],how='left',validate='many_to_one')
    if len(members)!=3068 or len(data)!=959347 or data.chain_group.isna().any() or not data.valid_rmsf.all():
        raise ValueError('Incomplete chain-group assignments')
    observed=data.groupby('sample_id').size().sort_index()
    expected=members.set_index('sample_id').n_residues.sort_index()
    if not np.array_equal(observed.to_numpy(),expected.to_numpy()):
        raise ValueError('Per-sample residue counts differ between inputs')
    for group,(n_samples,n_residues) in EXPECTED.items():
        if int(members.chain_group.eq(group).sum())!=n_samples or int(data.chain_group.eq(group).sum())!=n_residues:
            raise ValueError('Group coverage differs from the audited complete-cache chain counts')
    return data,members


def estimate_regions(data):
    plot_range=(np.log10(.2),np.log10(36))
    edges=np.linspace(plot_range[0]-.25,plot_range[1]+.25,257)
    centers=(edges[:-1]+edges[1:])/2
    sigma=.035/(edges[1]-edges[0])
    regions={}
    for group in GROUPS:
        subset=data.loc[data.chain_group.eq(group)]
        points=np.log10(subset[['md_rmsf_angstrom','pred_rmsf_angstrom']].to_numpy(float))
        counts,_,_=np.histogram2d(points[:,0],points[:,1],bins=(edges,edges))
        if int(counts.sum())!=len(subset): raise ValueError('Lost residues in density estimate')
        probability=gaussian_filter(counts,sigma=sigma,mode='constant')
        probability/=probability.sum()
        ordered=np.sort(probability.ravel())[::-1]
        threshold=ordered[np.searchsorted(np.cumsum(ordered),.95)]
        enclosed=probability>=threshold
        if enclosed[0].any() or enclosed[-1].any() or enclosed[:,0].any() or enclosed[:,-1].any():
            raise ValueError('Clipped density region')
        fig,ax=base.plt.subplots()
        contour=ax.contour(centers,centers,probability.T,levels=[threshold])
        paths=[p.copy() for p in contour.allsegs[0] if len(p)>=4]
        base.plt.close(fig)
        inside=np.zeros(len(points),dtype=bool)
        for p in paths:
            if not np.allclose(p[0],p[-1]) or p.min()<=plot_range[0] or p.max()>=plot_range[1]:
                raise ValueError('Open or visibly clipped contour')
            inside^=PlotPath(p).contains_points(points)
        if not .93<=inside.mean()<=.98: raise ValueError('Unexpected contour coverage')
        regions[group]=dict(paths=paths,counts=counts,probability=probability,
            n_samples=int(subset.sample_id.nunique()),n_residues=len(subset),
            estimated_probability_enclosed=float(probability[enclosed].sum()),
            observed_residue_fraction_enclosed=float(inside.mean()),observed_residues_enclosed=int(inside.sum()),
            threshold_probability_per_cell=float(threshold),contour_components=len(paths),
            color=COLORS[group],linestyle=STYLES[group],
            pooled_metrics=base.stats(subset.pred_rmsf_angstrom,subset.md_rmsf_angstrom))
    return regions,edges


def add_regions(ax,regions,enlarged=False):
    width=1.25 if enlarged else 1.05
    handles=[]
    for group in GROUPS:
        halo=[patheffects.Stroke(linewidth=width+1.1,foreground='white'),patheffects.Normal()]
        for p in regions[group]['paths']:
            line,=ax.plot(p[:,0],p[:,1],color=COLORS[group],lw=width,ls=STYLES[group],zorder=5)
            line.set_path_effects(halo)
        handles.append(Line2D([],[],color=COLORS[group],lw=width,ls=STYLES[group],
                              label=LABELS[group],path_effects=halo))
    legend=ax.legend(handles=handles,title='95% density regions',loc='upper left',
                     bbox_to_anchor=(.015,.79),borderaxespad=0,handlelength=2.9,
                     fontsize=7 if enlarged else 5.6,title_fontsize=7 if enlarged else 5.6,
                     labelspacing=.35,frameon=True,facecolor='white',edgecolor='none',framealpha=.9)
    legend.set_zorder(10)


def plot(data,members,regions):
    fig,ax=base.plt.subplots(figsize=(110*base.MM,112*base.MM))
    fig.subplots_adjust(left=.16,right=.89,bottom=.18,top=.82)
    base.heading(ax,'A','Absolute-value agreement','959,347 held-out residues \u00b7 both axes log-scaled')
    density=base.density(fig,ax,data.pred_rmsf_angstrom,data.md_rmsf_angstrom,'RMSF [\u00c5]',True)
    add_regions(ax,regions,enlarged=True)
    fig.text(.52,.055,'Samples: 1 chain, 2,193 \u00b7 2 chains, 424 \u00b7 \u22653 chains, 451\n'
             '95% regions by input chain count; short peptide chains are included.',
             fontsize=6.2,ha='center',va='center',linespacing=1.5)
    paths=base.export(fig,OUT,'fig_A_rmsf_chain_count_contours_enlarged')
    systems=pd.read_csv(PRIOR/'all_protein_rmsf_metrics_and_baseline.csv')
    folds=pd.read_csv(PRIOR/'paired_five_fold_comparison.csv')
    # Replace only the A overlay and output filename during this call. The
    # unchanged four-panel function preserves B, C, D and their dataset labels.
    overlay_hook,export_hook=previous.add_regions,previous.export
    try:
        previous.add_regions=add_regions
        previous.export=lambda fig,out,name:base.export(fig,out,'fig_v3_overall_rmsf_capability_chain_count_contours')
        full=previous.rmsf_figure(data,systems,folds,regions,OUT)
    finally:
        previous.add_regions,previous.export=overlay_hook,export_hook
    return dict(density=density,panel_a_paths=paths,full_figure=full)


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    protected=[PRIOR/'fig_v3_overall_rmsf_capability.png',
        ROOT/'output/flexppi_v3_overall_capability_dataset_contours/fig_A_rmsf_dataset_contours_enlarged.png']
    before={str(p):previous.sha256(p) for p in protected}
    data,members=load_groups()
    regions,edges=estimate_regions(data)
    audit=plot(data,members,regions)
    members.to_csv(OUT/'sample_chain_group_assignments.csv',index=False)
    crossing=members.groupby(['chain_group','dataset'],sort=False).agg(n_samples=('sample_id','size'),n_residues=('n_residues','sum')).reset_index()
    crossing.to_csv(OUT/'chain_groups_by_dataset.csv',index=False)
    summaries={group:{k:v for k,v in r.items() if k not in ['paths','counts','probability']} for group,r in regions.items()}
    pd.DataFrame([dict(chain_group=g,**{k:v for k,v in r.items() if k not in ['pooled_metrics','n_residues']},**r['pooled_metrics']) for g,r in summaries.items()]).to_csv(OUT/'chain_group_statistics.csv',index=False)
    coordinates=[]
    for group,region in regions.items():
        for component,path in enumerate(region['paths']):
            for vertex,(x,y) in enumerate(path):
                coordinates.append(dict(chain_group=group,component=component,vertex=vertex,md_log10=x,pred_log10=y,md_rmsf_angstrom=10**x,pred_rmsf_angstrom=10**y))
    pd.DataFrame(coordinates).to_csv(OUT/'chain_count_contour_coordinates.csv',index=False)
    np.savez_compressed(OUT/'chain_count_density_grids.npz',edges_log10=edges,
        **{g+'_'+k:r[k] for g,r in regions.items() for k in ['counts','probability']})
    after={str(p):previous.sha256(p) for p in protected}
    if before!=after: raise ValueError('Original figures were changed')
    from PIL import Image
    original=np.asarray(Image.open(PRIOR/'fig_v3_overall_rmsf_capability.png'))
    modified=np.asarray(Image.open(OUT/'fig_v3_overall_rmsf_capability_chain_count_contours.png'))
    y,x=np.where(np.any(original!=modified,axis=2))
    if x.max()>=1100 or y.max()>=900: raise ValueError('Changes extend outside panel A')
    audit.update(groups=summaries,coverage=dict(n_samples=len(members),n_residues=len(data)),
        method=dict(grouping='Number of input protein/peptide chains: 1, 2, >=3',
            target_probability=.95,coordinate_space='log10(MD RMSF), log10(predicted RMSF)',
            grid_bins=256,gaussian_sigma_log10=.035,normalization='independent for each chain-count group',
            unit_weight='one weight per residue; longer inputs contribute more residues',
            chain_assignment_source='Complete-coverage per_chain_metrics.csv; full-cache chain equivalence verified for this frozen run',
            caveat='Input chain count includes short peptides and does not determine biological oligomeric state. Chain count and data source are associated.'),
        originals_unchanged=True,original_sha256=before,
        only_panel_a_changed=True,changed_pixel_bbox_xyxy=[int(x.min()),int(y.min()),int(x.max()),int(y.max())],
        input_sha256={str(p):previous.sha256(p) for p in [ANALYSIS/'per_chain_metrics.csv',ANALYSIS/'all_residue_predictions.csv.gz',Path(__file__)]})
    (OUT/'chain_count_statistics_and_QA.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
    report='''# RMSF A panel grouped by input chain count

The blue background includes all 959,347 valid held-out residues from 3,068 samples. Each residue inherits its sample's full input chain count. Protein and short peptide chains both count. The groups are 1 chain, 2 chains, and 3 or more chains; they do not assert biological monomer or oligomer state.

The outlines enclose 95% of each group's independently normalized, Gaussian-smoothed density estimated in the displayed log10 coordinate plane. All original residues remain in the background, including those outside the outlines. This is not a confidence interval or a percentage of protein samples. The grid and bandwidth are identical to the previous dataset-contour figure.

Red dashed: 1 chain, 2,193 samples / 496,093 residues. Orange dashed: 2 chains, 424 samples / 149,016 residues. Purple dash-dot: 3 or more chains, 451 samples / 314,238 residues.

Source composition: the 1-chain group comprises ATLAS 1,929 and DynaRepo 264 samples. The 2-chain group contains DynaRepo only. The >=3-chain group comprises DynaRepo 207 and ITrace 244 samples. Group differences cannot be attributed to chain count alone.

The copied four-panel figure changes A only; B, C, and D are pixel-identical to the original. B/C remain grouped by dataset. Existing dataset-contour files are unchanged. Group metrics, membership, contour coordinates and source-count tables accompany the figures.
'''
    (OUT/'README_chain_count_contours.md').write_text(report,encoding='utf-8')
    print(json.dumps({'output':str(OUT),'groups':summaries,'only_A_changed':True},indent=2))


if __name__=='__main__': main()
