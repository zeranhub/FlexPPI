# Academic Figure Skill Asset Confirmation (MarginalDensity and manuscript Fig2 inspected)
# direct joint density + marginal histograms -> native hexbin, Fig2 style -> param inherit
# Template regression/inset tests are semantically incompatible and are not reused.
from pathlib import Path
from hashlib import sha256
import json
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','Helvetica','Liberation Sans'],
 'font.size':9,'axes.titlesize':10,'axes.labelsize':9,'xtick.labelsize':9,'ytick.labelsize':9,'legend.fontsize':9,
 'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.6,'xtick.direction':'out','ytick.direction':'out',
 'xtick.major.width':.6,'ytick.major.width':.6,'legend.frameon':False,'axes.unicode_minus':False,
 'pdf.fonttype':42,'svg.fonttype':'none','savefig.dpi':450,'savefig.bbox':None})
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap,LogNorm
from scipy.stats import spearmanr

ROOT=Path(__file__).resolve().parent
SOURCE=Path(r'E:\Project\FlexPPI\experiments\core_corrected_pipeline\neural_exports\evaluation\oof_per_residue.csv.gz')
N_EXPECTED=824988
CMAP=LinearSegmentedColormap.from_list('paper_blue_mass',['#F0F7FC','#9BC9DE','#2166AC'])
POINTS=ROOT/'same_residue_MD_and_Full_all824988.npz'
STATS=ROOT/'direct_relationship_statistics.json'

def prepare():
    columns=['sample_id','dataset','fold','variant','residue_index','amino_acid','target_rmsf','pred_rmsf','valid_rmsf','target_chi_entropy','pred_chi_entropy','valid_chi_entropy']
    frame=pd.read_csv(SOURCE,usecols=columns,float_precision='round_trip',low_memory=False)
    full=frame[frame.variant.eq('full')].copy()
    fields=['target_rmsf','target_chi_entropy','pred_rmsf','pred_chi_entropy']
    mask=full.valid_rmsf&full.valid_chi_entropy&np.isfinite(full[fields]).all(axis=1)
    a=full.loc[mask].copy();assert len(a)==N_EXPECTED and a.sample_id.nunique()==3063
    assert not a.duplicated(['sample_id','residue_index']).any()
    ids=sorted(a.sample_id.unique());codes=pd.Categorical(a.sample_id,categories=ids).codes
    np.savez_compressed(POINTS,MD_RMSF=a.target_rmsf.to_numpy(),MD_Hchi=a.target_chi_entropy.to_numpy(),
        Full_RMSF=a.pred_rmsf.to_numpy(),Full_Hchi=a.pred_chi_entropy.to_numpy(),record_code=codes,
        residue_index=a.residue_index.to_numpy(),amino_acid=a.amino_acid.to_numpy(str),fold=a.fold.to_numpy(),sample_ids=np.array(ids,str))
    a.groupby(['sample_id','dataset','fold'],as_index=False).size().rename(columns={'size':'valid_positions'}).to_csv(ROOT/'cohort_all3063_records.csv',index=False)
    statistics={'source':str(SOURCE),'source_sha256':sha256(SOURCE.read_bytes()).hexdigest(),
       'selection':'variant full; both target masks valid; all four MD/prediction columns finite; identical rows in both panels',
       'full_rows':len(full),'common_valid_positions':len(a),'records':len(ids),'masked_or_invalid_full_rows':int((~mask).sum()),
       'source_counts':a.groupby('dataset').agg(positions=('residue_index','size'),records=('sample_id','nunique')).to_dict('index'),
       'Pro_positions_retained':int(a.amino_acid.eq('P').sum()),'no_new_training_inference_or_MD':True,
       'prediction_protocol':'existing corresponding-heldout-fold OOF predictions, fixed epoch40; not training-fit or new external evaluation',
       'rho_unit':'pooled residue-level Spearman; descriptive, length-weighted; no independent-residue significance test',
       'density':'hexagonal cell counts with common log count-to-color mapping for MD/Full; no smoothing/subsampling',
       'central_view':'0–10 Angstrom as user reference; out-of-view tails retained in full-range companion, included in all statistics',
       'panels':{}}
    for label,xcol,ycol in [('MD','target_rmsf','target_chi_entropy'),('Full','pred_rmsf','pred_chi_entropy')]:
        x=a[xcol];y=a[ycol]
        centers=a.groupby('amino_acid')[[xcol,ycol]].transform('mean')
        rec_centers=a.groupby(['sample_id','amino_acid'])[[xcol,ycol]].transform('mean')
        statistics['panels'][label]={'rho':float(spearmanr(x,y).statistic),'AA_centered_rho':float(spearmanr(x-centers[xcol],y-centers[ycol]).statistic),
           'record_AA_centered_rho':float(spearmanr(x-rec_centers[xcol],y-rec_centers[ycol]).statistic),
           'RMSF_minmax':[float(x.min()),float(x.max())],'Hchi_minmax':[float(y.min()),float(y.max())],
           'above_10_A':int((x>10).sum()),'above_10_fraction':float((x>10).mean())}
    STATS.write_text(json.dumps(statistics,ensure_ascii=False,indent=2),encoding='utf-8')
    return statistics

def load_data():
    if not POINTS.exists() or not STATS.exists():prepare()
    with np.load(POINTS,allow_pickle=False) as z:arrays={k:z[k] for k in ['MD_RMSF','MD_Hchi','Full_RMSF','Full_Hchi']}
    return arrays,json.loads(STATS.read_text(encoding='utf-8'))

def draw_pair(fig,canvas_mm,boxes,xmax=10.,gridsize=(56,44),marginal_mm=5.,statistics_font=9):
    arrays,stats=load_data();W,H=canvas_mm;collections=[];audit=[]
    def axis(x,top,w,h):return fig.add_axes([x/W,1-(top+h)/H,w/W,h/H])
    # Both marginal panels are frequencies of the complete common denominator.
    hx_edges=np.arange(0,35.0001,.2);hy_edges=np.linspace(0,4.5,61)
    hx={k:np.histogram(arrays[k+'_RMSF'],bins=hx_edges)[0]/N_EXPECTED for k in ['MD','Full']}
    hy={k:np.histogram(arrays[k+'_Hchi'],bins=hy_edges)[0]/N_EXPECTED for k in ['MD','Full']}
    hxmax=max(v.max() for v in hx.values());hymax=max(v.max() for v in hy.values())
    for label,box in zip(['MD','Full'],boxes):
        x0,top,w,h=box;x=arrays[label+'_RMSF'];y=arrays[label+'_Hchi']
        a=axis(x0,top,w,h)
        visible=(x>=0)&(x<=xmax)&(y>=0)&(y<=4.5)
        hb=a.hexbin(x[visible],y[visible],gridsize=gridsize,extent=(0,xmax,0,4.5),mincnt=1,cmap=CMAP,linewidths=0,rasterized=False)
        assert int(hb.get_array().sum())==int(visible.sum()),(label,hb.get_array().sum(),visible.sum())
        collections.append(hb)
        a.set(xlim=(0,xmax),ylim=(0,4.5),xlabel=('MD' if label=='MD' else 'Predicted')+' RMSF (Å)',
            ylabel=('MD' if label=='MD' else 'Predicted')+' Hχ (nats)')
        a.set_xticks([0,2,4,6,8,10] if xmax==10 else [0,10,20,30,35]);a.set_yticks([0,1,2,3,4])
        a.tick_params(length=2.5,pad=2,width=.6);a.xaxis.labelpad=2.5;a.yaxis.labelpad=2.5
        a.text(.97,.035,f'Spearman ρ = {stats["panels"][label]["rho"]:.3f}',transform=a.transAxes,ha='right',va='bottom',fontsize=statistics_font,color='#262626',bbox=dict(fc='white',ec='none',pad=.6,alpha=.84))
        topax=axis(x0,top-marginal_mm-1,w,marginal_mm)
        topax.bar(hx_edges[:-1],hx[label],width=np.diff(hx_edges),align='edge',color='#9BC9DE',edgecolor='none')
        topax.set(xlim=(0,xmax),ylim=(0,hxmax*1.03));topax.axis('off')
        rightax=axis(x0+w+1,top,marginal_mm,h)
        rightax.barh(hy_edges[:-1],hy[label],height=np.diff(hy_edges),align='edge',color='#9BC9DE',edgecolor='none')
        rightax.set(ylim=(0,4.5),xlim=(0,hymax*1.03));rightax.axis('off')
        cells=pd.DataFrame(hb.get_offsets(),columns=['RMSF_bin_center_A','Hchi_bin_center_nats']);cells['positions']=hb.get_array().astype(int)
        cells.to_csv(ROOT/f'density_cells_{label}_x0_{xmax:g}.csv',index=False)
        audit.append({'panel':label,'all_input_positions':len(x),'visible_positions':int(visible.sum()),'out_of_view_positions':int((~visible).sum()),'density_count_sum':int(hb.get_array().sum()),'binned_nonempty_cells':len(cells)})
    shared_max=max(float(c.get_array().max()) for c in collections)
    for c in collections:c.set_norm(LogNorm(vmin=1,vmax=shared_max))
    return {'panels':audit,'shared_log_color_count_range':[1,shared_max],'common_physical_axis_limits':True,'both_same_common_valid_positions':True}

def standalone(xmax,stem):
    W,H=183.,83.
    fig=plt.figure(figsize=(W/25.4,H/25.4),facecolor='white')
    fig.text(3/W,1-2/H,'a',fontsize=11,fontweight='bold',va='top')
    fig.text(98/W,1-2/H,'b',fontsize=11,fontweight='bold',va='top')
    fig.text(47/W,1-2/H,'MD reference',fontsize=10,ha='center',va='top')
    fig.text(142/W,1-2/H,'Full prediction',fontsize=10,ha='center',va='top')
    audit=draw_pair(fig,(W,H),[(14,19,66,50),(109,19,66,50)],xmax=xmax,gridsize=(56,44) if xmax==10 else (100,44))
    fig.canvas.draw();renderer=fig.canvas.get_renderer();outside=[]
    for t in fig.findobj(mpl.text.Text):
        if t.get_visible() and t.get_text():
            b=t.get_window_extent(renderer)
            if b.x0<-.5 or b.y0<-.5 or b.x1>fig.bbox.width+.5 or b.y1>fig.bbox.height+.5:outside.append(t.get_text())
    assert not outside,outside
    for ext in ['png','pdf','svg']:fig.savefig(ROOT/f'{stem}.{ext}',dpi=450,bbox_inches=fig.bbox_inches)
    plt.close(fig)
    audit.update(size_mm=[W,H],outside_text=outside,visual_review='pending')
    (ROOT/f'{stem}_QA.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
    return audit

if __name__=='__main__':
    ROOT.mkdir(parents=True,exist_ok=True)
    stats=prepare()
    standalone(10,'RMSF_Hchi_direct_density')
    standalone(35,'RMSF_Hchi_direct_density_full_range')
    print(json.dumps({'positions':stats['common_valid_positions'],'records':stats['records'],'rho':{k:v['rho'] for k,v in stats['panels'].items()},'out_of_view_central':{k:v['above_10_A'] for k,v in stats['panels'].items()}},ensure_ascii=False,indent=2))
