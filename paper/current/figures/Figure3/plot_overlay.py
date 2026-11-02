# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# f -> MarginalDensity + existing v103 native hexbin -> param inherit
# Restore every nonzero MD density cell, including sparse right tails.
# User requests no contours: use only softer transparent gray/blue density layers.
# The two empirical counts are never smoothed, fitted or subsampled.
from pathlib import Path
from hashlib import sha256
import json, shutil
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
mpl.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','Helvetica','Liberation Sans'],
 'font.size':9,'axes.titlesize':10,'axes.labelsize':9,'xtick.labelsize':9,'ytick.labelsize':9,'legend.fontsize':9,
 'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.6,'xtick.direction':'out','ytick.direction':'out',
 'xtick.major.width':.6,'ytick.major.width':.6,'legend.frameon':False,'axes.unicode_minus':False,
 'pdf.fonttype':42,'svg.fonttype':'none','savefig.dpi':450,'savefig.bbox':None})
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb
from matplotlib.patches import Patch

ROOT=Path(__file__).resolve().parent
OLD=ROOT.parent/'revision_v103_direct_RMSF_Hchi'
N=824988
from palette import BLUE
COLORS={'MD':'#62676B','Full':BLUE}

def load_data():
    for name in ['same_residue_MD_and_Full_all824988.npz','direct_relationship_statistics.json','cohort_all3063_records.csv']:
        assert (ROOT/name).exists(), 'Run prepare_frozen_inputs.py first'
        frozen=json.loads((ROOT/'frozen_input_manifest.json').read_text(encoding='utf-8'))
        expected={r['relative_path']:r['sha256'] for r in frozen['copied_inputs']}
        assert sha256((ROOT/name).read_bytes()).hexdigest()==expected[name]
    with np.load(ROOT/'same_residue_MD_and_Full_all824988.npz',allow_pickle=False) as z:
        arrays={k:z[k] for k in ['MD_RMSF','MD_Hchi','Full_RMSF','Full_Hchi']}
    assert all(len(v)==N and np.isfinite(v).all() for v in arrays.values())
    return arrays,json.loads((ROOT/'direct_relationship_statistics.json').read_text(encoding='utf-8'))

def density_colors(counts):
    peak=max(v.max() for v in counts.values())
    rgba={}
    for label in ['MD','Full']:
        # Shared monotonic scale and identical nonzero floor keep both sparse tails visible.
        alpha=np.where(counts[label]>0,.10+.25*np.log1p(counts[label])/np.log1p(peak),0.)
        rgba[label]=np.column_stack([np.tile(to_rgb(COLORS[label]),(len(alpha),1)),alpha])
        assert np.all(rgba[label][counts[label]>0,3]>=.10)
        assert np.all(rgba[label][counts[label]==0,3]==0)
    return rgba,float(peak)

def draw_overlay(fig,canvas_mm,box,xmax=10.,gridsize=(65,44),marginal_mm=5.,legend=True):
    arrays,stats=load_data();W,H=canvas_mm
    def axis(x,top,w,h):return fig.add_axes([x/W,1-(top+h)/H,w/W,h/H])
    left,top,w,h=box;a=axis(*box)
    counts={};collections={};audit=[]
    for label in ['MD','Full']:
        x=arrays[label+'_RMSF'];y=arrays[label+'_Hchi']
        visible=(x>=0)&(x<=xmax)&(y>=0)&(y<=4.5)
        hb=a.hexbin(x[visible],y[visible],gridsize=gridsize,extent=(0,xmax,0,4.5),mincnt=None,linewidths=0)
        counts[label]=np.asarray(hb.get_array(),dtype=np.int64).copy();collections[label]=hb
        assert counts[label].sum()==visible.sum()
        audit.append({'population':label,'all_positions':N,'visible_positions':int(visible.sum()),
                      'out_of_view_positions':int((~visible).sum()),'bin_count_sum':int(counts[label].sum())})
    offsets=collections['MD'].get_offsets()
    assert np.array_equal(offsets,collections['Full'].get_offsets())
    rgba,peak=density_colors(counts)
    for layer,label in enumerate(['MD','Full']):
        collections[label].set_array(None);collections[label].set_facecolors(rgba[label])
        collections[label].set_edgecolors('none');collections[label].set_zorder(1+layer)
    assert all(np.isfinite(c).all() and c.min()>=0 and c.max()<=1 for c in rgba.values())
    a.set(xlim=(0,xmax),ylim=(0,4.5),xlabel='RMSF (Å)',ylabel='Hχ (nats)')
    a.set_xticks([0,2,4,6,8,10] if xmax==10 else [0,10,20,30,35]);a.set_yticks([0,1,2,3,4])
    a.tick_params(length=2.5,pad=2,width=.6);a.xaxis.labelpad=2.5;a.yaxis.labelpad=2.5
    hx_edges=np.arange(0,35.0001,.2);hy_edges=np.linspace(0,4.5,61)
    hx={k:np.histogram(arrays[k+'_RMSF'],bins=hx_edges)[0]/N for k in ['MD','Full']}
    hy={k:np.histogram(arrays[k+'_Hchi'],bins=hy_edges)[0]/N for k in ['MD','Full']}
    assert all(abs(v.sum()-1)<1e-12 for v in list(hx.values())+list(hy.values()))
    topax=axis(left,top-marginal_mm-1,w,marginal_mm)
    rightax=axis(left+w+1,top,marginal_mm,h)
    # Both population outlines remain visible in the shared marginals.
    for label in ['MD','Full']:
        style='--' if label=='MD' else '-'
        topax.stairs(hx[label],hx_edges,color=COLORS[label],lw=.85,ls=style,baseline=None)
        rightax.stairs(hy[label],hy_edges,orientation='horizontal',color=COLORS[label],lw=.85,ls=style,baseline=None)
    topax.set(xlim=(0,xmax),ylim=(0,max(v.max() for v in hx.values())*1.06));topax.axis('off')
    rightax.set(ylim=(0,4.5),xlim=(0,max(v.max() for v in hy.values())*1.06));rightax.axis('off')
    if legend:
        handles=[Patch(facecolor=COLORS[k],alpha=1.,edgecolor='none') for k in ['MD','Full']]
        labels=['MD','Full']
        fig.legend(handles=handles,labels=labels,loc='lower right',
                 bbox_to_anchor=((left+w+1+marginal_mm)/W,1-(top-marginal_mm-2)/H),
                 ncol=2,fontsize=8,handlelength=1.15,handletextpad=.35,labelspacing=.25,columnspacing=.7,
                 borderpad=0,borderaxespad=0,frameon=False)
    cells=pd.DataFrame(offsets,columns=['RMSF_bin_center_A','Hchi_bin_center_nats'])
    cells['MD_positions']=counts['MD'];cells['Full_positions']=counts['Full']
    for label in ['MD','Full']:
        for col,channel in enumerate(['R','G','B','opacity']):cells[label+'_'+channel]=rgba[label][:,col]
        assert int((counts[label]>0).sum())==int((rgba[label][:,3]>0).sum())
    cells.to_csv(ROOT/f'overlay_cells_x0_{xmax:g}_grid_{gridsize[0]}_{gridsize[1]}.csv',index=False)
    return {'populations':audit,'common_grid':True,'same_residue_cohort':True,'count_scale':'both: nonzero alpha=0.10+0.25*log1p(count)/log1p(shared_max); zero alpha=0',
            'shared_max_count':peak,'colors':COLORS,'overlay':'complete MD gray density + complete Full light blue density; subtle natural alpha blend',
            'all_nonzero_cells_visible':{k:int((rgba[k][:,3]>0).sum())==int((counts[k]>0).sum()) for k in ['MD','Full']},
            'no_density_threshold_for_filled_cells':True,
            'contours':False,
            'no_kernel_smoothing_or_subsampling':True,'marginals':'same bins, denominator824988, MD dashed, Full solid',
            'rho_all_positions':{k:stats['panels'][k]['rho'] for k in ['MD','Full']}}

def finish(fig,stem,audit):
    fig.canvas.draw();renderer=fig.canvas.get_renderer();outside=[]
    for text in fig.findobj(mpl.text.Text):
        if text.get_visible() and text.get_text():
            bb=text.get_window_extent(renderer)
            if bb.x0<-.5 or bb.y0<-.5 or bb.x1>fig.bbox.width+.5 or bb.y1>fig.bbox.height+.5:outside.append(text.get_text())
    assert not outside,outside
    for ext in ['png','pdf','svg']:fig.savefig(ROOT/f'{stem}.{ext}',dpi=450,bbox_inches=fig.bbox_inches)
    audit.update(outside_text=outside,visual_review='pending')
    (ROOT/f'{stem}_QA.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
    plt.close(fig)

def standalone(xmax=10):
    W,H=125.,82.
    fig=plt.figure(figsize=(W/25.4,H/25.4),facecolor='white')
    audit=draw_overlay(fig,(W,H),(14,20,97,49),xmax=xmax,gridsize=(65,44) if xmax==10 else (110,44),marginal_mm=5)
    finish(fig,'RMSF_Hchi_overlay'+('' if xmax==10 else '_full_range'),audit)

if __name__=='__main__':
    ROOT.mkdir(exist_ok=True);standalone(10);standalone(35)
    print('Soft gray/blue alpha overlay, no contours, all nonzero density cells retained.')
