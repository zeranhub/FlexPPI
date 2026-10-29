# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a,b) authentic molecular images -> cross-type inherit -> param inherit (PyMOL)
# (c,d) residue profiles -> LineTrend/plot_trend.py -> param inherit
# (e,f) distance and paired mutations -> cross-type inherit -> param inherit
# (S1,S2) fold/perturbation heatmaps -> heatmap/plot_comparison.py -> param inherit
# (S3) external energy / negative controls -> cross-type inherit -> param inherit
# RULE: "native run" = load pre-rendered PNG via Image.open().ax.imshow().
#       "param inherit" = drawing function below that copies Class A/B/C values.
#       Cumulative counts and per-column normalization are semantically incompatible;
#       use only clean axes, aligned grid and direct labels; retain all raw scales.
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
import json
import shutil
import numpy as np
import pandas as pd
from PIL import Image
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap,Normalize
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import spearmanr,rankdata

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'artifacts/pdz3_alpha3_focus_20260914'
OUT=ROOT/'output/flexppi_pdz3_alpha3_focus_20260914'
ASSETS=OUT/'structure_assets'
BG=['apo1bfe','strip5heb','strip1be9']
BG_LABEL=['1BFE apo backbone','5HEB stripped backbone','1BE9 stripped backbone']
BG_COL=['#222222','#762A83','#2166AC']
CM=LinearSegmentedColormap.from_list('delta',DIVERGING)
CM.set_bad('#D8DBDE')
REGIONS=[(318,323,'318-323'),(330,336,'330-336'),(341,356,'341-356')]

def title(fig,head,sub):
    fig.text(.055,.987,head,fontsize=11.5,weight='bold',ha='left',va='top')
    fig.text(.055,.960,sub,fontsize=7.2,color='#666666',ha='left',va='top')

def panel(ax,label,text):
    ax.text(-.095,1.06,label,transform=ax.transAxes,fontsize=9.5,weight='bold',va='bottom')
    ax.set_title(text,loc='left',pad=9,weight='bold')

def zero(ax):
    ax.axhline(0,color='#999999',lw=.65,ls=(0,(3,3)),zorder=0)

def picture(ax,name):
    im=np.asarray(Image.open(ASSETS/name).convert('RGB'))
    mask=np.any(im<245,axis=2);ys,xs=np.where(mask)
    x0,x1=max(0,xs.min()-45),min(im.shape[1],xs.max()+46)
    y0,y1=max(0,ys.min()-45),min(im.shape[0],ys.max()+46)
    ax.imshow(im[y0:y1,x0:x1],extent=(0,1,0,1))
    ax.set_aspect((y1-y0)/(x1-x0));ax.axis('off')
    return x0,y0,x1,y1

def annotate_structure(ax,viewname,box,labels):
    # PyMOL orthographic camera projection; labels use actual CA coordinates.
    rec=json.loads((ASSETS/'render_provenance.json').read_text())
    view=np.array(rec['views'][viewname]);rot=view[:9].reshape(3,3).T
    halfh=-view[11]*np.tan(np.deg2rad(20/2));halfw=halfh*1700/1300
    x0,y0,x1,y1=box
    for res,text,xytext in labels:
        q=rot@(np.array(rec['label_ca_coordinates'][str(res)])-view[12:15])+view[9:12]
        px=(.5+q[0]/(2*halfw))*1700;py=(.5-q[1]/(2*halfh))*1300
        xy=((px-x0)/(x1-x0),1-(py-y0)/(y1-y0))
        ax.annotate(text,xy=xy,xycoords='axes fraction',xytext=xytext,textcoords='axes fraction',fontsize=7,
            color=BLACK,ha='left',va='center',arrowprops=dict(arrowstyle='-',lw=.65,color='#666666'),
            bbox=dict(boxstyle='square,pad=.10',fc='white',ec='none',alpha=.82))

def residue_line(ax,data,metric,shade=True):
    value={'hchi':'hchi_mean','rmsf':'rmsf_mean'}[metric];spread={'hchi':'hchi_sd','rmsf':'rmsf_sd'}[metric]
    for bg,color,ls,label in zip(BG[:2],BG_COL[:2],['-','--'],BG_LABEL[:2]):
        g=data[data.background==bg].set_index('author_number').sort_index()
        if metric=='hchi':g=g[g.hchi_valid]
        x=g.index.to_numpy();y=g[value].to_numpy();sd=g[spread].to_numpy()
        if shade and bg=='apo1bfe':ax.fill_between(x,y-sd,y+sd,color=color,alpha=.10,lw=0)
        ax.plot(x,y,color=color,lw=1.1,ls=ls,label=label)
    for lo,hi,label in REGIONS:ax.axvspan(lo-.5,hi+.5,color='#2166AC',alpha=.055,lw=0)
    ax.axvspan(390.5,395.5,color='#F1A340',alpha=.17,lw=0)
    zero(ax);ax.set_xlim(305.5,395.5)
    ax.set_xticks([306,318,330,341,356,372,390,395])

def main_figure(s,d,g):
    f=plt.figure(figsize=(183/25.4,220/25.4))
    title(f,'PDZ3 | Predicted α3 constraints are residue-specific','No ligand in any input · fixed-backbone interventions · five frozen models')
    ax=f.add_axes([.055,.710,.43,.204]);box=picture(ax,'pdz3_whole.png');panel(ax,'a','α3-adjacent residues respond differently')
    annotate_structure(ax,'whole',box,[(337,'F337',(.09,.35)),(358,'Q358',(.01,.79)),(334,'E334',(.45,.05)),(399,'α3',(.87,.60))])
    ax=f.add_axes([.56,.710,.39,.204]);box=picture(ax,'pdz3_local.png');panel(ax,'b','Local aromatic packing near F337')
    annotate_structure(ax,'local',box,[(337,'F337',(.07,.36)),(397,'Y397',(.1,.9)),(400,'F400',(.86,.47)),(334,'E334',(.60,.08))])
    cax=f.add_axes([.12,.681,.76,.010]);f.colorbar(mpl.cm.ScalarMappable(norm=Normalize(-.55,.55),cmap=CM),cax=cax,orientation='horizontal',ticks=[-.55,0,.55])
    cax.tick_params(labelsize=6,length=2);cax.set_xlabel('ΔHχ (nats), Δ7ct − intact; orange: removed 396–402 shown as reference',fontsize=6.6,labelpad=2)
    wt=s[(s.kind=='deletion')&(s.context=='d7')&(s.mutant=='WT')]
    ax=f.add_axes([.105,.514,.837,.111]);panel(ax,'c','Side-chain response across every shared residue')
    residue_line(ax,wt,'hchi');ax.set_ylabel('ΔHχ (nats)');ax.tick_params(labelbottom=False)
    handles,labels=ax.get_legend_handles_labels()
    f.legend(handles,labels,loc='center',bbox_to_anchor=(.525,.066),ncol=2,fontsize=6.7,handlelength=2.6)
    apos=wt[(wt.background=='apo1bfe')].set_index('author_number')
    for r,dx,dy in [(334,-10,.1),(337,3,.08),(358,2,.08),(393,-8,.15)]:
        ax.annotate(str(apos.loc[r,'aa'])+str(r),xy=(r,apos.loc[r,'hchi_mean']),xytext=(r+dx,apos.loc[r,'hchi_mean']+dy),fontsize=6.5,arrowprops=dict(arrowstyle='-',lw=.5,color=GREY))
    ax.set_ylim(-.18,.70)
    ax=f.add_axes([.105,.378,.837,.095]);panel(ax,'d','Backbone response: local loop versus distant sites')
    residue_line(ax,wt,'rmsf');ax.set_ylabel('ΔRMSF (Å)');ax.set_xlabel('PSD-95 residue number',labelpad=3)
    ax.text(.01,.90,'Blue windows: MD-defined regions; amber: cut-adjacent 391–395',transform=ax.transAxes,fontsize=6.3,va='top',color='#666666')
    ax=f.add_axes([.105,.133,.345,.170]);panel(ax,'e','Response decreases with distance to α3')
    q=wt[(wt.background=='apo1bfe')&wt.hchi_valid&(wt.author_number<=390)].merge(g,on=['background','author_number','aa'])
    ax.scatter(q.alpha3_distance,q.hchi_mean,s=12,c='#999999',alpha=.65,edgecolor='white',linewidth=.3)
    for r,off in [(337,(4,-3)),(334,(9,8)),(358,(5,5)),(372,(-17,12))]:
        v=q[q.author_number==r].iloc[0];ax.scatter([v.alpha3_distance],[v.hchi_mean],s=20,c=ACCENT_RED,edgecolor='white',linewidth=.4,zorder=5)
        ax.annotate(v.aa+str(r),(v.alpha3_distance,v.hchi_mean),xytext=off,textcoords='offset points',fontsize=6.5)
    zero(ax);ax.set_xlabel('Closest α3 heavy atom (Å)');ax.set_ylabel('ΔHχ (nats)');ax.set_ylim(-.07,.61)
    ax.text(.98,.91,'Core ≤390; n = 61',transform=ax.transAxes,ha='right',fontsize=6.5,color='#666666')
    ax=f.add_axes([.585,.133,.357,.170]);panel(ax,'f','Y397 identity modulates the F337 response')
    muts=['Y397A','Y397F','F400A','F400Y','R399Q','E401Q']
    for bi,bg in enumerate(BG[:2]):
        for j,m in enumerate(muts):
            v=d[(d.background==bg)&(d.context=='full')&(d.kind=='mutation')&(d.mutant==m)&(d.author_number==337)].sort_values('fold').delta_chi_entropy.to_numpy()
            x=j+(bi-.5)*.23
            ax.scatter(x+np.linspace(-.04,.04,5),v,s=9,c=BG_COL[bi],alpha=.50,marker=['o','s'][bi],linewidth=0)
            ax.errorbar(x,v.mean(),yerr=v.std(ddof=1),fmt=['o','s'][bi],ms=3.3,color=BG_COL[bi],capsize=2,lw=.7)
    zero(ax);ax.set_xticks(range(len(muts)));ax.set_xticklabels(muts,rotation=55,ha='right',fontsize=6.5)
    ax.set_ylabel('F337 ΔHχ (nats)');ax.set_ylim(-.16,.30)
    f.text(.055,.038,'Δ7ct: delete 396–402 without relaxing the shared backbone. Positive values: predicted flexibility increases.',fontsize=6.6,color='#555555')
    f.text(.055,.023,'Lines/markers: five-model mean; shading/error bars: model SD. A/G/P omitted only from Hχ panels; all raw values retained.',fontsize=6.4,color='#555555')
    return f

def hm(ax,frame,rows,columns,val,lim,xlabels=None,ylabels=None):
    p=frame.pivot(index=rows,columns=columns,values=val)
    im=ax.imshow(np.ma.masked_invalid(p.to_numpy()),aspect='auto',cmap=CM,vmin=-lim,vmax=lim,interpolation='nearest')
    ax.set_yticks(range(len(p)));ax.set_yticklabels(ylabels or p.index,fontsize=5.8)
    if xlabels is None:
        ticks=[i for i,n in enumerate(p.columns) if int(n)%10==0]
        ax.set_xticks(ticks);ax.set_xticklabels([p.columns[i] for i in ticks],fontsize=6)
    else:ax.set_xticks(range(len(p.columns)));ax.set_xticklabels(xlabels,rotation=90,fontsize=5.7)
    ax.tick_params(length=0);[sp.set_visible(False) for sp in ax.spines.values()]
    return im,p

def robustness_page(d,s):
    f=plt.figure(figsize=(183/25.4,226/25.4));title(f,'PDZ3 | Fold, backbone and truncation-boundary checks','Every model retained; white is near zero, grey is excluded / unavailable')
    dd=d[(d.kind=='deletion')&(d.mutant=='WT')&d.context.isin(['d7','d8'])].copy()
    dd['row']=dd.apply(lambda r:f'{BG.index(r.background)+1} {r.context} F{r.fold+1}',axis=1)
    for metric,y,lim,label in [('chi_entropy',.555,.55,'a'),('rmsf',.19,.75,'b')]:
        q=dd.copy()
        if metric=='chi_entropy':q.loc[~q.hchi_valid,'delta_'+metric]=np.nan
        ax=f.add_axes([.13,y,.78,.29]);panel(ax,label,'All shared residues: '+('ΔHχ (nats)' if metric=='chi_entropy' else 'ΔRMSF (Å)'))
        im,p=hm(ax,q,'row','author_number','delta_'+metric,lim)
        for n in [5,10,15,20,25]:ax.axhline(n-.5,color='white',lw=1)
        cax=f.add_axes([.928,y,.013,.29]);f.colorbar(im,cax=cax,ticks=[-lim,0,lim]);cax.tick_params(labelsize=6)
        ax.set_xlabel('PSD-95 residue number')
    f.text(.10,.893,'1: 1BFE apo backbone    2: 5HEB stripped backbone    3: 1BE9 stripped backbone',fontsize=7)
    f.text(.10,.875,'d7 removes 396–402; d8 removes 395–402 (2026 experimental truncation boundary).',fontsize=7)
    f.text(.10,.10,'Sequence-only I328V normalization is used for 1BFE. No ligand is included in any prediction.',fontsize=7)
    f.text(.10,.077,'F1–F5 are trained models, not biological replicates. Core and mutation comparisons use exact author IDs.',fontsize=7)
    f.text(.10,.054,'The cut-adjacent region 391–395 is shown, not used to claim a distant allosteric response.',fontsize=7)
    return f

def scan_page(s):
    f=plt.figure(figsize=(183/25.4,240/25.4));title(f,'PDZ3 | Systematic mutation responses at unchanged sites','All non-alanine positions in 312–377 mutated to alanine; altered position is masked')
    # Four full matrices retain both metrics, both primary backgrounds, and ALL tested alanines.
    q=s[(s.kind=='mutation')&(s.context=='full')&s.mutant.str.match(r'^[A-Z](3[1-6][0-9]|37[0-7])A$')].copy()
    q=q[q.mutant.str[1:-1].astype(int).between(312,377)]
    q['mpos']=q.mutant.str[1:-1].astype(int);q['row']=q.mpos.astype(str)+' '+q.mutant
    for k,(bg,metric) in enumerate([(BG[0],'hchi'),(BG[1],'hchi'),(BG[0],'rmsf'),(BG[1],'rmsf')]):
        x=.08 if k%2==0 else .575;y=.525 if k<2 else .09
        ax=f.add_axes([x,y,.35,.36]);panel(ax,'abcd'[k],BG_LABEL[BG.index(bg)].replace(' backbone','')+' · '+('ΔHχ' if metric=='hchi' else 'ΔRMSF'))
        g=q[q.background==bg].copy();v=metric+'_mean'
        if metric=='hchi':g.loc[~g.hchi_valid,v]=np.nan
        im,p=hm(ax,g,'row','author_number',v,.15 if metric=='hchi' else .5,ylabels=None)
        ax.set_yticklabels([x.split()[1] for x in p.index],fontsize=5)
        cax=f.add_axes([x+.36,y,.010,.36]);f.colorbar(im,cax=cax,ticks=[-.15,0,.15] if metric=='hchi' else [-.5,0,.5]);cax.tick_params(labelsize=5.5)
        ax.set_xlabel('Responding residue')
    f.text(.055,.032,'Each cell is the mean of five frozen models; values are not row- or column-normalized. Units: nats / Å.',fontsize=6.5)
    f.text(.055,.016,'Color saturation improves visibility of small effects; untruncated numbers and all five models are in source_data.',fontsize=6.5)
    return f

def diagnostic_page(d,exp,es,sat,geometry):
    f=plt.figure(figsize=(183/25.4,220/25.4));title(f,'PDZ3 | Mechanistic tests and external evidence boundaries','Predicted residue flexibility and experimental binding energetics are distinct observables')
    ax=f.add_axes([.10,.66,.84,.23]);panel(ax,'a','V328A / V328I: the predicted local pattern does not reproduce the MD ordering')
    for j,(m,col,ls) in enumerate([('V328A','#B2182B','-'),('V328I','#2166AC','--')]):
        for bi,bg in enumerate(BG[:2]):
            g=d[(d.background==bg)&(d.kind=='mutation')&(d.context=='full')&(d.mutant==m)&d.author_number.between(329,336)]
            v=g.groupby('author_number').delta_rmsf.agg(['mean','std'])
            x=v.index.to_numpy()+(bi-.5)*.1
            ax.plot(x,v['mean'],color=col,ls=ls,lw=1.1,marker=['o','s'][bi],ms=3,label=m+' / '+['apo','stripped'][bi])
    zero(ax);ax.set_xticks(range(329,337));ax.set_xlabel('β2–β3-loop residue');ax.set_ylabel('ΔRMSF (Å)')
    ax.legend(loc='upper right',ncol=2,fontsize=6.3)
    ax=f.add_axes([.10,.365,.35,.21]);panel(ax,'b','F337: experiment (19 substitutions)')
    ee=exp[(exp.assay=='Binding')&(exp.Pos==337)&exp.energy_term.isin(['ddg_NpdzC','ddg_Npdz'])]
    p=ee.pivot(index='id',columns='energy_term',values='total_mean_kcal_mol').sort_index()
    for i,(_,r) in enumerate(p.iterrows()):
        jitter=(i-(len(p)-1)/2)*.010
        ax.plot([0+jitter,1+jitter],r[['ddg_NpdzC','ddg_Npdz']],c='#999999',alpha=.5,lw=.6)
        ax.scatter([jitter,1+jitter],r[['ddg_NpdzC','ddg_Npdz']],s=10,c=['#762A83','#666666'],zorder=3)
    zero(ax);ax.set_xticks([0,1]);ax.set_xticklabels(['With α3','Without α3']);ax.set_ylabel('Experimental ΔΔGbinding\n(kcal/mol)');ax.set_xlim(-.25,1.25)
    ax=f.add_axes([.59,.365,.35,.21]);panel(ax,'c','F337: model context dependence')
    q=sat[(sat.metric=='chi_entropy')].groupby(['background','context','mutant']).mean_abs_response.mean().reset_index()
    for bi,bg in enumerate(BG[:2]):
        p=q[q.background==bg].pivot(index='mutant',columns='context',values='mean_abs_response')
        for i,(_,r) in enumerate(p.iterrows()):
            j=(i-9)*.008+(bi-.5)*.15
            ax.plot([j,1+j],r[['full','d8']],c=BG_COL[bi],lw=.5,alpha=.35)
            ax.scatter([j,1+j],r[['full','d8']],c=BG_COL[bi],s=8,marker=['o','s'][bi],alpha=.6)
    ax.set_xticks([0,1]);ax.set_xticklabels(['With α3','Without α3']);ax.set_xlim(-.3,1.3);ax.set_ylabel('Mean |ΔHχ| at other\ncore sites (nats)')
    ax=f.add_axes([.10,.105,.35,.18]);panel(ax,'d','External site-level correspondence')
    q=es[es.background=='apo1bfe']
    ax.scatter(q.exp_abs_energy,q.hchi_mean.abs(),c='#999999',s=12,alpha=.7,edgecolor='white',lw=.3)
    for r in [334,337,355]:
        row=q[q.author_number==r].iloc[0];ax.annotate(row.aa+str(r),(row.exp_abs_energy,abs(row.hchi_mean)),xytext=(4,3),textcoords='offset points',fontsize=6.5)
    ax.set_xlabel('Median |experimental α3 coupling|\n(kcal/mol)',fontsize=7);ax.set_ylabel('|Predicted ΔHχ| (nats)')
    ax.text(.03,.94,'n = 46; ρ = 0.33',ha='left',transform=ax.transAxes,fontsize=6.5)
    ax=f.add_axes([.59,.105,.35,.18]);panel(ax,'e','Does this exceed a distance-only explanation?')
    paired=q.merge(geometry[['background','author_number','alpha3_distance']],on=['background','author_number'],validate='one_to_one')
    design=np.column_stack([np.ones(len(paired)),rankdata(paired.alpha3_distance)])
    residuals=[]
    for v in [paired.hchi_mean.abs(),paired.exp_abs_energy]:
        ranked=rankdata(v);residuals.append(ranked-design@np.linalg.lstsq(design,ranked,rcond=None)[0])
    vals=[float(spearmanr(paired.hchi_mean.abs(),paired.exp_abs_energy)[0]),
          float(spearmanr(-paired.alpha3_distance,paired.exp_abs_energy)[0]),float(np.corrcoef(*residuals)[0,1])]
    (BASE/'distance_adjusted_external_check.json').write_text(json.dumps(dict(n=len(paired),raw_spearman=vals[0],distance_only_spearman=vals[1],partial_rank_correlation=vals[2]),indent=2),encoding='utf-8')
    labels=['Model vs experiment','−Distance vs experiment','Distance-adjusted model']
    ax.scatter(vals,[2,1,0],s=26,c=['#222222','#999999','#762A83'])
    for y,v,label in zip([2,1,0],vals,labels):
        ax.plot([0,v],[y,y],c='#AAAAAA',lw=.8);ax.text(v+.015,y,f'{v:.2f}',va='center',fontsize=7)
        ax.text(.012,y+.27,label,fontsize=6.3)
    ax.set_yticks([]);ax.spines['left'].set_visible(False)
    ax.set_xlim(0,.57);ax.set_ylim(-.5,2.65);ax.set_xlabel('Rank correlation (descriptive)')
    f.text(.055,.032,'Experimental source: Hidalgo-Carcedo et al., Nature Communications (2026), Supplementary Data 2; no experimental data fitted.',fontsize=6.3)
    f.text(.055,.016,'Distance-adjusted: rank residual correlation after removing the linear effect of α3 distance. No claim of independent validation.',fontsize=6.3)
    return f

def helix_fold_page(d):
    f=plt.figure(figsize=(183/25.4,195/25.4))
    title(f,'PDZ3 | All tested α3 identities, all five models',
          'Unchanged responding residues · every predeclared single-site helix mutation · ΔHχ in nats')
    muts=json.loads((BASE/'protocol.json').read_text())['helix_mutations']
    for k,target in enumerate([337,358,334,355]):
        ax=f.add_axes([.13 if k%2==0 else .60,.535 if k<2 else .13,.305,.32])
        p=np.empty((len(muts),10))
        for i,m in enumerate(muts):
            for j,bg in enumerate(BG[:2]):
                v=d[(d.kind=='mutation')&(d.context=='full')&(d.background==bg)&
                    (d.mutant==m)&(d.author_number==target)].sort_values('fold')
                assert len(v)==5
                p[i,j*5:j*5+5]=v.delta_chi_entropy.to_numpy()
        im=ax.imshow(p,aspect='auto',cmap=CM,vmin=-.30,vmax=.30,interpolation='nearest')
        ax.axvline(4.5,color='white',lw=2)
        ax.set_xticks(range(10));ax.set_xticklabels(['F1','F2','F3','F4','F5']*2,fontsize=6)
        ax.set_yticks(range(len(muts)));ax.set_yticklabels(muts,fontsize=6.5)
        ax.tick_params(length=0);[sp.set_visible(False) for sp in ax.spines.values()]
        aa={337:'F',358:'Q',334:'E',355:'K'}[target]
        panel(ax,'abcd'[k],f'Response at {aa}{target}')
        ax.set_xlabel('1BFE apo           5HEB stripped',fontsize=6.7)
    cax=f.add_axes([.17,.065,.66,.012])
    f.colorbar(im,cax=cax,orientation='horizontal',ticks=[-.30,0,.30])
    cax.tick_params(labelsize=6,length=2)
    f.text(.055,.023,'Five models are not five experiments. Shared color scale; raw, unsaturated data remain in source_data.',fontsize=6.5)
    return f


def save_fixed_canvas(fig,path):
    # Keep the actual PDF/SVG page at the declared 183-mm width; no post-export scaling.
    for ext in ['pdf','png','svg']:
        fig.savefig(str(path)+'.'+ext,bbox_inches=fig.bbox_inches,pad_inches=0,dpi=300)


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    s=pd.read_csv(BASE/'residue_summary.csv');d=pd.read_csv(BASE/'residue_deltas.csv');g=pd.read_csv(BASE/'alpha3_residue_geometry.csv')
    exp=pd.read_csv(BASE/'experiment_2026_data2.csv');es=pd.read_csv(BASE/'experimental_site_comparison.csv');sat=pd.read_csv(BASE/'f337_saturation_summary.csv')
    assert s.n_models.eq(5).all() and len(d)>0 and np.isfinite(d.delta_rmsf).all()
    fig=main_figure(s,d,g);name=OUT/'fig01_pdz3_alpha3_residue_mechanism';save_fixed_canvas(fig,name);plt.close(fig)
    pages=[robustness_page(d,s),scan_page(s),diagnostic_page(d,exp,es,sat,g),helix_fold_page(d)]
    with PdfPages(OUT/'figS01_pdz3_alpha3_controls.pdf') as pdf:
        for i,fig in enumerate(pages,1):
            pdf.savefig(fig,bbox_inches=fig.bbox_inches,pad_inches=0,dpi=300)
            fig.savefig(OUT/f'figS01_page{i}.png',dpi=300,bbox_inches=fig.bbox_inches,pad_inches=0)
            fig.savefig(OUT/f'figS01_page{i}.svg',dpi=300,bbox_inches=fig.bbox_inches,pad_inches=0);plt.close(fig)
    sd=OUT/'source_data';sd.mkdir(exist_ok=True)
    for p in BASE.iterdir():
        if p.suffix in ('.csv','.json') and not p.name.startswith(('audit_existing','psd95_rcsb','pdz3_experimental')):
            shutil.copy2(p,sd/p.name)
    print('MAIN + FOUR SUPPLEMENTARY PAGES COMPLETE',flush=True)

if __name__=='__main__':main()
