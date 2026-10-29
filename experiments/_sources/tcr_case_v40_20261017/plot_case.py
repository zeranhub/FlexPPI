# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a) All66 frozen joint-Hchi profile -> LineTrend/plot_trend.py -> param inherit
# (b) All151 defined-chi scatter -> cross-type inherit -> param inherit
# (c) Original432 frozen MD Fstar curve rows -> LineTrend/plot_trend.py -> param inherit
# RULE: no native-run panels; retain all source observations without smoothing or selection.
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'completed_results_20261012'/'plot_runtime'))
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
import json,hashlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from PIL import Image
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'figure6';OUT.mkdir(exist_ok=True)
read=lambda p:pd.read_csv(p,float_precision='round_trip')
data=read(ROOT/'change_screen/pair_2ak4_4jry_all_positions.csv')
p=data[data.role.eq('peptide')].sort_values('position_a').copy()
rep=read(ROOT/'sb27_sb47/all13_peptide_per_replica_labels.csv')
ctx=read(ROOT/'sb27_sb47/context_inference.csv')
q=rep[rep.position.eq(7)]
assert len(p)==13 and p.hchi_valid.sum()==11 and len(rep)==78
assert ''.join(p.aa_a)=='LPEPLPQGQLTAY' and p.same_aa.all()
assert np.isclose(np.corrcoef(p.rmsf_md_delta.rank(),p.rmsf_full_delta.rank())[0,1],.945054945054945)
BLUE='#2166AC';MD='#555555';ORANGE='#CF8C27';RED='#B2182B'
fig=plt.figure(figsize=(183/25.4,158/25.4))
# Top: matched molecular context and an actual same-model input perturbation.
a=fig.add_axes([.025,.53,.575,.405]);a.axis('off')
a.text(0,1.09,'a',weight='bold',fontsize=10,transform=a.transAxes)
a.text(.08,1.09,'Same peptide–MHC, different TCR docking',fontsize=8,transform=a.transAxes)
ims=[Image.open(OUT/f'{pdb}_binding.png').convert('RGB') for pdb in ['2ak4','4jry']]
bounds=[]
for im in ims:
    pix=np.asarray(im);yy,xx=np.where(np.any(pix<235,axis=2));bounds.append((xx.min(),yy.min(),xx.max(),yy.max()))
crop=(min(b[0] for b in bounds)-35,min(b[1] for b in bounds)-30,max(b[2] for b in bounds)+35,max(b[3] for b in bounds)+30)
crop=tuple(int(v) for v in crop)
for i,(im,title) in enumerate(zip(ims,['SB27 · 2AK4','SB47 · 4JRY'])):
    ax=fig.add_axes([.038+i*.277,.567,.266,.353]);ax.imshow(im.crop(crop));ax.axis('off');ax.set_title(title,pad=1,fontsize=8)
a.text(.10,.005,'TCRα',color='#4D8FB5',transform=a.transAxes,fontsize=7)
a.text(.28,.005,'TCRβ',color='#6EAD91',transform=a.transAxes,fontsize=7)
a.text(.47,.005,'MHC',color='#777777',transform=a.transAxes,fontsize=7)
a.text(.64,.005,'Peptide',color=ORANGE,transform=a.transAxes,fontsize=7)
a.text(.88,.005,'Q7',color=RED,transform=a.transAxes,fontsize=7)

b=fig.add_axes([.708,.59,.255,.30])
b.text(-.21,1.31,'b',weight='bold',fontsize=10,transform=b.transAxes)
b.text(-.03,1.31,'Q7 depends on TCR context',fontsize=8,transform=b.transAxes)
for i,pdb in enumerate(['2ak4','4jry']):
    vals=q[q.pdb.eq(pdb)].Hchi_nats.to_numpy()
    b.scatter(i+np.array([-.065,0,.065]),vals,s=16,facecolors='white',edgecolors=MD,lw=.8,zorder=5)
mdmeans=q.groupby('pdb').Hchi_nats.mean().reindex(['2ak4','4jry']).to_numpy()
b.plot([0,1],mdmeans,color=MD,lw=.9,ls=':',zorder=2)
for input_name,color,marker,ls in [('full',BLUE,'s','-'),('pMHC_only',ORANGE,'D','--')]:
    vals=ctx[ctx.position.eq(7)&ctx.input.eq(input_name)].set_index('pdb').reindex(['2ak4','4jry']).Hchi.to_numpy()
    b.plot([0,1],vals,color=color,lw=1.2,marker=marker,ms=4,ls=ls,zorder=4)
b.set(xlim=(-.24,1.24),ylim=(1.35,3.4),xticks=[0,1],xticklabels=['SB27','SB47'],ylabel=r'Q7 $H_\chi$ (nats)')
b.set_yticks([1.5,2,2.5,3]);b.tick_params(length=3)
handles=[Line2D([],[],color=MD,marker='o',mfc='white',ms=3,lw=.9,ls=':',label='MD (3 replicas)'),Line2D([],[],color=BLUE,marker='s',ms=3,lw=1,label='Model · full complex'),Line2D([],[],color=ORANGE,marker='D',ms=3,lw=1,ls='--',label='Model · TCR removed')]

# Bottom: every peptide position and every valid side chain, including Q9.
for label,metric,rect,ylim,yticks,ylabel,title in [
    ('c','rmsf',[.085,.115,.389,.275],(-.24,.82),[-.2,0,.2,.4,.6,.8],r'$\Delta$ RMSF (Å)','Backbone change · all 13 positions'),
    ('d','hchi',[.594,.115,.379,.275],(-.25,1.36),[0,.4,.8,1.2],r'$\Delta H_\chi$ (nats)','Side-chain change · all 11 valid positions')]:
    ax=fig.add_axes(rect)
    ax.text(-.16,1.25,label,weight='bold',fontsize=10,transform=ax.transAxes)
    ax.text(-.01,1.25,title,fontsize=8,transform=ax.transAxes)
    valid=p[f'{metric}_valid'].to_numpy(bool)
    x=p.position_a.to_numpy();y=p[f'{metric}_md_delta'].to_numpy();pred=p[f'{metric}_full_delta'].to_numpy();sem=p[f'{metric}_delta_sem'].to_numpy()
    ax.axhline(0,color='#BBBBBB',lw=.6,zorder=0)
    ax.axvspan(6.55,7.45,color='#B2182B',alpha=.065,lw=0)
    ax.errorbar(x[valid]-.10,y[valid],yerr=sem[valid],fmt='o',mfc='white',mec=MD,color=MD,ms=3.3,lw=.75,capsize=1.8,zorder=3)
    ax.plot(x[valid]+.10,pred[valid],ls='none',marker='s',ms=3.3,color=BLUE,zorder=4)
    ax.set(xlim=(.5,13.5),ylim=ylim,ylabel=ylabel,xticks=x,yticks=yticks,xlabel='Peptide position · SB47 − SB27')
    ax.set_xticklabels([f'{aa}{i}' for aa,i in zip(p.aa_a,x)],rotation=55,ha='right',rotation_mode='anchor')
    for tick in ax.get_xticklabels():
        if tick.get_text()=='Q7':tick.set_color(RED);tick.set_fontweight('bold')
    if metric=='rmsf':
        rho=np.corrcoef(pd.Series(y).rank(),pd.Series(pred).rank())[0,1]
        mae=np.mean(abs(y-pred));text=f'ρ = {rho:.3f}   |   Change MAE = {mae:.3f} Å'
    else:
        text='G8 / A12: no χ angles'
        for xi in x[~valid]:ax.plot(xi,-.18,marker='x',color='#999999',ms=3,mew=.6)
    ax.text(0,1.065,text,fontsize=7,transform=ax.transAxes)
delta_handles=[Line2D([],[],color=MD,marker='o',mfc='white',ms=3,lw=.8,ls=':',label='MD reference'),handles[1],Line2D([],[],color=ORANGE,marker='D',ms=3,lw=1,ls='--',label='Model · TCR removed (b)')]
fig.legend(handles=delta_handles,loc='lower center',bbox_to_anchor=(.52,-.006),ncol=3,fontsize=7,handlelength=1.6)
fig.savefig(OUT/'fig6_tcr_case.pdf',dpi=450,bbox_inches='tight')
fig.savefig(OUT/'fig6_tcr_case.svg',dpi=450,bbox_inches='tight')
fig.savefig(OUT/'fig6_tcr_case.png',dpi=450,bbox_inches='tight')
visible=[t.get_text() for ax in fig.axes for t in [*ax.texts,ax.xaxis.label,ax.yaxis.label,ax.title,*ax.get_xticklabels(),*ax.get_yticklabels()]]
assert not any('\u4e00'<=c<='\u9fff' for s in visible for c in s)
p.to_csv(OUT/'all13_model_and_MD_changes.csv',index=False)
audit={'status':'passed','peptide_sequence':'LPEPLPQGQLTAY','n_RMSF':13,'n_Hchi':11,'no_chi':['G8','A12'],'all_positions_retained':True,'RMSF_delta_spearman':float(rho),'RMSF_change_MAE_A':float(mae),'RMSF_direction_agreement':int((np.sign(p.rmsf_md_delta)==np.sign(p.rmsf_full_delta)).sum()),'MD_error_bars':'independent-groups mean-difference SEM sqrt(s1^2/3+s2^2/3); n=3 per structure','Q7_MD':mdmeans.tolist(),'Q7_full':ctx[ctx.position.eq(7)&ctx.input.eq('full')].sort_values('pdb').Hchi.tolist(),'Q7_TCR_removed':ctx[ctx.position.eq(7)&ctx.input.eq('pMHC_only')].sort_values('pdb').Hchi.tolist(),'structure_crop_shared_pixels':crop,'no_significance_test':True,'scientific_message':'Test changes in the same peptide–MHC sequence under two TCR docking contexts; same-model removal tests input dependence, not physical apo dynamics.'}
(OUT/'figure6_data_audit.json').write_text(json.dumps(audit,indent=2),encoding='utf8')
(OUT/'visible_text.json').write_text(json.dumps(visible,ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps(audit,indent=2))
