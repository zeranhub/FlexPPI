# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a) Actual molecular views -> existing tcr_case_v40 render_structures.py -> param inherit
# (b) All10 paired RMSF-change profiles -> LineTrend/plot_trend.py -> param inherit
# RULE: no native-run panels; retain all ten observations and all six MD replicas.
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parent
EXP=ROOT.parent.parent
sys.path.insert(0,str(EXP/'completed_results_20261012'/'plot_runtime'))
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

import json,hashlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from PIL import Image

SOURCE=EXP/'itrace_case_search_only_20261017'
source_csv=SOURCE/'promising_cases_replica_validation/7ow5_7ow6_all_positions_three_replicas.csv'
allpos=pd.read_csv(source_csv,float_precision='round_trip')
p=allpos[allpos.role.eq('peptide')].sort_values('position_a').copy()
reps=pd.concat([pd.read_csv(SOURCE/'replica_recheck'/pdb/f'{pdb}_all_residue_per_replica_labels.csv',float_precision='round_trip') for pdb in ['7ow5','7ow6']],ignore_index=True)
pr=reps[reps.role.eq('peptide')].copy()
assert len(p)==10 and len(pr)==60 and p.rmsf_valid.all()
assert ''.join(p.aa_a)=='VVVGAGGVGK' and ''.join(p.aa_b)=='VVVGADGVGK'
assert p.loc[~p.same_aa,'position_a'].tolist()==[6]
for state,pdb in [('a','7ow5'),('b','7ow6')]:
    raw=pr[pr.pdb.eq(pdb)].groupby('position_in_chain').rmsf_A
    assert np.allclose(raw.mean(),p[f'rmsf_md_{state}'],atol=1e-9)
    assert np.allclose(raw.var(ddof=1),p[f'rmsf_var_{state}'],atol=1e-8)
assert np.allclose(p.rmsf_delta_sem,np.sqrt((p.rmsf_var_a+p.rmsf_var_b)/3),atol=1e-9)

def metrics(data):
    obs=data.rmsf_md_delta.to_numpy();pred=data.rmsf_full_delta.to_numpy()
    mae=float(np.mean(abs(obs-pred)));null=float(np.mean(abs(obs)))
    return {'n':len(data),'rho':float(np.corrcoef(pd.Series(obs).rank(),pd.Series(pred).rank())[0,1]),
            'change_MAE_A':mae,'zero_change_MAE_A':null,'MAE_reduction_fraction':1-mae/null,
            'direction_agreement':int((np.sign(obs)==np.sign(pred)).sum()),
            'mean_MD_change_A':float(obs.mean()),'mean_model_change_A':float(pred.mean())}
stats={'all10':metrics(p),'unchanged9':metrics(p[p.same_aa])}
assert np.isclose(stats['all10']['rho'],.806060606060606)
assert np.isclose(stats['unchanged9']['rho'],.7333333333333333)

def bounds(im):
    pix=np.asarray(im.convert('RGB')); yy,xx=np.where(np.any(pix<235,axis=2))
    return [int(xx.min()),int(yy.min()),int(xx.max()),int(yy.max())]
ims=[Image.open(ROOT/f'{pdb}_binding.png').convert('RGB') for pdb in ['7ow5','7ow6']]
bb=[bounds(im) for im in ims]
shared=(min(v[0] for v in bb)-20,min(v[1] for v in bb)-20,max(v[2] for v in bb)+20,max(v[3] for v in bb)+20)
overlay=Image.open(ROOT/'peptide_superposition.png').convert('RGB');pb=bounds(overlay)
pepcrop=(pb[0]-35,pb[1]-45,pb[2]+35,pb[3]+45)
fig=plt.figure(figsize=(183/25.4,130/25.4))
BLUE='#2166AC';MD='#555555';WT='#CF8C27';MUT='#762A83'
fig.text(.025,.96,'a',weight='bold',fontsize=10)
fig.text(.065,.96,'Single-residue antigen change in the same TCR–HLA complex',fontsize=8)
for i,(im,title,col) in enumerate(zip(ims,['WT · 7OW5','G12D · 7OW6'],[WT,MUT])):
    ax=fig.add_axes([.068+i*.265,.595,.255,.315]);ax.imshow(im.crop(shared));ax.axis('off')
    ax.set_title(title,color=col,pad=1,fontsize=8)
# The peptide is enlarged after the same MHC fit; no additional peptide fit is used.
axp=fig.add_axes([.64,.654,.33,.187]);axp.imshow(overlay.crop(pepcrop));axp.axis('off')
axp.set_title('Peptide superposition',fontsize=8,pad=8)
axp.annotate('G6 / D6',xy=(.57,.52),xytext=(.69,1.02),xycoords='axes fraction',textcoords='axes fraction',
             ha='center',fontsize=7,color=MUT,arrowprops={'arrowstyle':'-','lw':.65,'color':MUT})
fig.text(.806,.602,r'Peptide Cα RMSD = 0.61 Å',fontsize=7,ha='center',color='#444444')
fig.text(.806,.574,'MHC-aligned coordinates',fontsize=7,ha='center',color='#666666')
# WT and G12D sequence labels use precisely aligned letters; mutation is shown without shortening.
for y,name,seq,col,kd in [(.55,'WT','VVVGAGGVGK',WT,r'$K_D$ = 3 μM'),(.517,'G12D','VVVGADGVGK',MUT,r'$K_D$ = 743 pM')]:
    fig.text(.095,y,name,color=col,fontsize=7.5,ha='right')
    for j,aa in enumerate(seq):
        fig.text(.12+j*.025,y,aa,color=col if j==5 else '#444444',fontsize=8,weight='bold' if j==5 else 'normal',ha='center')
    fig.text(.4,y,kd,fontsize=7.5,color=col)
fig.text(.574,.533,'Binding measurements from ref. [28]',fontsize=7,color='#666666')

fig.text(.025,.455,'b',weight='bold',fontsize=10)
fig.text(.065,.455,'Model captures backbone-motion changes beyond the substituted residue',fontsize=8)
ax=fig.add_axes([.095,.13,.87,.245])
x=p.position_a.to_numpy();y=p.rmsf_md_delta.to_numpy();pred=p.rmsf_full_delta.to_numpy();sem=p.rmsf_delta_sem.to_numpy()
ax.axhline(0,color='#999999',lw=.65,zorder=0,ls=(0,(3,2)))
ax.axvspan(5.6,6.4,color=MUT,alpha=.065,lw=0,zorder=0)
ax.errorbar(x-.07,y,yerr=sem,fmt='o-',mfc='white',mec=MD,color=MD,ms=3.5,lw=.8,capsize=2,zorder=3,label='MD mean ± SEM')
ax.plot(x+.07,pred,marker='s',ms=3.5,color=BLUE,lw=1.05,zorder=4,label='Model prediction')
ax.set(xlim=(.6,10.4),ylim=(-1.42,.12),ylabel=r'Δ RMSF (Å)',xticks=x,yticks=[-1.2,-.8,-.4,0],xlabel='Peptide position · G12D − WT')
ax.set_xticklabels([f'{aa}{i}' if i!=6 else 'G6→D6' for aa,i in zip(p.aa_a,x)])
for tick in ax.get_xticklabels():
    if tick.get_text()=='G6→D6':tick.set_color(MUT);tick.set_fontweight('bold')
ax.tick_params(length=3,pad=3)
handles=[Line2D([],[],color=MD,marker='o',mfc='white',ms=3.5,lw=.8,label='MD mean ± SEM'),Line2D([],[],color=BLUE,marker='s',ms=3.5,lw=1.05,label='Model prediction')]
ax.legend(handles=handles,loc='lower right',bbox_to_anchor=(1.,1.015),ncol=2,fontsize=7,handlelength=1.4,columnspacing=1.6,borderaxespad=0)
fig.text(.095,.414,'All 10: ρ = 0.806, MAE = 0.229 Å',fontsize=7)
fig.text(.095,.391,'Unchanged 9: ρ = 0.733, MAE = 0.182 Å; zero-change MAE = 0.303 Å',fontsize=7)
fig.text(.095,.025,'MD: 3 independent trajectories per state. Error bars: SEM of the mean difference.',fontsize=7,color='#555555')
visible=[t.get_text() for t in fig.texts]+[t.get_text() for a in fig.axes for t in [*a.texts,a.xaxis.label,a.yaxis.label,a.title,*a.get_xticklabels(),*a.get_yticklabels()]]
assert not any('\u4e00'<=c<='\u9fff' for s in visible for c in s)
for ext in ['pdf','svg','png']:
    fig.savefig(ROOT/f'fig6_kras_case.{ext}',dpi=450,bbox_inches=None)
p.to_csv(ROOT/'source_all10_peptide_model_and_MD_changes.csv',index=False)
pr.to_csv(ROOT/'source_all10_peptide_six_MD_replicas.csv',index=False)
audit={'status':'passed','dimensions_mm':[183,130],'all_peptide_positions_retained':10,'unchanged_positions':9,'substituted_peptide_position':6,
       'WT_sequence':'VVVGAGGVGK','G12D_sequence':'VVVGADGVGK','MD_independent_replicas_per_state':3,
       'MD_error_bars':'sqrt(sample_variance_WT/3 + sample_variance_G12D/3), ddof=1; no pairing of trajectories',
       'MD_center':'mean of per-replica chainwise CA-aligned RMSF (Å)',
       'model':'family_full_no_ITrace, joint two-output model; only RMSF used in this figure',
       'no_old_strict_model_used':True,'no_significance_tests':True,'metrics':stats,
       'source_file':str(source_csv),'source_SHA256':hashlib.sha256(source_csv.read_bytes()).hexdigest(),
       'complex_shared_crop':shared,'peptide_crop':pepcrop,'structure_source':'actual model-input PDBs; rigid MHC alignment, no peptide refit',
       'binding_source':'ref. [28], doi:10.1038/s41467-022-32811-1; literature experiment, not model output',
       'visible_text':visible}
(ROOT/'figure6_data_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps(stats,indent=2))
