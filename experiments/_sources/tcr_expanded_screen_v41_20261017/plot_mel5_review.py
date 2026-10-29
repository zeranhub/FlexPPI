# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a) All-six-CDR residue change profile -> LineTrend/plot_trend.py -> param inherit
# (b) All natural pMHC interface positions -> LineTrend/plot_trend.py -> param inherit
# The monthly cumulative-area asset is structurally incompatible with residue changes.
# Cross-type inheritance keeps open spines, direct grouping, and restrained mark sizes.
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent/'completed_results_20261012'/'plot_runtime'))
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

import json
import hashlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parent
SRC = ROOT/'mel5_validation'/'plot_sources'
OUT = ROOT/'mel5_review_figure'
OUT.mkdir(exist_ok=True)
read = lambda p: pd.read_csv(p, float_precision='round_trip')
cdr = read(SRC/'mel5_all6_CDR_per_position_source.csv')
pmhc = read(SRC/'mel5_natural_pMHC_interface_source.csv')
context = read(ROOT/'mel5_context'/'all_natural_pMHC_context_comparison.csv')
context = context.set_index(['chain_a','id_a'])
pmhc['hchi_pMHC_only_delta'] = [float(context.loc[(ch,rid),'hchi_pMHC_only_delta']) for ch,rid in zip(pmhc.chain_a,pmhc.id_a)]
assert len(cdr) == 45 and (~cdr.same_aa).sum() == 17
assert cdr.region.str.startswith('CDR3').sum() == 22
assert pmhc.same_aa.all() and int(pmhc.hchi_valid.sum()) == 23
regions = ['CDR1α', 'CDR2α', 'CDR3α', 'CDR1β', 'CDR2β', 'CDR3β']
cdr['group'] = pd.Categorical(cdr.region, categories=regions, ordered=True)
cdr = cdr.sort_values(['group','position_a']).reset_index(drop=True)
pmhc['group'] = pd.Categorical(pmhc.role, categories=['mhc','peptide'], ordered=True)
pmhc = pmhc.sort_values(['group','position_a']).reset_index(drop=True)

def calculate(df, metric):
    d = df[df[f'{metric}_valid'].astype(bool)].copy()
    y = d[f'{metric}_md_delta'].to_numpy()
    p = d[f'{metric}_full_delta'].to_numpy()
    assert len(y) and np.isfinite(y).all() and np.isfinite(p).all()
    return dict(n=len(y), rho=float(np.corrcoef(pd.Series(y).rank(),pd.Series(p).rank())[0,1]),
                MAE=float(np.mean(abs(y-p))), zero_change_MAE=float(np.mean(abs(y))),
                MAE_gain=float(1-np.mean(abs(y-p))/np.mean(abs(y))))

stats = {'all6_CDR_RMSF':calculate(cdr,'rmsf'),
         'all_CDR3_RMSF':calculate(cdr[cdr.region.str.startswith('CDR3')],'rmsf'),
         'natural_pMHC_interface_Hchi':calculate(pmhc,'hchi')}
assert abs(stats['all_CDR3_RMSF']['rho'] - .7786561264822135) < 1e-12
assert abs(stats['natural_pMHC_interface_Hchi']['rho'] - .7549407114624506) < 1e-10

BLUE = CATEGORICAL[0]
MD = '#555555'
fig = plt.figure(figsize=(183/25.4,158/25.4))
fig.text(.105,.976,'MEL5 → α24β17  |  Same MART1 peptide–HLA-A2',fontsize=9,weight='bold')
fig.text(.105,.943,'Model and MD changes after affinity maturation',fontsize=8)

def profile(ax, df, metric, panel, title, ylabel, unit):
    df = df.copy()
    x = np.arange(len(df),dtype=float)
    groups = list(df.groupby('group',observed=True,sort=False))
    # Add visible gaps only between different CDRs or molecules, never within a loop.
    for k,(_,g) in enumerate(groups):
        x[g.index.to_numpy()] += .9*k
    valid = df[f'{metric}_valid'].astype(bool).to_numpy()
    y = df[f'{metric}_md_delta'].to_numpy()
    p = df[f'{metric}_full_delta'].to_numpy()
    sem = df[f'{metric}_delta_sem'].to_numpy()
    assert np.isfinite(sem[valid]).all()
    lo = min(float(np.min(y[valid]-sem[valid])),float(np.min(p[valid])))
    hi = max(float(np.max(y[valid]+sem[valid])),float(np.max(p[valid])))
    if metric=='hchi':
        removed=df.hchi_pMHC_only_delta.to_numpy()
        lo=min(lo,float(removed[valid].min()))
        hi=max(hi,float(removed[valid].max()))
    span = hi-lo
    ax.set_ylim(lo-.10*span,hi+.20*span)
    ax.axhline(0,color='#AAAAAA',lw=.6,zorder=1)
    ax.errorbar(x[valid]-.17,y[valid],yerr=sem[valid],fmt='o',mfc='white',mec=MD,
                color=MD,ms=3,lw=.75,capsize=1.5,zorder=3,label='MD mean ± SEM')
    ax.plot(x[valid]+.02,p[valid],ls='none',marker='s',color=BLUE,ms=3,zorder=4,label='Model')
    if metric=='hchi':
        ax.plot(x[valid]+.21,removed[valid],ls='none',marker='^',color='#CF8C27',ms=3,zorder=4,label='Model · TCR removed')
    ax.set(xlim=(x.min()-.7,x.max()+.7),ylabel=ylabel,xticks=x)
    if metric=='rmsf':
        labels=[f'{a}{i}' if same else f'{a}{i}{b}'
                for a,i,b,same in zip(df.aa_a,df.id_a,df.aa_b,df.same_aa)]
        xlabel='WT residue numbering; substitutions shown as WT–position–engineered'
    else:
        labels=[f'{a}{i}' for a,i in zip(df.aa_a,df.id_a)]
        xlabel='MHC or peptide residue numbering'
        for xi in x[~valid]:
            ax.plot(xi,lo-.045*span,marker='x',ms=3,mew=.7,color=GREY,zorder=3)
    ax.set_xticklabels(labels,rotation=90,ha='center',fontsize=6)
    for same,tick in zip(df.same_aa,ax.get_xticklabels()):
        if not same:
            tick.set_color(ACCENT_RED)
    ax.set_xlabel(xlabel,labelpad=5,fontsize=7)
    for name,g in groups:
        ids=g.index.to_numpy(); left=x[ids].min()-.36; right=x[ids].max()+.36
        ax.plot([left,right],[.968,.968],transform=ax.get_xaxis_transform(),color='#888888',lw=.7)
        name=str(name).replace('mhc','MHC').replace('peptide','Peptide')
        ax.text((left+right)/2,.985,name,transform=ax.get_xaxis_transform(),fontsize=7,ha='center',va='bottom')
    ax.text(-.095,1.17,panel,transform=ax.transAxes,fontsize=10,weight='bold')
    ax.text(0,1.17,title,transform=ax.transAxes,fontsize=8)
    st=calculate(df,metric)
    extra='  |  CDR3: ρ = 0.779 (n = 22)' if metric=='rmsf' else '  |  TCR removed: ρ = 0.057'
    ax.text(0,1.07,f"All {st['n']} valid positions: ρ = {st['rho']:.3f}  |  MAE = {st['MAE']:.3f} {unit}"+extra,
            transform=ax.transAxes,fontsize=7)
    return df

a=fig.add_axes([.105,.567,.868,.273])
cdr = profile(a,cdr,'rmsf','a','Backbone response across all six TCR recognition loops',r'$\Delta$ RMSF (Å)','Å')
b=fig.add_axes([.105,.190,.868,.220])
pmhc = profile(b,pmhc,'hchi','b','Side-chain response on the pMHC side of the interface',r'$\Delta H_\chi$ (nats)','nats')
handles=[Line2D([],[],marker='o',mfc='white',mec=MD,color=MD,lw=.75,ms=3,label='MD mean ± SEM (3 replicas per state)'),
         Line2D([],[],marker='s',color=BLUE,lw=0,ms=3,label='Model'),
         Line2D([],[],marker='^',color='#CF8C27',lw=0,ms=3,label='Model · TCR removed (b)'),
         Line2D([],[],marker='x',color=GREY,lw=0,ms=3,label='No χ angles (b)')]
fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.54,.030),ncol=2,fontsize=7,handlelength=1.4)
fig.text(.105,.010,'All changes: α24β17 − MEL5. Interface: Cα distance ≤ 8 Å in either structure.',fontsize=7)
for ext in ['pdf','svg','png']:
    fig.savefig(OUT/f'mel5_complete_region_review.{ext}',dpi=450,bbox_inches='tight')
visible=[t.get_text() for ax in fig.axes for t in [*ax.texts,ax.xaxis.label,ax.yaxis.label,*ax.get_xticklabels(),*ax.get_yticklabels()]]
visible += [t.get_text() for t in fig.texts]
assert not any('\u4e00' <= ch <= '\u9fff' for s in visible for ch in s)
cdr.to_csv(OUT/'panel_a_all45_CDR_positions.csv',index=False)
pmhc.to_csv(OUT/'panel_b_all_pMHC_interface_positions.csv',index=False)
audit={
    'statistics':stats,'status':'data_checks_passed_visual_review_pending',
    'all_six_CDR_positions':len(cdr),'all_CDR3_positions':int(cdr.region.str.startswith('CDR3').sum()),
    'engineered_CDR_positions':int((~cdr.same_aa).sum()),
    'pMHC_interface_all_positions':len(pmhc),'pMHC_interface_valid_Hchi_positions':int(pmhc.hchi_valid.sum()),
    'error_bars':'SEM of difference of independent three-replica means: sqrt(s_a^2/3+s_b^2/3)',
    'interface_definition':'Union of either-state TCR–pMHC C-alpha pairs within 8 Å; natural matching pMHC positions only. Independently recomputed from both input PDBs: zero annotation mismatches.',
    'replicate_unit':'A separate existing 200 ns MD trajectory; three per state. Residues and frames are not independent experimental replicas.',
    'model':'Joint RMSF/Hchi family fold 0, checkpoint ff2974406d5d1ee07e6dc7118e43e57c0acc9849cc8feba7f5290476fca15ec5; both complexes held out.',
    'metrics':'Spearman correlation and mean absolute error of engineered-minus-WT changes; baseline predicts zero change. Descriptive, no significance test.',
    'units':'RMSF Å; joint rotamer entropy Hchi nats.',
    'omission':'No CDR positions removed. No pMHC interface positions removed; positions without chi angles have x marks and are excluded only from entropy calculations.',
    'source_SHA256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [SRC/'mel5_all6_CDR_per_position_source.csv',SRC/'mel5_natural_pMHC_interface_source.csv']},
    'style':'Param inheritance: LineTrend companion preview inspected. Cumulative time-area data are incompatible; residue profile uses native vector error bars. Figure is for candidate review, not final manuscript layout.'
}
(OUT/'figure_statistics_and_QA.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps(audit,ensure_ascii=False,indent=2))
