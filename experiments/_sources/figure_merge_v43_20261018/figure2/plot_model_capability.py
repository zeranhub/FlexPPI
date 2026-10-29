# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a,b) ECDF -> prior verified fig2/redraw_quantitative_english.py -> param inherit
# (c,d) Categorical paired-MAE dot plot -> prior fig3 within-type gain -> param inherit
# (e,f) All-record chain-group scatter/box -> prior verified fig2 -> param inherit
# RULE: no native-run panels; all frozen source records and valid residue types retained.
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

import json,hashlib,shutil,xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter,NullLocator
from PIL import Image

# Explicit user target: source labels >=9pt so manuscript scaling to 167.6mm stays >=8.2pt.
mpl.rcParams.update({'font.size':9,'axes.titlesize':9.5,'axes.labelsize':9,
                    'xtick.labelsize':9,'ytick.labelsize':9,'legend.fontsize':9,
                    'axes.unicode_minus':False})
OLD=EXP/'reader_revision_v39_20261016'/'english_figures'
UPSTREAM=Path(r'C:\Users\DELL\.codex\worktrees\4666\FlexPPI\paper\revision_20261002\figures\corrected_main_results_preparation\completed\20261003T212155Z_0caa2efe')
ARCHIVE=ROOT/'source_archive';ARCHIVE.mkdir(exist_ok=True)
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
sources={}
def frozen(path):
    path=Path(path);target=ARCHIVE/path.name
    if not target.exists():shutil.copy2(path,target)
    assert sha(path)==sha(target)
    sources[str(path)]={'sha256':sha(path),'archived':str(target)}
    return pd.read_csv(target,float_precision='round_trip')

per=frozen(OLD/'fig2/fig2_all_3063_four_method_metrics.csv')
ecdf_old=frozen(OLD/'fig2/fig2_ecdf_all_24504_points.csv')
gains={t:frozen(OLD/f'fig2/fig2_{t}_all_3063_chain_gains.csv') for t in ['rmsf','chi_entropy']}
pairs_all=frozen(UPSTREAM/'amino_acid_paired_mae_support_and_macros.csv')
# Preserve the former panels' observations as an archive even though they leave the main figure.
for f in ['fig3_all_12252_scatter_points.csv','fig3_all_6126_full_metrics_and_reference.csv',
          'fig3_full_amino_acid_source.csv','fig3_amino_acid_paired_source.csv']:
    frozen(OLD/'fig3'/f)
M=['full','intrachain','static_sequence_complex','aa_only']
ML=dict(zip(M,['Full','Intra','Static','AA']))
MC=dict(zip(M,[CATEGORICAL[0],CATEGORICAL[3],CATEGORICAL[2],CATEGORICAL[5]]))
LS=dict(zip(M,['-','--','-.',':']))
TASKS=['rmsf','chi_entropy'];TN={'rmsf':'RMSF','chi_entropy':'Hχ'};UN={'rmsf':'Å','chi_entropy':'nats'}
AA=list('AGCSTVNDHILFPWYQEMRK')
assert len(per)==24504 and len(ecdf_old)==24504
assert not per.duplicated(['sample_id','metric','variant']).any()
assert per.groupby(['metric','variant']).size().eq(3063).all()
assert set(per.variant)==set(M)
for t in TASKS:
    assert len(gains[t])==3063 and not gains[t].sample_id.duplicated().any()
    assert gains[t].chain_group.astype(str).value_counts().to_dict()=={'1':2193,'3+':448,'2':422}
    assert np.allclose(gains[t].gain,gains[t].intrachain-gains[t].full)
paired=pairs_all[(pairs_all.dataset=='all')&pairs_all.comparator.isin(['intrachain','static_sequence_complex'])].copy()
assert len(paired)==80
paired['relative_MAE_reduction_percent']=100*paired.record_macro_improvement/paired.comparator_record_mean_common_support
assert np.allclose(paired.record_macro_improvement,
                   paired.comparator_record_mean_common_support-paired.full_record_mean_common_support,equal_nan=True)
assert paired.groupby('metric').relative_MAE_reduction_percent.count().to_dict()=={'chi_entropy':36,'rmsf':40}
paired.to_csv(ROOT/'source_type_paired_gains_all80_rows.csv',index=False)

fig=plt.figure(figsize=(183/25.4,184/25.4))
gs=fig.add_gridspec(3,2,left=.105,right=.973,bottom=.105,top=.918,wspace=.37,hspace=.63)
axs=[fig.add_subplot(gs[r,c]) for r in range(3) for c in range(2)]
def label(ax,letter,title):
    ax.set_title(title,loc='left',pad=9)
    ax.text(-.18,1.055,letter,transform=ax.transAxes,fontsize=11,fontweight='bold',ha='left',va='bottom')
def fmt(ax,axis):
    axis=getattr(ax,axis+'axis');axis.set_major_formatter(FuncFormatter(lambda v,p:format(v,'g')));axis.set_minor_locator(NullLocator())
handles=[Line2D([],[],color=MC[m],ls=LS[m],lw=1.5,label=ML[m]) for m in M]
fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.54,.998),ncol=4,columnspacing=1.65,handlelength=2.2)
stats=[];allpoints=[]
for j,t in enumerate(TASKS):
    ax=axs[j]
    for m in M:
        d=per[(per.metric==t)&(per.variant==m)].sort_values('mae',kind='mergesort')
        f=np.arange(1,len(d)+1)/len(d)
        ax.step(d.mae,f,where='post',color=MC[m],ls=LS[m],lw=1.45)
        check=ecdf_old[(ecdf_old.metric==t)&(ecdf_old.variant==m)].sort_values('x',kind='mergesort')
        assert np.allclose(check.x,d.mae) and np.allclose(check.y,f)
        for sid,x,y in zip(d.sample_id,d.mae,f):allpoints.append({'panel':'ab'[j],'sample_id':sid,'metric':t,'method':m,'x':x,'y':y})
        stats.append({'panel':'ab'[j],'metric':t,'comparator_or_method':m,'n':len(d),'mean':d.mae.mean(),'median':d.mae.median(),'min':d.mae.min(),'max':d.mae.max()})
    ax.set(ylim=(0,1.02),ylabel='Fraction of records',xlabel=f'MAE ({UN[t]})')
    if t=='rmsf':
        ax.set_xscale('symlog',linthresh=.1);ax.set_xticks([0,.1,1,10]);fmt(ax,'x')
    else:ax.set_xticks([0,.1,.2,.3])
    ax.set_xlim(0,per[per.metric==t].mae.max()*1.04)
    ax.set_yticks([0,.5,1]);label(ax,'ab'[j],TN[t]+': overall error')

    ax=axs[2+j]
    for k,m in enumerate(['intrachain','static_sequence_complex']):
        d=paired[(paired.metric==t)&(paired.comparator==m)].set_index('amino_acid').reindex(AA)
        valid=d.relative_MAE_reduction_percent.notna().to_numpy()
        x=np.arange(20)+(k-.5)*.38;y=d.relative_MAE_reduction_percent.to_numpy()
        ax.scatter(x[valid],y[valid],s=17,marker='o' if k==0 else '^',color=MC[m],edgecolors='none',zorder=3)
        for aa,xi,yi in zip(np.array(AA)[valid],x[valid],y[valid]):
            allpoints.append({'panel':'cd'[j],'metric':t,'method':m,'amino_acid':aa,'x':xi,'y':yi})
        for aa,r in d.iterrows():
            stats.append({'panel':'cd'[j],'metric':t,'comparator_or_method':m,'amino_acid':aa,
                          'n':r.paired_records_defined,'relative_MAE_reduction_percent':r.relative_MAE_reduction_percent,
                          'full_MAE_common_support':r.full_record_mean_common_support,
                          'comparator_MAE_common_support':r.comparator_record_mean_common_support,
                          'paired_residues':r.valid_residues})
    ax.axhline(0,color=GREY,lw=.7,ls=':');ax.set_xlim(-.65,19.65)
    ax.set_xticks(range(20),AA);ax.tick_params(axis='x',length=2.5,pad=4)
    ax.set_xlabel('Amino acid');ax.set_ylabel('Full MAE reduction (%)')
    if t=='rmsf':ax.set_ylim(-1,21);ax.set_yticks([0,5,10,15,20])
    else:
        ax.set_yscale('symlog',linthresh=10);ax.set_ylim(-5.5,130);ax.set_yticks([-5,0,5,10,30,100]);fmt(ax,'y')
        ax.text(.025,.90,'A/G: no χ',transform=ax.transAxes,ha='left',va='top',fontsize=9,color='#777777')
        for tick in ax.get_xticklabels()[:2]:tick.set_color('#999999')
    label(ax,'cd'[j],TN[t]+': within-type gain')
    if j==0:
        h=[Line2D([],[],marker='o',ls='none',ms=4,color=MC['intrachain'],label='vs Intra'),
           Line2D([],[],marker='^',ls='none',ms=4,color=MC['static_sequence_complex'],label='vs Static')]
        # Shared categorical-plot legend occupies the whitespace between the first and second rows.
        fig.legend(handles=h,loc='center',bbox_to_anchor=(.55,.653),ncol=2,columnspacing=1.75,handletextpad=.35)

    ax=axs[4+j];wide=gains[t].copy();bins=[]
    for i,g in enumerate(['1','2','3+']):
        d=wide[wide.chain_group.astype(str).eq(g)].sort_values('sample_id',kind='mergesort');y=d.gain.to_numpy();bins.append(y)
        # Deterministic offsets preserve all observations without random subsampling.
        x=i+((np.arange(len(d))%67)/66-.5)*.42
        ax.scatter(x,y,s=3,color=MC['full'],alpha=.32,linewidths=0,zorder=2)
        ax.scatter(i,y.mean(),s=35,color=BLACK,marker='D',edgecolor='white',linewidth=.6,zorder=5)
        for sid,xi,yi in zip(d.sample_id,x,y):allpoints.append({'panel':'ef'[j],'sample_id':sid,'metric':t,'method':'intrachain_minus_full','chain_group':g,'x':xi,'y':yi})
        stats.append({'panel':'ef'[j],'metric':t,'chain_group':g,'n':len(y),'mean':y.mean(),'median':np.median(y),
                      'positive':int((y>0).sum()),'zero':int((y==0).sum()),'negative':int((y<0).sum()),'min':y.min(),'max':y.max()})
    bp=ax.boxplot(bins,positions=[0,1,2],widths=.5,showfliers=False,patch_artist=True,
                  medianprops={'color':BLACK,'lw':1},whiskerprops={'color':GREY,'lw':.65},capprops={'color':GREY,'lw':.65})
    for box in bp['boxes']:box.set(facecolor='none',edgecolor=BLACK,linewidth=.8)
    ax.axhline(0,color=GREY,lw=.8,ls='--');ax.set_yscale('symlog',linthresh=.04 if j==0 else .003);fmt(ax,'y')
    ax.set_yticks([-1,-.1,0,.1,1] if j==0 else [-.01,0,.01,.1])
    ax.set_xticks([0,1,2],[f'{g}\nn = {len(v):,}' for g,v in zip(['1','2','3+'],bins)])
    ax.set_xlabel('Number of input chains');ax.set_ylabel('Intra − Full MAE ('+UN[t]+')')
    label(ax,'ef'[j],TN[t]+': chain context')
fig.text(.105,.018,'Points: all records; diamonds: mean; boxes: median and interquartile range.',fontsize=9,color='#555555')

allpoints=pd.DataFrame(allpoints);allpoints.to_csv(ROOT/'source_all_visible_points.csv',index=False)
pd.DataFrame(stats).to_csv(ROOT/'panel_statistics.csv',index=False)
assert allpoints.panel.value_counts().to_dict()=={'a':12252,'b':12252,'e':3063,'f':3063,'c':40,'d':36}
assert (allpoints.loc[allpoints.panel.eq('d'),'y']<0).sum()==2
fig.canvas.draw();renderer=fig.canvas.get_renderer();ignore=set()
for ax in fig.axes:
    for axis in [ax.xaxis,ax.yaxis]:
        lo,hi=sorted(axis.get_view_interval())
        for tick in axis.get_major_ticks()+axis.get_minor_ticks():
            if tick.get_loc()<lo-1e-10 or tick.get_loc()>hi+1e-10:ignore.update([tick.label1,tick.label2])
texts=[t for t in fig.findobj(mpl.text.Text) if t.get_visible() and t.get_text() and t not in ignore]
outside=[]
for text in texts:
    bb=text.get_window_extent(renderer)
    if bb.x0<-.5 or bb.y0<-.5 or bb.x1>fig.bbox.width+.5 or bb.y1>fig.bbox.height+.5:outside.append(text.get_text())
assert not outside,outside
assert min(t.get_fontsize() for t in texts)>=9
assert not any('\u4e00'<=c<='\u9fff' for t in texts for c in t.get_text())
for ext in ['png','pdf','svg']:fig.savefig(ROOT/f'fig2_model_capability.{ext}',bbox_inches=fig.bbox_inches,dpi=450)
svg=ET.parse(ROOT/'fig2_model_capability.svg');ns={'s':'http://www.w3.org/2000/svg'}
assert not svg.findall('.//s:image',ns)
qa={'status':'numeric_and_canvas_passed_visual_review_pending','width_mm':183,'height_mm':184,'dpi':450,
    'minimum_text_pt':min(t.get_fontsize() for t in texts),'minimum_text_pt_at_167_6mm':9*167.6/183,
    'text_outside_canvas':outside,'svg_raster_images':0,'all_visible_points':len(allpoints),
    'panel_points':allpoints.panel.value_counts().to_dict(),'visible_type_negative_gains':2,
    'Hchi_unavailable_types':['A','G'],'no_subsampling':True,'no_training':True,'no_new_inference':True,
    'type_gain_definition':'100*(comparator_record_mean_common_support-full_record_mean_common_support)/comparator_record_mean_common_support; ratio of record-equal means on identical paired support',
    'type_comparators':['intrachain','static_sequence_complex'],'all_comparator_upstream_rows_archived':len(pairs_all),
    'AA_type_points_in_current_main_panels':False,'source_files':sources,
    'panel_map':{'a':'RMSF MAE ECDF','b':'Hchi MAE ECDF','c':'RMSF within-type Full gain','d':'Hchi within-type Full gain','e':'RMSF chain-group paired gain','f':'Hchi chain-group paired gain'}}
(ROOT/'figure2_audit.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps({k:v for k,v in qa.items() if k!='source_files'},ensure_ascii=False,indent=2))
