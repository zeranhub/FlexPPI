# Academic Figure Skill Asset Confirmation
# MAE/AP point comparisons: parameter inheritance from BarComparison/plot_comparison_GeneRegulatory.py.
# Paired record gain: parameter inheritance from PairedBoxScatter/plot_PairedBoxScatter1.R.
# References were visually inspected. Data/statistical semantics differ from the templates:
# use frozen record means/AP intervals, and all paired records without new tests or simulated data.
# This independent candidate does not replace the v61 manuscript.

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
import hashlib
import json
import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from publication_style import BLUE, WARM, TEAL

ROOT=Path(__file__).resolve().parent
ARCHIVE=ROOT/'source_archive'
METHODS=['Full','Intra','Static','AA']
COLORS=dict(Full=BLUE,Intra=WARM,Static=TEAL,AA='#999999')
MARKERS=dict(Full='o',Intra='s',Static='D',AA='^')
summary=json.loads((ARCHIVE/'summary.json').read_text(encoding='utf-8'))
records=pd.read_csv(ARCHIVE/'per_record_metrics.csv',float_precision='round_trip')
within=pd.read_csv(ARCHIVE/'within_type_record_macro_metrics.csv',float_precision='round_trip')
assert len(records)==3063*2*4 and len(within)==3063*4
assert not records.duplicated(['record_code','target','method']).any()
assert not within.duplicated(['record_code','method']).any()
assert np.isfinite(records[['ap','auroc','mae','spearman']].to_numpy()).all()
assert np.isfinite(within[['ap','auroc','prevalence']].to_numpy()).all()
assert records.record_code.nunique()==within.record_code.nunique()==3063
assert records[records.method.eq('Full') & records.target.eq('rmsf')].dataset.value_counts().to_dict()==dict(atlas=1929,dynarepo=890,itrace=244)
TASKS=['rmsf','chi_entropy','within_type_hchi']
TITLES=['RMSF','Hχ','Hχ within type']
meta={t:summary['metrics'][t] for t in TASKS[:2]}
meta[TASKS[2]]=summary['within_type']
frames={t:records[records.target.eq(t)] for t in TASKS[:2]}
frames[TASKS[2]]=within
for task,frame in frames.items():
    assert frame.groupby('method').size().to_dict()==dict.fromkeys(METHODS,3063)
    for method in METHODS:
        assert np.isclose(frame[frame.method.eq(method)].ap.mean(),meta[task][method]['mean_ap'],rtol=0,atol=1e-14)
        if task in TASKS[:2]:
            assert np.isclose(frame[frame.method.eq(method)].mae.mean(),meta[task][method]['macro_mae'],rtol=0,atol=1e-14)

mpl.rcParams.update({'font.size':9,'axes.labelsize':9,'xtick.labelsize':9,'ytick.labelsize':9,
                     'legend.fontsize':9,'axes.unicode_minus':False})
W,H=183.,185.
fig=plt.figure(figsize=(W/25.4,H/25.4))
def axes_mm(x,y,w,h):return fig.add_axes([x/W,y/H,w/W,h/H])
def text_mm(x,y,s,**kw):return fig.text(x/W,y/H,s,**kw)
def heading(x,y,letter,title):
    text_mm(x-16,y,letter,fontsize=11,fontweight='bold',ha='left',va='bottom')
    text_mm(x,y,title,fontsize=10,ha='left',va='bottom')
def model_rows(ax,show_labels=True):
    # Direct labels identify methods; no legend is placed over the data.
    ax.set_ylim(-.55,3.55)
    ax.set_yticks([3,2,1,0],METHODS if show_labels else ['']*4)
    ax.tick_params(axis='y',length=0,pad=4)
    ax.spines['left'].set_visible(False)
    if show_labels:
        for tick,method in zip(ax.get_yticklabels(),METHODS):
            tick.set_color(COLORS[method])
            tick.set_fontweight('bold' if method=='Full' else 'normal')

text_mm(4,180,'Numerical prediction',fontsize=10,fontweight='bold',va='bottom')
text_mm(4,175,'Lower MAE is better',fontsize=9,va='bottom',color='#555555')
text_mm(179,180,'3,063 out-of-fold records',fontsize=9,ha='right',va='bottom',color='#555555')
means=[]
for j,task in enumerate(TASKS[:2]):
    x=[22,112][j]
    ax=axes_mm(x,135,59,29)
    model_rows(ax)
    heading(x,169,'ab'[j],TITLES[j])
    max_x=[1.08,.29][j]
    for k,method in enumerate(METHODS):
        y=3-k
        val=meta[task][method]['macro_mae']
        ax.hlines(y,0,val,color=COLORS[method],lw=1.2,alpha=.5)
        ax.plot(val,y,MARKERS[method],ms=5.5,color=COLORS[method],zorder=3)
        ax.text(val+max_x*.025,y,f'{val:.3f}',color=COLORS[method],va='center',fontsize=9)
        means.append(dict(task=task,method=method,MAE=val,spearman=meta[task][method]['mean_spearman']))
    ax.set_xlim(0,max_x)
    ax.set_xticks([0,.5,1] if j==0 else [0,.1,.2])
    ax.set_xlabel('Mean absolute error (Å)' if j==0 else 'Mean absolute error (nats)',labelpad=3)
    reduction=100*(1-meta[task]['Full']['macro_mae']/meta[task]['Static']['macro_mae'])
    text_mm(x,164.5,f'Full: {reduction:.1f}% lower MAE than Static',color=BLUE,fontsize=9,va='bottom')

text_mm(4,118,'High-dynamic site identification',fontsize=10,fontweight='bold',va='bottom')
text_mm(179,118,'Mean AP and 95% CI',fontsize=9,ha='right',va='bottom',color='#555555')
ap_data=[]
for j,task in enumerate(TASKS):
    x=[22,80,138][j]
    ax=axes_mm(x,77,33,29)
    model_rows(ax,show_labels=j==0)
    heading(x,109,'cde'[j],TITLES[j])
    for k,method in enumerate(METHODS):
        y=3-k
        m=meta[task][method]
        val=m['mean_ap'];lo,hi=m['ap_ci']
        assert 0<=lo<=val<=hi<=1
        ax.plot(val,y,MARKERS[method],ms=3.5,mfc='white',mec=COLORS[method],mew=.9,zorder=3)
        ax.errorbar(val,y,xerr=[[val-lo],[hi-val]],fmt='none',
                    color=COLORS[method],elinewidth=.85,capsize=2.3,capthick=.85,zorder=4)
        ax.text(1.045,y,f'{val:.3f}',color=COLORS[method],va='center',clip_on=False,fontsize=9)
        ap_data.append(dict(task=task,method=method,mean_AP=val,CI_low=lo,CI_high=hi))
    ax.set_xlim(0,1)
    ax.set_xticks([0,.5,1])
    ax.set_xlabel('Average precision',labelpad=3)

text_mm(4,61,'Paired gains across individual records',fontsize=10,fontweight='bold',va='bottom')
text_mm(179,61,'Right of zero: Full performs better',fontsize=9,ha='right',va='bottom',color='#555555')
gain_data=[]
gain_summary=[]
for j,task in enumerate(TASKS):
    x=[22,80,138][j]
    ax=axes_mm(x,26,38,22)
    heading(x,51,'fgh'[j],TITLES[j])
    paired=frames[task].pivot(index='record_code',columns='method',values='ap').sort_index()
    ax.axvline(0,color='#666666',lw=.7,ls=(0,(3,3)),zorder=0)
    for k,comp in enumerate(['Intra','Static']):
        y=1-k
        values=100*(paired.Full-paired[comp]).to_numpy()
        assert len(values)==3063 and np.isfinite(values).all() and np.abs(values).max()<80
        # Deterministic vertical offsets show every record; x remains the exact paired AP gain.
        offsets=.16*np.sin(paired.index.to_numpy()*2.399963229728653)
        ax.scatter(values,y+offsets,s=.65,alpha=.13,color=COLORS[comp],linewidths=0,zorder=1)
        q=np.quantile(values,[.05,.25,.5,.75,.95])
        ax.hlines(y,q[0],q[4],color=COLORS[comp],lw=.9,zorder=2)
        ax.hlines(y,q[1],q[3],color=COLORS[comp],lw=4,zorder=3)
        ax.vlines(q[2],y-.12,y+.12,color='white',lw=.9,zorder=4)
        ax.plot(values.mean(),y,'o',ms=3.5,mec=COLORS[comp],mfc='white',mew=1,zorder=5)
        wins=int((values>0).sum());losses=int((values<0).sum());ties=int((values==0).sum())
        gain_summary.append(dict(task=task,comparison='Full - '+comp,n_records=3063,
                                 mean_gain_percentage_points=float(values.mean()),
                                 median_gain_percentage_points=float(q[2]),
                                 positive_records=wins,negative_records=losses,ties=ties,
                                 positive_fraction=wins/3063,min_gain=float(values.min()),max_gain=float(values.max())))
        text_mm(x-9,12.5-k*5.5,f'{comp}: +{values.mean():.2f} pp · {100*wins/3063:.1f}% improved',
                color=COLORS[comp],fontsize=9,ha='left',va='bottom')
        gain_data.extend(dict(task=task,record_code=int(code),comparison='Full - '+comp,
                              AP_Full=float(full),AP_comparator=float(baseline),gain_percentage_points=float(gain))
                         for code,full,baseline,gain in zip(paired.index,paired.Full,paired[comp],values))
    ax.set(xlim=(-80,80),ylim=(-.5,1.5),xlabel='ΔAP (percentage points)')
    ax.set_xticks([-80,0,80])
    ax.set_yticks([1,0],['vs Intra','vs Static'] if j==0 else ['',''])
    ax.tick_params(axis='y',length=0,pad=4)
    ax.spines['left'].set_visible(False)
    ax.xaxis.labelpad=3

pd.DataFrame(means).to_csv(ROOT/'numerical_performance.csv',index=False)
pd.DataFrame(ap_data).to_csv(ROOT/'absolute_AP_and_frozen_CI.csv',index=False)
pd.DataFrame(gain_data).to_csv(ROOT/'all_paired_record_AP_gains.csv',index=False)
pd.DataFrame(gain_summary).to_csv(ROOT/'paired_gain_summary.csv',index=False)
fig.canvas.draw()
renderer=fig.canvas.get_renderer()
texts=[t for t in fig.findobj(mpl.text.Text) if t.get_visible() and t.get_text()]
outside=[]
for t in texts:
    bb=t.get_window_extent(renderer)
    if bb.x0<-.5 or bb.y0<-.5 or bb.x1>fig.bbox.width+.5 or bb.y1>fig.bbox.height+.5:
        outside.append(t.get_text())
assert not outside,outside
assert min(t.get_fontsize() for t in texts)>=9
assert not any('\u4e00'<=c<='\u9fff' for t in texts for c in t.get_text())
for ext in ['pdf','svg','png']:
    fig.savefig(ROOT/f'fig2_performance_candidate.{ext}',bbox_inches=fig.bbox_inches,dpi=450)
svg=ET.parse(ROOT/'fig2_performance_candidate.svg')
assert not svg.findall('.//{http://www.w3.org/2000/svg}image')
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
audit=dict(status='numeric_and_canvas_passed_visual_pending',width_mm=W,height_mm=H,dpi=450,
           min_font_pt=min(t.get_fontsize() for t in texts),all_3063_records_paired=True,
           plotted_paired_gain_points=len(gain_data),data_subsampling=False,smoothing=False,
           AP_CI='frozen 10000 source-stratified record bootstrap, seed 20261018',
           new_significance_tests=False,new_MAE_or_gain_CI=False,common_AP_axis=[0,1],
           common_gain_axis_percentage_points=[-80,80],all_gain_values_shown=True,
           jitter='vertical only; deterministic, does not change any AP value',
           source_sha256={p.name:sha(p) for p in ARCHIVE.iterdir() if p.is_file()},
           output_sha256={ext:sha(ROOT/f'fig2_performance_candidate.{ext}') for ext in ['pdf','svg','png']},
           outside_canvas_text=outside,main_manuscript='v61 preserved')
(ROOT/'data_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
(ROOT/'visible_text.json').write_text(json.dumps([t.get_text() for t in texts],ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(audit,ensure_ascii=False,indent=2))
print(pd.DataFrame(gain_summary).to_string(index=False))
