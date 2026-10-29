# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a,b) record-weighted complete density -> inherited_v61_plot.py / KernelDensity/plot_CorrelationDensity1.R -> param inherit
# (c,d) per-record paired AP scatter -> RMSF-net Fig2c/e composition / cross-type scatter -> param inherit
# (e,f,g) complete macro precision-recall -> inherited_v61_plot.py / AUROC/plot_AUROC2.py -> param inherit
# Companion source AP -> inherited_v61_plot.py / heatmap/plot_composition.py -> param inherit
# RULE: template unweighted KDE / ROC statistics are incompatible with frozen record-weighted density / PR.
# Keep complete data and native quantitative semantics; this is an independent trial, not a manuscript replacement.

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
from matplotlib.colors import LinearSegmentedColormap, LogNorm
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter, NullLocator
from publication_style import BLUE, WARM, TEAL, MD, unit_cmap, unit_norm, UNIT_TICKS

ROOT=Path(__file__).resolve().parent
ARCHIVE=ROOT/'source_archive'
METHODS=['Full','Intra','Static','AA']
TASKS=['rmsf','chi_entropy','within_type_hchi']
COLORS=dict(Full=BLUE,Intra=WARM,Static=TEAL,AA='#999999')
STYLES=dict(Full='-',Intra='--',Static='-.',AA=':')
WIDTHS=dict(Full=1.25,Intra=1.1,Static=1.1,AA=1.0)
ZORDER=dict(Full=5,Intra=6,Static=4,AA=3)
summary=json.loads((ARCHIVE/'summary.json').read_text(encoding='utf-8'))
pr=pd.read_csv(ARCHIVE/'macro_pr_curve.csv',float_precision='round_trip')
endpoints=pd.read_csv(ARCHIVE/'all_selected_endpoint.csv',float_precision='round_trip')
within_pr=pd.read_csv(ARCHIVE/'within_type_macro_pr.csv',float_precision='round_trip')
within_end=pd.read_csv(ARCHIVE/'within_type_all_selected_endpoint.csv',float_precision='round_trip')
trend=pd.read_csv(ARCHIVE/'four_method_binned_trend.csv',float_precision='round_trip')
source_ap=pd.read_csv(ARCHIVE/'source_ap_matrix.csv',float_precision='round_trip')
record_metrics=pd.read_csv(ARCHIVE/'per_record_metrics.csv',float_precision='round_trip')
within_record=pd.read_csv(ARCHIVE/'within_type_record_macro_metrics.csv',float_precision='round_trip')
points={t:np.load(ARCHIVE/f'prediction_reference_{t}.npz',allow_pickle=False) for t in TASKS[:2]}
assert len(record_metrics)==3063*2*4 and len(within_record)==3063*4
assert not record_metrics.duplicated(['record_code','target','method']).any()
assert not within_record.duplicated(['record_code','method']).any()
assert record_metrics.sample_id.nunique()==within_record.sample_id.nunique()==3063
assert np.isfinite(record_metrics[['ap','mae','spearman']].to_numpy()).all()
assert np.isfinite(within_record[['ap','auroc','prevalence']].to_numpy()).all()
nonzero=np.concatenate([points[t]['histogram_mass'][points[t]['histogram_mass']>0] for t in TASKS[:2]])
density_cmap=LinearSegmentedColormap.from_list('record_mass',['#EEEEEE','#C6C6C6','#5A5A5A'])
density_cmap.set_bad('white')
density_norm=LogNorm(vmin=nonzero.min(),vmax=nonzero.max())
mpl.rcParams.update({'font.size':9,'axes.labelsize':9,'axes.titlesize':10,
                     'xtick.labelsize':9,'ytick.labelsize':9,'legend.fontsize':9,'axes.unicode_minus':False})
W,H=183.,203.
fig=plt.figure(figsize=(W/25.4,H/25.4))

def axes_mm(canvas,x,y,w,h,height=H):
    return canvas.add_axes([x/W,y/height,w/W,h/height])

def text_mm(canvas,x,y,s,height=H,**kwargs):
    return canvas.text(x/W,y/height,s,**kwargs)

def heading(x,y,letter,title):
    text_mm(fig,x-12,y,letter,fontsize=11,fontweight='bold',va='bottom')
    text_mm(fig,x,y,title,fontsize=10,va='bottom')

def number_format(ax,axis):
    axis_obj=getattr(ax,axis+'axis')
    axis_obj.set_major_formatter(FuncFormatter(lambda value,pos:format(value,'g')))
    axis_obj.set_minor_locator(NullLocator())

panels={}
for j,task in enumerate(TASKS[:2]):
    x=[20,111][j]
    ax=axes_mm(fig,x,133,55,55)
    data=points[task]
    m=summary['metrics'][task]['Full']
    ref,pred=data['reference'],data['prediction']
    assert len(ref)==len(pred)==m['n_positions']
    assert np.isfinite(ref).all() and np.isfinite(pred).all()
    assert np.isclose(data['histogram_mass'].sum(),1)
    mesh=ax.pcolormesh(data['x_edges'],data['y_edges'],np.ma.masked_equal(data['histogram_mass'].T,0),
                       cmap=density_cmap,norm=density_norm,shading='flat',linewidth=0,rasterized=False)
    lo=min(float(ref.min()),float(pred.min()),float(data['x_edges'][0]),float(data['y_edges'][0]))
    hi=max(float(ref.max()),float(pred.max()),float(data['x_edges'][-1]),float(data['y_edges'][-1]))
    ax.plot([lo,hi],[lo,hi],color=MD,lw=.75,ls=(0,(3,3)))
    full_trend=trend[trend.target.eq(task)&trend.method.eq('Full')].sort_values('bin_index')
    assert len(full_trend)==20
    populated=full_trend.n_positions.gt(0)
    assert np.isfinite(full_trend.loc[populated,['prediction_mean','ci_low','ci_high']].to_numpy()).all()
    ax.fill_between(full_trend.x_mean,full_trend.ci_low,full_trend.ci_high,color=BLUE,alpha=.15,lw=0,zorder=3)
    ax.plot(full_trend.x_mean,full_trend.prediction_mean,color=BLUE,lw=1.25,zorder=4)
    units='Å' if j==0 else 'nats'
    ax.set(xlim=(lo,hi),ylim=(lo,hi),xlabel=f'MD reference ({units})',ylabel=f'Full prediction ({units})')
    ax.set_aspect('equal',adjustable='box')
    if j==0:
        ax.set_xscale('symlog',linthresh=1.0);ax.set_yscale('symlog',linthresh=1.0)
        ax.set_xticks([0,1,10,30]);ax.set_yticks([0,1,10,30])
        number_format(ax,'x');number_format(ax,'y')
    else:
        ax.set_xticks([0,1,2,3,4]);ax.set_yticks([0,1,2,3,4])
    ax.tick_params(pad=2)
    ax.xaxis.labelpad=3;ax.yaxis.labelpad=3
    heading(x,192,'ab'[j],['RMSF','Hχ'][j])
    ax.text(.035,.97,f'MAE {m["macro_mae"]:.3f} {units}\nρ {m["mean_spearman"]:.3f}',
            transform=ax.transAxes,ha='left',va='top',fontsize=9,linespacing=1.35,
            bbox=dict(facecolor='white',edgecolor='none',alpha=.9,pad=1.2))
    panels['ab'[j]]=dict(task=task,n_records=3063,n_positions=len(ref),density_mass=float(data['histogram_mass'].sum()),
                         limits=[lo,hi],trend_methods=['Full'],MAE=m['macro_mae'],spearman=m['mean_spearman'])

# One shared method key; only the Full trend appears in a/b, all four methods appear in e--g.
fig.legend(handles=[Line2D([],[],color=COLORS[m],ls=STYLES[m],lw=1.4,label=m) for m in METHODS],
           loc='center',bbox_to_anchor=(70/W,199/H),ncol=4,columnspacing=1.05,
           handlelength=2,handletextpad=.4,borderaxespad=0)
text_mm(fig,126,199,'Mass',ha='right',va='center',fontsize=9,color='#555555')
cbax=axes_mm(fig,131,198,39,1.6)
cb=fig.colorbar(mesh,cax=cbax,orientation='horizontal')
cb.solids.set_rasterized(False);cb.outline.set_visible(False)
cb.ax.xaxis.set_ticks_position('bottom');cb.set_ticks([1e-5,1e-3,1e-2])
cb.ax.xaxis.set_minor_locator(NullLocator());cb.ax.tick_params(labelsize=9,pad=.5,length=1.5)

paired_rows=[]
paired_summary=[]
for j,task in enumerate(['rmsf','within_type_hchi']):
    x=[22,113][j]
    ax=axes_mm(fig,x,73,40,40)
    frame=record_metrics[record_metrics.target.eq(task)] if task=='rmsf' else within_record
    paired=frame.pivot(index='record_code',columns='method',values='ap').sort_index()
    assert len(paired)==3063 and paired.notna().all().all()
    xx,yy=paired.Static.to_numpy(),paired.Full.to_numpy()
    assert ((xx>=0)&(xx<=1)&(yy>=0)&(yy<=1)).all()
    ax.plot([0,1],[0,1],color=MD,lw=.7,ls=(0,(3,3)),zorder=0)
    ax.scatter(xx,yy,s=2.0,color=BLUE,alpha=.25,linewidths=0,rasterized=False,zorder=1)
    ax.set(xlim=(0,1),ylim=(0,1),xlabel='Static AP',ylabel='Full AP')
    ax.xaxis.label.set_color(TEAL);ax.yaxis.label.set_color(BLUE)
    ax.set_aspect('equal',adjustable='box')
    ax.set_xticks([0,.5,1]);ax.set_yticks([0,.5,1]);number_format(ax,'x');number_format(ax,'y')
    ax.tick_params(pad=2);ax.xaxis.labelpad=3;ax.yaxis.labelpad=3
    heading(x-2,117.5,'cd'[j],['RMSF: record-level AP','Hχ within type: record-level AP'][j])
    wins=int((yy>xx).sum());ties=int((yy==xx).sum());losses=int((yy<xx).sum())
    text_mm(fig,x+45,96,'Full higher in',color=BLUE,fontsize=9,ha='left',va='bottom')
    text_mm(fig,x+45,91.4,f'{100*wins/3063:.1f}% of records',color=BLUE,fontsize=9,ha='left',va='bottom')
    paired_summary.append(dict(task=task,n_records=3063,Full_higher=wins,equal=ties,Full_lower=losses,
                               Full_higher_fraction=wins/3063,mean_AP_difference=float((yy-xx).mean())))
    paired_rows.extend(dict(task=task,record_code=int(code),Static_AP=float(a),Full_AP=float(b),difference=float(b-a))
                       for code,a,b in zip(paired.index,xx,yy))
    panels['cd'[j]]=paired_summary[-1]

vertices=[]
for j,task in enumerate(TASKS):
    x,width=[(15,43),(70,43),(125,51)][j]
    ax=axes_mm(fig,x,14,width,35)
    frame=within_pr if j==2 else pr[pr.target.eq(task)]
    ends=within_end if j==2 else endpoints[endpoints.target.eq(task)]
    metrics=summary['within_type'] if j==2 else summary['metrics'][task]
    for method in METHODS:
        current=frame[frame.method.eq(method)].sort_values('recall')
        end=ends[ends.method.eq(method)].iloc[0]
        assert len(current)==101 and end.recall==1
        expected=metrics[method]['mean_record_type_prevalence'] if j==2 else metrics[method]['prevalence_mean']
        assert np.isclose(end.precision,expected)
        xx=np.r_[current.recall.to_numpy(),end.recall]
        yy=np.r_[current.precision.to_numpy(),end.precision]
        ax.step(xx,yy,where='post',color=COLORS[method],ls=STYLES[method],lw=WIDTHS[method],zorder=ZORDER[method])
        vertices.extend(dict(task=task,method=method,recall=float(r),precision=float(p),final_endpoint=i==101)
                        for i,(r,p) in enumerate(zip(xx,yy)))
    prevalence=float(frame.no_skill_precision.iloc[0])
    assert np.isclose(prevalence,expected)
    ax.axhline(prevalence,color='#AAAAAA',ls=(0,(2,2)),lw=.65,zorder=0)
    ax.set(xlim=(0,1),ylim=(0,1.02),xlabel='Recall',ylabel='Precision' if j==0 else '')
    ax.set_xticks([0,.5,1]);ax.set_yticks([0,.5,1])
    if j>0:ax.set_yticklabels([])
    ax.tick_params(pad=2);ax.xaxis.labelpad=2.5
    heading(x,57.5,'efg'[j],['RMSF','Hχ','Hχ within type'][j])
    text_mm(fig,x-1,52.9,'AP',fontsize=9,ha='left',va='center',color='#444444')
    for k,method in enumerate(METHODS):
        text_mm(fig,x+11+k*10.9,52.9,f'{metrics[method]["mean_ap"]:.3f}',
                color=COLORS[method],fontsize=9,ha='center',va='center')
    panels['efg'[j]]=dict(task=task,n_records=3063,methods=METHODS,vertices_per_method=102,
                         all_selected_precision=prevalence,mean_AP={m:metrics[m]['mean_ap'] for m in METHODS})

pd.DataFrame(paired_rows).to_csv(ROOT/'all_paired_AP_points.csv',index=False)
pd.DataFrame(paired_summary).to_csv(ROOT/'paired_AP_summary.csv',index=False)
pd.DataFrame(vertices).to_csv(ROOT/'complete_PR_vertices.csv',index=False)

def export(canvas,stem,width,height):
    canvas.canvas.draw()
    renderer=canvas.canvas.get_renderer()
    ignored=set()
    for ax in canvas.axes:
        for axis in [ax.xaxis,ax.yaxis]:
            lo,hi=sorted(axis.get_view_interval())
            for tick in axis.get_major_ticks()+axis.get_minor_ticks():
                if tick.get_loc()<lo-1e-10 or tick.get_loc()>hi+1e-10:ignored.update([tick.label1,tick.label2])
    texts=[t for t in canvas.findobj(mpl.text.Text) if t.get_visible() and t.get_text() and t not in ignored]
    outside=[]
    for t in texts:
        box=t.get_window_extent(renderer)
        if box.x0<-.5 or box.y0<-.5 or box.x1>canvas.bbox.width+.5 or box.y1>canvas.bbox.height+.5:
            outside.append(t.get_text())
    assert not outside,outside
    assert min(t.get_fontsize() for t in texts)>=9
    assert not any('\u4e00'<=c<='\u9fff' for t in texts for c in t.get_text())
    for ext in ['png','pdf','svg']:canvas.savefig(ROOT/f'{stem}.{ext}',bbox_inches=canvas.bbox_inches,dpi=450)
    assert not ET.parse(ROOT/f'{stem}.svg').findall('.//{http://www.w3.org/2000/svg}image')
    return dict(width_mm=width,height_mm=height,minimum_font_pt=min(t.get_fontsize() for t in texts),outside_text=outside)

main_check=export(fig,'fig2_journal_trial',W,H)

# Preserve every source/task cell as a separate companion figure, with the original 0--1 blue scale.
sh=45.
source_fig=plt.figure(figsize=(W/25.4,sh/25.4))
ax=axes_mm(source_fig,24,4.5,145,20,height=sh)
matrix=np.empty((4,9))
sources=['atlas','dynarepo','itrace']
for i,method in enumerate(METHODS):
    for j,task in enumerate(TASKS):
        for k,source in enumerate(sources):
            cell=source_ap[source_ap.task.eq(task)&source_ap.dataset.eq(source)&source_ap.method.eq(method)]
            assert len(cell)==1
            matrix[i,3*j+k]=float(cell.mean_ap.iloc[0])
im=ax.pcolormesh(np.arange(10),np.arange(5),matrix,cmap=unit_cmap(),norm=unit_norm(),
                 linewidth=.8,edgecolor='white',rasterized=False)
ax.set(xlim=(0,9),ylim=(4,0));ax.set_xticks([]);ax.set_yticks(np.arange(4)+.5,METHODS)
ax.tick_params(axis='y',length=0,pad=5)
for label,method in zip(ax.get_yticklabels(),METHODS):label.set_color(COLORS[method])
for spine in ax.spines.values():spine.set_visible(False)
for i in range(4):
    for j in range(9):
        ax.text(j+.5,i+.5,f'{matrix[i,j]:.3f}',fontsize=9,ha='center',va='center',
                color='white' if matrix[i,j]>=.63 else '#222222')
for xx in [3,6]:ax.axvline(xx,color='white',lw=2)
text_mm(source_fig,5,41,'Mean AP by source',height=sh,fontsize=10,va='bottom')
for j,title in enumerate(['RMSF','Hχ','Hχ within type']):
    text_mm(source_fig,24+(3*j+1.5)*145/9,35.5,title,height=sh,fontsize=10,ha='center',va='bottom')
    for k,name in enumerate(['ATLAS','DynaRepo','ITrace']):
        count=source_ap[source_ap.task.eq(TASKS[j])&source_ap.dataset.eq(sources[k])].n_records.unique()
        assert len(count)==1
        text_mm(source_fig,24+(3*j+k+.5)*145/9,25.7,f'{name}\nn = {int(count[0]):,}',height=sh,
                fontsize=9,ha='center',va='bottom',linespacing=1.15)
cbax=axes_mm(source_fig,171,4.5,1.8,20,height=sh)
cb=source_fig.colorbar(im,cax=cbax,format='%g')
cb.solids.set_rasterized(False);cb.outline.set_visible(False)
cb.set_ticks(UNIT_TICKS);cb.ax.tick_params(length=1.5,pad=2,labelsize=9)
source_check=export(source_fig,'fig2_source_AP_companion',W,sh)

sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
audit=dict(status='numerical_and_canvas_passed_visual_pending',main=main_check,source_companion=source_check,
           n_records=3063,panel_info=panels,all_paired_records_plotted=True,paired_points=6126,
           common_AP_pair_limits=[0,1],all_PR_curves_complete=True,no_subsampling=True,no_smoothing=True,
           new_significance_tests=False,source_AP_cells=36,source_AP_limits=[0,1],
           source_sha256={p.name:sha(p) for p in ARCHIVE.iterdir() if p.is_file()},
           output_sha256={p.name:sha(p) for p in ROOT.iterdir() if p.suffix in ['.pdf','.svg','.png']},
           manuscript_replaced=False)
(ROOT/'data_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(audit,ensure_ascii=False,indent=2))
