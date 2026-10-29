# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a, dual x/y axes overlay) record-weighted complete density -> inherited_v61_plot.py / KernelDensity/plot_CorrelationDensity1.R -> param inherit
# (b,c,d) all-record AP violin -> GroupedViolin/plot_GroupedViolin.py -> param inherit
# (e,f,g) complete macro precision-recall -> inherited_v61_plot.py / AUROC/plot_AUROC2.py -> param inherit
# Companion source AP -> inherited_v61_plot.py / heatmap/plot_composition.py -> param inherit
# RULE: template unweighted KDE / ROC statistics are incompatible with frozen record-weighted density / PR.
# Selected C layout: both outputs overlaid in one plotting area, with independent x/y axes and original units.

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
import matplotlib.patheffects as path_effects
from matplotlib.transforms import Bbox
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

import fitz
W,H=183.,203.
all_metrics={**summary["metrics"],"within_type_hchi":summary["within_type"]}
frames={t:record_metrics[record_metrics.target.eq(t)] for t in TASKS[:2]}
frames["within_type_hchi"]=within_record

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


def number_format(ax,axis):
    axis_obj=getattr(ax,axis+'axis')
    axis_obj.set_major_formatter(FuncFormatter(lambda value,pos:format(value,'g')))
    axis_obj.set_minor_locator(NullLocator())

def canvas(height):
    return plt.figure(figsize=(W/25.4,height/25.4)),height

def ax_mm(fig,height,x,y,w,h):return fig.add_axes([x/W,y/height,w/W,h/height])
def txt(fig,height,x,y,s,**kw):return fig.text(x/W,y/height,s,**kw)
def title(fig,height,x,y,letter,label):
    txt(fig,height,x-11,y,letter,fontsize=11,fontweight='bold',va='bottom')
    txt(fig,height,x,y,label,fontsize=10,va='bottom')

def shared_key(fig,height):
    # The method-color key applies only to b--g; a colors encode outputs.
    txt(fig,height,15,119,'b–g:',fontsize=9,ha='left',va='center',color='#444444')
    fig.legend(bbox_to_anchor=(99/W,119/height),loc='center',
               handles=[Line2D([],[],color=COLORS[m],ls=STYLES[m],lw=1.4,label=m) for m in METHODS],
               ncol=4,columnspacing=1,handlelength=2,handletextpad=.4,borderaxespad=0)


def translucent_mass_cmap(color):
    # Hue identifies output; opacity shows the same frozen log mass for each output.
    rgb=mpl.colors.to_rgb(color)
    rgba=[(*rgb,.015+.335*t**1.15) for t in np.linspace(0,1,256)]
    cmap=mpl.colors.ListedColormap(rgba)
    cmap.set_bad((1,1,1,0))
    return cmap


def overlay_density(fig,height):
    axes={}
    colors=dict(rmsf=BLUE,chi_entropy='#9B5C1A')
    axes_position=(63,134,55,55)
    metadata={}
    for j,task in enumerate(TASKS[:2]):
        ax=ax_mm(fig,height,*axes_position)
        axes[task]=ax
        ax.set_zorder(j+1)
        if j:ax.patch.set_visible(False)
        data=points[task]
        lo=min(float(data['reference'].min()),float(data['prediction'].min()),float(data['x_edges'][0]),float(data['y_edges'][0]))
        hi=max(float(data['reference'].max()),float(data['prediction'].max()),float(data['x_edges'][-1]),float(data['y_edges'][-1]))
        cmap=translucent_mass_cmap(colors[task])
        ax.pcolormesh(data['x_edges'],data['y_edges'],np.ma.masked_equal(data['histogram_mass'].T,0),
                      cmap=cmap,norm=density_norm,shading='flat',linewidth=0,rasterized=False,zorder=1)
        ax.set(xlim=(lo,hi),ylim=(lo,hi))
        ax.set_aspect('equal',adjustable='box')
        if task=='rmsf':
            ax.set_xscale('symlog',linthresh=1);ax.set_yscale('symlog',linthresh=1)
            ax.set_xticks([0,1,10,30]);ax.set_yticks([0,1,10,30])
            number_format(ax,'x');number_format(ax,'y')
            ax.set(xlabel='MD RMSF (Å)',ylabel='Predicted RMSF (Å)')
            visible=('bottom','left')
        else:
            ax.set_xticks([0,1,2,3,4]);ax.set_yticks([0,1,2,3,4])
            ax.xaxis.set_ticks_position('top');ax.xaxis.set_label_position('top')
            ax.yaxis.set_ticks_position('right');ax.yaxis.set_label_position('right')
            ax.set(xlabel='MD Hχ (nats)',ylabel='Predicted Hχ (nats)')
            visible=('top','right')
        for name,spine in ax.spines.items():
            spine.set_visible(name in visible);spine.set_color(colors[task]);spine.set_linewidth(.6)
        ax.tick_params(axis='both',colors=colors[task],pad=2,length=2.5)
        ax.xaxis.label.set_color(colors[task]);ax.yaxis.label.set_color(colors[task])
        ax.xaxis.labelpad=3;ax.yaxis.labelpad=3
        metric=summary['metrics'][task]['Full']
        units='Å' if task=='rmsf' else 'nats'
        label='RMSF' if task=='rmsf' else 'Hχ'
        # Labels occupy the empty upper-left / lower-right corners inside the shared plot.
        position=(.035,.965) if task=='rmsf' else (.965,.035)
        alignment=dict(ha='left',va='top') if task=='rmsf' else dict(ha='right',va='bottom')
        ax.text(*position,f'{label}\nMAE {metric["macro_mae"]:.3f} {units}\nρ {metric["mean_spearman"]:.3f}',
                transform=ax.transAxes,fontsize=9,color=colors[task],linespacing=1.2,zorder=8,
                bbox=dict(facecolor='white',edgecolor='none',alpha=.88,pad=1.5),**alignment)
        # Vertical bars sit immediately beyond the corresponding y-axis title.
        cb_x=44 if task=='rmsf' else 132.5
        cbax=ax_mm(fig,height,cb_x,148,1.6,27)
        cb=fig.colorbar(plt.cm.ScalarMappable(norm=density_norm,cmap=cmap),cax=cbax,orientation='vertical')
        cb.solids.set_rasterized(False);cb.outline.set_visible(False)
        cb.set_ticks([1e-5,1e-2]);cb.ax.yaxis.set_minor_locator(NullLocator())
        if task=='rmsf':cb.ax.yaxis.set_ticks_position('left')
        cb.ax.tick_params(labelsize=9,pad=1.5,length=1.5,color=colors[task],labelcolor=colors[task])
        cb.ax.set_title('Mass',fontsize=9,color=colors[task],pad=4)
        metadata[task]=dict(limits=[lo,hi],scale='symlog, linthresh=1' if j==0 else 'linear',
                            n_positions=len(data['reference']),mass_sum=float(data['histogram_mass'].sum()),
                            color=colors[task],x_axis='bottom' if j==0 else 'top',y_axis='left' if j==0 else 'right')
    front=axes['chi_entropy']
    rmsf_text=axes['rmsf'].texts[0]
    annotation_text=rmsf_text.get_text()
    rmsf_text.remove()
    front.text(.035,.965,annotation_text,transform=front.transAxes,ha='left',va='top',
               fontsize=9,color=BLUE,linespacing=1.2,zorder=8,
               bbox=dict(facecolor='white',edgecolor='none',alpha=.88,pad=1.5))
    # Both symmetric coordinate transforms map y=x to the same geometric diagonal.
    front.plot([0,1],[0,1],transform=front.transAxes,color='#999999',lw=.55,ls=(0,(1.3,3.2)),zorder=3)
    # Draw both original trends in the front axes, each with its own data transform.
    # This prevents either density layer from painting over the other trend.
    for task,style in [('rmsf','-'),('chi_entropy','--')]:
        tr=trend[trend.target.eq(task)&trend.method.eq('Full')].sort_values('bin_index')
        transform=axes[task].transData
        front.fill_between(tr.x_mean,tr.ci_low,tr.ci_high,transform=transform,
                           color=colors[task],alpha=.10,lw=0,zorder=4)
        line,=front.plot(tr.x_mean,tr.prediction_mean,transform=transform,color=colors[task],
                         ls=style if task=='rmsf' else (0,(4,2.2)),lw=1.35 if task=='rmsf' else 1.9,zorder=5)
        if task=='chi_entropy':
            line.set_path_effects([path_effects.Stroke(linewidth=3.1,foreground='white',alpha=.9),path_effects.Normal()])
    # No data-limit autoscaling after cross-transform artists are added.
    for task,ax in axes.items():ax.set(xlim=metadata[task]['limits'],ylim=metadata[task]['limits'])
    fig.canvas.draw()
    assert np.allclose(axes['rmsf'].bbox.bounds,axes['chi_entropy'].bbox.bounds)
    for task,ax in axes.items():
        data=points[task]
        assert np.isclose(data['histogram_mass'].sum(),1)
        lo,hi=metadata[task]['limits']
        for value in [lo,(lo+hi)/2,hi]:
            px,py=ax.transData.transform((value,value))
            assert np.isclose((px-ax.bbox.x0)/ax.bbox.width,(py-ax.bbox.y0)/ax.bbox.height)
    return axes,metadata


def precision_recall(fig,height,x,y,w,h,task,letter,heading_y,ap_y,ylabel=True):
    ax=ax_mm(fig,height,x,y,w,h)
    frame=within_pr if task=='within_type_hchi' else pr[pr.target.eq(task)]
    ends=within_end if task=='within_type_hchi' else endpoints[endpoints.target.eq(task)]
    metrics=all_metrics[task]
    for method in METHODS:
        cur=frame[frame.method.eq(method)].sort_values('recall')
        end=ends[ends.method.eq(method)].iloc[0]
        assert len(cur)==101 and end.recall==1
        xx=np.r_[cur.recall.to_numpy(),end.recall];yy=np.r_[cur.precision.to_numpy(),end.precision]
        ax.step(xx,yy,where='post',color=COLORS[method],ls=STYLES[method],lw=WIDTHS[method],zorder=ZORDER[method])
    baseline=float(frame.no_skill_precision.iloc[0])
    ax.axhline(baseline,color='#AAAAAA',ls=(0,(2,2)),lw=.65,zorder=0)
    ax.set(xlim=(0,1),ylim=(0,1.02),xlabel='Recall',ylabel='Precision' if ylabel else '')
    ax.set_xticks([0,.5,1]);ax.set_yticks([0,.5,1])
    if not ylabel:ax.set_yticklabels([])
    ax.tick_params(pad=2);ax.xaxis.labelpad=2.5
    label={'rmsf':'RMSF','chi_entropy':'Hχ','within_type_hchi':'Hχ within type'}[task]
    title(fig,height,x,heading_y,letter,label)
    txt(fig,height,x-1,ap_y,'AP',fontsize=9,ha='left',va='center',color='#444444')
    # Same method ordering as the single shared legend, regardless of panel width.
    for k,method in enumerate(METHODS):
        txt(fig,height,x+11+k*10.9,ap_y,f'{metrics[method]["mean_ap"]:.3f}',
            color=COLORS[method],fontsize=9,ha='center',va='center')


# C: all 3,063 records for each method and task are visible as points and distribution shapes.
distribution_points=[]
distribution_summary=[]
fig,h=canvas(203.)
shared_key(fig,h)
# Two outputs share exactly one visible plotting rectangle, with two independent x/y coordinate systems.
title(fig,h,15,198,'a','Full vs MD')
overlay_axes,overlay_info=overlay_density(fig,h)
panel_a_axes=list(fig.axes)
for j,task in enumerate(TASKS):
    x,width=[(15,43),(70,43),(125,51)][j]
    ax=ax_mm(fig,h,x,68,width,37)
    frame=frames[task]
    arrays=[frame[frame.method.eq(method)].sort_values('record_code').ap.to_numpy() for method in METHODS]
    assert all(len(v)==3063 for v in arrays)
    parts=ax.violinplot(arrays,positions=[0,1,2,3],widths=.78,showmeans=False,showmedians=False,
                        showextrema=False,points=128,bw_method=.20)
    for k,(method,values,body) in enumerate(zip(METHODS,arrays,parts['bodies'])):
        body.set_facecolor(COLORS[method]);body.set_edgecolor(COLORS[method]);body.set_linewidth(.8);body.set_alpha(.25)
        codes=frame[frame.method.eq(method)].sort_values('record_code').record_code.to_numpy()
        horizontal_offset=.13*np.sin(codes*2.399963229728653)
        ax.scatter(k+horizontal_offset,values,s=.45,color=COLORS[method],alpha=.09,linewidths=0,zorder=2,rasterized=False)
        quartiles=np.quantile(values,[.25,.5,.75])
        ax.vlines(k,quartiles[0],quartiles[2],color=COLORS[method],lw=2.5,zorder=3)
        ax.plot([k-.085,k+.085],[quartiles[1],quartiles[1]],color='white',lw=1,zorder=4)
        distribution_summary.append(dict(task=task,method=method,n_records=3063,mean=float(values.mean()),
                                         q25=float(quartiles[0]),median=float(quartiles[1]),q75=float(quartiles[2]),
                                         minimum=float(values.min()),maximum=float(values.max())))
        distribution_points.extend(dict(task=task,method=method,record_code=int(code),AP=float(value)) for code,value in zip(codes,values))
    ax.set(xlim=(-.5,3.5),ylim=(-.02,1.02),ylabel='Record-level AP' if j==0 else '')
    ax.set_xticks([0,1,2,3],METHODS);ax.set_yticks([0,.5,1]);ax.tick_params(pad=2)
    for tick,method in zip(ax.get_xticklabels(),METHODS):tick.set_color(COLORS[method])
    if j>0:ax.set_yticklabels([])
    title(fig,h,x,113,'bcd'[j],['RMSF','Hχ','Hχ within type'][j])
for j,task in enumerate(TASKS):
    x,width=[(15,43),(70,43),(125,51)][j]
    precision_recall(fig,h,x,14,width,32,task,'efg'[j],55,50,j==0)
selected_check=export(fig,'fig2_C_selected',W,h)
# A tight PDF crop removes the former side-column whitespace from the detail preview.
renderer=fig.canvas.get_renderer()
panel_bounds=Bbox.union([ax.get_tightbbox(renderer) for ax in panel_a_axes]).padded(5)
panel_bounds_inches=panel_bounds.transformed(fig.dpi_scale_trans.inverted())
panel_clip_points=[panel_bounds_inches.x0*72,(h/25.4-panel_bounds_inches.y1)*72,
                   panel_bounds_inches.x1*72,(h/25.4-panel_bounds_inches.y0)*72]
plt.close(fig)
pd.DataFrame(distribution_points).to_csv(ROOT/'all_record_AP_distribution_points.csv',index=False)
pd.DataFrame(distribution_summary).to_csv(ROOT/'AP_distribution_summary.csv',index=False)


doc=fitz.open(ROOT/'fig2_C_selected.pdf')
factor=200/72*167.64/W
doc[0].get_pixmap(matrix=fitz.Matrix(factor,factor),alpha=False).save(ROOT/'print_preview.png')
fonts=[f[3] for f in doc[0].get_fonts()]
assert all('Arial' in name for name in fonts)
assert abs(doc[0].rect.width*25.4/72-W)<.01
assert abs(doc[0].rect.height*25.4/72-H)<.01
doc[0].get_pixmap(matrix=fitz.Matrix(2.5,2.5),clip=fitz.Rect(*panel_clip_points),alpha=False).save(ROOT/'panel_a_preview.png')
doc.close()
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
audit=dict(status='numeric_and_canvas_passed_visual_pending',selected_option='C',
           panel_count=7,internal_density_axes=2,same_plot_rectangle=True,dual_x_and_y_axes=True,overlay_info=overlay_info,
           density_units_separate=True,new_density_normalization=False,
           colorbars_vertical=True,metrics_inside_plot=True,hchi_trend_color="#9B5C1A",hchi_trend_linewidth_pt=1.9,
           panel_a_preview_width_mm=panel_bounds_inches.width*25.4,
           density_log_mass_range=[density_norm.vmin,density_norm.vmax],overlay_geometric_checks_passed=True,
           selected_check=selected_check,all_3063_records_in_each_distribution=True,
           all_original_PR_vertices_retained=True,no_new_statistics=True,
           source_sha256={p.name:sha(p) for p in ARCHIVE.iterdir() if p.is_file()},
           output_sha256={p.name:sha(p) for p in ROOT.iterdir() if p.name.startswith(('fig2_C_selected','print_preview','panel_a_preview'))},
           manuscript_replaced=False)
(ROOT/'layout_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(audit,ensure_ascii=False,indent=2))
