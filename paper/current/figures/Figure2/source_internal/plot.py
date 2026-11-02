# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a, dual x/y axes overlay) record-weighted complete density -> inherited_v61_plot.py / KernelDensity/plot_CorrelationDensity1.R -> param inherit
# Task panels / merged distribution: all-record AP violin -> GroupedViolin/plot_GroupedViolin.py -> param inherit
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
from matplotlib.path import Path as VectorPath
from matplotlib.patches import PathPatch
from shapely.geometry import LineString
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
FROZEN_AP=pd.read_csv(ROOT/'all_record_AP_distribution_points.csv',float_precision='round_trip')
FROZEN_PR=pd.read_csv(ROOT/'complete_PR_vertices.csv',float_precision='round_trip')
RENDER_COUNTS={}
ACTIVE_OPTION=''

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


def blend_stroke_intersections(fig,front,trend_lines,colors,width_pt):
    """Give both colors equal weight only where their actual stroke areas overlap."""
    fig.canvas.draw()
    stroke_shapes=[]
    radius_mm=width_pt*25.4/72/2
    for task in TASKS[:2]:
        xy=trend_lines[task].get_transform().transform(trend_lines[task].get_xydata())
        xy_mm=xy/fig.dpi*25.4
        segments=[]
        start=0
        finite=np.isfinite(xy_mm).all(axis=1)
        for end in range(len(xy_mm)+1):
            if end==len(xy_mm) or not finite[end]:
                if end-start>=2:segments.append(LineString(xy_mm[start:end]))
                start=end+1
        assert segments
        shape=segments[0].buffer(radius_mm,cap_style=2,join_style=2)
        for segment in segments[1:]:shape=shape.union(segment.buffer(radius_mm,cap_style=2,join_style=2))
        stroke_shapes.append(shape)
    overlap=stroke_shapes[0].intersection(stroke_shapes[1])
    color=np.mean([mpl.colors.to_rgb(colors[t]) for t in TASKS[:2]],axis=0)
    def polygons(geometry):
        if geometry.is_empty:return
        if geometry.geom_type=='Polygon':yield geometry
        elif hasattr(geometry,'geoms'):
            for part in geometry.geoms:yield from polygons(part)
    count=0
    for polygon in polygons(overlap):
        vertices=[];codes=[]
        for ring in [polygon.exterior,*polygon.interiors]:
            xy=np.asarray(ring.coords)/np.array([W,fig.get_figheight()*25.4])
            vertices.extend(xy)
            codes.extend([VectorPath.MOVETO]+[VectorPath.LINETO]*(len(xy)-2)+[VectorPath.CLOSEPOLY])
        patch=PathPatch(VectorPath(vertices,codes),transform=fig.transFigure,
                        facecolor=color,edgecolor='none',linewidth=0,zorder=6)
        patch.set_clip_path(front.patch)
        front.add_patch(patch)
        count+=1
    return dict(method='equal-color mixture within exact buffered stroke overlap',
                output_linewidth_pt=width_pt,overlap_regions=count,overlap_area_mm2=float(overlap.area),
                mixture_rgb=color.tolist(),trend_vertices_unchanged=True,no_line_gaps=True)


def overlay_density(fig,height,x=63,y=134,size=55):
    axes={}
    colors=dict(rmsf=BLUE,chi_entropy='#9B5C1A')
    axes_position=(x,y,size,size)
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
                **alignment)
        metadata[task]=dict(limits=[lo,hi],scale='symlog, linthresh=1' if j==0 else 'linear',
                            n_positions=len(data['reference']),mass_sum=float(data['histogram_mass'].sum()),
                            color=colors[task],x_axis='bottom' if j==0 else 'top',y_axis='left' if j==0 else 'right')
    # Both output colors share one physical mass scale and one set of ticks.
    # Two adjacent color strips are drawn inside the same colorbar Axes.
    cbax=ax_mm(fig,height,x-19,y+size/2-13.5,3.2,27)
    mass_edges=np.geomspace(density_norm.vmin,density_norm.vmax,257)
    mass_values=np.sqrt(mass_edges[:-1]*mass_edges[1:])[:,None]
    for k,task in enumerate(TASKS[:2]):
        cbax.pcolormesh([k,k+1],mass_edges,mass_values,cmap=translucent_mass_cmap(colors[task]),
                        norm=density_norm,shading='flat',linewidth=0,rasterized=False)
    cbax.set(xlim=(0,2),ylim=(density_norm.vmin,density_norm.vmax),yscale='log')
    cbax.set_xticks([]);cbax.set_yticks([1e-5,1e-2])
    cbax.yaxis.set_minor_locator(NullLocator());cbax.yaxis.set_ticks_position('left')
    cbax.tick_params(labelsize=9,pad=1.5,length=1.5,color='#555555',labelcolor='#555555')
    for spine in cbax.spines.values():spine.set_visible(False)
    cbax.set_title('Mass',fontsize=9,color='#555555',pad=4)
    front=axes['chi_entropy']
    rmsf_text=axes['rmsf'].texts[0]
    annotation_text=rmsf_text.get_text()
    rmsf_text.remove()
    front.text(.035,.965,annotation_text,transform=front.transAxes,ha='left',va='top',
               fontsize=9,color=BLUE,linespacing=1.2,zorder=8)
    assert all(text.get_bbox_patch() is None for text in front.texts)
    # Both symmetric coordinate transforms map y=x to the same geometric diagonal.
    front.plot([0,1],[0,1],transform=front.transAxes,color='#999999',lw=.55,ls=(0,(1.3,3.2)),zorder=3)
    # Draw both original trends in the front axes, each with its own data transform.
    # This prevents either density layer from painting over the other trend.
    trend_lines={}
    for task in TASKS[:2]:
        tr=trend[trend.target.eq(task)&trend.method.eq('Full')].sort_values('bin_index')
        transform=axes[task].transData
        front.fill_between(tr.x_mean,tr.ci_low,tr.ci_high,transform=transform,
                           color=colors[task],alpha=.10,lw=0,zorder=4)
        line,=front.plot(tr.x_mean,tr.prediction_mean,transform=transform,color=colors[task],
                         ls='-',lw=1.45,solid_capstyle='butt',solid_joinstyle='miter',zorder=5)
        trend_lines[task]=line
    # No data-limit autoscaling after cross-transform artists are added.
    for task,ax in axes.items():ax.set(xlim=metadata[task]['limits'],ylim=metadata[task]['limits'])
    metadata['stroke_intersections']=blend_stroke_intersections(fig,front,trend_lines,colors,1.45)
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


PR_Y_SCALES={
    'rmsf':dict(values=[0,.5,1],heights=[0,.25,1],ticks=[0,.5,1]),
    'chi_entropy':dict(values=[0,.5,.9,1],heights=[0,.10,.45,1],ticks=[0,.9,1]),
    'within_type_hchi':dict(values=[0,.45,.75,1],heights=[0,.12,.90,1],ticks=[0,.6,1]),
}


def map_precision(values,knots,heights):
    values=np.asarray(values,dtype=float)
    knots=np.asarray(knots,dtype=float)
    heights=np.asarray(heights,dtype=float)
    mapped=np.interp(values,knots,heights)
    lower=heights[0]+(values-knots[0])*(heights[1]-heights[0])/(knots[1]-knots[0])
    upper=heights[-1]+(values-knots[-1])*(heights[-1]-heights[-2])/(knots[-1]-knots[-2])
    return np.where(values<knots[0],lower,np.where(values>knots[-1],upper,mapped))


def precision_functions(task):
    scale=PR_Y_SCALES[task]
    def forward(values):
        return map_precision(values,scale['values'],scale['heights'])
    def inverse(values):
        return map_precision(values,scale['heights'],scale['values'])
    return forward,inverse


def precision_recall(fig,height,x,y,w,h,task,letter,heading_y,ap_y,ylabel=True,heading=True,ap_label=True):
    ax=ax_mm(fig,height,x,y,w,h)
    frame=within_pr if task=='within_type_hchi' else pr[pr.target.eq(task)]
    ends=within_end if task=='within_type_hchi' else endpoints[endpoints.target.eq(task)]
    metrics=all_metrics[task]
    for method in METHODS:
        cur=frame[frame.method.eq(method)].sort_values('recall')
        end=ends[ends.method.eq(method)].iloc[0]
        assert len(cur)==101 and end.recall==1
        xx=np.r_[cur.recall.to_numpy(),end.recall];yy=np.r_[cur.precision.to_numpy(),end.precision]
        frozen=FROZEN_PR[FROZEN_PR.task.eq(task)&FROZEN_PR.method.eq(method)]
        assert np.array_equal(frozen[['recall','precision']].to_numpy(),np.column_stack([xx,yy]))
        RENDER_COUNTS[ACTIVE_OPTION]['PR_vertices']+=len(xx)
        ax.step(xx,yy,where='post',color=COLORS[method],ls=STYLES[method],lw=WIDTHS[method],zorder=ZORDER[method])
    baseline=float(frame.no_skill_precision.iloc[0])
    ax.axhline(baseline,color='#AAAAAA',ls=(0,(2,2)),lw=.65,zorder=0)
    forward,inverse=precision_functions(task)
    ax.set_yscale('function',functions=(forward,inverse))
    # Use the same 2.5% display-space padding, rather than the same raw upper limit.
    ax.set(xlim=(0,1),ylim=(0,float(inverse(1.025))),xlabel='Recall',ylabel='Precision')
    ax.set_xticks([0,.5,1]);ax.set_yticks(PR_Y_SCALES[task]['ticks'])
    number_format(ax,'y')
    ax.tick_params(pad=2);ax.xaxis.labelpad=2.5
    label={'rmsf':'RMSF','chi_entropy':'Hχ','within_type_hchi':'Hχ within type'}[task]
    if heading:title(fig,height,x,heading_y,letter,label)
    if ap_label:txt(fig,height,x-1,ap_y,'AP',fontsize=9,ha='left',va='center',color='#444444')
    # Same method ordering as the single shared legend, regardless of panel width.
    for k,method in enumerate(METHODS):
        txt(fig,height,x+(k+.5)*w/4,ap_y,f'{metrics[method]["mean_ap"]:.3f}',
            color=COLORS[method],fontsize=9,ha='center',va='center')
    return ax


def method_key(fig,height,y):
    return fig.legend(bbox_to_anchor=(176/W,y/height),loc='center right',
               handles=[Line2D([],[],color=COLORS[m],ls=STYLES[m],lw=1.4,label=m) for m in METHODS],
               ncol=4,columnspacing=1,handlelength=2,handletextpad=.4,borderaxespad=0,borderpad=0)


def distributions(ax,task,positions=(0,1,2,3),labels=True,ylabel=False):
    frame=frames[task]
    arrays=[frame[frame.method.eq(method)].sort_values('record_code').ap.to_numpy() for method in METHODS]
    assert all(len(v)==3063 for v in arrays)
    parts=ax.violinplot(arrays,positions=positions,widths=.78,showmeans=False,showmedians=False,
                        showextrema=False,points=128,bw_method=.20)
    for k,(method,values,body) in enumerate(zip(METHODS,arrays,parts['bodies'])):
        center=positions[k]
        body.set_facecolor(COLORS[method]);body.set_edgecolor(COLORS[method]);body.set_linewidth(.8);body.set_alpha(.25)
        codes=frame[frame.method.eq(method)].sort_values('record_code').record_code.to_numpy()
        frozen=FROZEN_AP[FROZEN_AP.task.eq(task)&FROZEN_AP.method.eq(method)].sort_values('record_code')
        assert np.array_equal(codes,frozen.record_code.to_numpy())
        assert np.array_equal(values,frozen.AP.to_numpy())
        RENDER_COUNTS[ACTIVE_OPTION]['AP_values']+=len(values)
        horizontal_offset=.13*np.sin(codes*2.399963229728653)
        ax.scatter(center+horizontal_offset,values,s=.45,color=COLORS[method],alpha=.09,linewidths=0,zorder=2,rasterized=False)
        quartiles=np.quantile(values,[.25,.5,.75])
        ax.vlines(center,quartiles[0],quartiles[2],color=COLORS[method],lw=2.5,zorder=3)
        ax.plot([center-.085,center+.085],[quartiles[1],quartiles[1]],color='white',lw=1,zorder=4)
    ax.set(ylim=(-.02,1.02),ylabel='Record AP' if ylabel else '')
    ax.set_yticks([0,.5,1]);ax.tick_params(pad=2)
    if labels:
        ax.set_xticks(positions,METHODS)
        for tick,method in zip(ax.get_xticklabels(),METHODS):tick.set_color(COLORS[method])
    return ax


def save_option(fig,height,stem,expected_letters):
    check=export(fig,stem,W,height)
    doc=fitz.open(ROOT/f'{stem}.pdf')
    labels={span['text'] for block in doc[0].get_text('dict')['blocks'] if 'lines' in block
            for line in block['lines'] for span in line['spans']}
    assert set(expected_letters)<=labels
    assert all('Arial' in font[3] for font in doc[0].get_fonts())
    assert abs(doc[0].rect.width*25.4/72-W)<.01
    assert abs(doc[0].rect.height*25.4/72-height)<.01
    factor=200/72*167.64/W
    doc[0].get_pixmap(matrix=fitz.Matrix(factor,factor),alpha=False).save(ROOT/f'{stem}_preview.png')
    doc.close();plt.close(fig)
    return check


if __name__=='__main__':
    scale_grid=np.linspace(0,1,1001)
    for task,scale in PR_Y_SCALES.items():
        forward,inverse=precision_functions(task)
        assert np.all(np.diff(forward(scale_grid))>0)
        assert np.allclose(inverse(forward(scale_grid)),scale_grid,atol=1e-15)
        assert np.allclose(forward(scale['values']),scale['heights'])
        assert np.isclose(forward(inverse(1.025)),1.025)
    forward,_=precision_functions('rmsf')
    assert np.allclose(forward([0,.5,1]),[0,.25,1])
    reports={}
    ACTIVE_OPTION='merged_distributions'
    RENDER_COUNTS[ACTIVE_OPTION]=dict(AP_values=0,PR_vertices=0)
    fig,h=canvas(150.)
    title(fig,h,15,145,'a','Full vs MD')
    _,overlay_info=overlay_density(fig,h,30,82,48)
    title(fig,h,105,145,'b','Record AP distributions')
    ax=ax_mm(fig,h,105,82,71,48)
    centers=[]
    for j,task in enumerate(TASKS):
        pos=np.arange(4)+5*j
        distributions(ax,task,positions=pos,labels=False,ylabel=True)
        centers.append(float(pos.mean()))
    ax.set_xlim(-.7,13.7)
    ax.set_xticks(centers,['RMSF','Hχ','Within-type\nHχ'])
    ax.tick_params(axis='x',length=0,pad=4)
    method_legend=method_key(fig,h,69.5)
    txt(fig,h,4,60,'c',fontsize=11,fontweight='bold',va='bottom')
    txt(fig,h,5,55,'AP',fontsize=9,ha='left',va='center',color='#444444')
    for j,task in enumerate(TASKS):
        x,width=[(15,42),(76,42),(137,42)][j]
        txt(fig,h,x+width/2,60,['RMSF','Hχ','Hχ within type'][j],fontsize=10,ha='center',va='bottom')
        pr_axis=precision_recall(fig,h,x,11,width,40,task,'',0,55,j==0,heading=False,ap_label=False)
        pr_axis.set_xlabel('')
    txt(fig,h,95.5,1.3,'Recall',fontsize=9,ha='center',va='bottom')
    fig.canvas.draw()
    legend_right_mm=method_legend.get_window_extent(fig.canvas.get_renderer()).x1/fig.dpi*25.4
    assert np.isclose(legend_right_mm,176)
    reports['merged_distributions']=save_option(fig,h,'fig2_compact_merged','abc')

    assert all(v==dict(AP_values=36756,PR_vertices=1224) for v in RENDER_COUNTS.values())
    audit=dict(status='numeric_and_canvas_passed_visual_pending',options=reports,rendered_and_verified_counts=RENDER_COUNTS,
               original_width_mm=183,original_height_mm=203,
               full_density_data_preserved=True,no_metric_annotation_backgrounds=True,
               main_panels=["a","b","c"],shared_mass_scale=True,mass_colorbar_axes_count=1,
               output_trends_solid=True,distribution_width_mm=71,previous_distribution_width_mm=58,
               equal_output_linewidth_pt=1.45,stroke_intersections=overlay_info['stroke_intersections'],
               method_legend=dict(right_aligned=True,right_edge_mm=float(legend_right_mm),center_y_mm=69.5),
               PR_group_has_three_tasks=True,
               c_heading_removed=True,PR_plot_height_mm=40,
               PR_precision_axis=dict(type='task_specific_continuous_piecewise_linear',tasks=PR_Y_SCALES,
                                      display_top_padding=.025,shared_y_scale=False,
                                      y_ticks_and_label_on_each_plot=True,original_PR_vertices_preserved=True),
               all_36756_record_AP_values_plotted_in_each_option=True,
               all_1224_PR_vertices_plotted_in_each_option=True,
               all_frozen_metrics_unchanged=True,no_new_statistical_tests=True,
               manuscript_replaced=False)
    (ROOT/'layout_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(audit,ensure_ascii=False,indent=2))


