# Academic Figure Skill Asset Confirmation
# (a,b) original manuscript production functions, frozen all-data inputs -> visual adapt
# (c,d) original frozen complete PR and external curves -> visual adapt
# (e) original native bar PDF -> unchanged native crop; uniform letter/legend
# All physical font sizes retained; only axes sizes/layout and b's lower margin change.
from pathlib import Path
import sys,importlib.util,json,hashlib
import numpy as np,pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.transforms import Bbox
import fitz

OUT=Path(__file__).resolve().parent
sys.path.insert(0,str(OUT/'source_internal'))
def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m);return m
prod=module('original_internal',OUT/'source_internal'/'plot.py')
ext=module('original_external',OUT/'source_external.py')
ARCH=OUT/'source_archive'
ext.vertices=pd.read_csv(ARCH/'internal_PR_vertices.csv',float_precision='round_trip')
ext.summary=json.loads((ARCH/'internal_summary.json').read_text(encoding='utf-8'))
ext.curves=pd.read_csv(ARCH/'external_mean_curves.csv',float_precision='round_trip')
ext.curves=ext.curves[ext.curves.scope.eq('all')]
ext.metrics=pd.read_csv(ARCH/'external_summary_metrics.csv',float_precision='round_trip')
ext.metrics=ext.metrics[ext.metrics.scope.eq('all')]
ext.audit=json.loads((ARCH/'external_calculation_audit.json').read_text(encoding='utf-8'))
ext.prevalence=float(ext.audit['all_selected_precision'])
W,H=183.,167.
PT=72/25.4
mpl.rcParams.update({'font.size':9,'axes.labelsize':9,'axes.titlesize':10,
    'xtick.labelsize':9,'ytick.labelsize':9,'legend.fontsize':9,
    'pdf.fonttype':42,'svg.fonttype':'none'})

def ax_at(fig,x,top,w,h):return fig.add_axes([x/W,1-(top+h)/H,w/W,h/H])
def text(fig,x,baseline,s,fontsize=10,weight='normal',ha='left',color='black'):
    return fig.text(x/W,1-baseline/H,s,fontsize=fontsize,fontweight=weight,
        va='baseline',ha=ha,color=color)
def label(fig,x,baseline,s):return text(fig,x,baseline,s,fontsize=11,weight='bold')
fig=plt.figure(figsize=(W/25.4,H/25.4))
label(fig,4,5,'a');text(fig,15,5,'Full vs MD')
label(fig,94,5,'b');text(fig,105,5,'Record AP distributions')
# Original dense prediction-reference function uses all histogram cells/trends.
overlay,overlay_info=prod.overlay_density(fig,H,x=30,y=H-52,size=42)
prod.ACTIVE_OPTION='dense';prod.RENDER_COUNTS['dense']=dict(AP_values=0,PR_vertices=0)
bax=ax_at(fig,105,10,71,42)
centers=[]
for j,task in enumerate(prod.TASKS):
    positions=np.arange(4)+5*j
    prod.distributions(bax,task,positions=positions,labels=False,ylabel=True)
    centers.append(float(positions.mean()))
bax.set(xlim=(-.7,13.7),ylim=(-.02,1.02))
bax.set_xticks(centers,['RMSF','Hχ','Within-type\nHχ'])
bax.tick_params(axis='x',length=0,pad=4)
bhandles=[Line2D([],[],color=prod.COLORS[m],ls=prod.STYLES[m],lw=1.4,label=m) for m in prod.METHODS]
blegend=bax.legend(handles=bhandles,loc='center',bbox_to_anchor=(.5,.20),
    ncol=2,fontsize=9,columnspacing=.8,handlelength=1.2,handletextpad=.35,
    borderpad=0,borderaxespad=0,frameon=False)

ROW2_TOP,ROW2_H,ROW2_BASE=72.,34.,67.5
# Reuse original within-type PR function; no metric annotation is drawn.
cax=ext.original_pr(fig,H,'within_type_hchi',15,42,ROW2_TOP,ROW2_H,'',ROW2_BASE,None)
# Remove function-generated heading (its va='bottom' differs from the common baseline).
for t in list(fig.texts):
    if t.get_text()=='Hχ within type':t.remove()
label(fig,4,ROW2_BASE,'c');text(fig,36,ROW2_BASE,'Hχ within type',ha='center')
clegend=cax.legend(handles=bhandles,loc='upper right',bbox_to_anchor=(.99,1.0),
    fontsize=9,ncol=1,handlelength=2,handletextpad=.4,labelspacing=.1,
    borderpad=.15,borderaxespad=0,frameon=False)
label(fig,65,ROW2_BASE,'d')
curve_lines={}
for kind,x in [('PR',76),('ROC',134)]:
    ax=ax_at(fig,x,ROW2_TOP,42,ROW2_H)
    for method,name,color,ls,marker in reversed(ext.EXT):
        part=ext.curves[(ext.curves.method.eq(method))&(ext.curves.curve_type.eq(kind))].sort_values('x')
        xx,yy=part.x.to_numpy(),part['mean'].to_numpy()
        line,=ax.plot(xx,yy,color=color,ls=ls,lw=1.35 if method=='Full' else .95,
            marker=marker,markevery=[180,380,580,780] if marker else None,
            ms=2.25,mfc='white',mew=.55,zorder=10 if method=='Full' else 4)
        curve_lines[(kind,method)]=line
        if kind=='PR':ax.plot([1,1],[yy[-1],ext.prevalence],color=color,ls=ls,lw=.65,zorder=3,clip_on=False)
        else:ax.plot([0,0],[0,yy[0]],color=color,ls=ls,lw=.65,zorder=3,clip_on=False)
    if kind=='PR':
        ax.axhline(ext.prevalence,color='#AAAAAA',ls=(0,(2,2)),lw=.65,zorder=0)
        ax.set(xlabel='Recall',ylabel='Precision')
    else:
        ax.plot([0,1],[0,1],color='#AAAAAA',ls=(0,(2,2)),lw=.65,zorder=0)
        ax.set(xlabel='False-positive rate',ylabel='True-positive rate')
    ax.set(xlim=(0,1),ylim=(0,1.025));ax.set_xticks([0,.5,1]);ax.set_yticks([0,.5,1])
    ax.tick_params(pad=2,length=2.5,labelsize=8)
    ax.xaxis.labelpad=2.5;ax.yaxis.labelpad=2.5
    text(fig,x+21,ROW2_BASE,'RMSF '+kind,fontsize=9.5,ha='center')
    if kind=='ROC':
        handles=[Line2D([],[],color=c,ls=ls,lw=1.35 if m=='Full' else .95,
            marker=mark,ms=2.25,mfc='white',mew=.55,label=name) for m,name,c,ls,mark in ext.EXT]
        rlegend=ax.legend(handles=handles,loc='lower right',bbox_to_anchor=(.99,.005),
            fontsize=7.5,ncol=1,handlelength=1.35,handletextpad=.45,labelspacing=.1,
            borderpad=.1,borderaxespad=0,frameon=False)
        rlegend.get_texts()[0].set_fontweight('bold')

# Native original bottom crop: bars/axes/method labels retain their 29.5mm height.
BOTTOM_AXIS_TOP=124.;BOTTOM_SCALE=183/150
BOTTOM_SOURCE_AXIS_TOP=52*(1-.795)
BOTTOM_ORIGIN=BOTTOM_AXIS_TOP-BOTTOM_SOURCE_AXIS_TOP*BOTTOM_SCALE
label(fig,4,120.7,'e')
handles=[Patch(facecolor=c,edgecolor='#2B4A63',lw=.25*BOTTOM_SCALE,label=n)
    for n,c in [('1 interface','#B9CEE3'),('2–3 interfaces','#5E8FBF'),('≥4 interfaces','#2060A8')]]
elegend=fig.legend(handles=handles,loc='center',bbox_to_anchor=(.5,1-119.55/H),ncol=3,
    fontsize=7.442,handlelength=.85,handletextpad=.4,columnspacing=2.4,frameon=False)
fig.canvas.draw();renderer=fig.canvas.get_renderer()
legend_data_box=bax.transData.inverted().transform(blegend.get_window_extent(renderer).get_points())
minimum_ap=float(prod.FROZEN_AP.AP.min())
legend_physical=blegend.get_window_extent(renderer)
ap_points=[]
for collection in bax.collections:
    if isinstance(collection,mpl.collections.PathCollection):
        positions=np.asarray(collection.get_offsets())
        ap_points.append(collection.get_offset_transform().transform(positions))
all_display=np.concatenate(ap_points)
assert len(all_display)==36756
assert not np.any((all_display[:,0]>=legend_physical.x0)&(all_display[:,0]<=legend_physical.x1)
    &(all_display[:,1]>=legend_physical.y0)&(all_display[:,1]<=legend_physical.y1))
# The two-column key is under the Hchi distributions, whose lowest value is .375.
assert legend_data_box[1,1]<.375,legend_data_box
for collection in bax.collections:
    if isinstance(collection,mpl.collections.PolyCollection):
        for path in collection.get_paths():
            bb=path.get_extents(collection.get_transform())
            assert not bb.overlaps(legend_physical),('violin-legend bbox overlap',bb,legend_physical)
assert prod.RENDER_COUNTS['dense']['AP_values']==36756
for (kind,m),line in curve_lines.items():
    q=ext.curves[(ext.curves.method.eq(m))&(ext.curves.curve_type.eq(kind))].sort_values('x')
    assert np.array_equal(line.get_xydata(),q[['x','mean']].to_numpy())
box=Bbox.from_bounds(0,0,W/25.4,H/25.4)
fig.savefig(OUT/'dense_new_panels.pdf',bbox_inches=box,dpi=450)
plt.close(fig)

# Place the untouched native bar block on top of the new complete vector canvas.
new=fitz.open(OUT/'dense_new_panels.pdf');p=new[0]
bottom=fitz.open(ARCH/'original_fig2_panel_d.pdf')
# Isolated old d is removed; crop starts above the 1.0 ticks and every bar.
span=next(s for block in bottom[0].get_text('dict')['blocks'] if 'lines' in block
    for line in block['lines'] for s in line['spans'] if s['text']=='d')
r=fitz.Rect(span['bbox']);r.x0-=.25;r.x1+=.25;r.y1=span['origin'][1]+.75
bottom[0].add_redact_annot(r,fill=(1,1,1));bottom[0].apply_redactions(images=0,graphics=0,text=0)
source_clip=fitz.Rect(0,8*PT,150*PT,52*PT)
target=fitz.Rect(0,(BOTTOM_ORIGIN+8*BOTTOM_SCALE)*PT,183*PT,(BOTTOM_ORIGIN+52*BOTTOM_SCALE)*PT)
p.show_pdf_page(target,bottom,0,clip=source_clip)
new.subset_fonts();new.save(OUT/'fig2_dense.pdf',garbage=4,deflate=True)
new.close();bottom.close()
d=fitz.open(OUT/'fig2_dense.pdf');p=d[0]
p.get_pixmap(dpi=450,alpha=False).save(OUT/'fig2_dense.png')
p.get_pixmap(dpi=160,alpha=False).save(OUT/'fig2_dense_preview.png')
svg=p.get_svg_image(text_as_path=False);(OUT/'fig2_dense.svg').write_text(svg,encoding='utf-8')
assert len(p.get_images(full=True))==0 and '<image' not in svg
ss=[s for b in p.get_text('dict')['blocks'] if 'lines' in b for l in b['lines'] for s in l['spans']]
outside=[s['text'] for s in ss if s['bbox'][0]<-.5 or s['bbox'][1]<-.5 or s['bbox'][2]>p.rect.width+.5 or s['bbox'][3]>p.rect.height+.5]
assert not outside,outside
audit=dict(status='rendered_visual_review_pending',width_mm=W,height_mm=H,
    all_36756_AP_points_preserved=True,min_record_AP=minimum_ap,b_legend_bbox_data_coordinates=legend_data_box.tolist(),
    b_y_axis_limits=[-.02,1.02],b_y_ticks=[0,.5,1],b_legend_no_scatter_or_violin_overlap=True,
    density_all_data=overlay_info,external_mean_points_preserved=14014,
    internal_within_type_vertices_preserved=408,
    axes_sizes_mm=dict(a=[42,42],b=[71,42],c=[42,34],d_PR=[42,34],d_ROC=[42,34],e_native_bar_height=29.5),
    data_recomputed=False,statistics_changed=False,font_sizes_scaled=False,
    raster_images_pdf=0,raster_images_svg=0,all_text_inside_canvas=True)
(OUT/'dense_QA.json').write_text(json.dumps(audit,indent=2,ensure_ascii=False),encoding='utf-8')
print(json.dumps({k:audit[k] for k in ['status','width_mm','height_mm','min_record_AP','b_legend_bbox_data_coordinates','axes_sizes_mm']},indent=2))
