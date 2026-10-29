# Academic Figure Skill Asset Confirmation (verified against existing native production scripts)
# a -> original v96 PyMOL overview PNG/PSE -> native run, no image editing
# b/c -> v96 absolute residue profiles / LineTrend -> param inherit
# d/e/f -> v96 complete scatter and association plots -> param inherit
# All model and MD arrays are frozen; layout, typography, and annotations only.
# Style target: verified current manuscript Fig2; evidence-focused redesign with horizontal per-chi row.
# Academic Figure Skill Typography Baseline — COPY VERBATIM, place at TOP of script
import matplotlib as mpl
mpl.use("Agg")
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




# Reference-informed composition: MAVE-NN Fig4 / ENGINE Fig4 / RMSF-net Fig4.
# The published figures are design references, never source artwork for this figure.
from pathlib import Path
import json,hashlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from PIL import Image
ROOT=Path(__file__).resolve().parent.parent
OUT=Path(__file__).resolve().parent
P=OUT/'source_data';R=ROOT/'pymol'
W,H=183.,181.
BLUE='#2166AC';WARM='#B87936';TEAL='#3A8D8C';PURPLE='#8B72A5';GRAY='#BAC1C7';INK='#222222'
FS=9.0
mpl.rcParams.update({'savefig.bbox':None,'font.size':FS,'axes.labelsize':FS,'xtick.labelsize':FS,'ytick.labelsize':FS,'legend.fontsize':FS,'axes.unicode_minus':False,'axes.linewidth':.6})
cam=pd.read_csv(P/'cam_all25.csv').sort_values('CaMKI_native_position')
ub=pd.read_csv(P/'ubiquitin_72_position_profile_with_original_mask.csv').sort_values('residue_id')
gln=pd.read_csv(P/'ubiquitin_all6_Gln_Hchi_SASA_plot.csv').sort_values('residue_id')
corr=pd.read_csv(P/'retained_absolute_correlations_all498.csv')
stats=pd.read_csv(P/'retained_panel_statistics.csv')
manifest=json.loads((R/'render_manifest.json').read_text(encoding='utf-8'))
assert len(cam)==25 and ub.valid_hchi.sum()==66 and len(ub)==72 and len(gln)==6
fig=plt.figure(figsize=(W/25.4,H/25.4),facecolor='white');plotted=[]
def txt(x,y,s,size=FS,**kw):return fig.text(x/W,1-y/H,s,fontsize=size,va=kw.pop('va','top'),color=kw.pop('color',INK),**kw)
def ax(x,y,w,h):return fig.add_axes([x/W,1-(y+h)/H,w/W,h/H])
def panel(label,title,x,y):
    txt(x,y,label,11,fontweight='bold')
    if title:txt(x+6,y+.25,title,10,fontweight='normal')
def line(x1,y1,x2,y2,**kw):fig.add_artist(Line2D([x1/W,x2/W],[1-y1/H,1-y2/H],transform=fig.transFigure,**kw))
def values(panel,kind,ids,vals,**kw):
    for i,v in zip(ids,vals):
        if np.isfinite(v):plotted.append(dict(panel=panel,kind=kind,identity=str(i),value=float(v),**kw))
def molecule(key,x,y,w,h,labels):
    z=manifest['assets'][key];a=ax(x,y,w,h);a.imshow(Image.open(R/z['png']).convert('RGBA'));a.axis('off')
    # Display viewport only: keep original PNG/PSE and camera/coordinates intact.
    crop={'cam_overview':(200,110,1430,1315),'ub_overview':(90,205,1535,1195)}[key]
    a.set_xlim(crop[0],crop[2]);a.set_ylim(crop[3],crop[1]);a.set_aspect('equal')
    for name,location in labels.items():
        px=z['landmarks'][name]
        a.scatter([px[0]],[px[1]],s=11,facecolors='none',edgecolors=TEAL,linewidth=.85)
        a.annotate(name,xy=px,xytext=location,textcoords='axes fraction',color=INK,fontsize=FS,ha='center',va='center',arrowprops=dict(arrowstyle='-',color=TEAL,lw=.7),bbox=dict(facecolor='white',edgecolor='none',pad=.1,alpha=.92))
    return a

from publication_style import BLUE,WARM,TEAL,INK,MD,LIGHT
mpl.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','Helvetica','Liberation Sans'],'font.size':9,'axes.labelsize':9,'axes.titlesize':10,'xtick.labelsize':9,'ytick.labelsize':9,'legend.fontsize':9,'axes.linewidth':.6,'xtick.major.width':.6,'ytick.major.width':.6,'axes.unicode_minus':False,'pdf.fonttype':42,'svg.fonttype':'none'})
axis_records=[]
def style_axis(a,name):
    a.tick_params(length=2.5,pad=2,width=.6);a.xaxis.labelpad=2.5;a.yaxis.labelpad=2.5
    axis_records.append((name,a))
handles=[Line2D([],[],color=c,lw=lw,ls=ls,marker=mark,ms=2.2,mfc='white',mew=.6,label=name) for name,c,lw,ls,mark in [('MD',MD,.9,':','o'),('Full',BLUE,1.25,'-',None)]]

# a: retain the genuine molecular contexts together.
panel('a','',3,2)
txt(27,3,'CaM–CaMKI\n(1MXE)',ha='center',linespacing=1.12)
txt(70,3,'Ubiquitin–UBA\n(2QHO)',ha='center',linespacing=1.12)
molecule('cam_overview',7,14,39,37,{'W303':(.12,.68),'M316':(.91,.08)})
molecule('ub_overview',49,14,39,37,{'Q41':(.64,.96),'Q62':(.91,.10)})
for xx,label,col in [(13,'Analysed chain',PURPLE),(59,'Partner',GRAY)]:
    line(xx,55,xx+4,55,color=col,lw=2.2);txt(xx+5.5,53.4,label)

# b/c: both metrics per case; color = metric, line style = source.
import sys
sys.path.insert(0,str(ROOT))
from case_profiles import draw_case,method_handles,draw_sasa
# Original static context tracks and all their values are preserved.
# b: dual metrics plus original shielding track.
panel('b','CaMKI dynamics and solvent exposure',98,2)
fig.legend(handles=method_handles(),loc='center',bbox_to_anchor=(145/W,1-10/H),ncol=2,handlelength=1.8,columnspacing=1.2,frameon=False)
x=cam.CaMKI_native_position.to_numpy()
a,b,profile_values=draw_case(fig,(W,H),(111,18,60,22),'CaMKI',xlabel=False)
for (metric,series),z in profile_values.groupby(['metric','series'],sort=False):
    values('b','RMSF_A' if metric=='RMSF' else 'Hchi_nats',z.residue_id,z.value,series=series)
mae_full=float(np.mean(abs(cam.pred_rmsf-cam.target_rmsf)));mae_intra=float(np.mean(abs(cam.intra_rmsf-cam.target_rmsf)))
# Identify sites along the residue axis without implying equality across y axes.
for pos in [303,316]:
    a.axvline(pos,color=TEAL,lw=.4,ls=':',alpha=.6,zorder=0)
    a.text(pos,.98,('W' if pos==303 else 'M')+str(pos),transform=a.get_xaxis_transform(),ha='center',va='top',fontsize=9,color=TEAL)
a,sasa_values=draw_sasa(fig,(W,H),(111,44,60,9),'CaMKI')
values('b','static_sidechain_SASA_A2',sasa_values.residue_id,sasa_values.sidechain_SASA_A2,series='starting_structure')

# c/d: complete sequence accuracy, then same chemistry with different environments.
panel('c','Ubiquitin backbone and side chains',3,62)

a,b,profile_values=draw_case(fig,(W,H),(14,74,69,22),'Ubiquitin',xlabel=False)
for (metric,series),z in profile_values.groupby(['metric','series'],sort=False):
    values('c','RMSF_A' if metric=='RMSF' else 'Hchi_nats',z.residue_id,z.value,series=series)
for pos in [41,62]:a.axvspan(pos-.55,pos+.55,color=TEAL,alpha=.065,lw=0,zorder=0)
# Positional annotations sit on the shared x axis, outside the measured y values.
a.text(41,.99,'Q41',transform=a.get_xaxis_transform(),ha='center',va='top',fontsize=9,color=TEAL)
a.text(62,.99,'Q62',transform=a.get_xaxis_transform(),ha='center',va='top',fontsize=9,color=TEAL)
a,sasa_values=draw_sasa(fig,(W,H),(14,100,69,9),'Ubiquitin')
values('c','static_sidechain_SASA_A2',sasa_values.residue_id,sasa_values.sidechain_SASA_A2,series='starting_structure')
# d: compact single-axis RMSF/Hchi distribution, formerly panel f.
from plot_overlay import draw_overlay,load_data
panel('d','RMSF–Hχ relationship',98,62)
density_audit=draw_overlay(fig,(W,H),(111,80,60,31),xmax=10,gridsize=(50,35),marginal_mm=4.)
direct_arrays,direct_stats=load_data()

# e: all four angle comparisons in one horizontal row, with identical limits.
panel('e','Per-χ occupancy association',3,123)
eh=[Line2D([],[],ls='',marker=m,ms=3.6,color=c,label=name) for name,m,c in [('Other','o',BLUE),('Q41','D',TEAL),('Q2','s',TEAL)]]
fig.legend(handles=eh,loc='center',bbox_to_anchor=(145/W,1-126/H),ncol=3,handletextpad=.25,columnspacing=.7,borderpad=0,frameon=False)
defined=pd.read_csv(P/'defined_angles_long.csv',float_precision='round_trip')
chi_stats=pd.read_csv(P/'per_chi_correlations.csv',float_precision='round_trip');chi_stats=chi_stats[chi_stats.predictor.eq('predicted_joint_Hchi')&chi_stats.metric.eq('fine_entropy')]
for chi,left in [(1,14),(2,57.5),(3,101),(4,144.5)]:
    a=ax(left,140,35.5,25);style_axis(a,'e_chi'+str(chi))
    z=defined[defined.chi.eq('chi'+str(chi))];st=chi_stats[chi_stats.chi.eq('chi'+str(chi))].iloc[0]
    assert len(z)=={1:66,2:52,3:23,4:10}[chi]
    txt(left,134.5,'χ'+str(chi)+'  n='+str(len(z))+'  ρ='+format(st.rho_AA_centered,'.2f'))
    for name,marker,col in [('Other','o',BLUE),('Q41','D',TEAL),('Q2','s',TEAL)]:
        sub=z[~z.author_residue_id.isin([2,41])] if name=='Other' else z[z.author_residue_id.eq(int(name[1:]))]
        a.scatter(sub.predicted_Hchi_nats_AA_centered,sub.fine_entropy_four_origins_mean_nats_AA_centered,s=8.5 if name=='Other' else 19,marker=marker,color=col,alpha=.44 if name=='Other' else 1,linewidths=.25,edgecolors='white' if name=='Other' else col,zorder=3 if name=='Other' else 4)
    for row in z.itertuples():
        plotted.extend([dict(panel='e',kind='predicted_joint_Hchi_AA_centered',identity=str(int(row.author_residue_id)),chi=chi,series='Full',value=float(row.predicted_Hchi_nats_AA_centered)),dict(panel='e',kind='MD_single_chi_fine_entropy_AA_centered',identity=str(int(row.author_residue_id)),chi=chi,series='MD',value=float(row.fine_entropy_four_origins_mean_nats_AA_centered))])
    a.axhline(0,color=LIGHT,lw=.4,ls=':');a.axvline(0,color=LIGHT,lw=.4,ls=':')
    a.set(xlim=(-.7,.5),ylim=(-1,.75));a.set_xticks([-.5,0,.5]);a.set_yticks([-.5,0,.5])
    if chi!=1:a.tick_params(labelleft=False)
txt(3.1,152.5,'MD entropy deviation (nats)',rotation=90,ha='center',va='center')
txt(97,171,'Predicted joint Hχ deviation (nats)',ha='center')

df=pd.DataFrame(plotted);df.to_csv(OUT/'main_original_panels_a_to_e_values.csv',index=False)
assert df.groupby('panel').size().to_dict()=={'b':117,'c':348,'e':302}
fig.canvas.draw();renderer=fig.canvas.get_renderer();outside=[]
for t in fig.findobj(mpl.text.Text):
    if t.get_visible() and t.get_text():
        bb=t.get_window_extent(renderer)
        if bb.x0<-.5 or bb.y0<-.5 or bb.x1>fig.bbox.width+.5 or bb.y1>fig.bbox.height+.5:outside.append(t.get_text())
assert not outside,outside
for ext in ['png','pdf','svg']:fig.savefig(OUT/f'fig3_dual_case_entropy.{ext}',dpi=450,bbox_inches=fig.bbox_inches)
qa={'size_mm':[W,H],'font_pt':9,'a_b_c_e_data_unchanged_from_v109':True,'old_d_removed':True,'old_f_relocated_to_d':True,'density_axis_box_mm':[111,80,60,31],'b_c_static_SASA_tracks':True,'b_c_both_metrics_added':True,'e_four_horizontal':True,
 'd':'compact single-axis overlay of MD and Full RMSF-Hchi densities and marginals','common_d_positions':824988,
 'd_records':3063,'density':density_audit,'outside_text':outside,'model_and_MD_rerun':False,'word_changed':False,
 'visual_review':'pending','sha256':{e:hashlib.sha256((OUT/f'fig3_dual_case_entropy.{e}').read_bytes()).hexdigest() for e in ['png','pdf','svg']}}
(OUT/'figure_QA.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8')
plt.close(fig)
print(json.dumps(qa,ensure_ascii=False))
