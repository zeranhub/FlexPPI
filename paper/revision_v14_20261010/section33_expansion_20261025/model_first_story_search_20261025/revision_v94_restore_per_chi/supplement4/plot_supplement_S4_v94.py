# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a) Frozen ECDF -> original figure5 production + LineTrend -> param inherit
# (b) All83 paired correlations -> original figure5 production + PairedBoxScatter -> param inherit
# (c) All151 chi observations -> original figure5 production -> param inherit
# (d) All432 MD angular bins -> original figure5 production + LineTrend -> param inherit
# RULE: all native vector plots; retain every frozen point and original NaN gap.
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


# This module body is appended after the inherited verbatim style baselines.
from pathlib import Path
import hashlib
import json
import re
import shutil
import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
from PIL import Image

OUT = Path(__file__).resolve().parent
OLD = OUT.parent.parent/'revision_v80_redesign/figure5'
SOURCE = OUT/'supplement_source'
SOURCE.mkdir(exist_ok=True)
STEM = 'figS4_joint_readouts_angular_occupancy'
BLUE = '#2166AC'; WARM = '#B87936'; TEAL = '#3A8D8C'
BLACK = '#222222'; GREY = '#999999'
sha = lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
read = lambda p:pd.read_csv(p,float_precision='round_trip')
source_paths = [
    OLD/'source/original_fig4/fig4_joint_complementarity_all_ecdf_points.csv',
    OLD/'source/original_fig7/experiment/defined_angles_long.csv',
    OLD/'source/original_fig7/experiment/per_chi_correlations.csv',
    OLD/'source/original_fig7/source/histograms_501frames.csv',
    OLD/'source/original_fig7/source/all66_residue_summary.csv',
    OUT/'retained_absolute_correlations_all498.csv',
    OUT/'retained_panel_statistics.csv',
]
# Frozen input copies are bundled in supplement_source.

ecdf = read(SOURCE/source_paths[0].name)
defined = read(SOURCE/source_paths[1].name)
corr = read(SOURCE/source_paths[2].name)
corr = corr[corr.predictor.eq('predicted_joint_Hchi')&corr.metric.eq('fine_entropy')]
hist = pd.read_csv(SOURCE/source_paths[3].name) # match original float parser
ub = read(SOURCE/source_paths[4].name)
associations = read(SOURCE/source_paths[5].name)
statistics = read(SOURCE/source_paths[6].name)
statistics = statistics[statistics.panel_kind.eq('correlation')]
assert ecdf.groupby('role').size().to_dict()==dict(prediction=3063,reference=3063)
assert defined.groupby('chi').size().to_dict()=={'chi1':66,'chi2':52,'chi3':23,'chi4':10}
assert len(hist)==432 and hist.F_star.isna().sum()==211
assert len(associations)==498 and associations.sample_id.nunique()==83
font_manager.findfont('Arial',fallback_to_default=False)
mpl.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','Liberation Sans'],
    'font.size':8.5,'axes.labelsize':8.5,'axes.titlesize':8.5,'xtick.labelsize':8.5,
    'ytick.labelsize':8.5,'legend.fontsize':8.5,'axes.unicode_minus':False,
    'savefig.bbox':None})
W,H = 183.,168.
fig=plt.figure(figsize=(W/25.4,H/25.4),facecolor='white')
axis=lambda x,y,w,h:fig.add_axes([x/W,y/H,w/W,h/H])
text=lambda x,y,s,**kw:fig.text(x/W,y/H,s,**kw)
all_points=[]
axes={}
def panel(letter,x,y,title):
    text(x,y,letter,fontsize=11,weight='bold',va='baseline')
    text(x+6,y,title,fontsize=9,weight='bold',va='baseline')
def register(name,ax):
    axes[name]=ax
    ax.tick_params(length=2.5,pad=2)

# a: all 3063 ECDF vertices per curve, without averaging/sampling vertices.
panel('a',2,162,'RMSF–Hχ coupling')
ax=axis(15,117,45,36);register('a',ax)
for role,color,ls,name in [('reference',BLACK,'--','MD'),('prediction',BLUE,'-','Model')]:
    part=ecdf[ecdf.role.eq(role)].sort_values('y')
    ax.plot(part.x,part.y,color=color,lw=1.15,ls=ls)
    for r in part.itertuples():
        all_points.append(dict(panel='a',sample_id=r.sample_id,kind=role,x=r.x,y=r.y))
ax.set_xlim(-.25,1);ax.set_ylim(0,1.02)
ax.set_xticks([0,.5,1]);ax.set_yticks([0,.5,1])
ax.set_xlabel('RMSF–Hχ Spearman ρ',labelpad=3)
ax.set_ylabel('Fraction of records',labelpad=3)
text(15,156,'3,063 records',fontsize=8.5)
ax.text(.03,.96,'MD 0.391',transform=ax.transAxes,va='top',fontsize=8.5,color=BLACK)
ax.text(.03,.80,'Model 0.517',transform=ax.transAxes,va='top',fontsize=8.5,color=BLUE)

# b: all 83 paired records, three prediction readouts and two target descriptors.
panel('b',70,162,'Angular-occupancy association')
text(83,156,'83 records',fontsize=8.5)
settings=[('full_hchi',BLUE,'o',True),('intrachain_hchi',WARM,'D',True),('full_rmsf',BLUE,'^',False)]
for mi,(metric,left,title) in enumerate([('fine_entropy',83,'Fine-bin entropy'),('lowF_width',137,'Low-F* coverage')]):
    ax=axis(left,117,39,32);register('b_'+metric,ax)
    subset=associations[associations.metric_key.eq(metric)]
    paired=subset.pivot(index='sample_id',columns='predictor',values='rho_type_centered')
    assert paired.shape==(83,3)
    for sid,row in paired.iterrows():
        ax.plot([0,1,2],[row[p] for p,_,_,_ in settings],color='#C7C7C7',alpha=.45,lw=.35,zorder=1)
    for xi,(pred,color,marker,fill) in enumerate(settings):
        part=subset[subset.predictor.eq(pred)].sort_values('sample_id').reset_index(drop=True)
        jitter=((np.arange(83)*37)%83/82-.5)*.18
        ax.scatter(xi+jitter,part.rho_type_centered,s=8.5,color=color if fill else 'none',
            edgecolors=color,marker=marker,linewidths=.35,alpha=.55,zorder=2)
        st=statistics[statistics.metric_key.eq(metric)&statistics.predictor.eq(pred)].iloc[0]
        ax.errorbar(xi,st['mean'],yerr=[[st['mean']-st.bootstrap_95_low],[st.bootstrap_95_high-st['mean']]],
            fmt=marker,ms=4.5,mfc=color if fill else 'white',mec=BLACK,mew=.6,
            ecolor=BLACK,elinewidth=1.,capsize=2,zorder=5)
        for j,r in part.iterrows():
            all_points.append(dict(panel='b',metric=metric,sample_id=r.sample_id,kind=pred,
                x=float(xi+jitter[j]),y=r.rho_type_centered))
    ax.axhline(0,color=GREY,lw=.55,ls=':',zorder=0)
    ax.set_xticks([0,1,2],['Full Hχ','Intra Hχ','RMSF'],rotation=30,ha='right')
    ax.set_xlim(-.45,2.45);ax.set_ylim(-.06,.72);ax.set_yticks([0,.3,.6])
    ax.set_title(title,loc='left',pad=4,fontsize=8.5)
    if mi==0:ax.set_ylabel('Same-type Spearman ρ',labelpad=3)
    else:ax.tick_params(labelleft=False)

# c: prediction-selected actual MD metrics moved from the old main figure.
panel('c',2,100,'Actual MD occupancy in prediction-selected halves')
selection=read(OUT/'retained_prediction_groups_all498.csv')
selection_stats=read(OUT/'retained_panel_statistics.csv')
for metric,left,title,lim,ticks in [('fine_entropy',30,'Fine-bin entropy (nats/χ)',(0,3.4),[0,1.5,3]),('lowF_width',116,'Low-F* width (°/χ)',(0,60),[0,30,60])]:
    ax=axis(left,64,46,27);register('c_'+metric,ax)
    part=selection[selection.metric_key.eq(metric)&selection.predictor.eq('full_hchi')].sort_values('sample_id').reset_index(drop=True)
    assert len(part)==83
    for j,r in part.iterrows():
        jitter=((j*37)%83/82-.5)*.18
        ax.plot([0+jitter,1+jitter],[r.low_MD_mean,r.high_MD_mean],color='#AEB7BC',lw=.45,alpha=.5,zorder=1)
    for xi,field in enumerate(['low_MD_mean','high_MD_mean']):
        jitter=((np.arange(83)*37)%83/82-.5)*.18
        ax.scatter(xi+jitter,part[field],s=6,color=BLACK,alpha=.3,linewidths=0,zorder=2)
        st=selection_stats[selection_stats.panel_kind.eq('prediction_selection')&selection_stats.metric_key.eq(metric)&selection_stats.predictor.eq('full_hchi')&selection_stats.value_key.eq(field)].iloc[0]
        ax.errorbar(xi,st['mean'],yerr=[[st['mean']-st.bootstrap_95_low],[st.bootstrap_95_high-st['mean']]],fmt='D',ms=4.5,mfc='white',mec=BLACK,mew=.8,ecolor=BLACK,elinewidth=1,capsize=2,zorder=5)
        for j,r in part.iterrows():all_points.append(dict(panel='c',metric=metric,sample_id=r.sample_id,kind=field,x=float(xi+jitter[j]),y=r[field]))
    ax.set(xlim=(-.3,1.3),ylim=lim);ax.set_xticks([0,1],['Low','High']);ax.set_yticks(ticks)
    ax.set_title(title,fontsize=8.5,loc='left',pad=4)
text(99,55,'Groups selected by predicted joint Hχ',fontsize=8.5,ha='center')

# d: all 432 original histogram bins and original zero-count NaN gaps.
panel('d',2,44,'MD angular occupancy at prediction-selected Gln sites')
text(20,38,'Q41: lowest Hχ',fontsize=8.5,color=WARM)
text(89,38,'Q2: highest Hχ',fontsize=8.5,color=TEAL)
widths=[]
for chi,left in [(1,19),(2,75),(3,131)]:
    ax=axis(left,10,45,24);register('d_chi'+str(chi),ax)
    for site,color,ls in [('Q41',WARM,'-'),('Q2',TEAL,'--')]:
        part=hist[hist.residue.eq(site)&hist.chi.eq('chi'+str(chi))].sort_values('angle_center_degrees')
        assert len(part)==72
        ax.plot(part.angle_center_degrees,part.F_star,color=color,lw=1.1,ls=ls)
        widths.append(dict(residue=site,chi=chi,origin0_coverage_degrees=5*int((part.F_star.le(1)&part.F_star.notna()).sum())))
        for r in part.itertuples():
            all_points.append(dict(panel='d',residue_id=int(site[1:]),chi=chi,kind='MD_Fstar',
                x=r.angle_center_degrees,y=r.F_star))
    ax.axhline(1,color=GREY,lw=.6,ls=':')
    ax.text(.97,.94,'χ'+str(chi),transform=ax.transAxes,ha='right',va='top',fontsize=8.5,weight='bold')
    ax.set_xlim(-180,180);ax.set_ylim(-.12,5.9)
    ax.set_xticks([-180,-60,60,180]);ax.set_yticks([0,2,4])
    ax.set_xlabel('Side-chain angle (°)',labelpad=2)
    if chi==1:ax.set_ylabel('MD F*',labelpad=3)
    else:ax.tick_params(labelleft=False)

fig.canvas.draw()
texts=[t for t in fig.findobj(mpl.text.Text) if t.get_visible() and t.get_text()]
renderer=fig.canvas.get_renderer()
outside=[]
for t in texts:
    bbox=t.get_window_extent(renderer)
    if bbox.x0<-.5 or bbox.y0<-.5 or bbox.x1>fig.bbox.width+.5 or bbox.y1>fig.bbox.height+.5:
        outside.append(t.get_text())
assert not outside,outside
assert min(t.get_fontsize() for t in texts)>=8.5
assert not [t.get_text() for t in texts if re.search('[\u3400-\u9fff]',t.get_text())]
points=pd.DataFrame(all_points)
assert points.groupby('panel').size().to_dict()=={'a':6126,'b':498,'c':332,'d':432}
assert points[points.panel.eq('d')].y.isna().sum()==211
assert hist.F_star.max()<5.9
for ext in ['pdf','svg','png']:
    fig.savefig(OUT/(STEM+'.'+ext),dpi=450,bbox_inches=fig.bbox_inches)
assert not ET.parse(OUT/(STEM+'.svg')).findall('.//{http://www.w3.org/2000/svg}image')
points.to_csv(OUT/'figS4_all_plotted_points.csv',index=False,float_format='%.17g')
pd.DataFrame(widths).to_csv(OUT/'figS4_original_origin0_coverages.csv',index=False)
(OUT/'figS4_visible_text.json').write_text(json.dumps([t.get_text() for t in texts],ensure_ascii=False,indent=2),encoding='utf8')
audit=dict(status='numerical_passed_visual_pending',canvas_mm=[W,H],dpi=450,
    minimum_font_pt=min(t.get_fontsize() for t in texts),text_outside_canvas=outside,
    all_source_points_retained=True,no_smoothing=True,no_downsampling=True,
    plotted_rows=points.groupby('panel').size().to_dict(),NaN_Fstar_bins=211,
    all_native_vector=True,baseline_source=str(OLD/'compose_fig5.py'),
    style_asset_notes='Original scientific production inherited; generic LineTrend and PairedBoxScatter assets inspected. Their example time-series and before/after boxplot data do not match frozen ECDF/MD angular histograms, so parameter inheritance preserves scientific native vector data.',
    source_sha256={p.name:sha(p) for p in SOURCE.iterdir() if p.is_file()},
    output_sha256={ext:sha(OUT/(STEM+'.'+ext)) for ext in ['png','pdf','svg']},
    visual_QA='pending')
(OUT/'figS4_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps({k:audit[k] for k in ['status','canvas_mm','minimum_font_pt','plotted_rows']},ensure_ascii=False))
