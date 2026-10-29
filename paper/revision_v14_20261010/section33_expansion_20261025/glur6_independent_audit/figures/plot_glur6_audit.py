# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a) residue response line -> LineTrend/plot_sweep.py -> param inherit
# (b) grouped model comparison dots -> cross-type inherit BarComparison/plot_comparison_Trajectory.py -> param inherit
# (c) trajectory block comparison -> LineTrend/plot_sweep.py -> param inherit
# (d) largest-effect removal sweep -> LineTrend/plot_sweep.py -> param inherit
# RULE: "native run" = load pre-rendered PNG via Image.open().ax.imshow().
#       "param inherit" = drawing function below that copies Class A/B/C values.
#       If a panel says "native run" and you write a drawing function, you broke the contract.

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
import numpy as np
import pandas as pd
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator

# Explicit task styling after the unmodified shared baselines.
mpl.rcParams.update({'font.sans-serif':['Arial','DejaVu Sans'],
                     'xtick.labelsize':8, 'ytick.labelsize':8,
                     'axes.labelsize':8.5, 'axes.titlesize':8.5,
                     'text.color':BLACK, 'axes.labelcolor':BLACK,
                     'xtick.color':BLACK,'ytick.color':BLACK})
HERE=Path(__file__).resolve().parent
DATA=HERE.parent
FILES=['all_383_residue_deltas_and_windows.csv','delta_prediction_metrics.csv',
       'timeblock_prediction_metrics.csv']
res, metrics, blocks=[pd.read_csv(DATA/f) for f in FILES]
MODELS=['Full','Intra','AA_SASA_contacts','Static47_SASA','Zero_change']
LABELS={'Full':'Full','Intra':'Intra','AA_SASA_contacts':'AA + SASA + contacts',
        'Static47_SASA':'Static47 + SASA','Zero_change':'Zero change'}
STYLES={
    'Full':dict(color=CATEGORICAL[0],marker='o',linestyle='-',linewidth=1.55,zorder=5),
    'Intra':dict(color=GREY,marker='^',linestyle='--',linewidth=1.0,zorder=2),
    'AA_SASA_contacts':dict(color=CATEGORICAL[3],marker='s',linestyle='-',linewidth=1.0,zorder=3),
    'Static47_SASA':dict(color='#666666',marker='D',linestyle=':',linewidth=1.0,zorder=2),
    'Zero_change':dict(color=BLACK,marker='x',linestyle='--',linewidth=1.0,zorder=2),
}

# Scientific contract: selected-case review, not a mechanism claim. All 383
# exact-sequence positions are shown in (a). No smoothing or prediction-guided
# interval selection. Models use fold 2 predictions. No P values or error bars.
assert len(res)==3*383
assert set(res.dimer_target)=={'A','B','AB_mean'}
g=res[res.dimer_target=='AB_mean'].sort_values('sequence_position')
assert np.array_equal(g.sequence_position,np.arange(1,384))
assert g[['delta_md','delta_Full']].notna().all().all()
assert len(metrics[(metrics.scope=='all')])==15
assert len(blocks[(blocks.scope=='all')&(blocks.dimer_target=='AB_mean')])==25
source=[]

fig, axes=plt.subplots(2,2,figsize=(7.2,6.25))
fig.subplots_adjust(left=.095,right=.976,bottom=.218,top=.82,wspace=.32,hspace=.74)
fig.text(.095,.966,'3H6G: monomer-to-dimer changes across all 383 residues',
         fontsize=10,fontweight='bold',ha='left',va='top')
fig.text(.095,.936,'Same target sequence; fold 2 predictions; 500 ns per state',
         fontsize=8.5,ha='left',va='top',color='#555555')
for label,ax in zip('abcd',axes.flat):
    ax.text(-.16,1.16,label,transform=ax.transAxes,fontsize=10,fontweight='bold',va='top')
    ax.tick_params(length=3,pad=3)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=4))

ax=axes[0,0]
ax.set_title('All residues: mean of dimer chains A and B',loc='left',pad=21)
ax.axhline(0,color='#BBBBBB',linestyle='--',linewidth=.6,zorder=0)
ax.plot(g.sequence_position,g.delta_md,color='#B8B8B8',linewidth=1.0,label='MD',zorder=1)
ax.plot(g.sequence_position,g.delta_Full,color=CATEGORICAL[0],linewidth=1.0,label='Full',zorder=3)
ax.set(xlabel='Input sequence position',ylabel='Dimer − monomer ΔRMSF (Å)',xlim=(1,383))
ax.set_xticks([1,100,200,300,383])
ax.legend(loc='lower left',bbox_to_anchor=(0,1.005),ncol=2,borderaxespad=0,
          handlelength=1.8,columnspacing=1.2)
for _,r in g.iterrows():
    source.extend([dict(panel='a',series='MD',x=r.sequence_position,y=r.delta_md,n=383),
                   dict(panel='a',series='Full',x=r.sequence_position,y=r.delta_Full,n=383)])

ax=axes[0,1]
ax.set_title('Both dimer chains, then their mean',loc='left',pad=21)
targets=['A','B','AB_mean']
offsets=np.linspace(-.24,.24,len(MODELS))
for model,offset in zip(MODELS,offsets):
    d=metrics[(metrics.scope=='all')&(metrics.model==model)].set_index('dimer_target').loc[targets]
    style=STYLES[model]
    ax.plot(np.arange(3)+offset,d.delta_mae,linestyle='None',marker=style['marker'],
            color=style['color'],markersize=4.4,zorder=style['zorder'],markeredgewidth=.9)
    for x,(_,r) in zip(np.arange(3),d.iterrows()):
        source.append(dict(panel='b',series=model,x=targets[x],y=r.delta_mae,n=int(r.n)))
ax.set_xticks(np.arange(3),['Chain A','Chain B','A/B mean'])
ax.set(xlim=(-.5,2.5),ylim=(.16,.32),ylabel='ΔRMSF MAE (Å)',xlabel='Dimer target; same monomer reference')

ax=axes[1,0]
ax.set_title('Continuous time blocks, not independent repeats',loc='left',pad=16)
for model in MODELS:
    d=blocks[(blocks.scope=='all')&(blocks.dimer_target=='AB_mean')&(blocks.model==model)].sort_values('block')
    ax.plot(d.block,d.delta_mae,markersize=3.5,**STYLES[model])
    for _,r in d.iterrows():
        source.append(dict(panel='c',series=model,x=int(r.block),y=r.delta_mae,n=int(r.n)))
ax.set_xticks(range(1,6),['0–100','100–200','200–300','300–400','400–500'])
ax.tick_params(axis='x',labelsize=7.6)
ax.set(xlim=(.8,5.2),ylim=(.17,.43),xlabel='Time block (ns)',ylabel='ΔRMSF MAE (Å)')

ax=axes[1,1]
ax.set_title('Remove residues with the largest MD changes',loc='left',pad=16)
scopes=['all','exclude_top_abs_delta_0.05','exclude_top_abs_delta_0.1','exclude_top_abs_delta_0.2']
xs=[0,5,10,20]
for model in MODELS:
    d=metrics[(metrics.dimer_target=='AB_mean')&(metrics.model==model)].set_index('scope').loc[scopes]
    ax.plot(xs,d.delta_mae,markersize=3.5,**STYLES[model])
    for x,(_,r) in zip(xs,d.iterrows()):
        source.append(dict(panel='d',series=model,x=x,y=r.delta_mae,n=int(r.n)))
ax.set_xticks(xs)
ax.set(xlim=(-1,21),ylim=(.14,.31),xlabel='Largest |MD ΔRMSF| removed (%)',ylabel='ΔRMSF MAE (Å)')

handles=[Line2D([0],[0],markersize=4,**STYLES[m]) for m in MODELS]
fig.legend(handles,[LABELS[m] for m in MODELS],loc='lower center',bbox_to_anchor=(.54,.066),
           ncol=3,frameon=False,handlelength=2.1,columnspacing=1.7,labelspacing=.65,
           borderaxespad=0,fontsize=8.3)
fig.text(.095,.034,'MAE: mean absolute error across residues. A/B mean is calculated before the error.',
         ha='left',fontsize=8,color='#555555')
fig.text(.095,.012,'Panel d retains 383, 363, 344 and 306 residues, respectively. No P values or confidence intervals.',
         ha='left',fontsize=8,color='#555555')

fig.canvas.draw()
renderer=fig.canvas.get_renderer()
# Verify all text lies within the fixed canvas at the requested final width.
outside=[]
for text_artist in fig.findobj(mpl.text.Text):
    if not text_artist.get_visible() or not text_artist.get_text():
        continue
    bb=text_artist.get_window_extent(renderer=renderer)
    if bb.x0 < -1 or bb.y0 < -1 or bb.x1 > fig.bbox.width+1 or bb.y1 > fig.bbox.height+1:
        outside.append(text_artist.get_text())
assert not outside, f'Text outside fixed canvas: {outside}'
stem=HERE/'glur6_full_residue_review'
with mpl.rc_context({'savefig.bbox':None}):
    for ext in ['pdf','png','svg']:
        fig.savefig(stem.with_suffix('.'+ext),dpi=300,bbox_inches=None)
from PIL import Image
pixels=np.asarray(Image.open(stem.with_suffix('.png')).convert('RGB'))
densities={}
for label,ax in zip('abcd',axes.flat):
    bb=ax.get_window_extent(renderer=renderer)
    scale=300/fig.dpi
    x0,x1=int(bb.x0*scale),int(bb.x1*scale)
    y0,y1=pixels.shape[0]-int(bb.y1*scale),pixels.shape[0]-int(bb.y0*scale)
    crop=pixels[y0:y1,x0:x1]
    densities[label]=float(np.mean(np.any(crop<245,axis=2)))
assert min(densities.values())>.015
assert b'/FontFile2' in stem.with_suffix('.pdf').read_bytes()
assert '<text' in stem.with_suffix('.svg').read_text(encoding='utf-8')
pd.DataFrame(source).to_csv(HERE/'plotted_source_data.csv',index=False)
hashes={f:hashlib.sha256((DATA/f).read_bytes()).hexdigest() for f in FILES}
(HERE/'figure_data_audit.json').write_text(json.dumps(dict(
    source_sha256=hashes,width_inch=7.2,height_inch=6.25,dpi=300,
    panel_a_residues=383,panel_b_point_count=15,panel_c_point_count=25,panel_d_point_count=20,
    text_outside_canvas=outside,smoothed=False,confidence_intervals=False,
    panel_nonbackground_pixel_fraction=densities, pdf_truetype_font_embedded=True,
    svg_editable_text_present=True,
    scope='Selected case review. Continuous time blocks are not independent simulation repeats.',
    font='Arial with DejaVu Sans fallback',
    source='Existing audited per-residue and metric CSVs; raw trajectory not re-read.'
),indent=2),encoding='utf-8')
print(stem)
plt.close(fig)
