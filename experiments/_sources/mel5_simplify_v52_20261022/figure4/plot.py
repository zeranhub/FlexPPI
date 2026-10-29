# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a) Existing complex locator and compact paired model readouts -> existing v49 image / direct vector text -> param inherit
# (b) Complete CDR curves -> LineTrend/plot_trend.py, companion inspected in prior turn -> param inherit
# User-directed deletion: both v51 local molecular images removed. No replacement structure or new data panel.
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
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
    "savefig.bbox": "tight",
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
from PIL import Image
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.text import Text

ROOT=Path(__file__).resolve().parent
OLD=ROOT.parent.parent/'mel5_joint_case_v49_20261019'
cdr=pd.read_csv(OLD/'figure4/source_plotted_CDR_RMSF.csv',float_precision='round_trip')
key=json.loads((OLD/'case_sources/residue_story/trajectory_link/key_findings.json').read_text(encoding='utf-8'))
assert len(cdr)==45
WT='#CF8C27'; ENG=CATEGORICAL[4]; TEAL='#297B82'; GRAY='#6B6B6B'

fig=plt.figure(figsize=(183/25.4,100/25.4),facecolor='white')
fig.text(.02,.955,'a',fontsize=11,fontweight='bold')
fig.text(.062,.955,'Model-predicted responses at the TCR–HLA interface',fontsize=8.5)

im=Image.open(OLD/'figure4/4jff_overview.png').convert('RGBA')
arr=np.asarray(im);yy,xx=np.where((arr[:,:,3]>8)&np.any(arr[:,:,:3]<242,axis=2))
crop=[max(0,int(xx.min())-20),max(0,int(yy.min())-20),min(im.width,int(xx.max())+20),min(im.height,int(yy.max())+20)]
im=im.crop(crop)
oa=fig.add_axes([.045,.565,.16,.335]);oa.imshow(im);oa.axis('off')
oa.text(.00,.94,'TCR α',transform=oa.transAxes,fontsize=6.5,color=GRAY)
oa.text(.67,.94,'TCR β',transform=oa.transAxes,fontsize=6.5,color=ENG)
fig.text(.052,.544,'HLA-A*02:01 · ELA peptide',fontsize=6.5,color=TEAL)

fig.text(.255,.852,'MEL5',fontsize=8.5,fontweight='bold',color=WT)
fig.text(.345,.852,'→',fontsize=10,color=GRAY)
fig.text(.38,.852,'α24β17',fontsize=8.5,fontweight='bold',color=ENG)
fig.text(.255,.775,'SPR affinity [28]',fontsize=7.5,fontweight='bold')
fig.text(.255,.713,'18 μM → 0.6 nM',fontsize=9)
fig.text(.255,.628,'Same HLA–peptide complex',fontsize=6.8,color=GRAY)

# The two model quantities have distinct units. Direct paired numbers avoid
# a shared scale, effect-size encoding, fake variance or a new filler panel.
fig.text(.565,.852,'HLA Q72 · side-chain entropy (Hχ)',fontsize=7.6,fontweight='bold')
fig.text(.565,.786,f"{key['model']['Q72_Hchi_nats_WT']:.3f}",fontsize=10,color=WT)
fig.text(.657,.786,'→',fontsize=10,color=GRAY)
fig.text(.696,.786,f"{key['model']['Q72_Hchi_nats_engineered']:.3f}",fontsize=10,color=ENG)
fig.text(.798,.786,'nats',fontsize=7.4,color=GRAY)
fig.text(.565,.669,'TCR βGly54 · backbone RMSF',fontsize=7.6,fontweight='bold')
fig.text(.565,.603,f"{key['model']['beta_Gly_RMSF_A_WT']:.3f}",fontsize=10,color=WT)
fig.text(.657,.603,'→',fontsize=10,color=GRAY)
fig.text(.696,.603,f"{key['model']['beta_Gly_RMSF_A_engineered']:.3f}",fontsize=10,color=ENG)
fig.text(.798,.603,'Å',fontsize=7.4,color=GRAY)

fig.text(.02,.458,'b',fontsize=11,fontweight='bold')
fig.text(.062,.458,'Predicted backbone response across all six recognition loops',fontsize=8.5)
ax=fig.add_axes([.095,.13,.866,.25])
groups=['CDR1α','CDR2α','CDR3α','CDR1β','CDR2β','CDR3β']
ticks=[];ticklabels=[]
for region in groups:
    q=cdr[cdr.region==region].sort_values('plot_x');x=q.plot_x.to_numpy()
    if region=='CDR2β':ax.axvspan(x[0]-.6,x[-1]+.6,color=ENG,alpha=.055,lw=0)
    for field,color,marker,fill in [('MEL5_model_RMSF_A',WT,'o','white'),('alpha24beta17_model_RMSF_A',ENG,'s',ENG)]:
        ax.plot(x,q[field],color=color,marker=marker,ms=2.7,mfc=fill,mec=color,lw=.95,zorder=3)
    ax.text(x.mean(),1.035,region,transform=ax.get_xaxis_transform(),ha='center',fontsize=7.5)
    ax.plot([x[0]-.2,x[-1]+.2],[1.01,1.01],transform=ax.get_xaxis_transform(),color='#B5B5B5',lw=.6,clip_on=False)
    choose=[0,len(q)-1] if region!='CDR3β' else [0,6,len(q)-1]
    for j in choose:ticks.append(x[j]);ticklabels.append(str(int(q.iloc[j].common_TCR_position)))
    if region=='CDR2β':
        g=q[q.common_TCR_position==54].iloc[0]
        ax.annotate('Gly54',xy=(g.plot_x,g.alpha24beta17_model_RMSF_A),xytext=(g.plot_x+2.3,1.0),fontsize=7.3,color=ENG,arrowprops=dict(arrowstyle='-',lw=.6,color=ENG))
ax.set(xlim=(-1,cdr.plot_x.max()+1.5),ylim=(0,1.42),ylabel='RMSF (Å)',xticks=ticks,yticks=[0,.4,.8,1.2])
ax.set_xticklabels(ticklabels);ax.tick_params(length=2.5,pad=3)
ax.set_xlabel('Common TCR residue numbering',labelpad=3,color=GRAY,fontsize=7.5)
handles=[Line2D([],[],color=WT,marker='o',mfc='white',ms=3.5,lw=1,label='MEL5'),Line2D([],[],color=ENG,marker='s',mfc=ENG,ms=3.3,lw=1,label='α24β17')]
fig.legend(handles=handles,loc='upper right',bbox_to_anchor=(.969,.482),ncol=2,fontsize=7.6,columnspacing=1.2,handlelength=1.4,handletextpad=.4,borderaxespad=0)

fig.canvas.draw();renderer=fig.canvas.get_renderer();texts=[]
for item in fig.findobj(match=Text):
    if not item.get_visible() or not item.get_text():continue
    bb=item.get_window_extent(renderer).transformed(fig.transFigure.inverted())
    texts.append(dict(text=item.get_text(),font_pt=item.get_fontsize(),bounds=list(bb.extents),inside_canvas=bool(bb.x0>=-.001 and bb.y0>=-.001 and bb.x1<=1.001 and bb.y1<=1.001)))
assert all(t['inside_canvas'] for t in texts),[t for t in texts if not t['inside_canvas']]
assert min(t['font_pt'] for t in texts)>=6.5
for ext in ['png','pdf','svg']:
    with mpl.rc_context({'savefig.bbox':None}):fig.savefig(ROOT/f'fig4_mel5_compact.{ext}',dpi=450,bbox_inches=None)
plt.close(fig)
cdr.to_csv(ROOT/'source_CDR_RMSF_all45.csv',index=False,encoding='utf-8-sig')
(ROOT/'source_model_readouts.json').write_text(json.dumps(key['model'],indent=2),encoding='utf-8')
audit=dict(status='generated_pending_visual_QA',dimensions_mm=[183,100],raster_dpi=450,minimum_font_pt=min(t['font_pt'] for t in texts),visible_text=texts,
           removed=['Both v51 local molecular close-ups and every atom/distance/sidechain label attached to them'],
           plotted_CDR_positions=45,unchanged_model_and_curve_values=True,source_model_SHA256='ff2974406d5d1ee07e6dc7118e43e57c0acc9849cc8feba7f5290476fca15ec5',new_inference=0,new_training=0,new_MD=0)
(ROOT/'figure_data_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({k:audit[k] for k in ['status','dimensions_mm','plotted_CDR_positions']},ensure_ascii=False))

