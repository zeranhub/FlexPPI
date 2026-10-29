# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a) Contact map with paired model tracks -> PairedBoxScatter / LineTrend parameters -> param inherit
# (b) All45CDR paired curves -> LineTrend/plot_trend.py, companion inspected -> param inherit
# The contact grid has no matching production script; native data-dependent vector glyphs are used.
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
import json,hashlib
import numpy as np
import pandas as pd
from PIL import Image
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.text import Text

ROOT=Path(__file__).resolve().parent
BASE=ROOT.parent
OLD=BASE.parent/'mel5_joint_case_v49_20261019'
data=pd.read_csv(BASE/'analysis/local_HLA6_beta48_56_all54_pairs.csv',float_precision='round_trip')
hchi=pd.read_csv(BASE/'analysis/HLA6_Hchi_source.csv',float_precision='round_trip').sort_values('position_a')
cdr=pd.read_csv(OLD/'figure4/source_plotted_CDR_RMSF.csv',float_precision='round_trip')
beta=data.drop_duplicates('TCR_common_position').sort_values('TCR_common_position')
assert len(data)==54 and len(hchi)==6 and len(beta)==9 and len(cdr)==45
WT='#CF8C27'; ENG=CATEGORICAL[4]; TEAL='#297B82'; GRAY='#6B6B6B'
W,H=183,170
fig=plt.figure(figsize=(W/25.4,H/25.4),facecolor='white')
def axes(x,y,w,h):return fig.add_axes([x/W,1-(y+h)/H,w/W,h/H])
def label(x,y,text,**kw):return fig.text(x/W,1-y/H,text,**kw)
def heading(letter,title,y):
    label(3,y,letter,fontsize=11,fontweight='bold')
    label(11,y,title,fontsize=8.5)

heading('a','Residue contacts aligned with model dynamics',5)
state_handles=[Line2D([],[],color=WT,marker='o',mfc='white',lw=.9,ms=3,label='MEL5'),
               Line2D([],[],color=ENG,marker='s',mfc=ENG,lw=.9,ms=3,label='α24β17')]
fig.legend(handles=state_handles,ncol=2,loc='upper right',bbox_to_anchor=(.975,.994),fontsize=7.5,
           handlelength=1.35,handletextpad=.4,columnspacing=1.1,borderaxespad=0)

# The original whole-complex view is only a locator, occupying <5% canvas.
im=Image.open(OLD/'figure4/4jff_overview.png').convert('RGBA')
arr=np.asarray(im);yy,xx=np.where((arr[:,:,3]>8)&np.any(arr[:,:,:3]<242,axis=2))
im=im.crop((max(0,int(xx.min())-10),max(0,int(yy.min())-10),min(im.width,int(xx.max())+10),min(im.height,int(yy.max())+10)))
oa=axes(10,12,29,27);oa.imshow(im);oa.axis('off')
oa.text(.01,.93,'TCR',transform=oa.transAxes,fontsize=6.3,color=ENG)
oa.text(.04,.06,'HLA',transform=oa.transAxes,fontsize=6.3,color=TEAL)

row_positions=list(hchi.position_a.astype(int));column_positions=list(beta.TCR_common_position.astype(int))
row_labels=[f'{r.aa_a}{int(r.position_a)}' for r in hchi.itertuples()]
qrow=row_positions.index(72)

# Upper track: model RMSF aligned exactly with TCR columns below.
top=axes(60,16,115,23)
top.axvspan(1.5,6.5,color=ENG,alpha=.045,lw=0)
for field,color,marker,fill in [('TCR_RMSF_WT_A',WT,'o','white'),('TCR_RMSF_mature_A',ENG,'s',ENG)]:
    top.plot(range(9),beta[field],marker=marker,ms=3,color=color,mfc=fill,mec=color,lw=1.0)
top.set(xlim=(-.5,8.5),ylim=(.35,1.22),yticks=[.4,.8,1.2],xticks=[])
top.tick_params(length=2,pad=2)
top.set_ylabel('Pred. RMSF (Å)',fontsize=7,labelpad=4)
top.spines['bottom'].set_visible(False)
top.plot([1.55,6.45],[1.035,1.035],transform=top.get_xaxis_transform(),color=ENG,lw=.7,clip_on=False)
top.text(4,1.105,'CDR2β',transform=top.get_xaxis_transform(),ha='center',fontsize=7.2,color=ENG)

# Left track: all six local valid Hchi sites, in the same row order as contacts.
left=axes(12,45,30,56)
left.axhspan(qrow-.48,qrow+.48,color=TEAL,alpha=.075,lw=0)
for i,r in enumerate(hchi.itertuples()):
    left.plot([r.hchi_full_a,r.hchi_full_b],[i,i],color=TEAL if int(r.position_a)==72 else '#B8BEC2',lw=1.05)
    left.plot(r.hchi_full_a,i,'o',ms=3.2,mfc='white',mec=WT,mew=.9)
    left.plot(r.hchi_full_b,i,'s',ms=3,mfc=ENG,mec=ENG,mew=.8)
left.set(xlim=(0,3.2),ylim=(5.5,-.5),xticks=[0,1,2,3],yticks=[])
left.tick_params(length=2,pad=2)
left.spines['left'].set_visible(False)
left.set_xlabel('Pred. Hχ (nats)',fontsize=7,labelpad=3)
left.text(.5,1.055,'HLA side chain',transform=left.transAxes,ha='center',fontsize=7.3)

# Each cell is one real residue pair. The two markers are the two input states.
grid=axes(60,45,115,56)
grid.axhspan(qrow-.49,qrow+.49,color=TEAL,alpha=.075,lw=0,zorder=0)
grid.axvspan(1.5,6.5,color=ENG,alpha=.025,lw=0,zorder=0)
for boundary in np.arange(-.5,9,.999999999):grid.axvline(boundary,color='#E4E7E9',lw=.35,zorder=0)
for boundary in np.arange(-.5,6,.999999999):grid.axhline(boundary,color='#E4E7E9',lw=.35,zorder=0)
for r in data.itertuples():
    y=row_positions.index(int(r.HLA_common_position));x=column_positions.index(int(r.TCR_common_position))
    for state,dx,color,marker in [('WT',-.23,WT,'o'),('mature',.23,ENG,'s')]:
        d=float(getattr(r,state+'_minimum_heavy_distance_A'))
        if d<=4:ms,fc,ec,alpha=5,color,color,1
        elif d<=6:ms,fc,ec,alpha=4,'white',color,1
        else:ms,fc,ec,alpha=1.5,color,color,.42
        grid.plot(x+dx,y-.08,marker=marker,ms=ms,mfc=fc,mec=ec,mew=.8,alpha=alpha,zorder=3)
        if int(r.HLA_common_position)==72:
            grid.text(x+dx,y+.29,f'{d:.2f}',ha='center',va='center',fontsize=5.8,color=color,zorder=4)
for position in [53,54]:
    x=column_positions.index(position)
    grid.add_patch(Rectangle((x-.49,qrow-.48),.98,.96,facecolor='none',edgecolor=TEAL,lw=.8,zorder=2))
grid.set(xlim=(-.5,8.5),ylim=(5.5,-.5),yticks=range(6),yticklabels=row_labels,xticks=range(9))
column_labels=[]
for r in beta.itertuples():
    aas=r.TCR_WT_aa if r.TCR_WT_aa==r.TCR_mature_aa else r.TCR_WT_aa+'→'+r.TCR_mature_aa
    column_labels.append(f'{int(r.TCR_common_position)}\n{aas}')
grid.set_xticklabels(column_labels,fontsize=6.4)
grid.tick_params(length=0,pad=4)
for tick,pos in zip(grid.get_yticklabels(),row_positions):
    tick.set_fontsize(7.4)
    if pos==72:tick.set_color(TEAL);tick.set_fontweight('bold')
for spine in grid.spines.values():spine.set_visible(False)
grid.set_xlabel('TCR β residue · MEL5 → α24β17',fontsize=7.3,labelpad=4)

# Geometry legend: filled=contact; open=nearby; tiny marker=farther than6A.
leg=axes(17,115,156,7);leg.set(xlim=(0,1),ylim=(0,1));leg.axis('off')
leg.text(0,.5,'Nearest heavy-atom distance',va='center',fontsize=6.8)
for x,ms,fc,labeltext in [(.46,5,'#626A70','≤4 Å: contact'),(.67,4,'white','4–6 Å'),(.85,1.5,'#626A70','>6 Å')]:
    leg.plot(x,.5,'o',ms=ms,mfc=fc,mec='#626A70',mew=.8)
    leg.text(x+.025,.5,labeltext,va='center',fontsize=6.8)

heading('b','Model RMSF across all six recognition loops',130)
ax=axes(17,140,158,22)
groups=['CDR1α','CDR2α','CDR3α','CDR1β','CDR2β','CDR3β']
ticks=[];ticklabels=[]
for region in groups:
    q=cdr[cdr.region==region].sort_values('plot_x');x=q.plot_x.to_numpy()
    if region=='CDR2β':ax.axvspan(x[0]-.6,x[-1]+.6,color=ENG,alpha=.055,lw=0)
    for field,color,marker,fill in [('MEL5_model_RMSF_A',WT,'o','white'),('alpha24beta17_model_RMSF_A',ENG,'s',ENG)]:
        ax.plot(x,q[field],color=color,marker=marker,ms=2.5,mfc=fill,mec=color,lw=.95,zorder=3)
    ax.text(x.mean(),1.04,region,transform=ax.get_xaxis_transform(),ha='center',fontsize=7.2)
    ax.plot([x[0]-.2,x[-1]+.2],[1.01,1.01],transform=ax.get_xaxis_transform(),color='#B5B5B5',lw=.6,clip_on=False)
    choose=[0,len(q)-1] if region!='CDR3β' else [0,6,len(q)-1]
    for j in choose:ticks.append(x[j]);ticklabels.append(str(int(q.iloc[j].common_TCR_position)))
    if region=='CDR2β':
        g=q[q.common_TCR_position==54].iloc[0]
        ax.annotate('Gly54',xy=(g.plot_x,g.alpha24beta17_model_RMSF_A),xytext=(g.plot_x+2.3,1.0),fontsize=7,color=ENG,arrowprops=dict(arrowstyle='-',lw=.6,color=ENG))
ax.set(xlim=(-1,cdr.plot_x.max()+1.5),ylim=(0,1.42),ylabel='RMSF (Å)',xticks=ticks,yticks=[0,.4,.8,1.2])
ax.set_xticklabels(ticklabels);ax.tick_params(length=2,pad=2)
ax.set_xlabel('Common TCR residue numbering',fontsize=7.1,labelpad=2,color=GRAY)

fig.canvas.draw();renderer=fig.canvas.get_renderer();texts=[]
for item in fig.findobj(match=Text):
    if not item.get_visible() or not item.get_text():continue
    bb=item.get_window_extent(renderer).transformed(fig.transFigure.inverted())
    texts.append(dict(text=item.get_text(),font_pt=item.get_fontsize(),bounds=list(bb.extents),inside_canvas=bool(bb.x0>=-.001 and bb.y0>=-.001 and bb.x1<=1.001 and bb.y1<=1.001)))
assert min(t['font_pt'] for t in texts)>=5.8
assert all(t['inside_canvas'] for t in texts),[t for t in texts if not t['inside_canvas']]
for ext in ['png','pdf','svg']:
    with mpl.rc_context({'savefig.bbox':None}):fig.savefig(ROOT/f'fig4_mel5_contact_map.{ext}',dpi=450,bbox_inches=None)
plt.close(fig)
data.to_csv(ROOT/'source_all54_contacts.csv',index=False,encoding='utf-8-sig')
hchi.to_csv(ROOT/'source_all6_local_Hchi.csv',index=False,encoding='utf-8-sig')
beta.to_csv(ROOT/'source_all9_local_RMSF.csv',index=False,encoding='utf-8-sig')
cdr.to_csv(ROOT/'source_all45_CDR_RMSF.csv',index=False,encoding='utf-8-sig')
audit=dict(status='generated_pending_visual_QA',dimensions_mm=[W,H],raster_dpi=450,visible_text=texts,minimum_font_pt=min(t['font_pt'] for t in texts),local_pairs=54,HLA_positions=6,TCR_positions=9,complete_CDR_positions=45,geometry='minimum all-heavy atom residue-pair distance; contact<=4A; nearby4-6A; farther>6A',all_data_retained=True,causal_arrows=False,new_model_inference=0,new_training=0,new_MD=0)
(ROOT/'figure_data_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({k:audit[k] for k in ['status','dimensions_mm','local_pairs','minimum_font_pt']},ensure_ascii=False))

