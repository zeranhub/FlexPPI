# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a) solution overview + real profiles -> cross-type inherit -> param inherit
# (b) molecular geometry schematic -> cross-type inherit -> param inherit
# (c) geometric attention schematic -> cross-type inherit -> param inherit
# (d) supervision workflow -> cross-type inherit -> param inherit
# RULE: "native run" = load pre-rendered PNG via Image.open().ax.imshow().
#       "param inherit" = drawing function below that copies Class A/B/C values.
#       If a panel says "native run" and you write a drawing function, you broke the contract.
# No matching framework template exists in the skill asset catalogue. Real
# molecular images are reused byte-exact; all arrows/text/geometry are vectors.
from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parent
RUNTIME = ROOT.parents[1] / 'completed_results_20261012' / 'plot_runtime'
sys.path.insert(0, str(RUNTIME))

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

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle, FancyArrowPatch, FancyBboxPatch, Polygon
from matplotlib.collections import LineCollection
from PIL import Image
import json, hashlib, re

S = ROOT / 'source_archive'
W, H = 183., 153.
FS = 9.1
BLUE = '#3F7093'
GOLD = '#AE8535'
PURPLE = '#80709B'
GREEN = '#4F7D6D'
INK = '#28363E'
MUTED = '#6C7C86'
LIGHT = '#D8E0E4'
BG_BLUE = '#EDF3F7'
BG_GOLD = '#F8F2E5'
BG_PURPLE = '#F1EEF6'
BG_GREEN = '#EDF4F0'
fig = plt.figure(figsize=(W/25.4,H/25.4), facecolor='white')
ax = fig.add_axes([0,0,1,1]); ax.set(xlim=(0,W), ylim=(H,0)); ax.axis('off')
TEXT = []

def text(x,y,s,sz=FS,color=INK,ha='left',va='center',weight='normal',**kw):
    t=ax.text(x,y,s,fontsize=sz,color=color,ha=ha,va=va,fontweight=weight,**kw)
    TEXT.append(t);return t

def line(x1,y1,x2,y2,c=LIGHT,lw=.6,ls='-',z=1):
    return ax.plot([x1,x2],[y1,y2],color=c,lw=lw,ls=ls,zorder=z)[0]

def arrow(x1,y1,x2,y2,c=MUTED,lw=.8,style='-|>',ms=7,rad=0,z=3):
    p=FancyArrowPatch((x1,y1),(x2,y2),arrowstyle=style,mutation_scale=ms,
                      linewidth=lw,color=c,connectionstyle=f'arc3,rad={rad}',zorder=z)
    ax.add_patch(p); return p

def rect(x,y,w,h,fill='white',edge=LIGHT,lw=.6,r=.7,z=1):
    p=FancyBboxPatch((x,y),w,h,boxstyle=f'round,pad=0,rounding_size={r}',
                    facecolor=fill,edgecolor=edge,linewidth=lw,zorder=z)
    ax.add_patch(p);return p

def dot(x,y,r=.7,color=PURPLE,edge='white',lw=.5,z=5):
    p=Circle((x,y),r,facecolor=color,edgecolor=edge,linewidth=lw,zorder=z);ax.add_patch(p);return p

def image(path,x,y,w,h,alpha=1.):
    im=Image.open(path).convert('RGBA');ax.imshow(im,extent=(x,x+w,y+h,y),zorder=1,interpolation='lanczos',alpha=alpha)

def panel(letter,title,x,y):
    text(x,y,letter,sz=11.5,weight='bold');text(x+5,y,title,sz=10.1,weight='bold')

def bars(x,y,width=8,height=5,color=BLUE,n=6):
    for k in range(n):
        ax.add_patch(Rectangle((x+k*width/n,y),width/n-.15,height,
                              facecolor=color,alpha=.20+.10*(k%4),edgecolor='none',zorder=2))

def matrix(x,y,w=13,h=11,c=GOLD,n=5):
    # Schematic pair-feature glyph, not a measured numerical heatmap.
    for i in range(n):
        for j in range(n):
            ax.add_patch(Rectangle((x+j*w/n,y+i*h/n),w/n-.18,h/n-.18,
                                  facecolor=c,alpha=.13+.06*((i*3+j*2)%5),edgecolor='none'))

def axes3(o,vecs,scale=4,color=PURPLE,labels=False):
    for i,v in enumerate(np.asarray(vecs).T):
        p=np.asarray(o)+scale*v
        arrow(*o,*p,c=color,lw=.8,ms=5.5,z=5)
        if labels:
            p=p+.7*v/max(np.linalg.norm(v),1e-9)
            text(*p,'xyz'[i],color=color,ha='center')
    dot(*o,r=.6,color=color)

# a. The overview owns the full width. Fine details are unpacked in b-d.
panel('a','From static structure to residue dynamics',3,5)
for xx,tt in [(17,'Static structure'),(53,'Geometric encoding'),(92,'Shared encoder'),(124,'Joint readout'),(160,'Residue dynamics')]:
    text(xx,13,tt,ha='center',sz=9.3,weight='bold')
image(S/'input_structure.png',3,20,29,13)
text(17,37,'4g9f · 840 residues',ha='center',sz=9.1)
text(17,44,'Sequence + N–Cα–C',ha='center',sz=9.1,color=MUTED)
arrow(32.5,29,37.5,29)

# Two complementary feature streams plus a fixed geometric reference.
for p in [(40,23),(44,20),(48,25),(43,29),(51,28)]:
    dot(*p,r=.85,color=BLUE)
for p,q in [((40,23),(44,20)),((44,20),(48,25)),((40,23),(43,29)),((43,29),(48,25)),((48,25),(51,28))]:
    line(*p,*q,c=BLUE,lw=.6)
bars(59,20.3,width=8,height=4,color=BLUE)
text(61,27.6,'Residue',ha='center',color=BLUE)
matrix(40,33,w=14,h=10,c=GOLD,n=5)
text(60,36,'Pair',ha='center',color=GOLD)
text(60,43,'Frames',ha='center',color=PURPLE)
axes3((44.2,46.4),np.array([[1,-.6,0],[0,.5,-1]]),scale=3.1,color=PURPLE)
arrow(68.5,27,73,27,c=BLUE)
arrow(68.5,36,73,36,c=GOLD)
arrow(68.5,43,73,43,c=PURPLE)

# Repetition is communicated by a compact bracket, not decorative card piles.
rect(75,19,33,29,fill='#FAFBFC',edge=LIGHT)
rect(78,22,27,8,fill=BG_PURPLE,edge='#CABED7')
text(91.5,26,'Fixed-frame IPA',ha='center',color=PURPLE)
arrow(91.5,30.4,91.5,33.5,lw=.65,ms=5)
rect(78,34,27,7,fill=BG_BLUE,edge='#B6C9D6')
text(91.5,37.5,'Feed-forward',ha='center',color=BLUE)
text(91.5,45,'4 blocks · 8 heads',ha='center',sz=9.1,color=MUTED)
# Pre-normalisation and residual connections are defined in the method; this
# overview shows the two sequential transformations without an ambiguous loop.
arrow(108.7,32.5,112,32.5)
rect(113,22,22,23,fill=BG_GOLD,edge='#D1B775')
text(124,27,'Shared MLP',ha='center')
line(115,31,133,31,c='#D8C596',lw=.5)
text(124,35,'μ, log σ²',ha='center')
text(124,41,'RMSF · Hχ',ha='center',color=GOLD)
arrow(135.5,32.5,140,32.5)

# Real PDB rendering and all 840 raw predictions, never synthetic profiles.
image(S/'flexibility_envelope.png',143,17,36.5,16.1)
# Inherit the original numerical color mapping without recolouring the protein.
scale_source=json.loads((S/'rmsf_original_scale.json').read_text(encoding='utf8'))
for k,color in enumerate(scale_source['colors']):
    ax.add_patch(Rectangle((146+k*30.5/64,33.9),30.5/64+.005,1.3,facecolor=color,edgecolor='none'))
text(145,34.5,'0',sz=9.1,ha='right',color=MUTED)
text(178,34.5,'4 Å',sz=9.1,ha='left',color=MUTED)

df=pd.read_csv(S/'curve_data_all_840.csv')
local=pd.read_csv(S/'curve_local_means.csv')
assert len(df)==1680 and len(local)==840
assert set(df.global_index)==set(range(840))
assert sum((df.metric=='rmsf') & df.valid.astype(bool))==840
assert sum((df.metric=='hchi') & df.valid.astype(bool))==734
# Both physical quantities use the original raw units and their full fixed range.
pa=fig.add_axes([146/W,(H-51.5)/H,30.5/W,10.2/H])
for metric,color,maxy in [('rmsf',BLUE,4.),('hchi',GOLD,4.1)]:
    for chain in 'ABCDE':
        z=df[(df.metric==metric)&(df.chain==chain)].sort_values('global_index')
        vals=z.prediction.where(z.valid.astype(bool)).to_numpy(float)
        x=z.global_index.to_numpy()+1
        pa.plot(x,vals/maxy*4,color=color,lw=.32,alpha=.23)
        q=local[local.chain==chain]
        pa.plot(q.global_index+1,q[metric+'_local_mean']/maxy*4,color=color,lw=.9)
pa.set(xlim=(1,840),ylim=(0,4),xticks=[1,420,840],yticks=[0,4])
pa.tick_params(axis='x',labelsize=9.1,length=1.8,pad=1.4,width=.5)
pa.tick_params(axis='y',labelsize=9.1,length=1.4,pad=1.2,width=.5,colors=BLUE)
pa.spines['left'].set_visible(False);pa.spines['bottom'].set_color(LIGHT)
right=pa.twinx();right.set(ylim=(0,4.1),yticks=[0,4.1],yticklabels=['0','4.1']);right.tick_params(axis='y',labelsize=9.1,length=1.4,pad=1.2,width=.5,colors=GOLD)
right.spines['left'].set_visible(False);right.spines['bottom'].set_visible(False);right.spines['top'].set_visible(False);right.spines['right'].set_visible(False)
text(146,38.3,'RMSF (Å)',color=BLUE,sz=9.1)
text(176,38.3,'Hχ (nats)',color=GOLD,sz=9.1,ha='right')
# Profile line labels explicitly state their display, rather than imply MD data.
text(161.5,56.7,'Residue index',ha='center',color=MUTED,sz=9.1)
line(3,59,180,59,lw=.55)

# b. Two measured residue frames, one readable backbone construction.
panel('b','Backbone-anchored local frames',3,65)
text(3,71,'One chain, two residue frames',color=MUTED,sz=9.1)
protein=(3.5,78.,48.,36.35)
image(S/'whole_chain.png',*protein,alpha=.82)
projection=json.loads((S/'projection_receipt.json').read_text(encoding='utf8'))
allframes=json.loads((S/'local_frame_candidates.json').read_text(encoding='utf8'))['all_204_local_frames']
frame_by={x['residue_number']:x for x in allframes}
P=np.array(projection['camera_nominal_projection_px_per_A'])
camera=np.array(projection['PyMOL_rotation_convention']['camera_matrix'])
proj=camera[:2].copy();proj[1]*=-1
def screen(number):
    px=np.array(projection['pixels'][f'D{number}_origin'])
    return np.array(protein[:2])+px*np.array([protein[2]/1071,protein[3]/811])
for n,label in [(67,'j'),(166,'i')]:
    pos=screen(n);R=np.array(frame_by[n]['R_columns_xyz'])
    axes3(pos,proj@R,scale=4.6,color=PURPLE)
    import matplotlib.patheffects as pe
    if n==67:t=text(pos[0]-6,pos[1]+.5,label,sz=10,color=PURPLE,weight='bold')
    else:t=text(pos[0]+5,pos[1]-3,label,sz=10,color=PURPLE,weight='bold')
    t.set_path_effects([pe.withStroke(linewidth=2.5,foreground='white')])
# No global triad is superimposed on the protein. The shared coordinate mapping
# is unpacked on the right, at the same specific frame i = D:Met166.
pos=screen(166)
ax.add_patch(Circle(pos,2.7,fill=False,edgecolor=PURPLE,lw=.7,zorder=4))
line(pos[0]+2.7,pos[1],60,101,c='#BBB1CA',lw=.6,ls=(0,(2,2)))
line(60,101,60,80,c='#BBB1CA',lw=.6,ls=(0,(2,2)))
text(62,73,'Frame i',weight='bold',color=PURPLE)
# The local view lies in the real N–CA–C plane. These displayed atom coordinates
# use actual D:Met166 atom lengths/angle, not idealised bond geometry.
fi=frame_by[166];ca=np.array(fi['CA_world_A']);R=np.array(fi['R_columns_xyz'])
ca_view=np.array([70.5,94.]);sc=7.0
nv=R.T@(np.array(fi['N_world_A'])-ca)
cv=R.T@(np.array(fi['C_world_A'])-ca)
N=ca_view+sc*nv[:2]*np.array([1,-1]);C=ca_view+sc*cv[:2]*np.array([1,-1])
line(*N,*ca_view,c='#536672',lw=1.5);line(*ca_view,*C,c='#536672',lw=1.5)
for pt,label,dx,dy in [(N,'N',-2,-2),(ca_view,'Cα',-2,4),(C,'C',2,-2)]:
    dot(*pt,r=1.35,color='#CBD5DB',edge='white',lw=.6)
    text(pt[0]+dx,pt[1]+dy,label,ha='center')
arrow(*ca_view,88,94,c=PURPLE,lw=.75,ms=6)
arrow(*ca_view,70.5,78,c=PURPLE,lw=.75,ms=6)
text(89,94,'x',color=PURPLE)
text(71,77,'y',color=PURPLE)
dot(70.5,94,r=.45,color=PURPLE,edge='none')
# Dashed perpendicular projection identifies Gram–Schmidt rather than N-axis.
line(*N,70.5,N[1],c=PURPLE,lw=.55,ls=(0,(2,2)))
line(70.5,92,72.5,92,c=PURPLE,lw=.5);line(72.5,92,72.5,94,c=PURPLE,lw=.5)
# Out-of-plane z is a circled dot in the fixed atomic plane.
ax.add_patch(Circle((84.8,82.8),1.2,fill=False,edgecolor=PURPLE,lw=.7))
dot(84.8,82.8,r=.28,color=PURPLE,edge='none');text(87,82.8,'z',color=PURPLE)
q=np.array([79.5,85.5]);line(*ca_view,*q,c=PURPLE,lw=.65,ls=(0,(2,2)))
dot(*q,r=1.1,color=PURPLE);text(81.5,86,'qᵢ',color=PURPLE)
text(65,103.8,'Local → global',color=MUTED)
text(62,110.5,'pᵢ = Rᵢqᵢ + tᵢ',sz=9.5)
text(3,117,'Fixed from N, Cα and C',color=MUTED,sz=9.1)
line(96,63,96,118,lw=.5)

# c. Explicit mathematical contributions, rather than three empty text cards.
panel('c','Geometry-guided attention',101,65)
# Three glyphs show what is actually combined in the attention logit.
bars(103,74,width=7,height=2.4,color=BLUE,n=5)
bars(103,78.2,width=7,height=2.4,color=BLUE,n=5)
text(114,77.3,'q·k',color=BLUE,sz=10)
matrix(129,73.8,w=10,h=8.5,c=GOLD,n=4)
text(141,78,'bᵢⱼ',color=GOLD,sz=10)
# Learned Q/K points in one common global frame. Their coordinates are
# schematic latent points; these are not molecular atoms or predicted atoms.
line(155,81,177,81,c=LIGHT,lw=.5)
line(157,83,157,71,c=LIGHT,lw=.5)
dot(161,77.5,r=1,color=PURPLE);dot(173,73.5,r=1,color=PURPLE)
line(161,77.5,173,73.5,c=PURPLE,lw=.8,ls=(0,(2,2)))
text(160,72,'Qᵢ',color=PURPLE,sz=9.1)
text(173.5,78.5,'Kⱼ',color=PURPLE,sz=9.1)
for xx,tt,cc in [(110,'Scalar',BLUE),(136,'Pair',GOLD),(167,'Point distance',PURPLE)]:
    text(xx,85.8,tt,ha='center',color=cc,sz=9.1)
arrow(111,87.8,124,88.8,c=BLUE,lw=.7,ms=6)
arrow(136,87.8,136,88.8,c=GOLD,lw=.7,ms=6)
arrow(166,87.8,148,88.8,c=PURPLE,lw=.7,ms=6)
text(136,91,'Scalar + pair − distance penalty',ha='center',sz=9.1)
arrow(136,92.8,136,94.3,lw=.7,ms=5)
rect(122,94.4,28,6,fill=BG_GREEN,edge='#ADC6BB')
text(136,97.4,'Softmax weights',ha='center',color=GREEN,sz=9.1)
# One weight set is used for all three aggregation branches.
for xx,cc in [(110,BLUE),(136,GOLD),(164,PURPLE)]:
    arrow(136,100.5,xx,102.3,c=cc,lw=.65,ms=5)
bars(102,103,width=10,height=2.5,color=BLUE)
matrix(131,102.6,w=9,h=3.3,c=GOLD,n=3)
for xx,yy in [(158,104),(163,102.6),(166.5,105)]:dot(xx,yy,r=.65,color=PURPLE)
text(109,109,'Values',ha='center',color=BLUE,sz=9.1)
text(136,109,'Pair',ha='center',color=GOLD,sz=9.1)
text(165,109,'Points → i',ha='center',color=PURPLE,sz=9.1)
for xx,cc in [(110,BLUE),(136,GOLD),(164,PURPLE)]:
    arrow(xx,111,136,113,c=cc,lw=.55,ms=4)
rect(116,113,41,5.4,fill=BG_BLUE,edge='#B6C9D6')
text(136.5,115.7,'Merge + project → update',ha='center',color=BLUE,sz=9.1)
line(3,121,180,121,lw=.55)

# d. MD summaries and predictive parameters enter the same probabilistic loss.
panel('d','Learning from MD ensembles',3,127)
image(S/'d_snapshot_first.png',3,132,16,7.66)
image(S/'d_snapshot_middle.png',19.5,132,16,7.66)
image(S/'d_snapshot_last.png',36,132,16,7.66)
text(27.5,142,'MD snapshots',color=GREEN,ha='center',sz=9.1)
arrow(53,136,59,136,c=GREEN,lw=.8,ms=6)
text(65,133.5,'MD mean y',color=GREEN,sz=9.1)
text(65,139,'Variance v(MD)',color=GREEN,sz=9.1)
# Shared model parameterisation is already detailed in a; training shows only
# its two meaningful outputs and their role in the observation-aware NLL.
rect(3,145.3,49,5.8,fill=BG_GOLD,edge='#D1B775')
text(27.5,148.2,'Joint readout',ha='center',sz=9.1)
arrow(53,148,59,148,c=BLUE,lw=.8,ms=6)
text(65,148,'Model μ, variance σ²',color=BLUE,sz=9.1)
line(113,132,115,132,c=GREEN,lw=.65)
line(115,132,115,141,c=GREEN,lw=.65)
line(113,141,115,141,c=GREEN,lw=.65)
arrow(115,136.5,137,136.5,c=GREEN,lw=.8,ms=6)
arrow(119,148,137,148,c=BLUE,lw=.8,ms=6)
rect(138,131,42,20,fill=BG_GREEN,edge='#ADC6BB')
text(159,135.5,'Gaussian NLL',ha='center',weight='bold',sz=9.3)
text(159,141.5,'μ vs. MD mean',ha='center',sz=9.1)
text(159,147.5,'v = σ² + v(MD)',ha='center',sz=9.1)

# Render exact physical dimensions. No bbox trim changes the export size.
fig.canvas.draw()
renderer=fig.canvas.get_renderer()
alltexts=[o for o in fig.findobj(mpl.text.Text) if o.get_visible() and o.get_text()]
visible=[o.get_text() for o in alltexts]
assert not [s for s in visible if re.search('[\u3400-\u9fff]',s)]
outside=[]
for t in alltexts:
    bb=t.get_window_extent(renderer)
    if bb.x0<-.3 or bb.y0<-.3 or bb.x1>fig.bbox.width+.3 or bb.y1>fig.bbox.height+.3:
        outside.append(t.get_text())
assert not outside,outside
for ext in ['png','pdf','svg']:
    fig.savefig(ROOT/('fig1_framework.'+ext),dpi=450,bbox_inches=None,pad_inches=0)
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
(ROOT/'visible_text.json').write_text(json.dumps(visible,ensure_ascii=False,indent=2),encoding='utf8')
audit={
    'status':'rendered_pending_visual_QA','width_mm':W,'height_mm':H,
    'word_width_mm':167.64,'word_scale':167.64/W,
    'min_source_font_pt':min(t.get_fontsize() for t in alltexts),
    'min_word_font_pt':min(t.get_fontsize() for t in alltexts)*167.64/W,
    'out_of_canvas_text':outside,'visible_text_language':'English',
    'source_rows':1680,'positions':840,'rmsf_valid':840,'hchi_valid':734,
    'all_raw_predictions_plotted':True,'profile_smoothing':'11-site centered same-chain valid-position mean; invalid central Hchi remains absent',
    'profile_display':'full original fixed ranges RMSF 0..4 A, Hchi 0..4.1 nats; separate metric values linearly mapped to the shared compact plotting height; explicitly labelled as predictions',
    'b_frame_indices':[67,166],'b_chain':'D','b_true_backbone_atoms':True,
    'b_local_point_q':'schematic learned latent point, not an atom',
    'c_glyphs':'schematic operations, not observed attention weights',
    'c_frames_fixed':True,'new_training':0,'new_MD':0,'new_inference':0,
    'source_hashes':{p.name:sha(p) for p in sorted(S.glob('*')) if p.is_file()},
    'output_sha256':{ext:sha(ROOT/('fig1_framework.'+ext)) for ext in ['png','pdf','svg']}
}
(ROOT/'render_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps({k:audit[k] for k in ['status','min_source_font_pt','min_word_font_pt','width_mm','height_mm']},ensure_ascii=False))
