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




# Inherited audited full-profile plotting; native vector mechanism schematic.
from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle,Ellipse,FancyArrowPatch,PathPatch,Circle
from matplotlib.path import Path as MPath
from PIL import Image
OUT=Path(__file__).resolve().parent;ROOT=OUT.parent;P=OUT/'source_data';R=ROOT/'pymol_mechanism'
W,H=183.,178.;BLUE='#2166AC';WARM='#B87936';KEY='#3A8D8C';PARTNER='#8B72A5';NGRAY='#8A9299';CGRAY='#BAC1C7';INK='#222222'
MODELS=[('Full','pred',BLUE,'o','-'),('Intra','intra',WARM,'s','--')]
mpl.rcParams.update({'savefig.bbox':None,'axes.unicode_minus':False,'axes.labelsize':8.2,'xtick.labelsize':7.9,'ytick.labelsize':7.9})
df=pd.read_csv(P/'all25_dynamics_and_static_lobe_contacts.csv').sort_values('CaMKI_native_position');assert len(df)==25 and df.valid_chi_entropy.sum()==21
manifest=json.loads((R/'render_manifest.json').read_text(encoding='utf-8'))
fig=plt.figure(figsize=(W/25.4,H/25.4),facecolor='white');points=[];checks=[]
def txt(x,y,s,fs=8.2,**kw):return fig.text(x/W,1-y/H,s,fontsize=fs,va=kw.pop('va','top'),color=kw.pop('color',INK),**kw)
def ax(x,y,w,h):return fig.add_axes([x/W,1-(y+h)/H,w/W,h/H])
def panel(lab,title,x,y):txt(x,y,lab,11,fontweight='bold');txt(x+6,y,title,8.8,fontweight='bold')
def line(x1,y1,x2,y2,**kw):
 fig.add_artist(plt.Line2D([x1/W,x2/W],[1-y1/H,1-y2/H],transform=fig.transFigure,**kw))
def arr(x1,y1,x2,y2,**kw):
 patch=FancyArrowPatch((x1/W,1-y1/H),(x2/W,1-y2/H),transform=fig.transFigure,arrowstyle=kw.pop('arrowstyle','->'),mutation_scale=8,lw=kw.pop('lw',.8),color=kw.pop('color',INK),**kw);fig.add_artist(patch);return patch
panel('a','From inhibitory contacts to CaM capture',1,1.8)
txt(24,9.5,'Autoinhibited CaMKI',8.3,fontweight='bold',ha='center');txt(24,14,'1A06',7.9,ha='center')
txt(137,9.5,'CaM-bound regulatory peptide',8.3,fontweight='bold',ha='center');txt(137,14,'1MXE',7.9,ha='center')
# Images are observed crystal structures, each fitted into its frame. The rotation is identical.
a=ax(3,19,41,47);info=manifest['assets']['1A06_context'];a.imshow(Image.open(R/info['trim_file']));a.axis('off')
for key,loc in [('F298',(.96,.63)),('K296',(.09,.43))]:
 a.annotate(key,xy=info['trim_landmarks'][key],xytext=loc,textcoords='axes fraction',ha='center',va='center',fontsize=7.9,color=KEY,arrowprops=dict(arrowstyle='-',color=INK,lw=.55,shrinkA=2,shrinkB=1),bbox=dict(fc='white',ec='none',alpha=.9,pad=.3))
txt(24,66.1,'Catalytic core',7.9,ha='center')
# A local F298 pocket enlargement sits immediately right of the full kinase.
if 'F298_pocket' in manifest['assets']:
 info2=manifest['assets']['F298_pocket'];a2=ax(49,24,35,30);a2.imshow(Image.open(R/info2['trim_file']));a2.axis('off')
 for key,loc in [('F298',(.50,1.03)),('G214',(.96,.16)),('I210',(.92,.71))]:
  if key in info2['trim_landmarks']:a2.annotate(key,xy=info2['trim_landmarks'][key],xytext=loc,textcoords='axes fraction',ha='center',va='center',fontsize=7.9,color=KEY if key=='F298' else INK,arrowprops=dict(arrowstyle='-',lw=.5,color=INK),bbox=dict(fc='white',ec='none',alpha=.9,pad=.2))
else:
 a2=ax(49,24,35,30);a2.axis('off');a2.text(.5,.5,'F298 pocket',ha='center',va='center',fontsize=8)
line(40,38,48,31,color='#767F87',lw=.5);line(40,43,48,51,color='#767F87',lw=.5)
txt(67,57.4,'F298–kinase contacts',7.9,ha='center')
a=ax(102,19,71,48);info=manifest['assets']['1MXE_context'];a.imshow(Image.open(R/info['trim_file']));a.axis('off')
for key,loc in [('K296',(.35,1.03)),('F298',(.09,.80)),('W303',(-.01,.48)),('M316',(1.02,.27))]:
 a.annotate(key,xy=info['trim_landmarks'][key],xytext=loc,textcoords='axes fraction',ha='center',va='center',fontsize=7.9,color=KEY,arrowprops=dict(arrowstyle='-',color=INK,lw=.55,shrinkA=2,shrinkB=1),bbox=dict(fc='white',ec='none',alpha=.92,pad=.3))
txt(111,64.4,'C-lobe',7.9,ha='center');txt(164,64.4,'N-lobe',7.9,ha='center')
# Component colours remain distinct from Full/Intra throughout.
line(52,66.9,56,66.9,color=PARTNER,lw=3);txt(58,65.5,'CaMKI regulatory segment',7.6)
line(52,71,56,71,color=KEY,lw=3);txt(58,69.6,'Functional residues',7.6)
# Two complete absolute profiles with the same native residue axis and functional markings.
panel('b','Predicted backbone restraint and side-chain heterogeneity',1,77)
for x,(name,prefix,color,marker,style) in zip([139,162],MODELS):
 line(x-4,79,x+1,79,color=color,marker=marker,ms=3.3,lw=1.1,ls=style);txt(x+3,77.4,name,8)
positions=df.CaMKI_native_position.to_numpy();lo,hi=293.5,318.5
for metric,x0,pw,top,ticks in [('rmsf',17.,69.,8.2,[0,2,4,6,8]),('chi_entropy',109.,69.,4.3,[0,1,2,3,4])]:
 coreleft=x0+(302.5-lo)/(hi-lo)*pw;corewidth=14/(hi-lo)*pw
 # Two tracks: the specific overlapping residues, then the span of the experimentally defined anchors.
 txt(x0+2.5/25*pw,85,'K296/F298',7.9,ha='center',color=KEY)
 txt(coreleft+corewidth/2,85,'W303–M316 core',7.9,ha='center',color=KEY)
 line(coreleft,89,coreleft+corewidth,89,color=KEY,lw=1.8)
 for pos in [296,298]:
  xx=x0+(pos-lo)/(hi-lo)*pw;line(xx,88.2,xx,90.1,color=KEY,lw=.9)
 a=ax(x0,93,pw,29.5);a.axvspan(302.5,316.5,color='#EEF4F3',ec='none',zorder=0);ys=[]
 for name,prefix,color,marker,style in MODELS:
  values=df[prefix+'_'+metric].where(df['valid_'+metric]).to_numpy();ys.append(values)
  a.plot(positions,values,c=color,lw=1.1,ls=style,marker=marker,ms=2.7,zorder=2)
  for pos in [296,298,303,316]:
   q=df[df.CaMKI_native_position==pos].iloc[0]
   if q['valid_'+metric] and name=='Full':a.plot(pos,q[prefix+'_'+metric],marker=marker,ms=4,mfc=color,mec=INK,mew=.5,zorder=4)
  points.extend(dict(panel='b',metric=metric,position=int(pos),model=name,value=float(v)) for pos,v in zip(positions,values) if np.isfinite(v))
 a.set(xlim=(lo,hi),ylim=(0,top),ylabel='RMSF (Å)' if metric=='rmsf' else 'Hχ (nats)');a.set_yticks(ticks);a.tick_params(length=2,pad=2)
 a.set_xticks([294,298,303,308,313,318]);a.set_xlabel('CaMKI residue number',labelpad=2)
 if metric=='chi_entropy':
  for pos,loc in [(297,(299.6,2.30)),(308,(310.3,1.55))]:
   q=df[df.CaMKI_native_position==pos].iloc[0];a.plot(pos,q.pred_chi_entropy,'o',ms=4.2,mfc=BLUE,mec=INK,mew=.5,zorder=4);a.annotate('N'+str(pos),xy=(pos,q.pred_chi_entropy),xytext=loc,fontsize=7.9,color=KEY,ha='center',arrowprops=dict(arrowstyle='-',lw=.5,color=INK))
 z=np.concatenate([v[np.isfinite(v)] for v in ys]);assert min(z)>=0 and max(z)<=top;checks.append(dict(metric=metric,ylim=[0,top],range=[float(z.min()),float(z.max())]))
# A functional, falsifiable interpretation. This native vector drawing is a schematic, not a predicted structure.
panel('c','Hypothesis: dual-anchor restraint may limit inhibitory re-contact',1,137)
s=ax(2,146,121,28);s.set(xlim=(0,121),ylim=(0,28));s.axis('off')
# Catalytic core represented by two neutral lobes. White cleft marks a possible regulatory contact site.
s.add_patch(Ellipse((16,14),29,25,fc='#D3D7DB',ec='#8A9299',lw=.7));s.add_patch(Ellipse((35,14),25,23,fc='#C2C8CE',ec='#8A9299',lw=.7))
s.plot([31,37,43],[19,14,18],color='white',lw=6,solid_capstyle='round');s.text(18,1,'Catalytic core',fontsize=7.9,ha='center',va='bottom')
# CaM lobes capture the anchor span while the N-terminal overlap remains outside the catalytic core.
s.add_patch(Ellipse((80,17),34,18,fc=CGRAY,ec='#858D94',lw=.6));s.add_patch(Ellipse((107,17),23,18,fc=NGRAY,ec='#687078',lw=.6))
s.text(80,24,'C-lobe',fontsize=7.9,ha='center',va='center');s.text(107,24,'N-lobe',fontsize=7.9,ha='center',va='center')
s.text(80,19,'Tight binding',fontsize=7.9,ha='center',va='center');s.text(107,19,'Required for',fontsize=7.9,ha='center',va='center');s.text(107,15.3,'activation',fontsize=7.9,ha='center',va='center')
s.plot([58,67],[12,16],color='#858D94',lw=.8)
# The purple line is a sequence/constraint schematic, not an atomistic conformation.
verts=[(49,9),(55,14),(59,9),(64,14),(70,12),(79,12),(85,12),(96,12),(105,12)]
codes=[MPath.MOVETO]+[MPath.CURVE3]*8
s.add_patch(PathPatch(MPath(verts,codes),fc='none',ec=PARTNER,lw=3.5,capstyle='round'))
for x,t in [(80,'W303'),(105,'M316')]:
 s.add_patch(Circle((x,12),1.25,fc=KEY,ec='white',lw=.4));s.text(x,6,t,fontsize=7.9,ha='center',va='center',color=KEY)
s.text(53,25,'K296/F298',fontsize=7.9,ha='center',color=KEY);s.annotate('',xy=(56,13),xytext=(54,22),arrowprops=dict(arrowstyle='-',lw=.55,color=INK))
s.annotate('',xy=(38,15),xytext=(56,12),arrowprops=dict(arrowstyle='->',connectionstyle='arc3,rad=.35',ls=(0,(3,2)),lw=.8,color=INK))
s.text(51,.4,'Possible\nre-contact',fontsize=7.9,ha='center',va='bottom')
# A test that distinguishes activation from simply losing peptide/CaM binding.
txt(128,147,'Testable prediction',8.2,fontweight='bold')
txt(128,153,'Matched CaM occupancy',7.9)
txt(128,158,'Weaken N-lobe–M316 contact',7.9)
arr(128,163,128,166,color=KEY,lw=.8)
txt(132,164,'Core restraint ↓',7.9)
txt(132,169,'Kinase activity ↓',7.9)
# All plotted values preserved: 25*2 RMSF plus 21*2 Hchi.
pl=pd.DataFrame(points);assert len(pl)==92;pl.to_csv(P/'all_plotted_dynamics.csv',index=False)
for text_obj in fig.findobj(mpl.text.Text):
 if text_obj.get_visible() and text_obj.get_text() and text_obj.get_fontsize()<7.9:text_obj.set_fontsize(7.9)
fig.canvas.draw();ren=fig.canvas.get_renderer();vis=[t for t in fig.findobj(mpl.text.Text) if t.get_visible() and t.get_text()]
outside=[t.get_text() for t in vis if t.get_window_extent(ren).x0<-.6 or t.get_window_extent(ren).y0<-.6 or t.get_window_extent(ren).x1>fig.bbox.width+.6 or t.get_window_extent(ren).y1>fig.bbox.height+.6];assert not outside,outside
for ext,dpi in [('png',450),('pdf',600),('svg',600)]:fig.savefig(OUT/('fig3_mechanism_story.'+ext),dpi=dpi,bbox_inches=fig.bbox_inches)
qa=dict(canvas_mm=[W,H],source='family OOF fold1',panels=['a','b','c'],all25_RMSF_and21_Hchi_positions=True,plotted_profile_values=len(pl),all_quantitative_axes_start_at_zero=True,axes=checks,no_MD_in_main=True,no_numeric_tables=True,no_gray_explanatory_notes=True,mechanism_explicitly_hypothesis=True,activity_not_measured=True,Intra_not_unbound=True,outside_text=outside,min_font=min(t.get_fontsize() for t in vis),visual_review='pending',sha256={ext:hashlib.sha256((OUT/('fig3_mechanism_story.'+ext)).read_bytes()).hexdigest() for ext in ['png','pdf','svg']})
(OUT/'figure_QA.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(qa,ensure_ascii=False));plt.close(fig)
