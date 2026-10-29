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




# Inherited audited full-profile plotting; data-led layout with no mechanism schematic.
from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from PIL import Image
OUT=Path(__file__).resolve().parent;ROOT=OUT.parent;P=OUT/'source_data';R=ROOT/'pymol'
W,H=183.,154.;BLUE='#2166AC';WARM='#B87936';KEY='#3A8D8C';PARTNER='#8B72A5';NGRAY='#8A9299';CGRAY='#BAC1C7';INK='#222222'
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
panel('a','Structural context of the regulatory segment',1,1.8)
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
 a=ax(x0,93,pw,16.5 if metric=='rmsf' else 50.5);a.axvspan(302.5,316.5,color='#EEF4F3',ec='none',zorder=0);ys=[]
 for name,prefix,color,marker,style in MODELS:
  values=df[prefix+'_'+metric].where(df['valid_'+metric]).to_numpy();ys.append(values)
  a.plot(positions,values,c=color,lw=1.1,ls=style,marker=marker,ms=2.7,zorder=2)
  for pos in [296,298,303,316]:
   q=df[df.CaMKI_native_position==pos].iloc[0]
   if q['valid_'+metric] and name=='Full':a.plot(pos,q[prefix+'_'+metric],marker=marker,ms=4,mfc=color,mec=INK,mew=.5,zorder=4)
  points.extend(dict(panel='b',metric=metric,position=int(pos),model=name,value=float(v)) for pos,v in zip(positions,values) if np.isfinite(v))
 a.set(xlim=(lo,hi),ylim=(0,top),ylabel='RMSF (Å)' if metric=='rmsf' else 'Hχ (nats)');a.set_yticks(ticks);a.tick_params(length=2,pad=2)
 a.set_xticks([294,298,303,308,313,318])
 if metric=='rmsf':a.tick_params(axis='x',labelbottom=False);a.set_yticks([0,4,8])
 else:a.set_xlabel('CaMKI residue number',labelpad=2)
 if metric=='chi_entropy':
  for pos,loc in [(297,(299.6,2.30)),(308,(310.3,1.55))]:
   q=df[df.CaMKI_native_position==pos].iloc[0];a.plot(pos,q.pred_chi_entropy,'o',ms=4.2,mfc=BLUE,mec=INK,mew=.5,zorder=4);a.annotate('N'+str(pos),xy=(pos,q.pred_chi_entropy),xytext=loc,fontsize=7.9,color=KEY,ha='center',arrowprops=dict(arrowstyle='-',lw=.5,color=INK))
 z=np.concatenate([v[np.isfinite(v)] for v in ys]);assert min(z)>=0 and max(z)<=top;checks.append(dict(metric=metric,ylim=[0,top],range=[float(z.min()),float(z.max())]))

# Replot all25 Full RMSF values on a magnified absolute y scale; not independent evidence.
txt(17,113,'Full: expanded y scale',7.9,color=BLUE)
a=ax(17,119,69,24.5);a.axvspan(302.5,316.5,color='#EEF4F3',ec='none',zorder=0)
a.plot(positions,df.pred_rmsf,c=BLUE,lw=1.1,marker='o',ms=2.7,zorder=2)
for pos,loc in [(296,(298.4,1.75)),(303,(303.7,1.05)),(316,(313.3,1.3))]:
 q=df[df.CaMKI_native_position==pos].iloc[0];a.plot(pos,q.pred_rmsf,'o',ms=4.2,mfc=BLUE,mec=INK,mew=.5,zorder=4)
 a.annotate(q.amino_acid+str(pos),xy=(pos,q.pred_rmsf),xytext=loc,fontsize=7.9,color=KEY,ha='center',arrowprops=dict(arrowstyle='-',lw=.5,color=INK))
a.set(xlim=(lo,hi),ylim=(0,2.3),ylabel='RMSF (Å)');a.set_yticks([0,1,2]);a.set_xticks([294,298,303,308,313,318]);a.set_xlabel('CaMKI residue number',labelpad=2);a.tick_params(length=2,pad=2)
# Frozen original values remain unchanged, including the repeated Full track.
pl=pd.DataFrame(points);assert len(pl)==92;pl.to_csv(P/'all_plotted_dynamics.csv',index=False)
zoom=pd.DataFrame({'position':positions,'model':'Full','metric':'rmsf','value':df.pred_rmsf,'display':'expanded_y_scale'});zoom.to_csv(P/'repeated_Full_RMSF_expanded.csv',index=False);assert len(zoom)==25
for text_obj in fig.findobj(mpl.text.Text):
 if text_obj.get_visible() and text_obj.get_text() and text_obj.get_fontsize()<7.9:text_obj.set_fontsize(7.9)
fig.canvas.draw();ren=fig.canvas.get_renderer();vis=[t for t in fig.findobj(mpl.text.Text) if t.get_visible() and t.get_text()]
outside=[t.get_text() for t in vis if t.get_window_extent(ren).x0<-.6 or t.get_window_extent(ren).y0<-.6 or t.get_window_extent(ren).x1>fig.bbox.width+.6 or t.get_window_extent(ren).y1>fig.bbox.height+.6];assert not outside,outside
for ext,dpi in [('png',450),('pdf',600),('svg',600)]:fig.savefig(OUT/('fig3_data_led_story.'+ext),dpi=dpi,bbox_inches=fig.bbox_inches)
qa=dict(canvas_mm=[W,H],source='family OOF fold1',panels=['a','b'],all25_RMSF_and21_Hchi_positions=True,unique_plotted_profile_values=len(pl),repeated_Full_RMSF_values=len(zoom),repeated_data_not_independent_evidence=True,all_quantitative_axes_start_at_zero=True,original_axes=checks,expanded_RMSF_ylim=[0,2.3],no_MD_in_main=True,no_numeric_tables=True,no_gray_explanatory_notes=True,no_mechanism_schematic=True,activity_not_measured=True,Intra_not_unbound=True,outside_text=outside,min_font=min(t.get_fontsize() for t in vis),visual_review='pending',sha256={ext:hashlib.sha256((OUT/('fig3_data_led_story.'+ext)).read_bytes()).hexdigest() for ext in ['png','pdf','svg']})
(OUT/'figure_QA.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(qa,ensure_ascii=False));plt.close(fig)
