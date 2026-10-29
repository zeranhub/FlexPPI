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



# Two panels only: structural location and complete model profiles.
# Numeric profile logic inherited from the audited v87 native residue tracks.
# LineTrend/plot_trend.py inspected: calendar-count semantics do not apply here.
from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from PIL import Image
OUT=Path(__file__).resolve().parent;ROOT=OUT.parent;P=OUT/'source_data';R=ROOT/'pymol'
W,H=183.,100.;BLUE='#2166AC';WARM='#B87936';KEY='#3A8D8C';PARTNER='#8B72A5';NGRAY='#8A9299';CGRAY='#BAC1C7'
MODELS=[('Full','pred',BLUE,'o','-'),('Intra','intra',WARM,'s','--')]
mpl.rcParams.update({'savefig.bbox':None,'axes.unicode_minus':False,'axes.labelsize':8,'xtick.labelsize':7.2,'ytick.labelsize':7.2})
df=pd.read_csv(P/'all25_dynamics_and_static_lobe_contacts.csv').sort_values('CaMKI_native_position');assert len(df)==25 and df.valid_chi_entropy.sum()==21
manifest=json.loads((R/'render_manifest.json').read_text(encoding='utf-8'))
fig=plt.figure(figsize=(W/25.4,H/25.4),facecolor='white');points=[];checks=[]
def txt(x,y,s,fs=8,**kw):return fig.text(x/W,1-y/H,s,fontsize=fs,va=kw.pop('va','top'),color=kw.pop('color',BLACK),**kw)
def ax(x,y,w,h):return fig.add_axes([x/W,1-(y+h)/H,w/W,h/H])
def panel(lab,title,x,y):txt(x,y,lab,11,fontweight='bold');txt(x+6,y,title,8.7,fontweight='bold')
panel('a','Structural location',1,2);panel('b','Predicted flexibility',76,2)
info=manifest['assets']['overview']['trim'];a=ax(1,19,67,54);a.imshow(Image.open(R/info['file']).convert('RGBA'));a.axis('off')
for key,loc in [('N297',(.46,1.02)),('W303',(.08,.67)),('N308',(.93,.82)),('M316',(.96,.20))]:
 a.annotate(key,xy=info['landmarks'][key],xytext=loc,textcoords='axes fraction',ha='center',va='center',fontsize=7.8,color=KEY,arrowprops=dict(arrowstyle='-',color=BLACK,lw=.55,shrinkA=2,shrinkB=1),bbox=dict(fc='white',ec='none',alpha=.85,pad=.5))
for x,label,c in [(6,'CaM N-lobe',NGRAY),(34,'C-lobe',CGRAY),(57,'CaMKI',PARTNER)]:
 fig.add_artist(plt.Line2D([(x-4)/W,(x-1)/W],[1-82/H]*2,transform=fig.transFigure,color=c,lw=3));txt(x,80.6,label,7.2)
# One global numerical legend; colours always denote model variants below.
for x,(name,prefix,color,marker,style) in zip([132,157],MODELS):
 fig.add_artist(plt.Line2D([(x-4)/W,(x+1)/W],[1-4/H]*2,transform=fig.transFigure,color=color,marker=marker,ms=3.5,lw=1,ls=style));txt(x+3,2.5,name,8)
positions=df.CaMKI_native_position.to_numpy();lo,hi=293.5,318.5;x0,pw=84.,94.
coreleft=x0+(302.5-lo)/(hi-lo)*pw;corewidth=14/(hi-lo)*pw
fig.add_artist(Rectangle((coreleft/W,1-20/H),corewidth/W,.8/H,transform=fig.transFigure,fc=KEY,ec='none'))
txt(coreleft+corewidth/2,15.5,'W303–M316 recognition core',7.6,ha='center')
for _,q in df.iterrows():
 x=x0+(q.CaMKI_native_position-lo)/(hi-lo)*pw
 txt(x,23,q.amino_acid,7.5,ha='center',fontweight='bold' if q.CaMKI_native_position in [297,303,308,316] else 'normal',color=KEY if q.CaMKI_native_position in [297,303,308,316] else BLACK)
for metric,y,top,ticks in [('rmsf',31,8.2,[0,4,8]),('chi_entropy',64,4.3,[0,2,4])]:
 a=ax(x0,y,pw,24);a.axvspan(302.5,316.5,color='#F0F2F2',ec='none',zorder=0);ys=[]
 for name,prefix,color,marker,style in MODELS:
  values=df[prefix+'_'+metric].where(df['valid_'+metric]).to_numpy();ys.append(values);a.plot(positions,values,c=color,lw=1.1,ls=style,marker=marker,ms=2.7,zorder=2)
  for pos in [297,303,308,316]:
   q=df[df.CaMKI_native_position==pos].iloc[0]
   if q['valid_'+metric] and name=='Full':a.plot(pos,q[prefix+'_'+metric],marker=marker,ms=4,mfc=color,mec=BLACK,mew=.5,zorder=4)
  points.extend(dict(panel='b',metric=metric,position=int(pos),model=name,value=float(v)) for pos,v in zip(positions,values) if np.isfinite(v))
 a.set(xlim=(lo,hi),ylim=(0,top),ylabel='RMSF (Å)' if metric=='rmsf' else 'Hχ (nats)');a.set_yticks(ticks);a.tick_params(length=2,pad=2)
 if metric=='rmsf':a.tick_params(axis='x',labelbottom=False)
 else:
  a.set_xticks([294,297,303,308,312,316,318]);a.set_xlabel('CaMKI residue number',labelpad=2)
  for pos,loc in [(297,(299.0,2.65)),(308,(310.0,1.85))]:
   q=df[df.CaMKI_native_position==pos].iloc[0]
   a.annotate('N'+str(pos),xy=(pos,q.pred_chi_entropy),xytext=loc,fontsize=7.3,color=KEY,ha='center',arrowprops=dict(arrowstyle='-',lw=.45,color=BLACK))
 z=np.concatenate([v[np.isfinite(v)] for v in ys]);assert min(z)>=0 and max(z)<=top;checks.append(dict(metric=metric,ylim=[0,top],range=[float(z.min()),float(z.max())]))
pl=pd.DataFrame(points);assert len(pl)==92;pl.to_csv(P/'all_plotted_dynamics.csv',index=False)
fig.canvas.draw();ren=fig.canvas.get_renderer();vis=[t for t in fig.findobj(mpl.text.Text) if t.get_visible() and t.get_text()]
outside=[t.get_text() for t in vis if t.get_window_extent(ren).x0<-.6 or t.get_window_extent(ren).y0<-.6 or t.get_window_extent(ren).x1>fig.bbox.width+.6 or t.get_window_extent(ren).y1>fig.bbox.height+.6];assert not outside,outside
for ext,dpi in [('png',450),('pdf',600),('svg',600)]:fig.savefig(OUT/('fig3_two_panel_story.'+ext),dpi=dpi,bbox_inches=fig.bbox_inches)
qa=dict(canvas_mm=[W,H],source='family OOF fold1',panels=['a','b'],removed_independent_panels=['b_static_context','d','e','f'],all25_RMSF_and21_Hchi_positions=True,plotted_profile_values=len(pl),all_quantitative_axes_start_at_zero=True,axes=checks,no_MD_in_main=True,no_numeric_tables=True,no_gray_explanatory_notes=True,static_context_retired_to_source_files=True,outside_text=outside,min_font=min(t.get_fontsize() for t in vis),visual_review='pending',sha256={ext:hashlib.sha256((OUT/('fig3_two_panel_story.'+ext)).read_bytes()).hexdigest() for ext in ['png','pdf','svg']})
(OUT/'figure_QA.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(qa,ensure_ascii=False));plt.close(fig)
