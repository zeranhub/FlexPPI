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

# Figure contract: true-coordinate a/b; d/e/f above c; complete 25/21 profiles.
# The established residue-track and paired-point logic is inherited from v83.
# LineTrend/plot_trend.py is a timeline template, so no timeline data are reused.
from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle,Ellipse,FancyArrowPatch
from PIL import Image
OUT=Path(__file__).resolve().parent;ROOT=OUT.parent;P=OUT/'source_data';R=ROOT/'pymol'
W,H=183.,190.;BLUE='#2166AC';WARM='#B87936';KEY='#3A8D8C';PARTNER='#8B72A5';NGRAY='#8A9299';CGRAY='#BAC1C7'
MODELS=[('Full','pred',BLUE,'o','-'),('Intra','intra',WARM,'s','--')]
mpl.rcParams.update({'savefig.bbox':None,'axes.unicode_minus':False,'axes.labelsize':7.7,'xtick.labelsize':7.1,'ytick.labelsize':7.1})
df=pd.read_csv(P/'peptide25_model_MD_structure_all_values.csv').sort_values('CaMKI_native_position');assert len(df)==25 and df.valid_chi_entropy.sum()==21
manifest=json.loads((R/'render_manifest.json').read_text(encoding='utf-8'))
fig=plt.figure(figsize=(W/25.4,H/25.4),facecolor='white');axes=[];points=[];checks=[]
def txt(x,y,s,fs=8,**kw):return fig.text(x/W,1-y/H,s,fontsize=fs,va=kw.pop('va','top'),color=kw.pop('color',BLACK),**kw)
def ax(x,y,w,h):
 a=fig.add_axes([x/W,1-(y+h)/H,w/W,h/H]);axes.append(a);return a
def panel(lab,title,x,y):txt(x,y,lab,11,fontweight='bold');txt(x+6,y,title,8.5,fontweight='bold')
def photo(key,x,y,w,h,trim=False):
 info=manifest['assets'][key];info=info['trim'] if trim else info
 a=ax(x,y,w,h);a.imshow(Image.open(R/info['file']).convert('RGBA'));a.axis('off');return a,info
panel('a','CaM–CaMKI recognition',1,2);panel('b','Anchor environments',94,2)
a,info=photo('overview',1,10,88,55,True)
for key,loc in [('N297',(.47,1.02)),('W303',(.05,.67)),('N308',(.85,.80)),('M316',(.93,.24))]:
 a.annotate(key,xy=info['landmarks'][key],xytext=loc,textcoords='axes fraction',ha='center',va='center',fontsize=7.5,color=KEY,arrowprops=dict(arrowstyle='-',color=BLACK,lw=.55,shrinkA=2,shrinkB=1),bbox=dict(fc='white',ec='none',alpha=.85,pad=.5))
for x,label,c in [(6,'CaM N-lobe',NGRAY),(37,'C-lobe',CGRAY),(61,'CaMKI',PARTNER)]:
 fig.add_artist(plt.Line2D([(x-4)/W,(x-1)/W],[1-69/H]*2,transform=fig.transFigure,color=c,lw=3));txt(x,67.7,label,7.1)
for key,x,label,pos in [('local_303',94,'W303 · C-lobe',303),('local_316',139,'M316 · N-lobe',316)]:
 txt(x+21,10,label,7.8,ha='center',fontweight='bold');a,inf=photo(key,x,15,42,50)
 landmark=inf['landmarks']['W303' if pos==303 else 'M316'];a.plot(*landmark,marker='o',mfc='none',mec=KEY,ms=17,mew=.7)
# Shared numerical legend; structural component colours have separate meanings.
for x,(name,prefix,color,marker,style) in zip([135,159],MODELS):
 fig.add_artist(plt.Line2D([(x-4)/W,(x+1)/W],[1-74/H]*2,transform=fig.transFigure,color=color,marker=marker,ms=3,lw=.9,ls=style));txt(x+3,72.5,name,7.7)
panel('d','Backbone restraint',1,78);panel('e','Same-residue contrast',64,78);panel('f','Mechanistic hypothesis',119,78)
a=ax(14,90,42,25)
for j,(mask,label) in enumerate([(df.CaMKI_native_position.between(303,316),'303–316'),(~df.CaMKI_native_position.between(303,316),'Other 11')]):
 vals=[]
 for dx,(name,prefix,color,marker,style) in zip([-.13,.13],MODELS):
  v=float(df.loc[mask,prefix+'_rmsf'].mean());vals.append(v);a.scatter(j+dx,v,s=23,c=color,marker=marker,zorder=3)
  a.annotate(f'{v:.2f}',(j+dx,v),xytext=((-6,0) if v<.7 else (0,4 if name=='Intra' else -4)),textcoords='offset points',ha='right' if v<.7 else 'center',va='center' if v<.7 else ('bottom' if name=='Intra' else 'top'),fontsize=7)
  points.append(dict(panel='d',metric='rmsf',region=label,model=name,value=v,n=int(mask.sum())))
 a.plot([j-.13,j+.13],vals,c=BLACK,lw=.65,zorder=1)
a.set(xlim=(-.55,1.55),ylim=(0,4.6),ylabel='Mean RMSF (Å)');a.set_yticks([0,2,4]);a.set_xticks([0,1],['303–316','Other 11'],fontsize=7);a.tick_params(length=2,pad=2)
a=ax(76,90,34,25)
for j,pos in enumerate([297,308]):
 q=df[df.CaMKI_native_position==pos].iloc[0];vals=[]
 for dx,(name,prefix,color,marker,style) in zip([-.13,.13],MODELS):
  v=float(q[prefix+'_chi_entropy']);vals.append(v);above=v==max(q.pred_chi_entropy,q.intra_chi_entropy);a.scatter(j+dx,v,s=23,c=color,marker=marker,zorder=3)
  a.annotate(f'{v:.2f}',(j+dx,v),xytext=(0,4 if above else -4),textcoords='offset points',ha='center',va='bottom' if above else 'top',fontsize=7)
  points.append(dict(panel='e',metric='chi_entropy',position=pos,model=name,value=v))
 a.plot([j-.13,j+.13],vals,c=BLACK,lw=.65,zorder=1)
a.set(xlim=(-.6,1.6),ylim=(0,2.2),ylabel='Hχ (nats)');a.set_yticks([0,1,2]);a.set_xticks([0,1],['N297','N308'],fontsize=7);a.tick_params(length=2,pad=2)
# Dashed arrow = proposed consequence, not an observed activation trajectory.
a=ax(123,87,59,33);a.set(xlim=(0,100),ylim=(0,74));a.axis('off')
for x,label,c in [(23,'C-lobe binding',CGRAY),(77,'N-lobe restraint',NGRAY)]:
 a.add_patch(Ellipse((x,57),40,21,fc=c,ec='none',alpha=.8));a.text(x,59,label,ha='center',va='center',fontsize=7.3)
a.plot([22,78],[47,47],c=PARTNER,lw=2.5);a.scatter([22,78],[47,47],s=17,color=KEY,zorder=3)
a.text(22,42,'W303',ha='center',va='top',fontsize=6.8,color=KEY);a.text(78,42,'M316',ha='center',va='top',fontsize=6.8,color=KEY)
a.text(50,27,'Regulatory segment retained',ha='center',va='center',fontsize=7)
a.add_patch(FancyArrowPatch((50,21),(50,10),arrowstyle='-|>',mutation_scale=7,ls=(0,(3,2)),lw=.8,color=BLACK))
a.text(50,6,'Less return to inhibitory contacts',ha='center',va='top',fontsize=7)
panel('c','Predicted CaMKI profiles',1,128)
positions=df.CaMKI_native_position.to_numpy();lo,hi=293.5,318.5
corex=19+(309.5-lo)/(hi-lo)*158;txt(corex,136,'W303–M316 core',7.4,ha='center',fontweight='bold')
for metric,y,top,ticks in [('rmsf',143,8.2,[0,4,8]),('chi_entropy',165,4.3,[0,2,4])]:
 a=ax(19,y,158,16)
 a.add_patch(Rectangle((302.5,.02),14,top-.04,fc='#F7F7F7',ec=BLACK,lw=.6,zorder=0))
 ys=[]
 for name,prefix,color,marker,style in MODELS:
  values=df[prefix+'_'+metric].where(df['valid_'+metric]).to_numpy();ys.append(values);a.plot(positions,values,c=color,lw=1,ls=style,marker=marker,ms=2.3,zorder=2)
  for pos in [297,303,308,316]:
   q=df[df.CaMKI_native_position==pos].iloc[0]
   if q['valid_'+metric] and name=='Full':a.plot(pos,q[prefix+'_'+metric],marker=marker,ms=3.5,mfc=color,mec=BLACK,mew=.5,zorder=4)
  points.extend(dict(panel='c',metric=metric,position=int(pos),model=name,value=float(v)) for pos,v in zip(positions,values) if np.isfinite(v))
 a.set(xlim=(lo,hi),ylim=(0,top),ylabel='RMSF (Å)' if metric=='rmsf' else 'Hχ (nats)');a.set_yticks(ticks);a.tick_params(length=2,pad=2)
 if metric=='rmsf':a.tick_params(axis='x',labelbottom=False)
 else:a.set_xticks([294,298,303,308,312,316,318]);a.set_xlabel('CaMKI residue number',labelpad=2)
 z=np.concatenate([v[np.isfinite(v)] for v in ys]);assert min(z)>=0 and max(z)<=top
 checks.append(dict(metric=metric,ylim=[0,top],range=[float(z.min()),float(z.max())]))
pl=pd.DataFrame(points);assert len(pl[pl.panel=='c'])==92 and len(pl)==100
pl.to_csv(P/'all_plotted_quantities.csv',index=False)
fig.canvas.draw();ren=fig.canvas.get_renderer();vis=[t for t in fig.findobj(mpl.text.Text) if t.get_visible() and t.get_text()];outside=[t.get_text() for t in vis if t.get_window_extent(ren).x0<-.6 or t.get_window_extent(ren).y0<-.6 or t.get_window_extent(ren).x1>fig.bbox.width+.6 or t.get_window_extent(ren).y1>fig.bbox.height+.6];assert not outside,outside
for ext,dpi in [('png',450),('pdf',600),('svg',600)]:fig.savefig(OUT/('fig3_1MXE_mechanistic_hypothesis.'+ext),dpi=dpi,bbox_inches=fig.bbox_inches)
qa=dict(canvas_mm=[W,H],source='family OOF fold1',panel_layout='a,b top; d,e,f middle; c bottom',no_MD_drawn_in_main=True,no_numeric_tables=True,no_gray_notes=True,all_25_RMSF_and_21_Hchi_positions=True,core_defined_by_structural_anchors=[303,316],plotted_counts=pl.groupby('panel').size().to_dict(),all_quantitative_axes_start_at_zero=True,axes=checks,full_intra_are_model_variants_not_biological_states=True,model_palette={'Full':BLUE,'Intra':WARM},outside_text=outside,min_font=min(t.get_fontsize() for t in vis),visual_review='pending',sha256={ext:hashlib.sha256((OUT/('fig3_1MXE_mechanistic_hypothesis.'+ext)).read_bytes()).hexdigest() for ext in ['png','pdf','svg']})
(OUT/'figure_QA.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(qa,ensure_ascii=False));plt.close(fig)