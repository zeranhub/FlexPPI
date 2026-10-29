# Academic Figure Skill Asset Confirmation (existing scientific profile and paired-point logic)
# Full sequence -> LineTrend inherited; region/site absolute pairs -> existing paired-point logic inherited.
# All frozen-model residue values retained. Fixed 9-residue regions are descriptive model-selected examples.
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
import json,hashlib
import numpy as np,pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
OUT=Path(__file__).resolve().parent;P=OUT/'source_data';W,H=183.,70.;A='#B87936';B='#2166AC'
STATES=[('NHE1','a',A,'o'),('IQCG','b',B,'s')]
mpl.rcParams.update({'font.size':8,'axes.labelsize':7.5,'xtick.labelsize':7.2,'ytick.labelsize':7.2,'savefig.bbox':None,'axes.unicode_minus':False})
fig=plt.figure(figsize=(W/25.4,H/25.4),facecolor='white');points=[];axes=[]
def txt(x,y,s,fs=8,**kw):return fig.text(x/W,1-y/H,s,fontsize=fs,va=kw.pop('va','top'),color=kw.pop('color',BLACK),**kw)
def ax(x,y,w,h):
    a=fig.add_axes([x/W,1-(y+h)/H,w/W,h/H]);axes.append(a);return a
df=pd.read_csv(P/'cam_md_state_contrast_all_common_residues.csv');assert len(df[df.metric=='rmsf'])==145 and len(df[df.metric=='hchi'])==122
rank=[]
for metric,g in df.groupby('metric'):
    g=g.copy();g['absolute_change_rank']=g.full0_delta.abs().rank(ascending=False,method='min').astype(int);g['relative_change_percent']=100*g.full0_delta/g.full0_a;rank.append(g)
pd.concat(rank).to_csv(P/'all_full_prediction_effect_ranks.csv',index=False)
regions=[(37,45),(111,119)];windows=[]
for l,h in regions:
    for metric,g in df.groupby('metric'):
        q=g[g.mature_cam_number.between(l,h)];windows.append(dict(first=l,last=h,metric=metric,n=len(q),full_NHE1=float(q.full0_a.mean()),full_IQCG=float(q.full0_b.mean()),full_delta=float(q.full0_delta.mean()),MD_NHE1=float(q.md_a.mean()),MD_IQCG=float(q.md_b.mean()),MD_delta=float(q.md_delta.mean()),Full_increased=int((q.full0_delta>0).sum()),Full_decreased=int((q.full0_delta<0).sum())))
win=pd.DataFrame(windows);win.to_csv(P/'regional_effects_full_and_reference.csv',index=False)
txt(1,1,'c',11,fontweight='bold');txt(7,1,'Partner-dependent flexibility',8.5,fontweight='bold')
for x,(state,s,c,m) in zip([134,160],STATES):
    fig.add_artist(plt.Line2D([x/W],[1-3/H],transform=fig.transFigure,color=c,marker=m,ms=3.6,lw=0));txt(x+3,1.5,state,7.5)
lo,hi=3.5,148.5
for metric,y in [('rmsf',15),('hchi',43)]:
    a=ax(14,y,108,15);q=df[df.metric==metric].set_index('mature_cam_number').reindex(range(4,149));top=4 if metric=='rmsf' else 4.1
    for l,h in regions:a.add_patch(Rectangle((l-.5,.04),h-l+1,top-.08,facecolor='#F5F5F5',edgecolor=BLACK,lw=.65,zorder=0))
    for state,s,c,m in STATES:
        a.plot(q.index,q['full0_'+s],color=c,lw=.95)
        if metric=='hchi':
            for pos in [37,42,60]:a.plot(pos,q.loc[pos,'full0_'+s],marker=m,ms=3.1,mfc=c,mec=BLACK,mew=.45,zorder=4)
        points.extend(dict(panel='full',metric=metric,position=int(p),state=state,value=float(v)) for p,v in q['full0_'+s].dropna().items())
    a.set(xlim=(lo,hi),ylim=(0,top),ylabel='RMSF (Å)' if metric=='rmsf' else 'Hχ (nats)');a.set_yticks([0,2,4]);a.tick_params(length=2,pad=2)
    if metric=='rmsf':a.tick_params(axis='x',labelbottom=False)
    else:a.set_xticks([4,20,40,60,80,100,120,140,148]);a.set_xlabel('Mature CaM residue number',labelpad=2)
for l,h in regions:txt(14+(((l+h)/2)-lo)/(hi-lo)*108,10,f'{l}–{h}',7.2,ha='center',fontweight='bold')

txt(142,10,'Region means',7.8,ha='center',fontweight='bold');a=ax(140,15,39,15)
for j,(l,h) in enumerate(regions):
    q=win[(win.metric=='rmsf')&(win['first']==l)].iloc[0];vals=[]
    for dx,(state,s,c,m) in zip([-.11,.11],STATES):
        v=q['full_NHE1'] if s=='a' else q['full_IQCG'];vals.append(v);a.scatter(j+dx,v,color=c,marker=m,s=21,zorder=3);above=v==max(q.full_NHE1,q.full_IQCG)
        a.annotate(f'{v:.2f}',(j+dx,v),xytext=(0,3 if above else -4),textcoords='offset points',ha='center',va='bottom' if above else 'top',fontsize=6.9,color=BLACK)
        points.append(dict(panel='region',metric='rmsf',region=f'{l}-{h}',state=state,value=float(v)))
    a.plot([j-.11,j+.11],vals,color=BLACK,lw=.6,zorder=1)
a.set(xlim=(-.55,1.55),ylim=(0,2.6),ylabel='RMSF (Å)');a.set_yticks([0,1,2]);a.set_xticks([0,1],['37–45','111–119'],fontsize=7);a.tick_params(length=2,pad=2);a.yaxis.labelpad=2

txt(157,37,'Side-chain sites',7.8,ha='center',fontweight='bold');a=ax(140,41,39,21)
for j,pos in enumerate([37,42,60]):
    q=df[(df.metric=='hchi')&(df.mature_cam_number==pos)].iloc[0];vals=[]
    for dx,(state,s,c,m) in zip([-.11,.11],STATES):
        v=q['full0_'+s];vals.append(v);a.scatter(j+dx,v,color=c,marker=m,s=21,zorder=3)
        above=v==max(q.full0_a,q.full0_b)
        a.annotate(f'{v:.2f}',(j+dx,v),xytext=(0,3 if above else -4),textcoords='offset points',ha='center',va='bottom' if above else 'top',fontsize=6.9,color=BLACK)
        points.append(dict(panel='site',metric='hchi',position=pos,state=state,value=float(v)))
    a.plot([j-.11,j+.11],vals,color=BLACK,lw=.6,zorder=1)
a.set(xlim=(-.55,2.55),ylim=(0,4.1),ylabel='Hχ (nats)');a.set_yticks([0,2,4]);a.set_xticks([0,1,2],['R37','N42','N60'],fontsize=7);a.tick_params(length=2,pad=2);a.yaxis.labelpad=2
fig.canvas.draw();ren=fig.canvas.get_renderer();vis=[t for t in fig.findobj(mpl.text.Text) if t.get_visible() and t.get_text()];outside=[t.get_text() for t in vis if t.get_window_extent(ren).x0<-.6 or t.get_window_extent(ren).y0<-.6 or t.get_window_extent(ren).x1>fig.bbox.width+.6 or t.get_window_extent(ren).y1>fig.bbox.height+.6];assert not outside,outside
assert len([q for q in points if q['panel']=='full'])==534
for ext,dpi in [('png',450),('pdf',600),('svg',600)]:fig.savefig(OUT/('fig3c_region_diagnostic.'+ext),dpi=dpi,bbox_inches=fig.bbox_inches)
pd.DataFrame(points).to_csv(P/'plotted_values.csv',index=False)
qa={'all_full_sequence_values':534,'all_axes_start_at_zero':True,'no_MD_drawn':True,'no_axis_cropping':True,'no_tables_or_gray_notes':True,'full_sequence_limits':{'rmsf':[0,4],'hchi':[0,4.1]},'outside_text':outside,'visual_review':'pending','output_sha256':{ext:hashlib.sha256((OUT/('fig3c_region_diagnostic.'+ext)).read_bytes()).hexdigest() for ext in ['png','pdf','svg']}}
(OUT/'figure_QA.json').write_text(json.dumps(qa,indent=2),encoding='utf-8');print(json.dumps(qa));plt.close(fig)