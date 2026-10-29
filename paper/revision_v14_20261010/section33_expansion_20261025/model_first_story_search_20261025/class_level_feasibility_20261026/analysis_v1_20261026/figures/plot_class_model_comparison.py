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

# Class cohort comparison: all included PDBs, no performance-based exemplar selection.
from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
OUT=Path(__file__).resolve().parent;R=OUT.parent
W,H=183.,130.
mpl.rcParams.update({'savefig.bbox':None,'font.size':8.5,'axes.labelsize':8.5,'xtick.labelsize':8.3,'ytick.labelsize':8.3})
df=pd.read_csv(R/'formal_results/per_record_all_methods_references.csv')
df=df[df.core_functional_peptide & df.reference.eq('uniform_rep1') & df.scope.eq('all')]
METHODS=['full','intra','static47_sasa'];NAMES=['Full','Intra','Structural']
COLORS=['#2166AC','#B87936','#767F87'];MARKERS=['o','s','D']
points=[];means=[];fig=plt.figure(figsize=(W/25.4,H/25.4),facecolor='white')
def txt(x,y,s,size=9,**kw):return fig.text(x/W,1-y/H,s,fontsize=size,va='top',**kw)
def ax(x,y,w,h):return fig.add_axes([x/W,1-(y+h)/H,w/W,h/H])
for k,(role,metric,title,x,y,h,lim,ticks) in enumerate([
 ('peptide','rmsf','Peptide · RMSF',18,13,43,4.2,[0,1,2,3,4]),
 ('peptide','hchi','Peptide · Hχ',110,13,43,.45,[0,.15,.30,.45]),
 ('receptor','rmsf','Receptor · RMSF',18,79,40,4.2,[0,1,2,3,4]),
 ('receptor','hchi','Receptor · Hχ',110,79,40,.45,[0,.15,.30,.45])
]):
 txt(x-16,y-10,chr(97+k),11,fontweight='bold')
 txt(x-10,y-10,title,9,fontweight='bold')
 a=ax(x,y,65,h);g=df[df.role.eq(role)&df.metric.eq(metric)&df.method.isin(METHODS)]
 # One value per PDB and method; duplicates average before plotting.
 z=g.groupby(['pdb_id','component_id','method']).mae.mean().unstack('method').reindex(columns=METHODS)
 assert len(z)==40 and z.notna().all().all()
 for key,row in z.iterrows():
  a.plot([0,1,2],row.to_numpy(),c='#CDD2D6',lw=.45,alpha=.65,zorder=1)
  for j,m in enumerate(METHODS):
   a.scatter(j,row[m],s=7,c=COLORS[j],alpha=.35,marker=MARKERS[j],linewidths=0,zorder=2)
   points.append(dict(panel=chr(97+k),pdb_id=key[0],component_id=key[1],role=role,metric=metric,method=m,mae=float(row[m])))
 fam=z.groupby('component_id').agg(['sum','count'])
 rng=np.random.default_rng(20261026);ix=rng.integers(0,len(fam),size=(10000,len(fam)))
 for j,m in enumerate(METHODS):
  b=fam[(m,'sum')].to_numpy()[ix].sum(1)/fam[(m,'count')].to_numpy()[ix].sum(1)
  mean=float(z[m].mean());lo,hi=np.quantile(b,[.025,.975])
  a.errorbar(j,mean,yerr=[[mean-lo],[hi-mean]],fmt=MARKERS[j],ms=5.4,mec='white',mew=.55,color=COLORS[j],elinewidth=1.35,capsize=3,capthick=1.0,zorder=6)
  means.append(dict(role=role,metric=metric,method=m,PDB_MAE=mean,CI95_lo=float(lo),CI95_hi=float(hi),n_PDB=40,n_split_groups=36))
 assert z.to_numpy().min()>=0 and z.to_numpy().max()<=lim
 a.set(xlim=(-.4,2.4),ylim=(0,lim),ylabel='MAE (Å)' if metric=='rmsf' else 'MAE (nats)')
 a.set_xticks([0,1,2],NAMES);a.set_yticks(ticks);a.tick_params(length=2,pad=2)
 a.spines['bottom'].set_bounds(-.25,2.25)
pd.DataFrame(points).to_csv(OUT/'figure_all480_values.csv',index=False)
pd.DataFrame(means).to_csv(OUT/'figure_macro_means_CI.csv',index=False)
assert len(points)==480
fig.canvas.draw();ren=fig.canvas.get_renderer()
texts=[t for t in fig.findobj(mpl.text.Text) if t.get_visible() and t.get_text()]
for t in texts:
 if t.get_fontsize()<8.2:t.set_fontsize(8.2)
fig.canvas.draw();outside=[t.get_text() for t in texts if t.get_window_extent(ren).x0<-.5 or t.get_window_extent(ren).y0<-.5 or t.get_window_extent(ren).x1>fig.bbox.width+.5 or t.get_window_extent(ren).y1>fig.bbox.height+.5]
assert not outside,outside
for ext in ['png','pdf','svg']:fig.savefig(OUT/('class_core40_absolute_MAE.'+ext),dpi=450,bbox_inches=fig.bbox_inches)
qa=dict(width_mm=W,height_mm=H,points=480,all_PDBs_per_panel=40,split_groups=36,methods=METHODS,reference='0-100 ns, 501 frames, same frozen predictions',axes_start_at_zero=True,all_individual_points=True,duplicate_PDBs_averaged_before_plotting=True,means='PDB equal; family connected-component cluster percentile interval 10000 draws',no_tables=True,no_gray_notes=True,no_experimental_claim=True,min_font=min(t.get_fontsize() for t in texts),outside_text=outside,visual_QA='pending',sha256={e:hashlib.sha256((OUT/('class_core40_absolute_MAE.'+e)).read_bytes()).hexdigest() for e in ['png','pdf','svg']})
(OUT/'figure_QA.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(qa,ensure_ascii=False));plt.close(fig)





