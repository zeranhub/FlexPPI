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

# Full 1MXE supplement: 25/21 peptide profiles and 169/145 complex comparisons.
from pathlib import Path
import json,hashlib,shutil
import numpy as np,pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
OUT=Path(__file__).resolve().parent;ROOT=OUT.parent;P=OUT/'source_data';W,H=183.,147.;BLUE='#2166AC';WARM='#B87936'
mpl.rcParams.update({'savefig.bbox':None,'axes.unicode_minus':False,'axes.labelsize':7.7,'xtick.labelsize':7.1,'ytick.labelsize':7.1})
full=pd.read_csv(P/'three_CaM_family_fold1_per_residue.csv');full=full[full.sample_id=='A001K'].sort_values('residue_index');assert len(full)==169 and int(full.valid_chi_entropy.sum())==145
full.to_csv(P/'full_complex_169_family_fold1.csv',index=False)
df=pd.read_csv(P/'peptide25_model_MD_structure_all_values.csv').sort_values('CaMKI_native_position')
wr=pd.read_csv(P/'MD_rmsf_five_windows_25_residues.csv');wh=pd.read_csv(P/'MD_hchi_five_windows_21_sidechains.csv')
# Frozen MD window inputs are already included in source_data.
fig=plt.figure(figsize=(W/25.4,H/25.4),facecolor='white');points=[];summaries=[];axes=[]
def txt(x,y,s,fs=8,**kw):return fig.text(x/W,1-y/H,s,fontsize=fs,va=kw.pop('va','top'),color=kw.pop('color',BLACK),**kw)
def ax(x,y,w,h):
 a=fig.add_axes([x/W,1-(y+h)/H,w/W,h/H]);axes.append(a);return a
def panel(lab,title,x,y):txt(x,y,lab,11,fontweight='bold');txt(x+6,y,title,8.5,fontweight='bold')
MODELS=[('Full','pred',BLUE,'o','-'),('Intra','intra',WARM,'s','--')]
for metric,letter,title,y,top,ticks,windows,wcol in [('rmsf','a','CaMKI backbone',18,8.2,[0,4,8],wr,'rmsf'),('chi_entropy','b','CaMKI side chains',88,4.4,[0,2,4],wh,'hchi')]:
 panel(letter,title,1,y-16);a=ax(15,y,75,35);x=df.CaMKI_native_position.to_numpy();valid=df['valid_'+metric]
 w=windows.groupby('CaMKI_native_position')[wcol].agg(['min','max','mean']).reindex(x)
 ref=df['target_'+metric].where(valid).to_numpy();np.testing.assert_allclose(w['mean'].to_numpy()[valid],ref[valid],atol=1e-6)
 a.add_patch(Rectangle((302.5,.02),14,top-.04,fc='none',ec=BLACK,lw=.6,zorder=0))
 a.fill_between(x,w['min'],w['max'],color='#C7C7C7',alpha=.38,lw=0,zorder=0);a.plot(x,ref,c=BLACK,lw=1.1,ls=':',marker='.',ms=2.2,label='MD')
 points.extend(dict(panel=letter,metric=metric,position=int(pos),series='MD',value=float(v)) for pos,v in zip(x,ref) if np.isfinite(v))
 for name,prefix,color,marker,style in MODELS:
  vals=df[prefix+'_'+metric].where(valid).to_numpy();a.plot(x,vals,c=color,lw=.95,ls=style,marker=marker,ms=2.1)
  points.extend(dict(panel=letter,metric=metric,position=int(pos),series=name,value=float(v)) for pos,v in zip(x,vals) if np.isfinite(v))
 a.set(xlim=(293.5,318.5),ylim=(0,top),ylabel='RMSF (Å)' if metric=='rmsf' else 'Hχ (nats)',xlabel='CaMKI residue number');a.set_xticks([294,303,308,316,318]);a.set_yticks(ticks);a.tick_params(length=2,pad=2);a.xaxis.labelpad=2
 assert np.nanmax(w['max'])<=top
for metric,letter,title,y,xmax,ymax in [('rmsf','c','All 169 backbone sites',18,2.6,8.2),('chi_entropy','d','All 145 side-chain sites',88,4.4,4.4)]:
 panel(letter,title,102,y-16);a=ax(115,y,62,35);good=full['valid_'+metric];ref=full.loc[good,'target_'+metric].to_numpy();n=len(ref)
 for name,prefix,color,marker,style in MODELS:
  vals=full.loc[good,prefix+'_'+metric].to_numpy();a.scatter(ref,vals,c=color,marker=marker,s=8,alpha=.65,edgecolors='none')
  rho=float(pd.Series(ref).rank().corr(pd.Series(vals).rank()));mae=float(abs(ref-vals).mean());summaries.append(dict(metric=metric,model=name,n=n,MAE=mae,Spearman=rho))
  points.extend(dict(panel=letter,metric=metric,position=int(i),series=name,x=float(t),value=float(v)) for i,t,v in zip(full.loc[good,'residue_index'],ref,vals))
  assert max(vals)<=ymax and min(vals)>=0 and max(ref)<=xmax
 a.plot([0,min(xmax,ymax)],[0,min(xmax,ymax)],ls='--',c=BLACK,lw=.65)
 unit='RMSF (Å)' if metric=='rmsf' else 'Hχ (nats)';a.set(xlim=(0,xmax),ylim=(0,ymax),xlabel='MD '+unit,ylabel='Predicted '+unit);a.set_xticks([0,1,2] if metric=='rmsf' else [0,2,4]);a.set_yticks([0,4,8] if metric=='rmsf' else [0,2,4]);a.tick_params(length=2,pad=2);a.xaxis.labelpad=2;a.yaxis.labelpad=2
for x,label,c,ls,m in [(16,'MD',BLACK,':','.'),(45,'Full',BLUE,'-','o'),(76,'Intra',WARM,'--','s')]:
 fig.add_artist(plt.Line2D([(x-6)/W,(x-1)/W],[1-139/H]*2,transform=fig.transFigure,color=c,lw=1,ls=ls,marker=m,ms=3));txt(x+1,137.5,label,7.8)
pl=pd.DataFrame(points);assert pl.groupby('panel').size().to_dict()=={'a':75,'b':63,'c':338,'d':290};pl.to_csv(P/'all_plotted_validation_values.csv',index=False);pd.DataFrame(summaries).to_csv(P/'full_complex_validation_metrics.csv',index=False)
fig.canvas.draw();ren=fig.canvas.get_renderer();vis=[t for t in fig.findobj(mpl.text.Text) if t.get_visible() and t.get_text()];outside=[t.get_text() for t in vis if t.get_window_extent(ren).x0<-.6 or t.get_window_extent(ren).y0<-.6 or t.get_window_extent(ren).x1>fig.bbox.width+.6 or t.get_window_extent(ren).y1>fig.bbox.height+.6];assert not outside,outside
for ext,dpi in [('png',450),('pdf',600),('svg',600)]:fig.savefig(OUT/('figS3_1MXE_validation.'+ext),dpi=dpi,bbox_inches=fig.bbox_inches)
qa=dict(canvas_mm=[W,H],source='family OOF fold1',complete_peptide_positions=[25,21],complete_complex_positions=[169,145],MD_window_ranges='five continuous 100ns windows, not independent replicates',outside_text=outside,plotted_counts=pl.groupby('panel').size().to_dict(),all_axes_zero=True,summary=summaries,visual_review='pending',sha256={ext:hashlib.sha256((OUT/('figS3_1MXE_validation.'+ext)).read_bytes()).hexdigest() for ext in ['png','pdf','svg']})
(OUT/'figure_QA.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(qa,ensure_ascii=False));plt.close(fig)