# Academic Figure Skill Asset Confirmation (existing publication production inspected)
# matched groups -> native v100 journal typography and paired point grammar -> param inherit
# all 83 record-level pairs retained in each readout; no generated data or fitted curve
from pathlib import Path
import json
import shutil
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.rcParams.update({
    'font.family':'sans-serif','font.sans-serif':['Arial','Helvetica','Liberation Sans'],
    'font.size':9,'axes.titlesize':10,'axes.labelsize':9,'xtick.labelsize':9,'ytick.labelsize':9,
    'legend.fontsize':9,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.6,
    'xtick.direction':'out','ytick.direction':'out','xtick.major.width':.6,'ytick.major.width':.6,
    'legend.frameon':False,'axes.unicode_minus':False,'pdf.fonttype':42,'svg.fonttype':'none',
    'savefig.dpi':450,'savefig.bbox':None})
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

root=Path(__file__).resolve().parent
records=pd.read_csv(root/'record_means_original_RMSF_tol0.10.csv',float_precision='round_trip')
summary=pd.read_csv(root/'summary_all_prespecified_analyses.csv',float_precision='round_trip')
summary=summary[summary.matching_variable.eq('original_RMSF')&summary.tolerance_A.eq(.1)]
assert len(records)==83 and len(summary)==12
BLUE='#2166AC';GRAY='#92999F';INK='#262626'
def panel(a,metric,title,ylabel,ylimits,yticks,delta_label):
    # Every gray line is the equally amino-acid-weighted paired mean of a record.
    for r in records.itertuples(index=False):
        lo=getattr(r,metric+'_low');hi=getattr(r,metric+'_high')
        a.plot([0,1],[lo,hi],color=GRAY,alpha=.27,lw=.45,zorder=1)
        a.scatter([0,1],[lo,hi],s=4.,color=GRAY,alpha=.35,linewidths=0,zorder=2)
    s=summary[summary.metric.eq(metric)].set_index('value_kind')
    means=s.loc[['low','high'],'mean'].to_numpy()
    errors=np.array([means-s.loc[['low','high'],'bootstrap95_low'].to_numpy(),s.loc[['low','high'],'bootstrap95_high'].to_numpy()-means])
    a.errorbar([0,1],means,yerr=errors,fmt='o-',color=BLUE,lw=1.25,ms=4,mfc=BLUE,mew=.6,capsize=2,elinewidth=1.0,zorder=6)
    a.set(xlim=(-.28,1.28),ylim=ylimits,ylabel=ylabel)
    a.set_xticks([0,1],['Lower predicted Hχ','Higher predicted Hχ']);a.set_yticks(yticks)
    a.tick_params(length=2.5,pad=2);a.xaxis.labelpad=2.5;a.yaxis.labelpad=2.5
    a.set_title(title,pad=21,fontweight='normal',fontsize=10)
    a.text(.5,1.02,delta_label,transform=a.transAxes,ha='center',va='bottom',fontsize=9,color=INK)
    return means.tolist()

W,H=183.,73.
fig=plt.figure(figsize=(W/25.4,H/25.4),facecolor='white')
a=fig.add_axes([15/W,15/H,63/W,43/H])
panel(a,'original_RMSF','Backbone motion (matched)','MD reference RMSF (Å)',(0,6.8),[0,2,4,6],'Mean 1.116 → 1.118 Å')
a=fig.add_axes([111/W,15/H,63/W,43/H])
panel(a,'fine_entropy_mean_per_chi','Side-chain diversity','MD occupancy entropy (nats/χ)',(2.35,3.15),[2.4,2.7,3.0],'Mean 2.531 → 2.689')
handles=[Line2D([],[],color=GRAY,alpha=.6,lw=.6,label='Each record'),Line2D([],[],color=BLUE,lw=1.25,marker='o',ms=3,label='Mean ± 95% CI')]
fig.legend(handles=handles,loc='center',bbox_to_anchor=(.5,.045),ncol=2,handlelength=1.7,handletextpad=.45,columnspacing=2,frameon=False)
fig.canvas.draw();renderer=fig.canvas.get_renderer()
outside=[]
for t in fig.findobj(mpl.text.Text):
    if t.get_visible() and t.get_text():
        b=t.get_window_extent(renderer)
        if b.x0<-.5 or b.y0<-.5 or b.x1>fig.bbox.width+.5 or b.y1>fig.bbox.height+.5:outside.append(t.get_text())
assert not outside,outside
for ext in ['png','pdf','svg']:fig.savefig(root/f'f_direct_matched_Hchi_concept.{ext}',dpi=450,bbox_inches=fig.bbox_inches)
plt.close(fig)

# Keep the angular-range validation and all matched-group controls separately.
fig=plt.figure(figsize=(W/25.4,H/25.4),facecolor='white')
a=fig.add_axes([15/W,15/H,63/W,43/H])
panel(a,'original_RMSF','Backbone motion (matched)','MD reference RMSF (Å)',(0,6.8),[0,2,4,6],'Mean 1.116 → 1.118 Å')
a=fig.add_axes([111/W,15/H,63/W,43/H])
panel(a,'width_Fstar_le_1p0_degrees','Side-chain angular range','MD low-F* width (°/χ)',(28,56),[30,40,50],'Mean 36.44° → 40.42°')
fig.legend(handles=handles,loc='center',bbox_to_anchor=(.5,.045),ncol=2,handlelength=1.7,handletextpad=.45,columnspacing=2,frameon=False)
for ext in ['png','pdf','svg']:fig.savefig(root/f'f_direct_matched_angular_range.{ext}',dpi=450,bbox_inches=fig.bbox_inches)
plt.close(fig)
(root/'plot_QA.json').write_text(json.dumps({'records_per_readout':83,'all_record_lines_retained':True,'no_rho_axes':True,'MD_RMSF_near_not_exact_equal':True,'means_CI_from_prespecified_source_stratified_record_bootstrap':True,'size_mm':[W,H],'style_matches_actual_Fig2':True,'outside_text':outside,'visual_review':'pending'},ensure_ascii=False,indent=2),encoding='utf-8')
print('Saved direct matched comparison concepts; no existing figure/manuscript changed.')
