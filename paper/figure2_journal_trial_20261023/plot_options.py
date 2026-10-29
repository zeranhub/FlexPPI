# Academic Figure Skill Asset Confirmation
# A: inherited_v61_plot.py -> frozen record-weighted density and complete PR; plot.py -> paired AP.
# B: same frozen data -> task-led asymmetric composition; no new metric or selected record.
# C: GroupedViolin/plot_GroupedViolin.py -> parameter inheritance for all-record AP distributions.
# The production violin template's unpaired significance tests are not applicable and are not run.
# plot.py supplies the verbatim typography, palette and export baseline, shared by all three options.
from pathlib import Path
import runpy
import shutil
import json
import hashlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import fitz

ROOT=Path(__file__).resolve().parent
state=runpy.run_path(str(ROOT/'plot.py'))
for key in ['METHODS','TASKS','COLORS','STYLES','WIDTHS','ZORDER','summary','pr','endpoints',
            'within_pr','within_end','trend','record_metrics','within_record','points','BLUE','TEAL','MD',
            'density_cmap','density_norm','export','number_format']:
    globals()[key]=state[key]
plt.close('all')
W=183.
all_metrics={**summary['metrics'],'within_type_hchi':summary['within_type']}
frames={t:record_metrics[record_metrics.target.eq(t)] for t in TASKS[:2]}
frames['within_type_hchi']=within_record
for extension in ['pdf','svg','png']:
    shutil.copyfile(ROOT/f'fig2_journal_trial.{extension}',ROOT/f'option_A_balanced.{extension}')

def canvas(height):
    return plt.figure(figsize=(W/25.4,height/25.4)),height

def ax_mm(fig,height,x,y,w,h):return fig.add_axes([x/W,y/height,w/W,h/height])
def txt(fig,height,x,y,s,**kw):return fig.text(x/W,y/height,s,**kw)
def title(fig,height,x,y,letter,label):
    txt(fig,height,x-11,y,letter,fontsize=11,fontweight='bold',va='bottom')
    txt(fig,height,x,y,label,fontsize=10,va='bottom')

def shared_key(fig,height):
    fig.legend(handles=[Line2D([],[],color=COLORS[m],ls=STYLES[m],lw=1.4,label=m) for m in METHODS],
               loc='center',bbox_to_anchor=(69/W,(height-4)/height),ncol=4,columnspacing=1,
               handlelength=2,handletextpad=.4,borderaxespad=0)
    txt(fig,height,126,height-4,'Mass',ha='right',va='center',fontsize=9,color='#555555')
    axis=ax_mm(fig,height,131,height-5,39,1.6)
    scalar=plt.cm.ScalarMappable(cmap=density_cmap,norm=density_norm)
    cb=fig.colorbar(scalar,cax=axis,orientation='horizontal')
    cb.solids.set_rasterized(False);cb.outline.set_visible(False)
    cb.ax.xaxis.set_ticks_position('bottom');cb.set_ticks([1e-5,1e-3,1e-2])
    cb.ax.tick_params(labelsize=9,pad=.5,length=1.5)
    from matplotlib.ticker import NullLocator
    cb.ax.xaxis.set_minor_locator(NullLocator())

def density(fig,height,x,y,size,task,letter,heading_y):
    ax=ax_mm(fig,height,x,y,size,size)
    z=points[task];m=summary['metrics'][task]['Full']
    ax.pcolormesh(z['x_edges'],z['y_edges'],np.ma.masked_equal(z['histogram_mass'].T,0),
                  cmap=density_cmap,norm=density_norm,shading='flat',linewidth=0,rasterized=False)
    lo=min(float(z['reference'].min()),float(z['prediction'].min()),float(z['x_edges'][0]),float(z['y_edges'][0]))
    hi=max(float(z['reference'].max()),float(z['prediction'].max()),float(z['x_edges'][-1]),float(z['y_edges'][-1]))
    ax.plot([lo,hi],[lo,hi],color=MD,lw=.75,ls=(0,(3,3)))
    tr=trend[trend.target.eq(task)&trend.method.eq('Full')].sort_values('bin_index')
    ax.fill_between(tr.x_mean,tr.ci_low,tr.ci_high,color=BLUE,alpha=.15,lw=0,zorder=3)
    ax.plot(tr.x_mean,tr.prediction_mean,color=BLUE,lw=1.25,zorder=4)
    units='Å' if task=='rmsf' else 'nats'
    ax.set(xlim=(lo,hi),ylim=(lo,hi),xlabel=f'MD reference ({units})',ylabel=f'Full prediction ({units})')
    ax.set_aspect('equal',adjustable='box')
    if task=='rmsf':
        ax.set_xscale('symlog',linthresh=1);ax.set_yscale('symlog',linthresh=1)
        ax.set_xticks([0,1,10,30]);ax.set_yticks([0,1,10,30]);number_format(ax,'x');number_format(ax,'y')
    else:ax.set_xticks([0,1,2,3,4]);ax.set_yticks([0,1,2,3,4])
    ax.tick_params(pad=2);ax.xaxis.labelpad=3;ax.yaxis.labelpad=3
    ax.text(.035,.97,f'MAE {m["macro_mae"]:.3f} {units}\nρ {m["mean_spearman"]:.3f}',
            transform=ax.transAxes,ha='left',va='top',fontsize=9,linespacing=1.35,
            bbox=dict(facecolor='white',edgecolor='none',alpha=.9,pad=1.2))
    title(fig,height,x,heading_y,letter,'RMSF' if task=='rmsf' else 'Hχ')

def precision_recall(fig,height,x,y,w,h,task,letter,heading_y,ap_y,ylabel=True):
    ax=ax_mm(fig,height,x,y,w,h)
    frame=within_pr if task=='within_type_hchi' else pr[pr.target.eq(task)]
    ends=within_end if task=='within_type_hchi' else endpoints[endpoints.target.eq(task)]
    metrics=all_metrics[task]
    for method in METHODS:
        cur=frame[frame.method.eq(method)].sort_values('recall')
        end=ends[ends.method.eq(method)].iloc[0]
        assert len(cur)==101 and end.recall==1
        xx=np.r_[cur.recall.to_numpy(),end.recall];yy=np.r_[cur.precision.to_numpy(),end.precision]
        ax.step(xx,yy,where='post',color=COLORS[method],ls=STYLES[method],lw=WIDTHS[method],zorder=ZORDER[method])
    baseline=float(frame.no_skill_precision.iloc[0])
    ax.axhline(baseline,color='#AAAAAA',ls=(0,(2,2)),lw=.65,zorder=0)
    ax.set(xlim=(0,1),ylim=(0,1.02),xlabel='Recall',ylabel='Precision' if ylabel else '')
    ax.set_xticks([0,.5,1]);ax.set_yticks([0,.5,1])
    if not ylabel:ax.set_yticklabels([])
    ax.tick_params(pad=2);ax.xaxis.labelpad=2.5
    label={'rmsf':'RMSF','chi_entropy':'Hχ','within_type_hchi':'Hχ within type'}[task]
    title(fig,height,x,heading_y,letter,label)
    txt(fig,height,x-1,ap_y,'AP',fontsize=9,ha='left',va='center',color='#444444')
    # Same method ordering as the single shared legend, regardless of panel width.
    for k,method in enumerate(METHODS):
        txt(fig,height,x+11+k*10.9,ap_y,f'{metrics[method]["mean_ap"]:.3f}',
            color=COLORS[method],fontsize=9,ha='center',va='center')

def pair(fig,height,x,y,size,task,letter,heading_y,annotation=None):
    ax=ax_mm(fig,height,x,y,size,size)
    p=frames[task].pivot(index='record_code',columns='method',values='ap').sort_index()
    assert len(p)==3063
    ax.plot([0,1],[0,1],color=MD,lw=.7,ls=(0,(3,3)),zorder=0)
    ax.scatter(p.Static,p.Full,s=2,color=BLUE,alpha=.25,linewidths=0,rasterized=False,zorder=1)
    ax.set(xlim=(0,1),ylim=(0,1),xlabel='Static AP',ylabel='Full AP')
    ax.set_aspect('equal',adjustable='box');ax.set_xticks([0,.5,1]);ax.set_yticks([0,.5,1])
    number_format(ax,'x');number_format(ax,'y')
    ax.tick_params(pad=2);ax.xaxis.labelpad=3;ax.yaxis.labelpad=3
    ax.xaxis.label.set_color(TEAL);ax.yaxis.label.set_color(BLUE)
    title(fig,height,x-2,heading_y,letter,'Hχ within type: record-level AP')
    if annotation:
        xx,yy=annotation
        wins=(p.Full>p.Static).sum()
        txt(fig,height,xx,yy,'Full higher in',color=BLUE,fontsize=9,ha='left',va='bottom')
        txt(fig,height,xx,yy-4.6,f'{100*wins/3063:.1f}% of records',color=BLUE,fontsize=9,ha='left',va='bottom')
    return ax

# B: two compact numerical panels, two global PR panels, then a larger within-type comparison.
fig,h=canvas(204.)
shared_key(fig,h)
density(fig,h,20,134,52,'rmsf','a',190)
density(fig,h,111,134,52,'chi_entropy','b',190)
precision_recall(fig,h,20,82,55,30,'rmsf','c',120,115,True)
precision_recall(fig,h,111,82,55,30,'chi_entropy','d',120,115,True)
precision_recall(fig,h,20,14,73,44,'within_type_hchi','e',66,61,True)
paired_axis=pair(fig,h,125,14,43,'within_type_hchi','f',66)
txt(fig,h,125,61,'Full higher: 81.2%',color=BLUE,fontsize=9,ha='left',va='center')
option_b=export(fig,'option_B_within_type',W,h)
plt.close(fig)

# C: all 3,063 records for each method and task are visible as points and distribution shapes.
distribution_points=[]
distribution_summary=[]
fig,h=canvas(203.)
shared_key(fig,h)
density(fig,h,20,139,49,'rmsf','a',192)
density(fig,h,111,139,49,'chi_entropy','b',192)
for j,task in enumerate(TASKS):
    x,width=[(15,43),(70,43),(125,51)][j]
    ax=ax_mm(fig,h,x,75,width,38)
    frame=frames[task]
    arrays=[frame[frame.method.eq(method)].sort_values('record_code').ap.to_numpy() for method in METHODS]
    assert all(len(v)==3063 for v in arrays)
    parts=ax.violinplot(arrays,positions=[0,1,2,3],widths=.78,showmeans=False,showmedians=False,
                        showextrema=False,points=128,bw_method=.20)
    for k,(method,values,body) in enumerate(zip(METHODS,arrays,parts['bodies'])):
        body.set_facecolor(COLORS[method]);body.set_edgecolor(COLORS[method]);body.set_linewidth(.8);body.set_alpha(.25)
        codes=frame[frame.method.eq(method)].sort_values('record_code').record_code.to_numpy()
        horizontal_offset=.13*np.sin(codes*2.399963229728653)
        ax.scatter(k+horizontal_offset,values,s=.45,color=COLORS[method],alpha=.09,linewidths=0,zorder=2,rasterized=False)
        quartiles=np.quantile(values,[.25,.5,.75])
        ax.vlines(k,quartiles[0],quartiles[2],color=COLORS[method],lw=2.5,zorder=3)
        ax.plot([k-.085,k+.085],[quartiles[1],quartiles[1]],color='white',lw=1,zorder=4)
        distribution_summary.append(dict(task=task,method=method,n_records=3063,mean=float(values.mean()),
                                         q25=float(quartiles[0]),median=float(quartiles[1]),q75=float(quartiles[2]),
                                         minimum=float(values.min()),maximum=float(values.max())))
        distribution_points.extend(dict(task=task,method=method,record_code=int(code),AP=float(value)) for code,value in zip(codes,values))
    ax.set(xlim=(-.5,3.5),ylim=(-.02,1.02),ylabel='Record-level AP' if j==0 else '')
    ax.set_xticks([0,1,2,3],METHODS);ax.set_yticks([0,.5,1]);ax.tick_params(pad=2)
    for tick,method in zip(ax.get_xticklabels(),METHODS):tick.set_color(COLORS[method])
    if j>0:ax.set_yticklabels([])
    title(fig,h,x,122,'cde'[j],['RMSF','Hχ','Hχ within type'][j])
for j,task in enumerate(TASKS):
    x,width=[(15,43),(70,43),(125,51)][j]
    precision_recall(fig,h,x,14,width,36,task,'fgh'[j],58,53,j==0)
option_c=export(fig,'option_C_distributions',W,h)
plt.close(fig)
pd.DataFrame(distribution_points).to_csv(ROOT/'all_record_AP_distribution_points.csv',index=False)
pd.DataFrame(distribution_summary).to_csv(ROOT/'AP_distribution_summary.csv',index=False)

# Print-sized previews are PDF renders, not raster edits of source figures.
options=[('A','option_A_balanced',203.),('B','option_B_within_type',204.),('C','option_C_distributions',203.)]
for letter,stem,height in options:
    doc=fitz.open(ROOT/f'{stem}.pdf')
    factor=200/72*167.64/183
    doc[0].get_pixmap(matrix=fitz.Matrix(factor,factor),alpha=False).save(ROOT/f'{letter}_print_preview.png')
    doc.close()

sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
audit=dict(status='numeric_canvas_passed_visual_pending',options=dict(A=state['main_check'],B=option_b,C=option_c),
           same_density_data=True,same_density_axes=True,same_density_normalization=True,
           all_three_tasks_complete_PR_in_every_option=True,frozen_AP_values_unchanged=True,
           all_3063_records_in_each_distribution=True,distribution_points=len(distribution_points),
           violin_density='Gaussian KDE for shape only; bw_method=0.20, 128 evaluation points, clipped to empirical min/max.',
           violin_raw_points='Every raw AP preserved; deterministic horizontal offset only.',
           distribution_intervals='25th-75th percentiles and median, not confidence intervals.',
           no_new_significance_tests=True,manuscript_replaced=False,
           output_sha256={p.name:sha(p) for p in ROOT.iterdir() if p.name.startswith(('option_','A_print','B_print','C_print'))})
(ROOT/'options_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(audit,ensure_ascii=False,indent=2))
