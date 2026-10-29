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
from collections import defaultdict
from datetime import datetime,timezone
import csv,json,hashlib
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'result_figures';OUT.mkdir(exist_ok=True)
read=lambda p:list(csv.DictReader(p.open(encoding='utf-8-sig',newline='')))
PANELS=[f'MSKCC_TCR{i}' for i in range(1,8)]
TICKS=['TCR1 · CMV*','TCR2 · CMV','TCR3 · CMV','TCR4 · gp100','TCR5 · gp100','TCR6 · gp100','TCR7 · GRL']
source_tables=[]

def export(fig,name):
    fig.canvas.draw()
    save_cns_figure(fig,OUT/name);fig.savefig(OUT/(name+'.svg'),bbox_inches='tight')
    plt.close(fig)

def overview():
    candidate=read(ROOT/'evaluation/results_stronger_baselines/all_candidate_budget_metrics.csv')
    sites=read(ROOT/'evaluation/results_primary/all_site_nomination_metrics.csv')
    specs=[('a','Site nomination · 38 assays · gain vs contact-count rule',sites,'mode1',38,
        ['static_contact_degree','static_relative_SASA','parent_R_percentile','parent_H_percentile','parent_joint_percentile','parent_joint_low','parent_partner_context_response'],
        ['Contact\ncount','Sidechain\nexposure','RMSF\nhigh','Hχ\nhigh','Joint\nhigh','Joint\nlow','TCR\ncontext'],4),
      ('b','One-shot selection · 18 + 30 assays · gain vs sequence-only',candidate,'mode2',30,
        ['random','sequence_only_control','BLOSUM_only_control','static','baseline','baseline_R','baseline_H','baseline_RH','published_HERMES_fixed','published_ProteinMPNN'],
        ['Random','Sequence','BLOSUM','Static\nscores','Sequence\n+ static','+ RMSF','+ Hχ','+ Both','HERMES\nzero-shot','MPNN\nzero-shot'],7),
      ('c','Feedback selection · 20 + 40 assays · gain vs sequence-only',candidate,'mode3',40,
        ['random','sequence_only_control','BLOSUM_only_control','baseline','baseline_R','baseline_H','baseline_RH'],
        ['Random','Sequence','BLOSUM','Sequence\n+ static','+ RMSF','+ Hχ','+ Both'],6)]
    grids=[];sources=[]
    for letter,title,rows,mode,budget,methods,names,jointcol in specs:
        grid=[]
        for pid in PANELS:
            pp=[r for r in rows if r['panel_id']==pid and (r.get('mode')==mode)]
            pp=[r for r in pp if int(r['assay_budget'] if mode=='mode1' else r['extra_assays'])==budget]
            reference='static_contact_degree' if mode=='mode1' else 'sequence_only_control'
            rand=np.mean([float(r['new_hits']) for r in pp if r['strategy']==reference])
            values=[]
            for method in methods:
                v=np.mean([float(r['new_hits']) for r in pp if r['strategy']==method]);assert np.isfinite(v)
                values.append(v-rand)
                sources.append(dict(figure='three_modes_yield',panel=letter,mode=mode,panel_id=pid,strategy=method,
                    assay_budget=budget,mean_new_hits=v,reference_strategy=reference,mean_reference_hits=float(rand),gain_vs_reference=v-rand))
            grid.append(values)
        grids.append(np.asarray(grid))
    limit=max(1,max(float(np.max(abs(x))) for x in grids));limit=np.ceil(limit)
    fig,axes=plt.subplots(3,1,figsize=(183/25.4,198/25.4),gridspec_kw={'height_ratios':[1,1,1]})
    fig.subplots_adjust(left=.16,right=.93,top=.96,bottom=.095,hspace=.57)
    cmap=LinearSegmentedColormap.from_list('yield_gain',DIVERGING)
    for ax,spec,data in zip(axes,specs,grids):
        letter,title,rows,mode,budget,methods,names,jointcol=spec
        im=ax.imshow(data,cmap=cmap,vmin=-limit,vmax=limit,aspect='auto',interpolation='nearest')
        ax.set_yticks(range(7),TICKS);ax.set_xticks(range(len(names)),names)
        ax.tick_params(length=0,pad=4,labelsize=6.5)
        ax.set_title(title,loc='left',pad=9,fontsize=8)
        ax.text(-.12,1.12,letter,transform=ax.transAxes,weight='bold',fontsize=10)
        for j in (2.5,5.5):ax.axhline(j,color='white',lw=1.1)
        for i in range(7):
            for j in range(len(methods)):
                v=data[i,j];ax.text(j,i,f'{v:+.1f}',ha='center',va='center',fontsize=6.4,
                    color='white' if abs(v)>.62*limit else BLACK)
        ax.add_patch(Rectangle((jointcol-.5,-.5),1,7,fill=False,ec=BLACK,lw=1.1,clip_on=False))
        for sp in ax.spines.values():sp.set_visible(False)
    cax=fig.add_axes([.951,.225,.018,.53]);cb=fig.colorbar(im,cax=cax)
    cb.ax.tick_params(labelsize=6);cb.set_label('Difference in recognition-retaining variants vs stated reference',fontsize=7)
    fig.text(.16,.032,'Highest-dose normalized response ≥ 0.3. All 1,197 variants; 20 paired sampling seeds.\nBlack outline: prespecified joint strategy. *TCR1 source-sequence discrepancy; see report.',fontsize=6.5)
    export(fig,'three_modes_yield')
    source_tables.extend(sources)

def feedback():
    rows=read(ROOT/'evaluation/results_stronger_baselines/all_candidate_budget_metrics.csv')
    methods=['random','sequence_only_control','BLOSUM_only_control','baseline','baseline_H','baseline_RH']
    names=['Random','Sequence only','BLOSUM only','Sequence + static','+ Hχ','+ Both']
    colors=[GREY,CATEGORICAL[5],BLACK,CATEGORICAL[0],CATEGORICAL[4],CATEGORICAL[1]]
    fig,axes=plt.subplots(2,4,figsize=(183/25.4,113/25.4),sharex=True,sharey=True)
    fig.subplots_adjust(left=.078,right=.99,top=.91,bottom=.16,wspace=.23,hspace=.48)
    for j,pid in enumerate(PANELS):
        ax=axes.flat[j];ax.set_title(TICKS[j],loc='left',fontsize=8)
        for method,name,color in zip(methods,names,colors):
            values=[]
            for budget in range(0,41,5):
                vv=[float(r['new_hits']) for r in rows if r['panel_id']==pid and r['mode']=='mode3' and r['strategy']==method and int(r['extra_assays'])==budget]
                assert len(vv)==20
                values.append(vv)
                source_tables.append(dict(figure='feedback_all_contexts',panel_id=pid,strategy=method,extra_assays=budget,
                    mean_new_hits=float(np.mean(vv)),seed_percentile10=float(np.quantile(vv,.1)),seed_percentile90=float(np.quantile(vv,.9))))
            values=np.asarray(values);x=np.arange(0,41,5);mu=values.mean(axis=1)
            ax.plot(x,mu,color=color,lw=1.4 if method=='baseline_RH' else 1.,ls='--' if method in {'random','BLOSUM_only_control'} else ':' if method=='sequence_only_control' else '-',label=name)
            if method in {'baseline','baseline_RH'}:
                ax.fill_between(x,np.quantile(values,.1,axis=1),np.quantile(values,.9,axis=1),color=color,alpha=.1,lw=0)
        ax.set_xlim(0,40);ax.set_ylim(0,40);ax.set_xticks([0,20,40]);ax.set_yticks([0,20,40]);ax.tick_params(labelsize=6)
        if j%4==0:ax.set_ylabel('New retained-response hits',fontsize=7)
        if j>=4:ax.set_xlabel('Additional assays',fontsize=7)
    ax=axes.flat[7];ax.axis('off')
    handles,labels=axes.flat[0].get_legend_handles_labels();ax.legend(handles,labels,loc='upper left',bbox_to_anchor=(-.05,1),fontsize=7)
    ax.text(0,.01,'Same 20 initial results\nEight batches × 5 candidates\nBands: 10–90% across seeds',transform=ax.transAxes,fontsize=6.5,va='bottom')
    fig.text(.078,.035,'All seven contexts shown. Seeds describe sampling variability, not independent biological replicates.',fontsize=6.5)
    export(fig,'feedback_all_contexts')

def controls():
    rows=read(ROOT/'evaluation/results_controls/all_candidate_budget_metrics.csv')
    stronger=read(ROOT/'evaluation/results_stronger_baselines/all_candidate_budget_metrics.csv')
    fig,axes=plt.subplots(1,2,figsize=(183/25.4,78/25.4),sharey=True)
    fig.subplots_adjust(left=.09,right=.98,bottom=.26,top=.85,wspace=.2)
    results=[]
    for ax,mode,budget,title in zip(axes,['mode2','mode3'],[30,40],['One-shot selection','Feedback selection']):
        ax.axhline(0,color=GREY,lw=.7,ls='--')
        for j,pid in enumerate(PANELS):
            subset=[r for r in rows if r['panel_id']==pid and r['mode']==mode and int(r['extra_assays'])==budget]
            by={(r['strategy'],int(r['seed'])):float(r['new_hits']) for r in subset}
            strong={(r['strategy'],int(r['seed'])):float(r['new_hits']) for r in stronger if r['panel_id']==pid and r['mode']==mode and int(r['extra_assays'])==budget}
            diffs=[[],[],[]]
            for seed in range(20261022,20261042):
                diffs[0].append(by['baseline_RH',seed]-np.mean([by[f'baseline_RH_shuffled{k:02d}',seed] for k in range(10)]))
                diffs[1].append(by['BLOSUM_static_RH',seed]-by['BLOSUM_static',seed])
                diffs[2].append(by['baseline_RH',seed]-strong['BLOSUM_only_control',seed])
            for k,(values,color,label) in enumerate(zip(diffs,[CATEGORICAL[1],CATEGORICAL[0],BLACK],['Joint vs shuffled profiles','Joint added to BLOSUM + static','Joint vs BLOSUM only'])):
                mean=float(np.mean(values));lo,hi=np.quantile(values,[.1,.9]);x=j+(k-1)*.2
                ax.vlines(x,lo,hi,color=color,lw=.8)
                ax.plot(x,mean,marker=['o','s','D'][k],ms=3.2,color=color,ls='',label=label if j==0 else None)
                results.append(dict(figure='profile_specificity_controls',panel_id=pid,mode=mode,contrast=label,
                    mean_extra_hits=mean,seed_percentile10=float(lo),seed_percentile90=float(hi)))
        ax.set_title(title,loc='left');ax.set_xticks(range(7),[f'TCR{i}' for i in range(1,8)],rotation=45,ha='right')
        ax.tick_params(labelsize=6.5);ax.set_xlim(-.6,6.6)
    axes[0].set_ylabel('Paired difference in new hits')
    handles,labels=axes[0].get_legend_handles_labels();fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.5,1.02),ncol=3,fontsize=6.5)
    fig.text(.09,.035,'Bars: 10–90% of 20 shared sampling seeds. Shuffled control averages 10 fixed profile permutations.',fontsize=6.5)
    export(fig,'profile_specificity_controls');source_tables.extend(results)

if __name__=='__main__':
    overview();feedback();controls()
    keys=list(dict.fromkeys(k for r in source_tables for k in r))
    with (OUT/'figure_source_values.csv').open('w',encoding='utf-8-sig',newline='') as h:
        w=csv.DictWriter(h,fieldnames=keys);w.writeheader();w.writerows(source_tables)
    report=dict(status='rendered_pending_personal_visual_QA',timestamp_utc=datetime.now(timezone.utc).isoformat(),
        all7contexts_and1197primary_variants_used=True,all_plotted_aggregates_exported=True,no_biological_pvalues_from_seeds=True,
        no_outcome_based_method_or_panel_selection=True,language='English',formats=['PDF embedded TrueType','SVG editable text','PNG300dpi'],
        archetype='quantitative_grid; three independent operational uses, plus all-context budget curves and profile-specificity controls',
        asset_scan='LineTrend/plot_trend.py and heatmap/plot_comparison.py inspected. Source charts are cumulative publication counts/benchmark summary examples; native statistical plot adapted to finite-pool selection metrics, preserving simple cumulative curves and signed comparison color logic.',
        source_values_sha256=hashlib.sha256((OUT/'figure_source_values.csv').read_bytes()).hexdigest())
    (OUT/'figure_QA.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf8')
    print(json.dumps(report))
