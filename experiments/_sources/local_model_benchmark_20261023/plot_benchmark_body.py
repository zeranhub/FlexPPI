"""Plot every full-cohort completed method, with fixed rows across strata.

The journal reference is AlphaMissense Fig2A-C; this adapts its aligned method
rows and independent benchmark columns to points/intervals, not truncated bars.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'figures'
OUT.mkdir(exist_ok=True)
GROUPS=['single_chain','other_multichain','tcr_pmhc_itrace']
TITLES=['Single-chain proteins','Other multichain systems','TCR–pMHC complexes']
COUNTS=[2193,626,244]
ORDER=['Full','Intra','backflip-2.1 (assembly)','backflip-2.1 (chains)',
    'backflip-1.0 (assembly)','backflip-1.0 (chains)',
    'backflip-1.0-seq (assembly)','backflip-1.0-seq (chains)',
    'DynaProt-M','PLANET-MD','PLANET-MD-seq','PLANET-MD-mini','FastProtFlex','ANM','GNM','Static','AA']
DISPLAY={'Full':'FlexPPI','Intra':'FlexPPI · intra',
    'backflip-2.1 (assembly)':'BackFlip-2.1 · assembly','backflip-2.1 (chains)':'BackFlip-2.1 · chains',
    'backflip-1.0 (assembly)':'BackFlip-1 · assembly','backflip-1.0 (chains)':'BackFlip-1 · chains',
    'backflip-1.0-seq (assembly)':'BackFlip-1 + seq · assembly',
    'backflip-1.0-seq (chains)':'BackFlip-1 + seq · chains','Static':'Static features','AA':'Residue-type mean'}
# User's manuscript palette overrides the generic palette's semantic roles.
BLUE='#31688E'; OCHRE='#B78038'; TEAL='#348980'; INK='#253846'

def visual(name,meta):
    context=meta[name]['input_context']
    color=BLUE if name in ['Full','Intra'] else TEAL if name in ['ANM','GNM','Static','AA'] else OCHRE
    if name=='Full': return color,'s',3.9
    if context in ['independent_chain','within_chain','independent_chain_predictions_reassembled_per_record']:
        return color,'D',3.0
    if 'sequence' in context: return color,'^',3.2
    return color,'o',3.1

def main():
    summary=pd.read_csv(ROOT/'evaluation/group_summary.csv')
    audit=json.loads((ROOT/'evaluation/scoring_audit.json').read_text())
    meta=audit['methods']
    complete=set(audit['complete_methods'])
    methods=[x for x in ORDER if x in complete]
    unknown=complete-set(methods)
    assert not unknown, f'Add methods explicitly to fixed order: {unknown}'
    assert len(methods)>=7
    # Equal quantitative panels. No nonlinear axis, statistical zero retained.
    height=max(128,36+len(methods)*6.6)
    fig,axes=plt.subplots(3,3,figsize=(183/25.4,height/25.4),sharey=True)
    fig.subplots_adjust(left=.257,right=.986,bottom=.092,top=.914,wspace=.21,hspace=.34)
    metric_specs=[('spearman','Spearman ρ ↑',(-.1,1.0),[0,.5,1]),
        ('top10_ap','High-RMSF average precision ↑',(0,1.0),[0,.5,1]),
        ('mae','RMSF MAE (Å) ↓',None,None)]
    abs_values=summary[(summary.method.isin(methods))&(summary.metric=='mae')&(summary.group.isin(GROUPS))]
    max_mae=max(1,float(abs_values[['mean','ci_high']].max().max())*1.12)
    max_mae=float(np.ceil(max_mae*2)/2)
    source=[]
    for ri,(metric,label,lim,ticks) in enumerate(metric_specs):
        for ci,group in enumerate(GROUPS):
            ax=axes[ri,ci]
            if ri==0:
                ax.set_title(TITLES[ci]+f'\nn = {COUNTS[ci]:,}',fontsize=7.6,pad=9)
            for yi,name in enumerate(methods):
                part=summary[(summary.method==name)&(summary.group==group)&(summary.metric==metric)]
                if part.empty:
                    ax.text(.98,yi,'relative score',transform=ax.get_yaxis_transform(),ha='right',va='center',fontsize=5.4,color='#7C8589')
                    continue
                assert len(part)==1
                row=part.iloc[0]
                assert row.n_records==COUNTS[ci], (name,group,metric,row.n_records)
                color,marker,size=visual(name,meta)
                if np.isfinite(row.ci_low) and np.isfinite(row.ci_high):
                    ax.plot([row.ci_low,row.ci_high],[yi,yi],color=color,lw=.85,solid_capstyle='round',zorder=2)
                ax.plot(row['mean'],yi,marker=marker,ms=size,color=color,mec='white',mew=.35,ls='',zorder=3)
                source.append(row.to_dict())
            ax.axhspan(-.45,.45,color=BLUE,alpha=.065,lw=0,zorder=0)
            ax.set_ylim(len(methods)-.4,-.65)
            ax.set_yticks(range(len(methods)),[DISPLAY.get(n,n) for n in methods],fontsize=6.3)
            ax.tick_params(axis='y',length=0,pad=5,labelleft=ci==0)
            ax.tick_params(axis='x',labelsize=6.1,length=2.5,pad=2)
            ax.spines['left'].set_visible(False)
            ax.spines['bottom'].set_color('#80909A')
            if lim is None:
                ax.set_xlim(0,max_mae)
                ax.set_xticks([0,max_mae/2,max_mae])
            else:
                ax.set_xlim(*lim);ax.set_xticks(ticks)
            ax.set_xlabel(label,fontsize=6.8,labelpad=4)
            if ci==0:
                ax.text(-.345,1.06,chr(97+ri),transform=ax.transAxes,fontsize=9,fontweight='bold',va='bottom')
    handles=[Line2D([],[],color=BLUE,marker='o',lw=0,ms=3.3,label='FlexPPI configurations'),
        Line2D([],[],color=OCHRE,marker='o',lw=0,ms=3.3,label='Published learned models'),
        Line2D([],[],color=TEAL,marker='o',lw=0,ms=3.3,label='Graph / simple controls')]
    fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.59,.017),ncol=1,fontsize=6.1,handletextpad=.5,labelspacing=.35)
    fig.text(.012,.016,f'{len(methods)} completed model/configuration entries · same 3,063 records',fontsize=6.1,color=INK)
    save_cns_figure(fig,OUT/'benchmark_current')
    fig.savefig(OUT/'benchmark_current.svg',bbox_inches='tight',dpi=300)
    pd.DataFrame(source).to_csv(OUT/'benchmark_current_source.csv',index=False)
    # Record-distribution companion: all observations in empirical CDF, no downsampling.
    per=pd.read_csv(ROOT/'evaluation/per_record_metrics.csv')
    fig2,ax2=plt.subplots(1,3,figsize=(183/25.4,67/25.4))
    fig2.subplots_adjust(left=.08,right=.99,bottom=.25,top=.8,wspace=.25)
    selected=[x for x in ['Full','backflip-2.1 (assembly)','backflip-2.1 (chains)','DynaProt-M','FastProtFlex','ANM'] if x in complete]
    colors=[visual(name,meta)[0] for name in selected]
    styles=['-','-','--',':','-.','-']
    for ci,group in enumerate(GROUPS):
        ax=ax2[ci]
        for name,color,style in zip(selected,colors,styles):
            values=np.sort(per.loc[(per.method==name)&(per.group==group),'spearman'].to_numpy())
            assert len(values)==COUNTS[ci]
            ax.step(values,np.arange(1,len(values)+1)/len(values),where='post',color=color,lw=1,ls=style,label=DISPLAY.get(name,name))
        ax.set(xlim=(-.2,1),ylim=(0,1),xticks=[0,.5,1],yticks=[0,.5,1],xlabel='Per-record Spearman ρ',title=TITLES[ci])
        if ci==0:ax.set_ylabel('Fraction of records')
        ax.tick_params(labelsize=6.2)
        ax.title.set_fontsize(7)
        ax.xaxis.label.set_fontsize(7);ax.yaxis.label.set_fontsize(7)
    h,l=ax2[0].get_legend_handles_labels()
    fig2.legend(h,l,loc='lower center',bbox_to_anchor=(.5,-.02),ncol=3,fontsize=6)
    save_cns_figure(fig2,OUT/'benchmark_distributions_current')
    fig2.savefig(OUT/'benchmark_distributions_current.svg',bbox_inches='tight',dpi=300)
    (OUT/'plot_audit.json').write_text(json.dumps(dict(methods=methods,records=3063,primary_groups=dict(zip(GROUPS,COUNTS)),
        style_reference='Science AlphaMissense Fig2A-C: fixed method rows and independent benchmark columns; adapted to dot intervals.',
        all_observations_in_ecdf=True,no_score_based_selection=True,
        confidence_intervals='Record-weighted family-cluster bootstrap; TCR-pMHC single family descriptive points only.',
        absolute_mae_missing='Relative outputs are explicitly labeled; no test-fitted scaling.',
        progress_snapshot=True,unfinished_methods_in_primary_figure=False),indent=2),encoding='utf-8')
    print('Rendered',len(methods),'complete model/configuration entries',flush=True)

if __name__=='__main__': main()
