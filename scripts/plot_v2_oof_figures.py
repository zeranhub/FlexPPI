# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a) joint density -> MarginalDensity/plot_VariableCorrelation.py -> param inherit
# (b) fold comparison -> BarComparison/plot_comparison_Trajectory.py -> param inherit
# The marginal asset requires genomic quality groups and log-log regression;
# the comparison asset has embedded three-dimensional synthetic benchmark data.
# Neither data structure matches these full-coverage OOF predictions. Preserve
# the inspected marginal layout and separated metric panels, not their filters.
# RULE: "native run" = load pre-rendered PNG via Image.open().ax.imshow().
#       "param inherit" = drawing function below that copies Class A/B/C values.
#       If a panel says "native run" and you write a drawing function, you broke the contract.

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

import argparse
import json
from pathlib import Path

mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT=Path(__file__).resolve().parents[1]
CMAP=LinearSegmentedColormap.from_list('v2_density',SEQUENTIAL)
MM=1/25.4


def export(fig, folder, stem):
    save_cns_figure(fig,str(folder/stem))
    fig.savefig(folder/(stem+'.svg'),dpi=300,bbox_inches='tight')
    plt.close(fig)


def label(ax, text):
    ax.text(-0.15,1.07,text,transform=ax.transAxes,fontsize=9,fontweight='bold')


def joint(fig, cell, x, y, xlabel, ylabel, title, equal, letter):
    x,y=np.asarray(x,float),np.asarray(y,float)
    if not len(x) or not np.isfinite(x).all() or not np.isfinite(y).all():
        raise ValueError('Invalid density data; do not silently drop rows')
    gs=cell.subgridspec(2,2,width_ratios=[1,.12],height_ratios=[.12,1],hspace=.04,wspace=.04)
    ax=fig.add_subplot(gs[1,0]); top=fig.add_subplot(gs[0,0],sharex=ax)
    right=fig.add_subplot(gs[1,1],sharey=ax)
    xmax=max(float(x.max()),1e-6)*1.035;ymax=max(float(y.max()),1e-6)*1.035
    xmin=min(0,float(x.min())*1.035);ymin=min(0,float(y.min())*1.035)
    if equal:
        xmax=ymax=max(xmax,ymax);xmin=ymin=min(xmin,ymin)
    artist=ax.hexbin(x,y,gridsize=65,bins='log',mincnt=1,cmap=CMAP,
                     linewidths=0,extent=(xmin,xmax,ymin,ymax),rasterized=True)
    top.hist(x,bins=70,color=SEQUENTIAL[1],linewidth=0)
    right.hist(y,bins=70,color=SEQUENTIAL[1],linewidth=0,orientation='horizontal')
    top.axis('off');right.axis('off')
    ax.set(xlim=(xmin,xmax),ylim=(ymin,ymax),xlabel=xlabel,ylabel=ylabel)
    if equal:
        ax.plot([xmin,xmax],[ymin,ymax],color=GREY,lw=.6,ls='--',zorder=0)
    rho=float(spearmanr(x,y).statistic)
    top.set_title(title,loc='left',pad=7)
    ax.text(.97,.04,f'Pooled ρ = {rho:.3f}\nn = {len(x):,} residues',
            ha='right',va='bottom',transform=ax.transAxes,fontsize=6.5,
            bbox=dict(facecolor='white',edgecolor='none',alpha=.8))
    label(ax,letter)
    return dict(n_residues=len(x),pooled_spearman=rho,x_range=[xmin,xmax],y_range=[ymin,ymax],
                outside_axes=int(((x<xmin)|(x>xmax)|(y<ymin)|(y>ymax)).sum()))


def draw_densities(data, out):
    d=data.loc[data.dataset=='itrace'];audit={}
    fig=plt.figure(figsize=(183*MM,88*MM));gs=fig.add_gridspec(1,2,wspace=.5)
    for i,(metric,name,unit) in enumerate([('chi_entropy',r'$H_\chi$','nats'),('rmsf','RMSF','Å')]):
        s=d.loc[d['valid_'+metric]]
        audit[metric]=joint(fig,gs[i],s['target_'+metric],s['pred_'+metric],
                            f'MD {name} ({unit})',f'FlexPPI-v2 {name} ({unit})',
                            f'ITrace OOF · {s.sample_id.nunique()} complexes',True,'ab'[i])
    export(fig,out,'fig_v2_itrace_prediction_truth')
    s=d.loc[d.valid_rmsf & d.valid_chi_entropy]
    fig=plt.figure(figsize=(183*MM,88*MM));gs=fig.add_gridspec(1,2,wspace=.5)
    for i,(prefix,title) in enumerate([('target','MD truth'),('pred','FlexPPI-v2 OOF')]):
        audit['coupling_'+prefix]=joint(fig,gs[i],s[prefix+'_rmsf'],s[prefix+'_chi_entropy'],
                                        f'{title} RMSF (Å)',f'{title} '+r'$H_\chi$ (nats)',
                                        title,False,'ab'[i])
    export(fig,out,'fig_v2_itrace_rmsf_hchi_relation')
    return audit


def draw_folds(folds, out):
    fig,axes=plt.subplots(2,2,figsize=(183*MM,125*MM))
    for row,metric in enumerate(['rmsf','chi_entropy']):
        subset=folds.loc[folds.metric==metric]
        datasets=[ds for ds in ['atlas','dynarepo','itrace'] if ds in set(subset.dataset)]
        for col,stat in enumerate(['spearman','mae']):
            ax=axes[row,col]
            for i,ds in enumerate(datasets):
                points=subset.loc[subset.dataset==ds].sort_values('fold')[stat].to_numpy()
                if len(points)!=5:raise ValueError('Expected all five folds')
                ax.scatter(i+np.linspace(-.14,.14,5),points,s=18,color=CATEGORICAL[i],zorder=3)
                ax.errorbar(i,points.mean(),yerr=points.std(ddof=1),fmt='_',markersize=15,
                            color=BLACK,capsize=4,lw=.8,zorder=4)
            ax.set_xticks(range(len(datasets)),[s.upper() if s=='atlas' else s.title() for s in datasets])
            ax.set_ylabel('System-macro Spearman ρ' if stat=='spearman' else (
                'System-macro MAE (Å)' if metric=='rmsf' else 'System-macro MAE (nats)'))
            ax.set_title('RMSF' if metric=='rmsf' else r'$H_\chi$',loc='left')
            label(ax,'abcd'[row*2+col])
    fig.subplots_adjust(hspace=.55,wspace=.45,bottom=.13)
    fig.text(.5,.025,'Each point: one held-out fold · black: mean ± SD (ddof = 1) · not a confidence interval',
             ha='center',fontsize=6.5)
    export(fig,out,'fig_v2_five_fold_stability')


def draw_hchi_correction(folds,out):
    fig,axes=plt.subplots(1,2,figsize=(183*MM,80*MM))
    metrics=['chi_entropy','chi_type_prior','chi_residual']
    for idx,ds in enumerate(['dynarepo','itrace']):
        ax=axes[idx]
        for i,m in enumerate(metrics):
            points=folds.loc[(folds.dataset==ds)&(folds.metric==m)].sort_values('fold').spearman.to_numpy()
            if len(points)!=5:raise ValueError('Missing Hchi prior/correction folds')
            ax.scatter(i+np.linspace(-.12,.12,5),points,s=18,color=CATEGORICAL[i])
            ax.errorbar(i,points.mean(),yerr=points.std(ddof=1),fmt='_',color=BLACK,capsize=4,markersize=14)
        ax.set_xticks(range(3),['Raw prediction','AA-type prior','Type-subtracted\nresidual'])
        ax.set(title=ds.title(),ylabel='System-macro Spearman ρ')
        label(ax,'ab'[idx])
    fig.subplots_adjust(wspace=.4,bottom=.25)
    fig.text(.5,.045,'AA prior fitted only on training folds · residuals use a different target · mean ± fold SD',ha='center',fontsize=6.5)
    export(fig,out,'fig_v2_hchi_type_correction')


def draw_bridge_table(analysis,out):
    table=pd.read_csv(analysis/'bridge_macro_metrics.csv')
    metrics=['spearman','pearson','z_mae','top20_auroc','top20_average_precision','top20_precision_at_k']
    for stem,cohorts,title in [
        ('static',['dynarepo_complex','dynarepo_single_chain','atlas_single_chain'],'Legacy-cohort bridge · not full 3068-system benchmark'),
        ('itrace',['itrace_oof_bridge'],'ITrace 244 complexes · cached baselines vs v2 OOF')]:
        selected=table.loc[table.cohort.isin(cohorts)]
        rows=[];colors=[]
        for cohort in cohorts:
            part=selected.loc[selected.cohort==cohort].sort_values('spearman',ascending=False)
            for j,(_,r) in enumerate(part.iterrows()):
                rows.append([cohort.replace('_',' ') if j==0 else '',r.method]+[f'{r[k]:.4f}' for k in metrics])
                colors.append(['#F0F0F0' if r.method=='FlexPPI-v2' else 'white']*8)
        fig,ax=plt.subplots(figsize=(183*MM,(26+len(rows)*5)*MM));ax.axis('off')
        tab=ax.table(cellText=rows,colLabels=['Cohort','Method','ρ ↑','r ↑','z-MAE ↓','AUROC ↑','AP ↑','P@20% ↑'],
                     cellColours=colors,cellLoc='center',colWidths=[.24,.20]+[.0933]*6,bbox=[0,0,1,.92])
        tab.auto_set_font_size(False);tab.set_fontsize(6)
        for (r,c),cell in tab.get_celld().items():
            cell.set_linewidth(.25);cell.set_edgecolor('#DDDDDD')
            if c in [0,1]:cell.set_text_props(ha='left')
            if r==0:cell.set_text_props(weight='bold')
        ax.set_title(title,loc='left',fontsize=8,pad=8)
        export(fig,out,'table_v2_'+stem+'_baseline_bridge')


def main():
    p=argparse.ArgumentParser();p.add_argument('--analysis-dir',type=Path,default=ROOT/'artifacts/evaluation_v2/analysis')
    p.add_argument('--output-dir',type=Path,default=ROOT/'output/flexppi_v2_figures')
    args=p.parse_args();a=args.analysis_dir;o=args.output_dir;o.mkdir(parents=True,exist_ok=True)
    if json.loads((a/'analysis_summary.json').read_text())['status']!='complete':raise ValueError('Incomplete analysis')
    data=pd.read_csv(a/'oof_per_residue.csv.gz',dtype={'author_residue_id':str},keep_default_na=True)
    folds=pd.read_csv(a/'fold_metrics.csv')
    audit=draw_densities(data,o)
    draw_folds(folds,o);draw_hchi_correction(folds,o);draw_bridge_table(a,o)
    (o/'data_range_audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
    print('Generated six figures/tables in',o,flush=True)


if __name__=='__main__':main()
