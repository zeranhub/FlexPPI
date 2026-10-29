# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (A) All-residue density -> MarginalDensity/plot_VariableCorrelation.py -> param inherit
# (B,C) Full protein distributions -> Violin/plot_Violin.py -> param inherit
# (D) Paired fold baseline comparison -> cross-type inherit -> param inherit
# Assets inspected with companion PNGs. The genomic filtering/log regression and
# categorical significance tests do not answer held-out prediction accuracy.
# Inherit clean axes, bounded violins and explicit count labels; no data filters.

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
import hashlib
import json
from pathlib import Path

mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LogNorm
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd
from scipy.stats import pearsonr

ROOT = Path(__file__).resolve().parents[1]
MM = 1/25.4
SOURCES = ['atlas', 'dynarepo', 'itrace']
NAMES = {'atlas':'ATLAS', 'dynarepo':'DynaRepo', 'itrace':'ITrace'}
COLORS = {'atlas':'#2166AC', 'dynarepo':'#C78124', 'itrace':'#762A83'}
CMAP = LinearSegmentedColormap.from_list('residue_count', SEQUENTIAL)


def stats(pred, target):
    p, y = np.asarray(pred, float), np.asarray(target, float)
    if len(y)<2 or not np.isfinite(p).all() or not np.isfinite(y).all():
        raise ValueError('Invalid paired regression input')
    return dict(n_residues=len(y), pearson=float(pearsonr(p,y).statistic),
                mae=float(np.mean(np.abs(p-y))), rmse=float(np.sqrt(np.mean((p-y)**2))),
                bias=float(np.mean(p-y)))


def export(fig, out, name):
    paths=[]
    with mpl.rc_context({'savefig.bbox':None}):
        for extension in ['png','pdf','svg']:
            path=out/f'{name}.{extension}'
            fig.savefig(path,dpi=300,bbox_inches=None,facecolor='white')
            paths.append(str(path))
    plt.close(fig)
    return paths


def heading(ax, letter, title, subtitle):
    ax.set_title(title,loc='left',pad=27,fontsize=8.5)
    ax.text(0,1.035,subtitle,transform=ax.transAxes,fontsize=6.6,va='bottom')
    ax.text(-.17,1.16,letter,transform=ax.transAxes,fontsize=11,fontweight='bold')


def density(fig, ax, pred, target, label, logarithmic=False, short=False):
    result=stats(pred,target)
    p,y=np.asarray(pred,float),np.asarray(target,float)
    if logarithmic:
        if np.min(p)<=0 or np.min(y)<=0:
            raise ValueError('Log density requires positive RMSF')
        p,y=np.log10(p),np.log10(y)
        lo=min(np.log10(.2),float(min(p.min(),y.min())))
        hi=max(np.log10(36),float(max(p.max(),y.max())))
        ticks=np.array([.25,.5,1,2,4,8,16,32])
        ax.set_xticks(np.log10(ticks),[f'{v:g}' for v in ticks],fontsize=6)
        ax.set_yticks(np.log10(ticks),[f'{v:g}' for v in ticks],fontsize=6)
    else:
        lo=min(0,float(min(p.min(),y.min())))
        hi=float(max(p.max(),y.max()))
        pad=(hi-lo)*.025
        lo-=pad;hi+=pad
        ax.xaxis.set_major_locator(MaxNLocator(4))
        ax.yaxis.set_major_locator(MaxNLocator(4))
    edges=np.linspace(lo,hi,91)
    counts,_,_=np.histogram2d(y,p,bins=(edges,edges))
    if int(counts.sum()) != len(y):
        raise ValueError('Density bins lost residues')
    artist=ax.pcolormesh(edges,edges,np.ma.masked_where(counts.T==0,counts.T),
                         cmap=CMAP,norm=LogNorm(vmin=1,vmax=counts.max()),
                         shading='flat',rasterized=True)
    ax.plot([lo,hi],[lo,hi],color='#5F5F5F',lw=.7,ls='--')
    ax.set(xlim=(lo,hi),ylim=(lo,hi),xlabel=f'MD {label}',ylabel=f'Predicted {label}')
    ax.text(.025,.97,f"r = {result['pearson']:.3f}\nMAE = {result['mae']:.3f}",
            va='top',transform=ax.transAxes,fontsize=7,
            bbox=dict(facecolor='white',edgecolor='none',alpha=.85,pad=2))
    bar=fig.colorbar(artist,ax=ax,fraction=.035,pad=.025,shrink=.82)
    bar.set_label('Residues per bin',fontsize=6.5)
    bar.ax.tick_params(labelsize=6,length=2)
    bar.outline.set_visible(False)
    ax.set_aspect('equal',adjustable='box')
    result.update(histogram_count=int(counts.sum()),log_axis=logarithmic,
                  range_transformed=[lo,hi],plotted_residue_fraction=1.0)
    return result


def distribution(ax, groups, names, colors, logarithmic=False, limits=None):
    labels=[]; summaries=[]
    for i,(values,name,color) in enumerate(zip(groups,names,colors)):
        values=np.asarray(values,float)
        if not np.isfinite(values).all() or (logarithmic and np.min(values)<=0):
            raise ValueError('Invalid per-system metric')
        v=np.log10(values) if logarithmic else values
        violin=ax.violinplot([v],positions=[i],widths=.78,showextrema=False,
                             points=180,bw_method=.20)
        body=violin['bodies'][0]
        body.set_facecolor(color);body.set_edgecolor(color);body.set_alpha(.18)
        body.set_linewidth(.7)
        # Horizontal offsets only spread overlapping dots; every sample is drawn.
        jitter=((np.arange(len(v))*.6180339887498949)%1-.5)*.44
        ax.scatter(i+jitter,v,s=3.4,alpha=.20,color=color,lw=0,rasterized=True,zorder=2)
        q=np.quantile(v,[.25,.5,.75])
        ax.vlines(i,q[0],q[2],color=BLACK,lw=2.8,zorder=3)
        ax.scatter([i],[q[1]],s=12,facecolor='white',edgecolor=BLACK,lw=.65,zorder=4)
        labels.append(f'{name}\nn = {len(v):,}')
        summaries.append(dict(group=name,n_samples=len(values),median=float(np.median(values)),
                               mean=float(values.mean()),minimum=float(values.min()),maximum=float(values.max())))
    ax.set_xticks(range(len(groups)),labels,fontsize=7)
    ax.set_xlim(-.62,len(groups)-.38)
    if logarithmic:
        ticks=np.array([.05,.1,.2,.5,1,2,5,10,20])
        data=np.concatenate(groups);low=float(np.log10(data.min()))-.07;high=float(np.log10(data.max()))+.07
        ticks=ticks[(np.log10(ticks)>=low)&(np.log10(ticks)<=high)]
        ax.set_yticks(np.log10(ticks),[f'{v:g}' for v in ticks])
        ax.set_ylim(low,high)
    elif limits:
        ax.set_ylim(*limits)
    ax.tick_params(axis='x',length=0,pad=5)
    return summaries


def load_and_audit(analysis,out):
    summary=json.loads((analysis/'test_summary.json').read_text())
    if summary['status']!='complete' or summary['systems']!=3068:
        raise ValueError('Need completed full evaluation')
    columns=['sample_id','dataset','fold','chain','author_residue_id','amino_acid',
             'valid_rmsf','md_rmsf_angstrom','pred_rmsf_angstrom',
             'valid_chi_entropy','md_chi_entropy_nats','pred_chi_entropy_nats']
    data=pd.read_csv(analysis/'all_residue_predictions.csv.gz',usecols=columns,
                     dtype={'sample_id':str,'chain':str,'author_residue_id':str,'amino_acid':str})
    if len(data)!=959347 or data.duplicated(['sample_id','chain','author_residue_id']).any():
        raise ValueError('Unexpected coverage or duplicate residue identity')
    if int(data.valid_rmsf.sum())!=959347 or int(data.valid_chi_entropy.sum())!=431023:
        raise ValueError('Valid-residue count mismatch')
    systems=pd.read_csv(analysis/'per_system_metrics.csv')
    rmsf=systems.loc[systems.metric=='rmsf'].copy()
    if len(rmsf)!=3068 or rmsf.sample_id.duplicated().any():
        raise ValueError('Expected one RMSF result per sample')
    data['baseline_abs_error']=np.nan
    priors=[]
    for fold in range(5):
        train=data.fold!=fold
        value=float(data.loc[train,'md_rmsf_angstrom'].mean())
        server=json.loads((ROOT/f'artifacts/performance_first_v3_difficulty_balanced/evaluation/fold_{fold}_test.json').read_text())
        expected=server['metrics']['rmsf']['training_mean_baseline']['fitted_mean']
        if abs(value-expected)>1e-6:
            raise ValueError('Training-fold mean mismatch')
        data.loc[~train,'baseline_abs_error']=np.abs(data.loc[~train,'md_rmsf_angstrom']-value)
        priors.append(dict(fold=fold,training_only_rmsf_mean=value,n_train_residues=int(train.sum())))
    baseline=data.groupby('sample_id').baseline_abs_error.mean().rename('baseline_mae')
    rmsf=rmsf.join(baseline,on='sample_id',validate='one_to_one')
    if rmsf.baseline_mae.isna().any():
        raise ValueError('Missing baseline')
    rmsf['relative_mae_reduction']=1-rmsf.mae/rmsf.baseline_mae
    rmsf.to_csv(out/'all_protein_rmsf_metrics_and_baseline.csv',index=False)
    pd.DataFrame(priors).to_csv(out/'training_only_rmsf_priors.csv',index=False)
    folds=rmsf.groupby('fold')[['mae','baseline_mae','pearson']].mean().reset_index()
    folds.to_csv(out/'paired_five_fold_comparison.csv',index=False)
    hc=data.loc[data.valid_chi_entropy].copy()
    aa_priors=pd.read_csv(analysis/'hchi_training_only_type_priors.csv')
    # The mapping is the mean of that amino-acid type in the OTHER four folds.
    hc=hc.merge(aa_priors[['fold','amino_acid','mean_training_hchi']],
                on=['fold','amino_acid'],validate='many_to_one',how='left')
    if len(hc)!=431023 or hc.mean_training_hchi.isna().any():
        raise ValueError('Hchi training prior mapping lost residues')
    hc['md_corrected']=hc.md_chi_entropy_nats-hc.mean_training_hchi
    hc['pred_corrected']=hc.pred_chi_entropy_nats-hc.mean_training_hchi
    hrows=[]
    for sid,g in hc.groupby('sample_id',sort=True):
        result=stats(g.pred_corrected,g.md_corrected)
        hrows.append(dict(sample_id=sid,dataset=g.dataset.iloc[0],fold=int(g.fold.iloc[0]),
                          **result))
    hsystems=pd.DataFrame(hrows)
    hsystems.to_csv(out/'all_protein_hchi_type_corrected_metrics.csv',index=False)
    reference=pd.read_csv(analysis/'hchi_type_corrected_metrics.csv')
    for fold in range(5):
        value=hsystems.loc[hsystems.fold==fold].pearson.mean()
        target=reference.loc[(reference.fold==fold)&(reference.scope=='macro_by_system'),'pearson'].iloc[0]
        if abs(value-target)>1e-6:
            raise ValueError('Hchi corrected sample statistics mismatch')
    bins=np.array([0,.5,1,2,3,5,10,np.inf])
    rows=[]
    for low,high in zip(bins[:-1],bins[1:]):
        g=data.loc[data.md_rmsf_angstrom.ge(low)&data.md_rmsf_angstrom.lt(high)]
        rows.append(dict(md_min=low,md_max=high,n_residues=len(g),
                         mean_md=float(g.md_rmsf_angstrom.mean()),mean_prediction=float(g.pred_rmsf_angstrom.mean()),
                         bias=float((g.pred_rmsf_angstrom-g.md_rmsf_angstrom).mean())))
    pd.DataFrame(rows).to_csv(out/'rmsf_amplitude_bias_by_md_range.csv',index=False)
    return data,rmsf,folds,hc,hsystems,systems


def rmsf_figure(data,systems,folds,out):
    fig=plt.figure(figsize=(183*MM,163*MM))
    gs=fig.add_gridspec(2,2,left=.095,right=.98,bottom=.095,top=.83,hspace=.66,wspace=.40)
    axes=[fig.add_subplot(gs[i,j]) for i in range(2) for j in range(2)]
    fig.text(.53,.984,'FlexPPI-v3 · overall RMSF prediction',ha='center',va='top',fontsize=11)
    fig.text(.53,.949,'3,068 held-out protein samples · 959,347 residues · fixed epoch-40 models',
             ha='center',va='top',fontsize=7.8)
    audit={}
    heading(axes[0],'A','Absolute-value agreement','All valid residues · both axes log-scaled')
    audit['density']=density(fig,axes[0],data.pred_rmsf_angstrom,data.md_rmsf_angstrom,'RMSF [Å]',True)
    heading(axes[1],'B','Profile agreement by protein','One dot per sample; white: median, black: IQR')
    groups=[systems.loc[systems.dataset==s,'pearson'].to_numpy() for s in SOURCES]
    audit['per_protein_r']=distribution(axes[1],groups,[NAMES[s] for s in SOURCES],
                                      [COLORS[s] for s in SOURCES],limits=(0,1.04))
    axes[1].set_ylabel('Per-protein Pearson r')
    heading(axes[2],'C','Absolute error by protein','All samples shown · logarithmic error axis')
    groups=[systems.loc[systems.dataset==s,'mae'].to_numpy() for s in SOURCES]
    audit['per_protein_mae']=distribution(axes[2],groups,[NAMES[s] for s in SOURCES],
                                        [COLORS[s] for s in SOURCES],logarithmic=True)
    axes[2].set_ylabel('Per-protein MAE [Å]')
    heading(axes[3],'D','Improvement over a simple baseline','Each connected pair: one held-out fold')
    for _,row in folds.iterrows():
        offset=(int(row.fold)-2)*.045
        axes[3].plot([offset,1+offset],[row.baseline_mae,row.mae],color='#A5A9AF',lw=.75,zorder=1)
        axes[3].scatter([offset],[row.baseline_mae],color='#7A7F87',s=24,zorder=3)
        axes[3].scatter([1+offset],[row.mae],color='#1B9474',s=24,zorder=3)
    means=folds[['baseline_mae','mae']].mean().to_numpy()
    for x,value in enumerate(means):
        axes[3].plot([x-.16,x+.16],[value,value],color=BLACK,lw=1.5,zorder=4)
        axes[3].text(x,value+.065,f'{value:.3f} Å',ha='center',fontsize=7.5)
    improvement=float(1-means[1]/means[0])
    axes[3].text(.98,.98,f'{improvement:.1%} lower MAE',ha='right',va='top',
                 transform=axes[3].transAxes,color='#15785E',fontsize=9)
    axes[3].set_xticks([0,1],['Training-fold\nmean predictor','FlexPPI-v3'],fontsize=7.5)
    axes[3].set(xlim=(-.45,1.45),ylim=(0,1.15),ylabel='Mean per-protein MAE [Å]')
    axes[3].tick_params(axis='x',length=0)
    audit['baseline']=dict(fold_mean_baseline_mae=float(means[0]),fold_mean_model_mae=float(means[1]),
                           relative_reduction=improvement,folds=folds.to_dict('records'))
    fig.text(.53,.025,'Every sample is predicted only by its held-out fold model. No illustrative-case selection.',
             fontsize=7,ha='center')
    audit['paths']=export(fig,out,'fig_v3_overall_rmsf_capability')
    return audit


def hchi_figure(hc,hsystems,systems,out):
    fig=plt.figure(figsize=(183*MM,87*MM))
    gs=fig.add_gridspec(1,3,left=.065,right=.985,bottom=.24,top=.69,wspace=.64)
    axes=[fig.add_subplot(gs[0,i]) for i in range(3)]
    fig.text(.53,.983,'Overall side-chain entropy prediction',ha='center',va='top',fontsize=10.5)
    fig.text(.53,.92,'1,139 protein samples with Hχ labels · 431,023 residues · ATLAS has no Hχ labels in this run',
             ha='center',va='top',fontsize=7.4)
    heading(axes[0],'A','Raw Hχ agreement','All valid Hχ residues')
    raw=density(fig,axes[0],hc.pred_chi_entropy_nats,hc.md_chi_entropy_nats,'Hχ [nat]',short=True)
    heading(axes[1],'B','After amino-acid correction','Training-fold type mean removed')
    corrected=density(fig,axes[1],hc.pred_corrected,hc.md_corrected,'Hχ residual [nat]',short=True)
    heading(axes[2],'C','Corrected profile agreement','Each dot: one protein sample')
    groups=[hsystems.loc[hsystems.dataset==s,'pearson'].to_numpy() for s in SOURCES[1:]]
    summary=distribution(axes[2],groups,[NAMES[s] for s in SOURCES[1:]],
                         [COLORS[s] for s in SOURCES[1:]],limits=(min(0,hsystems.pearson.min()-.03),1.03))
    axes[2].set_ylabel('Per-protein corrected Pearson r',fontsize=7)
    fig.text(.53,.055,'Correction tests differences beyond the average entropy of each amino-acid type; no test-fold labels fit the mean.',
             fontsize=6.5,ha='center')
    return dict(raw=raw,corrected=corrected,protein_corrected=summary,
                 paths=export(fig,out,'fig_v3_overall_hchi_capability'))


def report(data,rmsf,folds,hc,hsystems,audit,out,analysis):
    thresholds=[dict(threshold=value,count=int((rmsf.pearson>=value).sum()),
                     fraction=float((rmsf.pearson>=value).mean())) for value in [.7,.8,.9]]
    baseline_improvement=1-rmsf.mae.mean()/rmsf.baseline_mae.mean()
    wins=int((rmsf.mae<rmsf.baseline_mae).sum())
    metrics=dict(rmsf_pooled=stats(data.pred_rmsf_angstrom,data.md_rmsf_angstrom),
                 rmsf_macro_pearson=float(rmsf.pearson.mean()),rmsf_macro_mae=float(rmsf.mae.mean()),
                 rmsf_median_pearson=float(rmsf.pearson.median()),rmsf_median_mae=float(rmsf.mae.median()),
                 rmsf_thresholds=thresholds,baseline_mae=float(rmsf.baseline_mae.mean()),
                 baseline_relative_mae_reduction=float(baseline_improvement),baseline_win_samples=wins,
                 baseline_win_fraction=wins/len(rmsf),
                 hchi_corrected_macro_pearson=float(hsystems.pearson.mean()))
    audit['summary']=metrics
    audit['coverage']=dict(rmsf_samples=3068,rmsf_residues=959347,hchi_samples=1139,hchi_residues=431023,
                           unique_biological_systems=int(rmsf.biological_system.nunique()),
                           rmsf_samples_dropped=0,hchi_missing_label_residues=528324)
    audit['source_files']={}
    for name in ['all_residue_predictions.csv.gz','per_system_metrics.csv','test_summary.json',
                 'hchi_training_only_type_priors.csv']:
        path=analysis/name
        audit['source_files'][name]=dict(path=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest())
    audit['visual_qa']='Pending final rendered inspection'
    (out/'overall_statistics_and_QA.json').write_text(json.dumps(audit,indent=2,ensure_ascii=False),encoding='utf-8')
    rows=[]
    for dataset,g in rmsf.groupby('dataset'):
        rows.append(dict(dataset=dataset,n_samples=len(g),mean_pearson=float(g.pearson.mean()),
                          median_pearson=float(g.pearson.median()),mean_mae=float(g.mae.mean()),
                          median_mae=float(g.mae.median()),baseline_mae=float(g.baseline_mae.mean())))
    pd.DataFrame(rows).to_csv(out/'dataset_summary.csv',index=False)
    text=f'''# FlexPPI 整体能力评价

这组图使用所有符合已冻结标签掩码的测试数据，没有按拟合效果挑选案例。RMSF 共 3,068 个蛋白样本、959,347 个残基，Hχ 共 1,139 个有标签的样本、431,023 个残基。3,068 是 sample_id 数，对应 {rmsf.biological_system.nunique():,} 个不同 biological_system，不代表 3,068 个彼此独立的蛋白实体。每个样本只使用所属测试折的 epoch_040.pt；未重新训练、归一化目标或更改数据划分。

## 怎样读主图

A 把所有残基的 MD RMSF 与预测 RMSF 放在一起。横坐标为 MD 标签，纵坐标为模型预测；虚线表示二者相等，颜色越深表示该区域的残基越多。两轴用对数间隔展开小数值、压缩大数值，但刻度仍为真实 Å，所有残基都纳入计数，指标使用未经变换的原始值。此处合并全体预测后 r={metrics['rmsf_pooled']['pearson']:.3f}，MAE={metrics['rmsf_pooled']['mae']:.3f} Å；它们不是五折 r 的算术平均。

B 每个点是一整个蛋白样本的相关系数，先在每个蛋白内比较各残基的峰谷，再展示所有蛋白的分布。三种颜色区分 ATLAS、DynaRepo 和 ITrace。形状宽处表示蛋白较多，黑线为中间 50% 的范围，白点为中位数；所有样本散点均绘出。全部蛋白等权平均 r={rmsf.pearson.mean():.3f}，中位数 r={rmsf.pearson.median():.3f}；{thresholds[0]['count']:,}/{len(rmsf):,}（{thresholds[0]['fraction']:.1%}）的蛋白达到 r≥0.7。

C 每个点是一个蛋白的平均绝对误差，越低越好。这样既能看到大多数蛋白的表现，也保留误差较大的长尾；纵轴为对数间隔。整体每个蛋白等权平均 MAE={rmsf.mae.mean():.3f} Å，中位数为 {rmsf.mae.median():.3f} Å。每个蛋白等权与每个残基等权是两个不同的问题：后者让大蛋白贡献更多权重。

D 与简单基线比较。基线在每个测试折都预测同一个数，该数只取另外四折训练标签的平均 RMSF；没有用测试标签拟合。每条线连接同一测试折的基线和模型误差，各折内部每个蛋白等权。五折均值从 {folds.baseline_mae.mean():.3f} Å 降到 {folds.mae.mean():.3f} Å，相对降低 {1-folds.mae.mean()/folds.baseline_mae.mean():.1%}。这是相对于训练均值基线的改进，不是优于其它已发表模型的证据。全体样本中有 {wins:,}/{len(rmsf):,}（{wins/len(rmsf):.1%}）的模型误差低于该基线。

## 模型目前的能力与不足

模型能够在多数蛋白上识别残基波动的相对高低，但绝对幅度仍有系统性低估。全体残基的平均偏差（预测减 MD）为 {metrics['rmsf_pooled']['bias']:.3f} Å；高波动残基尤其容易被低估。因此论文中可以支持“模型捕捉多数蛋白的柔性分布并优于训练均值基线”，不能仅凭相关性声称“准确复现全部运动幅度”。不同来源的性能不同，需要连同样本数一起展示。

Hχ 附图同时给出原始与氨基酸类别校正后的结果。不同氨基酸的侧链天然具有不同自由度，原始高相关性一部分来自这种差异。减去该氨基酸在训练折中的平均值后，全体残基 r={audit['hchi']['corrected']['pearson']:.3f}，每个蛋白等权平均 r={hsystems.pearson.mean():.3f}。校正均值仅来自各自训练折，MAE 不因两边减去相同均值而改变。本轮 ATLAS 没有 Hχ 标签，不能把缺失当作零加入评价。

## 图注与复现

主图：全测试集 RMSF 能力。A，原始 RMSF 值的残基计数密度与相等参考线；B，按来源展示全部蛋白的 Pearson r；C，按来源展示全部蛋白的 MAE；D，五个留出折上的训练均值基线与 FlexPPI 平均 MAE 对照。B/C 中白点为中位数、粗线为四分位范围，所有散点均保留；D 中黑色短线为五折均值，连线不表示独立重复实验。无显著性检验、无置信区间。

Hχ 附图：A，全部有效残基的原始 Hχ 密度；B，扣除训练集氨基酸类型均值后的密度；C，各蛋白的校正 Pearson r。所有密度图按格子累计真实残基数，不抽样；PDF/SVG 的文字、坐标及统计标记为矢量，密度和散点层栅格化以控制体积。两图宽 183 mm，PNG 300 dpi。

运行 `.venv-md/Scripts/python.exe scripts/plot_v3_overall_capability.py` 即可复现。图中逐蛋白统计、训练均值基线、五折比较、Hχ 校正后结果和按 MD 幅度分组的偏差均保存为 CSV。完整逐残基数据来自 `{analysis/'all_residue_predictions.csv.gz'}`，文件 SHA-256 记于 overall_statistics_and_QA.json。

本评价限于已冻结的 v3 五折协议，属于复合物组合层面的交叉验证，不等同于完全未见同源家族的外部泛化评价。v3 划分是在查看 v2 结果后建立，需作为独立版本说明。若要论证优于其它方法，还需在相同测试样本、标签定义、输入条件和掩码下做配对基线比较。
'''
    (out/'HOW_TO_READ_OVERALL_RESULTS.md').write_text(text,encoding='utf-8')
    print(json.dumps(metrics,ensure_ascii=False,indent=2))


def main():
    parser=argparse.ArgumentParser(description='Plot full v3 model capability from all held-out predictions')
    parser.add_argument('--analysis-dir',type=Path,default=ROOT/'artifacts/evaluation_v3/analysis')
    parser.add_argument('--output-dir',type=Path,default=ROOT/'output/flexppi_v3_overall_capability')
    args=parser.parse_args();args.output_dir.mkdir(parents=True,exist_ok=True)
    data,rmsf,folds,hc,hsystems,systems=load_and_audit(args.analysis_dir,args.output_dir)
    audit={'rmsf':rmsf_figure(data,rmsf,folds,args.output_dir),
           'hchi':hchi_figure(hc,hsystems,systems,args.output_dir)}
    report(data,rmsf,folds,hc,hsystems,audit,args.output_dir,args.analysis_dir)


if __name__=='__main__':
    main()
