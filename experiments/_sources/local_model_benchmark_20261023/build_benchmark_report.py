"""Produce current factual report and paired differences from completed cohorts."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'evaluation'
audit=json.loads((OUT/'scoring_audit.json').read_text())
summary=pd.read_csv(OUT/'group_summary.csv')
per=pd.read_csv(OUT/'per_record_metrics.csv')
complete=audit['complete_methods']
groups=['single_chain','other_multichain','tcr_pmhc_itrace','all_multichain','all']
labels=dict(zip(groups,['单链蛋白','其他多链体系','TCR–pMHC复合物','全部多链（汇总）','全部记录（汇总）']))
numbers=dict(zip(groups,[2193,626,244,870,3063]))
exposure_path=ROOT/'references/training_exposure_protocol_sources_sha256.json'
exposure=json.loads(exposure_path.read_text(encoding='utf-8')) if exposure_path.exists() else None
if exposure:
    membership_path=ROOT/'references/training_exposure_membership_frozen.csv'
    assert hashlib.sha256(membership_path.read_bytes()).hexdigest()==exposure['output_hashes'][membership_path.name]
    own_oof=exposure['own_Full_Intra_OOF']
    assert own_oof['per_variant']['full']['sample_in_corresponding_train']==0
rng=np.random.default_rng(20261023)
def mean_value(method,group,metric):
    rows=summary[(summary.method==method)&(summary.group==group)&(summary.metric==metric)]
    assert len(rows)==1 and int(rows.iloc[0].n_records)==numbers[group]
    return float(rows.iloc[0]['mean'])
paired=[]
for group in groups:
    source=per[per.method.isin(complete)]
    if group=='all_multichain':source=source[source.group!='single_chain']
    elif group!='all':source=source[source.group==group]
    for comparator in complete:
        if comparator=='Full':continue
        for metric in ['spearman','top10_ap','mae','z_mae']:
            left=source[source.method=='Full'][['sample_id','family',metric]]
            right=source[source.method==comparator][['sample_id',metric]]
            match=left.merge(right,on='sample_id',validate='one_to_one',suffixes=('_full','_external'))
            match=match[np.isfinite(match[metric+'_full'])&np.isfinite(match[metric+'_external'])].copy()
            if match.empty:continue
            match['delta']=match[metric+'_full']-match[metric+'_external']
            fam=match.groupby('family').delta.agg(['sum','count'])
            low=high=np.nan
            if len(fam)>1:
                index=rng.integers(0,len(fam),(2000,len(fam)))
                samples=fam['sum'].to_numpy()[index].sum(1)/fam['count'].to_numpy()[index].sum(1)
                low,high=np.quantile(samples,[.025,.975])
            paired.append(dict(group=group,comparator=comparator,metric=metric,n_records=len(match),n_families=len(fam),
                full_minus_comparator=float(match.delta.mean()),ci_low=low,ci_high=high,
                direction='positive favors Full' if metric in ['spearman','top10_ap'] else 'negative favors Full'))
pd.DataFrame(paired).to_csv(OUT/'paired_differences.csv',index=False)
lines=['# 本机模型横评：当前完成结果','',
    f'目前有 {len(complete)} 个方法/配置在相同的全部3,063条记录上达到完整覆盖。这个数包含同一模型的输入配置和本文消融，不等于独立外部模型家族数量。未完成模型见覆盖率表与控制器状态，不进入主图。','',
    '用户确认三组：单链2,193、其他多链626、TCR–pMHC复合物244。其他多链组也包含生物学复合物，不能称为“非复合物”。','',
    '每条记录先在相同冻结有效残基上计分，再对记录等权平均。Spearman衡量柔性排序，AP衡量MD参考最高10%柔性残基的筛选（含阈值并列），MAE只用于原作者明确输出物理幅度的模型。相对分数不给Å误差，也不在测试集拟合缩放。','']
if 'backflip-2.1 (chains)' in complete:
    comparator='backflip-2.1 (chains)'
    sg='single_chain';mg='other_multichain';tg='tcr_pmhc_itrace'
    difference=next(x for x in paired if x['group']==mg and x['comparator']==comparator and x['metric']=='spearman')
    lines += [f"单链组中，BackFlip-2.1的平均Spearman为{mean_value(comparator,sg,'spearman'):.3f}，FlexPPI为{mean_value('Full',sg,'spearman'):.3f}。其他多链组中两者排序成绩接近（FlexPPI {mean_value('Full',mg,'spearman'):.3f}，BackFlip-2.1逐链 {mean_value(comparator,mg,'spearman'):.3f}；配对差值95%区间{difference['ci_low']:.3f}至{difference['ci_high']:.3f}）；AP分别为{mean_value('Full',mg,'top10_ap'):.3f}与{mean_value(comparator,mg,'top10_ap'):.3f}，MAE分别为{mean_value('Full',mg,'mae'):.3f}与{mean_value(comparator,mg,'mae'):.3f} Å。TCR–pMHC组的描述性Spearman分别为{mean_value('Full',tg,'spearman'):.3f}与{mean_value(comparator,tg,'spearman'):.3f}。这些比较仅对当前冻结参考和已完成方法成立。",'']
if exposure:
    lines += [f"训练条件审计：Full/Intra的测试记录本身均不在该折训练集，且生物学组不跨折；但{own_oof['families_crossing_original_folds']}个序列家族跨原折，共{own_oof['per_variant']['full']['records_with_same_family_in_train']}条测试记录有同家族其他样本在训练中。244条TCR–pMHC分散于原五折，不能将它们的成绩解释为全新序列家族留出测试。家族bootstrap不会改变这一训练条件。",'',
        '外部模型也存在已文档化的训练名单ID重叠：DynaProt-M与本队列ATLAS记录有1,262条精确ID匹配；BF1/seq的Flexpert官方候选名单各有1,108条（原checkpoint CSV本地缺失）。BF2分组清单未能绑定，保持未知。这些是名单ID审计，不证明实际训练轨迹或同源暴露完全相同。详见references/training_exposure_audit_zh.md及已冻结membership；公开ID组的同记录补充比较见TRAINING_EXPOSURE_SENSITIVITY_zh.md。','']
for group in groups:
    lines += [f'## {labels[group]}：n={numbers[group]:,}','',
        '| 方法/配置 | Spearman ρ ↑ | AP ↑ | MAE (Å) ↓ | 标准化MAE ↓ |','|---|---:|---:|---:|---:|']
    rows=summary[(summary.group==group)&summary.method.isin(complete)]
    pivot=rows.pivot(index='method',columns='metric',values='mean')
    for name in pivot.sort_values('spearman',ascending=False).index:
        def fmt(metric):
            value=pivot.loc[name,metric] if metric in pivot else np.nan
            return f'{value:.3f}' if np.isfinite(value) else '不适用'
        lines.append('| '+name+' | '+' | '.join(fmt(m) for m in ['spearman','top10_ap','mae','z_mae'])+' |')
    lines.append('')
lines += ['## 解释边界','',
    'ANM/GNM采用本机SciPy实现的官方ProDy方程（15 Å、gamma=1、20个非零模态）。GNM按官方默认CKDTree的float32接触判定修正了一条恰在截断边界的接触；原始本机结果与修正覆盖文件均保留。另9条记录的第20/21模态存在简并，保留本机本征基，未按分数择值；数值与指标影响见rmsfnet_nma/gnm_ckdtree_v2/。','',
    'Full/Intra/Static/AA是已冻结的五折OOF预测复用；其他方法的来源逐项记录在scoring_audit.json。对外部发布权重不能保证与ATLAS/mdCATH等完全无训练交叠，所以这不是所有模型采用相同训练集的竞赛。参考RMSF按链独立去除整体平移/转动；整体对齐目标训练的模型在此属于标签迁移评价。','',
    '区间采用家族成组重抽样，保留记录等权平均的目标。244条ITrace记录仅属于一个依赖家族，图上只显示描述性点和分布，不给独立家族置信区间。配对差值见paired_differences.csv，未使用显著性星号。','',
    '外部柔性模型没有与本文相同定义的侧链联合状态熵Hχ输出，因此不能把SHP、方差或其他熵直接充当Hχ横评。','',
    '主图、分布图与图注见figures/。完整逐记录分数、分组结果、覆盖率与原始预测身份审计见evaluation/及各模型目录。未完成的未来图仍需最后视觉复核，当前图不自动写入论文。','']
(ROOT/'RESULTS_zh.md').write_text('\n'.join(lines),encoding='utf-8')
caption='''Figure. Protein flexibility prediction across single-chain proteins, other multichain systems and TCR–pMHC complexes.
(a) Record-mean Spearman rank correlation; (b) average precision for high-RMSF residues, defined by reference RMSF at or above each record's 90th percentile, retaining ties; (c) RMSF mean absolute error in Angstroms for methods with physical amplitude outputs. All completed entries use the same 3,063 records and frozen valid residue identities; each record contributes equally. Colors distinguish FlexPPI configurations (blue), published learned models (ochre) and graph/simple controls (teal). Shapes identify complete-structure (circles), independent/within-chain (diamonds) and sequence (triangles) configurations, with FlexPPI highlighted by a square. Assemblies and independent-chain predictions are listed separately when tested. Intervals are 95% record-weighted family-cluster bootstrap intervals (2,000 resamples); the single-family ITrace stratum has descriptive points without an independent-family interval. Relative flexibility scores are not calibrated to Angstroms. Full, Intra, Static and the residue-type mean reuse frozen five-fold out-of-fold predictions, whereas external pretrained weights follow the documented inference protocols; training exposure is not equalized and overlap is unresolved. The reference removes chain-wise rigid-body motion, which may differ from external training targets. Only full-cohort completed entries are plotted; missing and unfinished outputs are disclosed in the coverage audit. The companion empirical cumulative distributions include all records without downsampling.
'''
if exposure:
    caption=caption.replace('training exposure is not equalized and overlap is unresolved.',
        'training exposure differs and documented training-list ID matches exist. Full/Intra use record/biological-group OOF, with related-family training exposure; this is not sequence-family-held-out evaluation. The exposure audit and published-ID secondary analysis are reported separately.')
(ROOT/'figures/CAPTION_en.txt').write_text(caption,encoding='utf-8')
draft=[f'Completed-cohort results ({len(complete)} method/configuration entries; snapshot scope defined in the scoring audit).',
    'All entries were evaluated on the same 3,063 records, comprising 2,193 single-chain proteins, 626 other multichain systems and 244 TCR–pMHC complexes. Performance was first calculated within each record using identical frozen valid residues and then averaged with equal record weights. The reference RMSF removes rigid-body motion separately for each chain.']
comparators=[x for x in ['DynaProt-M','backflip-2.1 (assembly)','backflip-2.1 (chains)',
    'backflip-1.0 (assembly)','backflip-1.0 (chains)','backflip-1.0-seq (assembly)',
    'backflip-1.0-seq (chains)','FastProtFlex'] if x in complete]
for group in groups[:3]:
    values=[f"FlexPPI {mean_value('Full',group,'spearman'):.3f}"]
    values += [f'{name} {mean_value(name,group,"spearman"):.3f}' for name in comparators]
    english_group={'single_chain':'single-chain proteins','other_multichain':'other multichain systems','tcr_pmhc_itrace':'TCR–pMHC complexes'}[group]
    draft.append(f'Mean Spearman correlation for {english_group} was '+', '.join(values)+'.')
if all(x in complete for x in ['backflip-2.1 (assembly)','backflip-2.1 (chains)']):
    changes=[f'{group}: {mean_value("backflip-2.1 (assembly)",group,"spearman")-mean_value("backflip-2.1 (chains)",group,"spearman"):+.3f}' for group in groups[1:3]]
    draft.append('For BackFlip-2.1, the change in mean Spearman correlation from independent-chain to complete-structure input was '+', '.join(changes)+'. These comparisons retain interchain geometry but do not establish that the model explicitly encodes chain identities.')
    group='other_multichain';comparator='backflip-2.1 (chains)'
    difference=next(x for x in paired if x['group']==group and x['comparator']==comparator and x['metric']=='spearman')
    draft.append(f"In other multichain systems, the paired Spearman difference between FlexPPI and chain-wise BackFlip-2.1 was {difference['full_minus_comparator']:.3f} (95% family-cluster bootstrap interval {difference['ci_low']:.3f} to {difference['ci_high']:.3f}). Mean high-RMSF average precision was {mean_value('Full',group,'top10_ap'):.3f} and {mean_value(comparator,group,'top10_ap'):.3f}, respectively; mean absolute RMSF error was {mean_value('Full',group,'mae'):.3f} and {mean_value(comparator,group,'mae'):.3f} Angstroms, respectively.")
draft.append('These are evaluations of released weights under a common reference, rather than equal-training comparisons; overlap between external training data and this cohort has not been excluded. Family-cluster bootstrap intervals are reported for strata with multiple families. The TCR–pMHC stratum belongs to a single family component and is described without an independent-family confidence interval. Relative flexibility scores are not assigned Angstrom errors or substituted for the joint side-chain state entropy Hχ.')
if exposure:
    draft.append(f"The published training-ID audit identifies 1,262 ATLAS matches for DynaProt-M and 1,108 candidate Flexpert-split matches for each BackFlip-1 variant; the original BackFlip-1 checkpoint CSVs were unavailable and BackFlip-2 membership remains unresolved. Full/Intra exclude each evaluated record and its biological group from the corresponding training fold, but {own_oof['per_variant']['full']['records_with_same_family_in_train']} records have related-family members in training. Thus the primary TCR–pMHC comparisons describe record-level OOF performance under the frozen reference, rather than generalization to unseen sequence families. All published training/validation/test/unlisted ID categories are retained in the separate matched-record sensitivity analysis.")
(ROOT/'figures/RESULTS_en_draft.txt').write_text('\n\n'.join(draft)+'\n',encoding='utf-8')
print('Report and paired differences written for',len(complete),'complete entries')
