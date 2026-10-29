"""Prepare small frozen-source plotting tables and Chinese interpretation."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import sys
import pandas as pd

BASE = Path(__file__).resolve().parent
REMOTE = '/root/flexppi_completion_20261010/energy_landscape_v21_20261012/oof_analysis'

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def main():
    sys.stdout.reconfigure(encoding='utf-8')
    audit = json.loads((BASE/'run_audit.json').read_text())
    local_checks = {name:sha(BASE/name)==digest for name,digest in audit['files'].items() if name!='paired_common_residues.csv.gz'}
    assert all(local_checks.values())
    validation = json.loads((BASE/'validation/independent_validation.json').read_text())
    assert validation['status'].startswith('passed') and validation['run_audit_sha256']==sha(BASE/'run_audit.json')
    for name,digest in validation['validation_output_sha256'].items():
        assert sha(BASE/'validation'/name)==digest
    coupling = pd.read_csv(BASE/'correlation_record_macro.csv',float_precision='round_trip')
    classes = pd.read_csv(BASE/'reference_class_record_macro.csv',float_precision='round_trip')
    per = pd.read_csv(BASE/'per_record_reference_classes.csv',float_precision='round_trip')
    counts = pd.read_csv(BASE/'raw_source_counts.csv').set_index('dataset')
    confusion = pd.read_csv(BASE/'class_confusion_residue_counts.csv')
    summary_rows = []
    class_labels = {'low_RMSF_low_Hchi':'低RMSF／低Hχ','low_RMSF_high_Hchi':'低RMSF／高Hχ','high_RMSF_low_Hchi':'高RMSF／低Hχ','high_RMSF_high_Hchi':'高RMSF／高Hχ','unclassified':'未分类'}
    for r in classes.to_dict(orient='records'):
        r['class_label_zh'] = class_labels[r['reference_class']]
        common = 824988 if r['dataset']=='all' else int(counts.loc[r['dataset'],'common_valid'])
        selected = per[per.reference_class.eq(r['reference_class'])]
        if r['dataset']!='all':
            selected = selected[selected.dataset.eq(r['dataset'])]
        r['cohort_records'] = len(selected)
        r['records_with_reference_class'] = int(selected['count'].gt(0).sum())
        r['fraction_micro'] = r['common_residues']/common
        r['rank_class_recall_micro'] = float(selected.intersection_count.sum()/selected['count'].sum()) if r['reference_class']!='unclassified' and selected['count'].sum() else None
        r['rank_class_precision_micro'] = float(selected.intersection_count.sum()/selected.predicted_class_count.sum()) if r['reference_class']!='unclassified' and selected.predicted_class_count.sum() else None
        summary_rows.append(r)
    compact = pd.DataFrame(summary_rows)
    output = BASE/'plot_data'
    output.mkdir(exist_ok=True)
    columns = ['dataset','mode','unique_records','common_residues','reference_spearman_record_macro','prediction_spearman_record_macro','reference_pearson_record_macro','prediction_pearson_record_macro','spearman_absolute_discrepancy_record_macro','pearson_absolute_discrepancy_record_macro','reference_spearman_n_records','prediction_spearman_n_records']
    coupling[columns].to_csv(output/'coupling_summary_compact.csv',index=False)
    columns = ['dataset','reference_class','class_label_zh','cohort_records','records_with_reference_class','common_residues','fraction_micro','reference_fraction_record_macro','rank_class_recall_record_macro','rank_class_recall_n_records','rank_class_recall_micro','rank_class_precision_record_macro','rank_class_precision_micro']+[label+'_'+measure for label in ['rmsf','hchi'] for measure in ['reference_mean_record_macro','prediction_mean_record_macro','mae_record_macro','mae_n_records','spearman_record_macro','spearman_n_records','pearson_record_macro','pearson_n_records']]
    compact[columns].to_csv(output/'quartile_class_summary_compact.csv',index=False)
    coupling_all = coupling[coupling.dataset.eq('all')].set_index('mode')
    class_all = compact[compact.dataset.eq('all')].set_index('reference_class')
    discordant = ['low_RMSF_high_Hchi','high_RMSF_low_Hchi']
    discordant_count = int(class_all.loc[discordant,'common_residues'].sum())
    discordant_fraction_macro = float(class_all.loc[discordant,'reference_fraction_record_macro'].sum())
    discordant_fraction_micro = discordant_count/824988
    rows = []
    for dataset in ['all','atlas','dynarepo','itrace']:
        selected = per[per.reference_class.isin(discordant)]
        if dataset!='all':
            selected = selected[selected.dataset.eq(dataset)]
        summed = selected.groupby('sample_id')[['count','intersection_count','reference_fraction']].sum()
        valid = summed['count'].gt(0)
        rows.append(dict(dataset=dataset,cohort_records=len(summed),records_with_discordant_positions=int(valid.sum()),discordant_positions=int(summed['count'].sum()),reference_discordant_fraction_record_macro=float(summed.reference_fraction.mean()),exact_same_class_recall_record_macro=float((summed.loc[valid,'intersection_count']/summed.loc[valid,'count']).mean()),exact_same_class_recall_micro=float(summed.intersection_count.sum()/summed['count'].sum())))
    pd.DataFrame(rows).to_csv(output/'discordant_combinations_summary.csv',index=False)
    raw = coupling_all.loc['cross_type_raw']
    centered = coupling_all.loc['within_record_AA_centered']
    method = '在已冻结的3,063条完整模型折外预测中，按原RMSF与Hχ有效掩码取交集，共824,988个位点。每条记录的预测来自其对应留出折，采用种子7、第40轮权重。分别计算真实标签和预测输出之间的Spearman与Pearson相关；类型控制版本在每条记录内，对同一种氨基酸分别减去两个指标各自的均值，再计算相关。相应类型内指标仅在至少四个位点且相关变量均有变化时计算，所有不能计算的原因另行保留。四类局部组合按每条记录内同类型真实标签的第25和第75百分位确定；中间位置、同类型不足四个位点及两个百分位重合的情况保留为未分类。预测类别由预测自身的同类型百分位确定。先逐记录计算相关、误差、类别比例和召回率，再对记录等权汇总；原始位置数与位置等权结果另外报告。'
    results = f'两类动态指标呈相关但不同的残基排序。跨类型原始数值的真实和预测Spearman相关分别为{raw.reference_spearman_record_macro:.3f}与{raw.prediction_spearman_record_macro:.3f}；减去记录内同类型均值后分别为{centered.reference_spearman_record_macro:.3f}与{centered.prediction_spearman_record_macro:.3f}。真实标签中有{discordant_count:,}个位点落入低RMSF／高Hχ或高RMSF／低Hχ组合，占共同有效位置的{100*discordant_fraction_micro:.3f}%，逐记录等权比例为{100*discordant_fraction_macro:.3f}%。这些位置显示主链活动幅度与侧链状态多样性可以提供不同的局部读数。预测总体也保留两类组合，但同类型去均值后的相关高于参考值，表明两个预测读数的共同变化更强。'
    reliability = f'以各输出自身的同类型四分位排序检验组合恢复时，低RMSF／高Hχ与高RMSF／低Hχ参考组的记录等权召回率分别为{100*class_all.loc[discordant[0],"rank_class_recall_record_macro"]:.1f}%和{100*class_all.loc[discordant[1],"rank_class_recall_record_macro"]:.1f}%；同向的低／低和高／高组分别为{100*class_all.loc["low_RMSF_low_Hchi","rank_class_recall_record_macro"]:.1f}%和{100*class_all.loc["high_RMSF_high_Hchi","rank_class_recall_record_macro"]:.1f}%。低RMSF／高Hχ组的两项MAE为{class_all.loc[discordant[0],"rmsf_mae_record_macro"]:.3f} Å和{class_all.loc[discordant[0],"hchi_mae_record_macro"]:.3f} nats；高RMSF／低Hχ组为{class_all.loc[discordant[1],"rmsf_mae_record_macro"]:.3f} Å和{class_all.loc[discordant[1],"hchi_mae_record_macro"]:.3f} nats。该组合指标比逐目标排序更严格：参考属于某类而预测落在自身任一指标的中间范围时，也计为未恢复。因此，新结果支持联合读数的互补用途，同时表明不同步尾部组合的筛选可靠性仍低于同向组合。'
    limitations = '这里的Hχ为固定角度分箱和现有轨迹采样得到的联合状态熵。跨蛋白的RMSF–Hχ位置图及其出现频率描述的是不同局部动态读数的分布；本分析没有用两指标拟合或构造自由能。四分位组的“低”和“高”是同一记录同一种残基内的相对位置，因此表中的组均值与误差不是跨蛋白通用的物理阈值。来源及生物体系中的记录依赖仍按主文原体系留出定义处理，本分析为描述性总体计算，没有把3,063条记录当作独立生物学重复进行显著性推断。'
    caption = '建议新增图使用：（a）真实与预测的逐记录RMSF–Hχ Spearman相关分布，分别显示原值和记录内同类型去均值版本；（b）四类参考组合及未分类的记录等权比例，另在图注给出原始计数；（c,d）同一参考组合中的真实与预测RMSF均值（Å）及Hχ均值（nats），分开坐标；如展示类别恢复，使用完整5×5混淆矩阵保留未分类列。点或条形为记录等权汇总，误差范围可直接使用逐记录分布分位数，不以五折SD或独立重复显著性表述。不要把经验频率画为物理自由能或给不存在的动力学速率。'
    table = '\n'.join('| '+class_labels[name]+' | '+str(int(class_all.loc[name,'common_residues']))+' | '+f'{100*class_all.loc[name,"reference_fraction_record_macro"]:.3f}%'+' | '+f'{class_all.loc[name,"rmsf_mae_record_macro"]:.3f}'+' | '+f'{class_all.loc[name,"hchi_mae_record_macro"]:.3f}'+' | '+str(int(class_all.loc[name,'records_with_reference_class']))+' |' for name in ['low_RMSF_low_Hchi','low_RMSF_high_Hchi','high_RMSF_low_Hchi','high_RMSF_high_Hchi','unclassified'])
    markdown = f'{method}\n\n{results}\n\n{reliability}\n\n{limitations}\n\n| 真实同类型相对组合 | 位置数 | 记录等权比例 | RMSF MAE（Å） | Hχ MAE（nats） | 有该组的记录数 |\n|---|---:|---:|---:|---:|---:|\n{table}\n\n{caption}\n\n来源：`run_audit.json`、`input_freeze.json`、`correlation_record_macro.csv`、`reference_class_record_macro.csv`、`class_confusion_residue_counts.csv`；原始逐残基配对表仅保留在远程 `{REMOTE}/paired_common_residues.csv.gz`（SHA256 {audit["files"]["paired_common_residues.csv.gz"]}）。整体RMSF原有效位置为961,321，本配对分析移除Hχ不适用或无效的136,333个位点，不能将配对RMSF误差直接代替主文的全RMSF分数。\n'
    (BASE/'manuscript_summary_zh.md').write_text(markdown,encoding='utf-8')
    numbers = dict(protocol=audit['protocol'],records=3063,common_valid=824988,run_audit_sha256=sha(BASE/'run_audit.json'),input_freeze_sha256=sha(BASE/'input_freeze.json'),analysis_code_sha256=audit['analysis_code_sha256'],current_v20_document_sha256='4b4597b42fdc8a1300ab8f4155e98a4960872a2368fc887397d133af7d264d79',coupling_summary=coupling.to_dict(orient='records'),reference_class_counts=audit['reference_class_counts'],prediction_class_counts=audit['prediction_class_counts'],discordant_reference_count=discordant_count,discordant_reference_fraction_micro=discordant_fraction_micro,discordant_reference_fraction_record_macro=discordant_fraction_macro,method_zh=method,results_zh=results,reliability_zh=reliability,scope_zh=limitations,local_certified_file_hashes_passed=len(local_checks),independent_validation_checks=len(validation['checks']))
    (BASE/'manuscript_ready_numbers.json').write_text(json.dumps(numbers,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    provenance = dict(status='local_sync_and_derivative_summary_passed',utc=datetime.now(timezone.utc).isoformat(),run_audit_sha256=sha(BASE/'run_audit.json'),independent_validation_sha256=sha(BASE/'validation/independent_validation.json'),rare_nonAG_mask_inventory_sha256=sha(BASE/'validation/rare_nonAG_mask_exclusions.json'),remote=REMOTE,local_certified_file_checks=local_checks,large_paired_residue_gzip_downloaded=False,large_paired_residue_gzip_remote_sha256=audit['files']['paired_common_residues.csv.gz'],command=['env','OMP_NUM_THREADS=2','OPENBLAS_NUM_THREADS=2','MKL_NUM_THREADS=2','/root/flexppi_completion_20261010/cpu_env/bin/python',REMOTE+'/analyze_oof_complementarity.py','--training-root','/root/flexppi_revision_20261002/corrected_training_v1','--cohort','/root/flexppi_completion_20261010/baselines/cohort','--output',REMOTE],summary_code_sha256=sha(__file__),derivative_files={p.relative_to(BASE).as_posix():sha(p) for p in [BASE/'manuscript_summary_zh.md',BASE/'manuscript_ready_numbers.json']+sorted(output.glob('*.csv'))})
    (BASE/'local_sync_and_summary_audit.json').write_text(json.dumps(provenance,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status=provenance['status'],correlation_record_macro=coupling[coupling.dataset.eq('all')][['mode','reference_spearman_record_macro','prediction_spearman_record_macro']].to_dict(orient='records'),discordant_reference_count=discordant_count,discordant_record_macro_percent=100*discordant_fraction_macro,reference_class_summary=compact[compact.dataset.eq('all')][['reference_class','records_with_reference_class','common_residues','rank_class_recall_record_macro','rmsf_mae_record_macro','hchi_mae_record_macro']].to_dict(orient='records'),markdown=str(BASE/'manuscript_summary_zh.md')),ensure_ascii=False,indent=2))

if __name__=='__main__':
    main()
