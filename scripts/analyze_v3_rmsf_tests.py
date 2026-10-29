"""Summarize every held-out v3 prediction and select explicitly labelled cases."""
from __future__ import annotations

import argparse
import csv
import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd

from export_v2_oof_predictions import ROOT, METRICS, load_cache, atomic_json
from flexppi.evaluation import _correlation, _spearman


def regression(p, y):
    p, y = np.asarray(p, np.float64), np.asarray(y, np.float64)
    if not len(y) or not np.isfinite(p).all() or not np.isfinite(y).all():
        raise ValueError("Empty/nonfinite regression input")
    error = p - y
    return dict(mae=float(np.abs(error).mean()), rmse=float(np.sqrt((error**2).mean())),
                bias=float(error.mean()), pearson=_correlation(p, y), spearman=_spearman(p, y))


def choose_cases(chains, manifest, out):
    """Select once after full testing; this is descriptive case selection, not CV."""
    eligible = chains.loc[(chains.metric == 'rmsf') & chains.n_residues.between(80, 350)
                          & (chains.valid_fraction >= .95) & chains.pearson.notna()].copy()
    chosen, used_systems, selection_audit = [], set(), {}
    # Represent the middle of the eligible distribution in BOTH correlation and error.
    for dataset in ['atlas', 'dynarepo', 'itrace']:
        group = eligible.loc[eligible.dataset == dataset].copy()
        if group.empty:
            raise ValueError('No eligible representative chain in ' + dataset)
        group['r_percentile'] = group.pearson.rank(pct=True, method='average')
        group['mae_percentile'] = group.mae.rank(pct=True, method='average')
        group['selection_distance'] = (group.r_percentile-.5).abs() + (group.mae_percentile-.5).abs()
        row = group.sort_values(['selection_distance', 'sample_id', 'chain']).iloc[0]
        chosen.append((row, 'representative'))
        used_systems.add(row.biological_system)
        selection_audit[dataset] = dict(eligible_chains=len(group),
            median_pearson=float(group.pearson.median()), median_mae=float(group.mae.median()))
    # A publication-style gallery of favourable fits. Explicitly NOT representative.
    for dataset, count in [('atlas', 3), ('dynarepo', 3), ('itrace', 2)]:
        group = eligible.loc[eligible.dataset == dataset].copy()
        group['r_percentile'] = group.pearson.rank(pct=True, method='average')
        group['mae_percentile'] = group.mae.rank(pct=True, method='average')
        group['selection_distance'] = (1-group.r_percentile) + group.mae_percentile
        number = 0
        for _, row in group.sort_values(['selection_distance', 'sample_id', 'chain']).iterrows():
            if row.biological_system in used_systems:
                continue
            chosen.append((row, 'showcase'))
            used_systems.add(row.biological_system)
            number += 1
            if number == count:
                break
        if number != count:
            raise ValueError('Insufficient gallery cases for ' + dataset)
    lookup = {r['sample_id']: r for r in manifest}
    cases = []
    for number, (row, kind) in enumerate(chosen, 1):
        item = lookup[row.sample_id]
        with np.load(item['cache_path'], allow_pickle=False) as z:
            meta = json.loads(str(z['metadata_json'].item()))
        case = {key: row[key].item() if hasattr(row[key], 'item') else row[key] for key in
                ['sample_id', 'dataset', 'fold', 'chain', 'n_residues', 'pearson', 'spearman',
                 'mae', 'rmse', 'r_percentile', 'mae_percentile', 'selection_distance',
                 'biological_system']}
        case.update(case_id=f'{kind}_{number:02d}', selection_group=kind,
                    cache_path=item['cache_path'], structure_path=meta['structure_path'],
                    trajectory_paths=item['record'].get('trajectory_paths', []))
        pdb_id = str(item['record'].get('metadata', {}).get('pdb_id') or row.sample_id)
        case['pdb_id'] = pdb_id
        case['display_name'] = (str(row.sample_id) if str(row.chain) == '_' else
                                f'{pdb_id.lower()} · chain {row.chain}')
        cases.append(case)
    atomic_json(out/'selected_cases.json', cases)
    pd.DataFrame(cases).drop(columns=['trajectory_paths']).to_csv(out/'selected_cases.csv', index=False)
    atomic_json(out/'selection_protocol.json', dict(
        eligibility='All RMSF chains with 80–350 valid residues, >=95% label/geometry coverage and finite Pearson.',
        representative='One chain per source minimizing |Pearson percentile-0.5| + |MAE percentile-0.5|; ties by sample ID, chain.',
        showcase='3 ATLAS, 3 DynaRepo, 2 ITrace chains minimizing (1-Pearson percentile)+MAE percentile; unique biological systems, excluding representatives.',
        caveat='Cases were selected after testing for illustration. The favourable gallery is not an estimate of typical performance. Full test statistics include all valid residues.',
        structure_colour='Same input structure colored by MD target or prediction, common linear 0–5 Å scale, values above 5 Å saturate; curves remain unclipped.',
        sources=selection_audit))
    return cases


def run(args):
    folder, out = Path(args.oof_dir), Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    summary = json.loads((folder/'run_summary.json').read_text())
    if summary['status'] != 'complete' or summary['partial'] or summary['processed_samples'] != 3068:
        raise ValueError('Requires complete 3068-system OOF inference')
    manifest = [json.loads(s) for s in (folder/'prediction_manifest.jsonl').read_text().splitlines()]
    states, chains, arrays = [], [], {metric: [] for metric in METRICS}
    all_residues, missing = 0, dict.fromkeys(METRICS, 0)
    fields = ['sample_id', 'dataset', 'fold', 'chain', 'author_residue_id', 'residue_index', 'amino_acid',
              'valid_rmsf', 'md_rmsf_angstrom', 'pred_rmsf_angstrom', 'md_rmsf_variance', 'model_rmsf_variance',
              'valid_chi_entropy', 'md_chi_entropy_nats', 'pred_chi_entropy_nats']
    with gzip.open(out/'all_residue_predictions.csv.gz', 'wt', newline='', encoding='utf-8') as handle:
        writer = csv.writer(handle)
        writer.writerow(fields)
        for counter, item in enumerate(manifest, 1):
            if item['status'] != 'complete':
                raise ValueError('Incomplete sample: ' + item['sample_id'])
            a = load_cache(Path(item['cache_path']), item['identity'])
            raw = item['record']
            base = dict(sample_id=item['sample_id'], dataset=raw['dataset'], fold=item['fold'],
                        biological_system=raw['metadata']['biological_system'])
            n = len(a['chain'])
            all_residues += n
            for j in range(n):
                writer.writerow([base['sample_id'], base['dataset'], base['fold'], a['chain'][j],
                    a['author_residue_id'][j], int(a['residue_index'][j]), a['amino_acid'][j],
                    bool(a['valid'][j, 0]), float(a['target'][j, 0]), float(a['mean'][j, 0]),
                    float(a['label_variance'][j, 0]), float(a['model_variance'][j, 0]),
                    bool(a['valid'][j, 1]), float(a['target'][j, 1]), float(a['mean'][j, 1])])
            for i, metric in enumerate(METRICS):
                valid = a['valid'][:, i]
                missing[metric] += int((~valid).sum())
                if not valid.any():
                    continue
                p, y = a['mean'][valid, i], a['target'][valid, i]
                states.append(dict(**base, metric=metric, n_residues=int(valid.sum()), **regression(p, y)))
                arrays[metric].append(dict(fold=item['fold'], dataset=raw['dataset'], p=p, y=y,
                                           aa=a['amino_acid'][valid], sid=item['sample_id']))
                for chain in dict.fromkeys(a['chain']):
                    selected = valid & (a['chain'] == chain)
                    if not selected.any():
                        continue
                    chains.append(dict(**base, chain=str(chain), metric=metric,
                        n_residues=int(selected.sum()), chain_length=int((a['chain'] == chain).sum()),
                        valid_fraction=float(selected.sum()/(a['chain'] == chain).sum()),
                        **regression(a['mean'][selected, i], a['target'][selected, i])))
            if counter % 500 == 0:
                print(f'Analyzed {counter}/{len(manifest)} systems', flush=True)
    systems, chain_data = pd.DataFrame(states), pd.DataFrame(chains)
    systems.to_csv(out/'per_system_metrics.csv', index=False)
    chain_data.to_csv(out/'per_chain_metrics.csv', index=False)
    pooled_rows, macro_rows, checks = [], [], []
    statistics = ['mae', 'rmse', 'bias', 'pearson', 'spearman']
    for metric in METRICS:
        s = systems.loc[systems.metric == metric]
        for fold in list(range(5))+['all']:
            for dataset in ['all', 'atlas', 'dynarepo', 'itrace']:
                parts = [v for v in arrays[metric] if (fold == 'all' or v['fold'] == fold)
                         and (dataset == 'all' or v['dataset'] == dataset)]
                group = s.loc[((s.fold == fold) if fold != 'all' else np.ones(len(s), bool)) &
                              ((s.dataset == dataset) if dataset != 'all' else np.ones(len(s), bool))]
                if not parts:
                    continue
                p, y = np.concatenate([v['p'] for v in parts]), np.concatenate([v['y'] for v in parts])
                result = regression(p, y)
                pooled_rows.append(dict(metric=metric, fold=fold, dataset=dataset,
                                        n_residues=len(y), n_systems=len(group), **result))
                macro = group[statistics].mean().to_dict()
                macro_rows.append(dict(metric=metric, fold=fold, dataset=dataset,
                                       n_residues=len(y), n_systems=len(group), **macro))
                if fold != 'all':
                    ref = json.loads((ROOT/f'artifacts/performance_first_v3_difficulty_balanced/evaluation/fold_{fold}_test.json').read_text())['metrics'][metric]
                    ref = ref['overall'] if dataset == 'all' else ref['by_dataset'][dataset]
                    for scope, values, target in [('pooled', result, ref), ('macro_by_system', macro, ref['macro_by_system'])]:
                        for key in ['mae', 'rmse', 'pearson', 'spearman']:
                            delta = float(values[key]-target[key])
                            checks.append(dict(metric=metric, fold=fold, dataset=dataset,
                                scope=scope, statistic=key, local_value=values[key], server_value=target[key],
                                absolute_difference=abs(delta), passed=abs(delta)<2e-4))
    pooled, macros = pd.DataFrame(pooled_rows), pd.DataFrame(macro_rows)
    pooled.to_csv(out/'pooled_metrics.csv', index=False)
    macros.to_csv(out/'macro_system_metrics.csv', index=False)
    # Remove amino-acid-type averages using labels from only the four training folds.
    corrected_rows, prior_rows = [], []
    for fold in range(5):
        training = [v for v in arrays['chi_entropy'] if v['fold'] != fold]
        aa = np.concatenate([v['aa'] for v in training])
        y = np.concatenate([v['y'].astype(np.float64) for v in training])
        priors = {str(k): float(y[aa == k].mean()) for k in np.unique(aa)}
        reference = json.loads((ROOT/f'artifacts/performance_first_v3_difficulty_balanced/evaluation/fold_{fold}_test.json').read_text())['metrics']['chi_entropy']['amino_acid_corrected']
        for k, value in priors.items():
            delta = abs(value-reference['training_priors'][k])
            if delta > 1e-6:
                raise ValueError('Training-only Hchi prior disagrees with server')
            prior_rows.append(dict(fold=fold, amino_acid=k, mean_training_hchi=value,
                                   n_training_residues=int((aa == k).sum())))
        local_p, local_y, local_stats = [], [], []
        for part in arrays['chi_entropy']:
            if part['fold'] != fold:
                continue
            prior = np.asarray([priors[str(k)] for k in part['aa']])
            p, y = part['p'].astype(float)-prior, part['y'].astype(float)-prior
            local_p.append(p)
            local_y.append(y)
            local_stats.append(regression(p, y))
        corrected = regression(np.concatenate(local_p), np.concatenate(local_y))
        macro = pd.DataFrame(local_stats)[statistics].mean().to_dict()
        for scope, result, target in [('pooled', corrected, reference['residual_overall']),
                                      ('macro_by_system', macro, reference['residual_overall']['macro_by_system'])]:
            corrected_rows.append(dict(fold=fold, scope=scope, **result))
            for key in ['mae', 'rmse', 'pearson', 'spearman']:
                difference = abs(result[key]-target[key])
                checks.append(dict(metric='chi_entropy_type_corrected', fold=fold, dataset='all',
                    scope=scope, statistic=key, local_value=result[key], server_value=target[key],
                    absolute_difference=difference, passed=difference < 2e-4))
    corrected_table = pd.DataFrame(corrected_rows)
    corrected_table.to_csv(out/'hchi_type_corrected_metrics.csv', index=False)
    pd.DataFrame(prior_rows).to_csv(out/'hchi_training_only_type_priors.csv', index=False)
    pd.DataFrame(checks).to_csv(out/'server_reproduction_checks.csv', index=False)
    if not all(c['passed'] for c in checks):
        raise ValueError('Local predictions disagree with server metrics; inspect reproduction CSV')
    five_fold = []
    for scope, table in [('pooled', pooled), ('macro_by_system', macros)]:
        table = table.loc[table.fold != 'all']
        for (metric, dataset), group in table.groupby(['metric', 'dataset']):
            if len(group) != 5:
                raise ValueError('Incomplete five-fold metric aggregation')
            for stat in statistics:
                five_fold.append(dict(scope=scope, metric=metric, dataset=dataset, statistic=stat,
                                      mean=float(group[stat].mean()), sd=float(group[stat].std(ddof=0)),
                                      sample_sd_ddof1=float(group[stat].std(ddof=1))))
    five = pd.DataFrame(five_fold)
    five.to_csv(out/'five_fold_mean_sd.csv', index=False)
    counts = {metric: int(sum(len(x['y']) for x in arrays[metric])) for metric in METRICS}
    if counts != {'rmsf':959347, 'chi_entropy':431023}:
        raise ValueError('Valid label counts mismatch')
    cases = choose_cases(chain_data, manifest, out)
    audit = dict(status='complete', systems=len(manifest), all_structural_residues=all_residues,
                 valid_residues=counts, masked_residues=missing, silently_dropped_samples=0,
                 server_checks=len(checks), all_server_checks_passed=True,
                 maximum_server_difference=max(c['absolute_difference'] for c in checks),
                 selected_cases=len(cases), oof_run_summary=summary)
    atomic_json(out/'test_summary.json', audit)
    lines = ['# FlexPPI v3：五折测试与残基级展示', '',
        '本次对 3,068 个样本逐一重新推理，每个样本仅使用其所属测试折的 epoch_040.pt。使用完整结构上下文；训练权重、数据划分与标签均未修改。', '',
        f'有效 RMSF 标签 {counts["rmsf"]:,} 个；有效 Hχ 标签 {counts["chi_entropy"]:,} 个。全部 {len(checks)} 项本机/服务器指标核对通过，最大绝对差为 {audit["maximum_server_difference"]:.3g}。', '',
        '|指标|统计方式|Pearson r（五折均值 ± SD）|MAE（五折均值 ± SD）|', '|---|---|---|---|']
    for metric, title, unit in [('rmsf', 'RMSF', 'Å'), ('chi_entropy', 'Hχ', 'nat')]:
        for scope, label in [('macro_by_system', '每个蛋白分别计算后平均'), ('pooled', '每折合并所有残基')]:
            rows = five.loc[(five.metric == metric)&(five.scope == scope)&(five.dataset == 'all')].set_index('statistic')
            lines.append(f'|{title}|{label}|{rows.loc["pearson","mean"]:.3f} ± {rows.loc["pearson","sd"]:.3f}|{rows.loc["mae","mean"]:.3f} ± {rows.loc["mae","sd"]:.3f} {unit}|')
    lines += ['', '表中 SD 沿用已冻结服务器汇总方式（ddof=0，即描述这五折的总体标准差），不是置信区间；CSV 另列 ddof=1 的样本标准差供查阅。每个蛋白样本等权按 sample_id 计算，不代表同源簇等权。', '',
        'RMSF 是每条链独立对齐后计算的 Cα 内部波动幅度，单位 Å，去除了链之间整体平移和旋转的影响；与参考论文的局部对齐窗口定义未必相同。r 表示峰谷变化趋势是否一致，MAE 表示预测值平均偏离 MD 标签多少；较高 r 不保证幅度没有系统偏差。', '',
        '主图每个来源选一个接近其可展示链群体中位表现的案例。补充图选取高相关、低误差案例以展示拟合细节，不能据此声称所有蛋白都达到该水平。选择规则、百分位与所有案例完整数据见 selection_protocol.json、selected_cases.csv；全测试结果使用所有有效标签。', '',
        '结构图将相同输入结构分别按 MD 标签和 FlexPPI 预测着色。它们不是模型生成的新构象或模型预测的轨迹。若显示阴影，阴影来自标签文件所存方差的平方根，表示 MD 标签离散程度，不是预测置信区间。颜色范围固定为 0–5 Å，超过上限处饱和，曲线保留完整取值。', '',
        'Hχ 为侧链构象熵。未标注 Hχ 的位置（包括本轮 ATLAS 标签）不参与 Hχ 评价；缺失标签没有被补成零。两项任务的有效位置掩码分别记录于完整 CSV。不同氨基酸天然具有不同的侧链自由度，会抬高未校正的总体相关性，因此另用各训练折按氨基酸类别的均值作校正，详细结果见 hchi_type_corrected_metrics.csv。', '',
        '本结果属于已冻结 v3 五折划分的留出评价，不能称为完全未见同源家族的外部测试；v3 在查看 v2 结果后建立，须作为独立版本报告。', '',
        '复现流程：先运行 scripts/export_v3_oof_predictions.py，再运行 scripts/analyze_v3_rmsf_tests.py、scripts/render_v3_rmsf_structures.py 和 scripts/plot_v3_rmsf_reference_figures.py。所有预测缓存保留检查点、结构、标签、清单与实现的 SHA-256。']
    corrected_macro = corrected_table.loc[corrected_table.scope == 'macro_by_system', 'pearson']
    rmsf_bias = five.loc[(five.metric == 'rmsf') & (five.scope == 'pooled') &
                        (five.dataset == 'all') & (five.statistic == 'bias'), 'mean'].iloc[0]
    lines += ['', f'补充结果：RMSF 的折均值预测偏差（预测减 MD）为 {rmsf_bias:.3f} Å；负值说明模型平均低估波动幅度。Hχ 扣除氨基酸类别均值后，样本等权 Pearson 为 {corrected_macro.mean():.3f} ± {corrected_macro.std(ddof=0):.3f}，该结果更直接反映同类氨基酸在不同环境中的构象差异。']
    (out/'TEST_REPORT.md').write_text('\n'.join(lines)+'\n', encoding='utf-8')
    print(json.dumps({k:v for k,v in audit.items() if k!='oof_run_summary'}, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--oof-dir', default=str(ROOT/'artifacts/evaluation_v3/oof'))
    parser.add_argument('--output-dir', default=str(ROOT/'artifacts/evaluation_v3/analysis'))
    run(parser.parse_args())
