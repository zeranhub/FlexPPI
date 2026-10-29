# -*- coding: utf-8 -*-
"""Export six previously computed rows for plotting; no new statistical analysis."""
from pathlib import Path
import json
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent


def main():
    ledger = json.loads((HERE/'main_numbers.json').read_text(encoding='utf-8'))
    base = pd.read_csv(HERE/'matched_summary_family_cluster_bootstrap.csv')
    absolute = pd.read_csv(HERE/'sasa_absolute_sensitivity/summary_family_cluster_bootstrap.csv')
    base_record = pd.read_csv(HERE/'matched_per_record.csv')
    absolute_record = pd.read_csv(HERE/'sasa_absolute_sensitivity/matched_per_record.csv')
    cases = [
        ('primary_rmsf_to_hchi', 'primary_fixed', 'rmsf_to_hchi', 'physical_0p10', 'all_within_record', 'RMSF差≤0.10 Å，考察Hχ'),
        ('primary_hchi_to_rmsf', 'primary_fixed', 'hchi_to_rmsf', 'physical_0p10', 'all_within_record', 'Hχ差≤0.10 nats，考察RMSF'),
        ('sasa10_record_rmsf_to_hchi', 'post_exploration_robustness', 'rmsf_to_hchi', 'sasa_abs_10A2_and_condition_abs_0p10', 'all_within_record', 'SASA差≤10 Å²、RMSF差≤0.10 Å，考察Hχ（探索后检查）'),
        ('sasa10_record_hchi_to_rmsf', 'post_exploration_robustness', 'hchi_to_rmsf', 'sasa_abs_10A2_and_condition_abs_0p10', 'all_within_record', 'SASA差≤10 Å²、Hχ差≤0.10 nats，考察RMSF（探索后检查）'),
        ('sasa10_chain_rmsf_to_hchi', 'post_exploration_robustness', 'rmsf_to_hchi', 'sasa_abs_10A2_and_condition_abs_0p10', 'same_chain_only', None),
        ('sasa10_chain_hchi_to_rmsf', 'post_exploration_robustness', 'hchi_to_rmsf', 'sasa_abs_10A2_and_condition_abs_0p10', 'same_chain_only', None),
    ]
    rows = []
    for row_id, tier, direction, rule, scope, ledger_key in cases:
        isbase = tier == 'primary_fixed'
        frame, records = (base, base_record) if isbase else (absolute, absolute_record)
        filt = ((frame.label_mode == 'uniform_rep1')&(frame.target_subset == 'all_close')&
                (frame.dataset_scope == 'all')&(frame.direction == direction)&(frame.match_rule == rule))
        rec_filt = ((records.label_mode == 'uniform_rep1')&(records.target_subset == 'all_close')&
                    (records.direction == direction)&(records.match_rule == rule))
        if not isbase:
            filt &= frame.pair_scope == scope
            rec_filt &= records.pair_scope == scope
        block, individual = frame[filt].set_index('metric'), records[rec_filt]
        assert len(individual) == 83 and block.index.is_unique
        if ledger_key:
            for metric, expected in ledger[ledger_key].items():
                np.testing.assert_allclose(block.loc[metric, 'record_macro_mean'], expected['record_macro_mean'], atol=1e-11, rtol=0)
        f = block.loc['full_direction_accuracy']
        condition, target = direction.split('_to_')
        row = dict(row_id=row_id, evidence_tier=tier, label_mode='uniform_rep1',
            target_subset='all_close', dataset_scope='all', direction=direction,
            match_rule=rule, pair_scope=scope, condition_metric=condition,
            condition_abs_caliper=.10, condition_unit='angstrom' if condition == 'rmsf' else 'nats',
            target_metric=target, target_unit='angstrom' if target == 'rmsf' else 'nats',
            sasa_abs_caliper_A2=np.nan if isbase else 10,
            cohort_records=83, defined_records=int(f.defined_records),
            undefined_records=83-int(f.defined_records),
            matched_pairs=int(f.matched_pairs), matched_positions=int(f.matched_positions),
            available_family_clusters=int(f.available_family_clusters),
            defined_family_clusters=int(f.defined_family_clusters),
            true_target_abs_delta_macro=float(block.loc['target_true_abs_delta_mean', 'record_macro_mean']),
            full_above_chance_records=int((individual.full_direction_accuracy > .5).sum()),
            full_below_chance_records=int((individual.full_direction_accuracy < .5).sum()),
            full_equal_chance_records=int((individual.full_direction_accuracy == .5).sum()),
            full_at_or_below_chance_records=int((individual.full_direction_accuracy <= .5).sum()),
            no_pair_record_ids=';'.join(individual.loc[individual.matched_pairs == 0, 'sample_id']))
        for name in ('full', 'intra', 'surrogate'):
            r = block.loc[name+'_direction_accuracy']
            assert int(r.defined_records) == row['defined_records']
            row[name+'_accuracy_pct'] = 100*r.record_macro_mean
            row[name+'_accuracy_ci95_low_pct'] = 100*r.family_cluster_ci95_low
            row[name+'_accuracy_ci95_high_pct'] = 100*r.family_cluster_ci95_high
        for name in ('full_minus_surrogate', 'full_minus_intra'):
            r = block.loc[name+'_accuracy']
            row[name+'_pp'] = 100*r.record_macro_mean
            row[name+'_ci95_low_pp'] = 100*r.family_cluster_ci95_low
            row[name+'_ci95_high_pp'] = 100*r.family_cluster_ci95_high
        rows.append(row)
    result = pd.DataFrame(rows)
    expected = result[result.row_id == 'sasa10_record_rmsf_to_hchi'].iloc[0]
    assert expected.defined_records == 82 and expected.undefined_records == 1
    assert expected.no_pair_record_ids == '2gib_A'
    assert expected.full_at_or_below_chance_records == 11
    assert expected.full_below_chance_records == 9 and expected.full_equal_chance_records == 2
    assert result.label_mode.eq('uniform_rep1').all() and result.target_subset.eq('all_close').all()
    result.to_csv(HERE/'figure3_matching_contract.csv', index=False, float_format='%.12g', encoding='utf-8')
    prose = r'''# 图3匹配数据与措辞合同

本CSV仅提取已经完成的结果，没有重新筛选位置或重新计算置信区间。六行全部固定为`label_mode=uniform_rep1`、`target_subset=all_close`、`dataset_scope=all`。真实指标来自每条记录第一条轨迹的0—100 ns共同窗口，200 ps间隔、501帧；预测为不变的原OOF。禁止把`separated_target_p50`高对比分离子集的82%或93%一致率放入这些主结果。原标签与可用重复均值不是本表的真值模式。

前两行是原先固定的物理匹配：在同记录、同氨基酸类型、组内至少4个位点中穷尽全部无序残基对。`rmsf_to_hchi`要求真实RMSF绝对差≤0.10 Å，检验真实Hχ的排序；`hchi_to_rmsf`要求真实Hχ绝对差≤0.10 nats，检验真实RMSF的排序。这两行均有83条合格记录，不额外限制SASA。

中间两行保留同记录同类型规则，并要求起始完整复合物中的侧链SASA绝对差≤10 Å²，同时应用相同的0.10物理条件。它们是原SASA百分位匹配完成后的探索后稳健性检验，不得称为原始预注册主分析，也不应覆盖掉原百分位结果。最后两行再要求每一对的链标识相同；≥4位点的资格与目标百分位仍按记录/类型定义，并非重新要求每条链有≥4位点。

`all_within_record`中的记录可为多链装配体，残基对可能跨链。应写“同一复合物记录内、同类型残基”，不能写“同一蛋白内”或“同一条链内”。只有`same_chain_only`可以写“进一步限制在同一条蛋白链的残基对”。同链限制是稳健性检查，不改变最初记录级匹配的定义。

方向一致率先对同记录同类型内的合格对等权，再对类型等权，最后对记录等权。真实目标差为零的对仅从方向分母剔除；预测平局记0.5。`*_accuracy_pct`及其CI已经乘100，是百分数；`full_minus_surrogate_pp`及其CI是百分点，不能再乘100。`surrogate`不是另一个模型：目标为Hχ时用Full预测RMSF代排，目标为RMSF时用Full预测Hχ代排。专属增量是同一批记录的配对“Full专属输出－Full另一输出”，不是两个独立均值CI相减。`full_minus_intra_*`同样是配对差；本表这些区间均跨零，不能宣称伙伴信息产生了显著额外优势。

置信区间是10,000次可得家族连通簇整簇bootstrap的百分位区间。83条记录只有42个簇，最大簇包含38条；记录数、残基对数都不应标作独立蛋白数。`matched_pairs`和`matched_positions`是该规则下的描述性计数，不是统计检验样本数。`defined_family_clusters`是有定义记录对应的簇数；`available_family_clusters`还可能包含最终没有匹配对的记录所在簇。家族簇重抽不能将原OOF变成家族留出预测。

样本量必须逐行标注：SASA≤10 Å²且真实RMSF差≤0.10 Å的Hχ组只有82条可评分记录，`2gib_A`没有合格对，保留为1条无定义记录；同链限制下也是82条。反向RMSF组及不限制SASA的两行均为83条。因此应写“原83记录队列中82条满足该匹配条件”，不能写“全部83条均成立”。SASA10的记录内Hχ组仍有11条记录一致率≤50%，其中9条低于50%、2条等于50%；严格说是11条不优于该基准的记录，不能把平局也说成严格负效应。同链Hχ组相应有12条（10条低于、2条等于），不要误用11条。CSV各行均保留对应计数。

所有匹配依赖真实MD指标，是解释性分析。结论可写“两种专属输出在另一读数近似相同的位点中仍提供互补排序信息，且对SASA与同链控制稳健”；不能写成可部署时已知真实RMSF的预测器、接触或结合的因果效应、四种动态组合均可靠恢复、联合训练优于两个独立模型，或每个体系都成功。原全队列异向组合LH/HL召回只有15.79%/16.48%，不应被本表的条件排序结果遮蔽。

来源为`main_numbers.json`（前四个对应条目的点估计已逐项一致断言）、`matched_summary_family_cluster_bootstrap.csv`以及`sasa_absolute_sensitivity/summary_family_cluster_bootstrap.csv`。同链两行来自后者，因为旧main_numbers.json没有同链条目。负例及无匹配记录数来自对应的逐记录源表。导出代码为`export_figure3_matching_contract.py`。
'''
    (HERE/'figure3_matching_contract.md').write_text(prose, encoding='utf-8')
    print(result[['row_id', 'defined_records', 'matched_pairs', 'full_accuracy_pct',
                  'intra_accuracy_pct', 'surrogate_accuracy_pct', 'full_minus_surrogate_pp',
                  'full_minus_surrogate_ci95_low_pp', 'full_minus_surrogate_ci95_high_pp',
                  'full_at_or_below_chance_records']].to_string(index=False))


if __name__ == '__main__':
    main()
