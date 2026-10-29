# -*- coding: utf-8 -*-
"""Generate a concise evidence ledger; never edits the manuscript or raw inputs."""
from pathlib import Path
import json
import pandas as pd

HERE = Path(__file__).resolve().parent


def select(frame, **filters):
    keep = pd.Series(True, index=frame.index)
    for key, value in filters.items():
        keep &= frame[key] == value
    found = frame[keep]
    assert len(found) == 1, filters
    return found.iloc[0].to_dict()


def main():
    base = pd.read_csv(HERE/'matched_summary_family_cluster_bootstrap.csv')
    sasa = pd.read_csv(HERE/'sasa_extension/matched_summary_family_cluster_bootstrap.csv')
    absolute = pd.read_csv(HERE/'sasa_absolute_sensitivity/summary_family_cluster_bootstrap.csv')
    classes = pd.read_csv(HERE/'four_combinations_summary_family_cluster_bootstrap.csv')
    ledger = {}
    cases = [
        ('RMSF差≤0.10 Å，考察Hχ', base, dict(direction='rmsf_to_hchi', match_rule='physical_0p10')),
        ('Hχ差≤0.10 nats，考察RMSF', base, dict(direction='hchi_to_rmsf', match_rule='physical_0p10')),
        ('SASA百分位差≤0.25、RMSF差≤0.10 Å，考察Hχ', sasa,
         dict(direction='rmsf_to_hchi', match_rule='sasa_q25_and_rmsf_physical_0p10')),
        ('SASA差≤10 Å²、RMSF差≤0.10 Å，考察Hχ（探索后检查）', absolute,
         dict(direction='rmsf_to_hchi', match_rule='sasa_abs_10A2_and_condition_abs_0p10', pair_scope='all_within_record')),
        ('SASA差≤10 Å²、Hχ差≤0.10 nats，考察RMSF（探索后检查）', absolute,
         dict(direction='hchi_to_rmsf', match_rule='sasa_abs_10A2_and_condition_abs_0p10', pair_scope='all_within_record')),
    ]
    lines = ['# 可用于3.3扩充的互补证据与边界', '',
        '此次分析得到一个可以跨位置、跨记录检验的结论：同一种氨基酸的主链活动接近时，侧链联合角度占据的多样性仍可明显不同；反向条件也成立，而且对应的专属预测输出比借用另一输出排序更准确。它说明两种读数提供互补的位点信息，超出了挑选一个环或两个谷氨酰胺位置作展示的证据范围。这里“同记录”在多链数据中指同一装配体；另给同链限制的稳健性检查。', '',
        '所有主分析使用原先选定83条记录、0—100 ns共同窗口和不变的OOF预测，先穷尽每记录同类型位置对，再对类型、记录依次等权。原33,493个共同位置中，1,296个同类型至少4位点的组保留33,142个位点；所有潜在同类型无序对共624,602个，匹配时不依据预测成败选择。原标签和可用独立重复均值作为敏感性对照。', '',
        '| 条件与目标 | 合格记录/残基对 | 真实目标绝对差宏均值 | 专属输出方向一致率 [整簇95%CI] | 另一输出代排 | 专属增量 [整簇95%CI] |',
        '|---|---:|---:|---:|---:|---:|']
    for title, frame, filt in cases:
        filt = dict(label_mode='uniform_rep1', target_subset='all_close', dataset_scope='all', **filt)
        rows = {metric: select(frame, metric=metric, **filt) for metric in
                ['target_true_abs_delta_mean', 'full_direction_accuracy', 'surrogate_direction_accuracy',
                 'full_minus_surrogate_accuracy', 'full_minus_intra_accuracy']}
        ledger[title] = rows
        f, s, delta, gap = (rows[k] for k in ['full_direction_accuracy', 'surrogate_direction_accuracy',
                                            'full_minus_surrogate_accuracy', 'target_true_abs_delta_mean'])
        lines.append(f"| {title} | {int(f['defined_records'])} / {int(f['matched_pairs']):,} | {gap['record_macro_mean']:.4f} | {100*f['record_macro_mean']:.2f}% [{100*f['family_cluster_ci95_low']:.2f}, {100*f['family_cluster_ci95_high']:.2f}] | {100*s['record_macro_mean']:.2f}% | +{100*delta['record_macro_mean']:.2f}点 [{100*delta['family_cluster_ci95_low']:.2f}, {100*delta['family_cluster_ci95_high']:.2f}] |")
    lines += ['', 'SASA百分位相近并不等于绝对面积近似相等：原百分位匹配的物理侧链面积差仍为记录宏平均14.00 Å²。因此另行、明确标注为探索后稳健性检查，固定绝对SASA差5、10、20 Å²并执行两个方向；主结果未被替换。在中心10 Å²阈值下，再限制匹配对来自同一条链，Hχ与RMSF专属输出相对另一读数的优势仍分别为14.29与22.79个百分点，家族簇区间分别为[2.60,19.41]与[13.00,27.00]个百分点。真实RMSF近似且面积近似并不消除接触几何等全部混杂，是否超出简单结构基线，仍应结合其他分支的强基线结果。', '',
        '结果支持“两个读数不能互相替代”，但不支持“完整伙伴信息在这里产生额外优势”：本分支所有主要匹配条件中，Full相对Intra的方向一致率差的家族簇区间均跨越零，且存在Intra点估计较高的条件。因此不能将互补性证据转写为伙伴导致的动态变化，也不能宣称联合训练本身优于两个单独训练的模型；后者需要相应训练消融。', '',
        '四组合全队列复核逐位点重现原3063条记录、824,988个位置的分类，原逐记录计数和召回全部一致。同向LL与HH的记录宏召回分别为54.62%与61.37%，异向LH与HL仅15.79%与16.48%。异向位置真实占比分别为3.18%与2.92%，70.05%位置是中间或未分类。故模型能够恢复部分条件性排序差异，但不能称为可靠恢复所有主链/侧链组合状态；分类采用记录内同类型分位，亦不能解释为绝对物理高低。', '',
        '不确定性以可得家族连通簇整簇重抽计算。83条记录实际上仅42个簇，其中最大簇含38条；30条ITrace全部位于同一个连通簇，因此ITrace内部的家族簇置信区间保留为空，不伪造独立重复。全体结果同时给出家族等权与去最大簇敏感性；后者仍保留专属输出增量，但幅度通常下降。家族簇重抽只处理评价样本之间的依赖，不能把原始OOF变成家族留出预测。', '',
        '这些匹配依赖真实MD RMSF或Hχ，属于解释性分析，不是部署时可用的输入。真实标签与预测的方向一致性不建立埋藏、接触或伙伴结合的因果机制；Hχ也不等同完整热力学熵。预先固定的目标百分位差至少0.50子集结果更高，但必须与全部匹配对主结果一起呈现，不能只选高值。所有原始文件、预测与Word均未修改。', '',
        '独立核验见independent_validation.json，代码与完整逐类型、逐记录、逐匹配对源表均保存在本目录。主口径详见protocol.json；SASA质量与匹配扩展见sasa_extension/sasa_extension_protocol.json；探索后绝对面积检查见sasa_absolute_sensitivity/protocol_post_exploration.json。']
    (HERE/'conclusions_for_section33_zh.md').write_text('\n'.join(lines)+'\n', encoding='utf-8')
    (HERE/'main_numbers.json').write_text(json.dumps(ledger, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    # Correct a prose shorthand in the earlier diagnostic report; all numbers unchanged.
    path = HERE/'sasa_extension/reviewer_interpretation_zh.md'
    path.write_text(path.read_text(encoding='utf-8').replace('在同一蛋白中', '在同一复合物记录中'), encoding='utf-8')
    readme = r'''# 复现方法与文件说明

本目录仅写派生分析，不改原始数据、OOF预测或Word。使用Python 3.9.13、NumPy 2.0.2、pandas 2.3.1；读取parquet另需当前环境已有的pyarrow。全部命令在本目录执行。

```powershell
python .\analyze_matched_complementarity.py --freeze-only
python .\analyze_matched_complementarity.py
python .\run_sasa_extension.py
python .\run_absolute_sasa_sensitivity.py
python .\verify_analysis_independently.py
python .\write_conclusions.py
```

先读conclusions_for_section33_zh.md和main_numbers.json。protocol.json是计算前固定的主规则，source_manifest.json记录337个只读来源hash。absolute SASA是探索后检查，独立目录和protocol明确标记，不替换最初主分析。

common_positions_with_uniform_truth_and_predictions.csv.gz给出共同窗口真值、原标签和不变预测。matched_per_record_AA.csv.gz、matched_per_record.csv分别为类型和记录层源表；uniform_rep1_matched_pairs.csv.gz保存满足任一主匹配口径的对，residue_index是从0开始的OOF顺序，delta均为j减i。所有量先在类型内、再在记录内等权；matched_pairs只是描述计数，不能作独立样本数。

matched_summary_family_cluster_bootstrap.csv含来源分层、全部阈值与三个真值模式；primary_matched_summary.csv为核心规则。weak_or_negative_records_all_rules.csv保留弱结果及负结果。four_combinations_*重现原全队列分类，并提供83共同窗口对照。sasa_extension是质量合格SASA百分位匹配；sasa_absolute_sensitivity是5/10/20 Å²绝对面积和same_chain_only检查。后者同时覆盖两个方向。

independent_validation.json由独立标量残基对枚举核验，另逐项验证原始来源未变、类型/记录聚合和原3063分类/召回。completion_audit.json是主分析完成时的快照；后续扩展各有自己的完成审计。gzip文件的压缩时间戳可在重跑时变化，不代表CSV内容变化。
'''
    (HERE/'README_zh.md').write_text(readme, encoding='utf-8')


if __name__ == '__main__':
    main()
