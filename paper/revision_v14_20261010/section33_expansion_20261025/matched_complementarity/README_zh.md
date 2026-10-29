# 复现方法与文件说明

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
