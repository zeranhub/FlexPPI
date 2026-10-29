# 同类型匹配与主链/侧链互补：审稿式结果说明

本分析锁定原先83条共同窗口记录，使用不变的OOF预测，在同记录同氨基酸类型内穷尽全部符合规则的残基对。主要真实指标来自0—100 ns、每200 ps保存一次、共501帧的统一窗口；可用独立重复均值与原标签作为敏感性分析。规则在结果计算前写入protocol.json；没有依据模型预测或成功程度筛选位置。

这是一项以真实MD读数分层的解释分析。它不能被称为用可部署输入构造的新预测器，不能从相关或匹配对直接得出结构约束的因果效应，也不能单独证明结合导致的动态变化。

83条记录共有33,493个共同位置，其中同类型至少4个位点的组覆盖33,142个位置。可用家族映射仅包含42个簇，最大簇含38条记录；必须以整簇重抽评价不确定性。家族簇置信区间不能消除原始OOF训练/测试之间可能存在的家族相关。

| 匹配方向 | 条件 | 指标 | 记录宏平均 [家族簇95%区间] | 去最大簇 |
|---|---|---|---:|---:|
| hchi_to_rmsf | percentile_0p25 | target_true_abs_delta_mean | 0.6846 [0.5928, 0.9146] | 0.7949 |
| hchi_to_rmsf | percentile_0p25 | full_direction_accuracy | 0.7983 [0.7390, 0.8266] | 0.7555 |
| hchi_to_rmsf | percentile_0p25 | surrogate_direction_accuracy | 0.5949 [0.5861, 0.6165] | 0.6049 |
| hchi_to_rmsf | percentile_0p25 | intra_direction_accuracy | 0.7984 [0.7396, 0.8239] | 0.7602 |
| hchi_to_rmsf | percentile_0p25 | full_minus_surrogate_accuracy | 0.2034 [0.1333, 0.2378] | 0.1506 |
| hchi_to_rmsf | percentile_0p25 | full_minus_intra_accuracy | -0.0001 [-0.0208, 0.0122] | -0.0048 |
| hchi_to_rmsf | physical_0p10 | target_true_abs_delta_mean | 0.7085 [0.6142, 0.9503] | 0.8244 |
| hchi_to_rmsf | physical_0p10 | full_direction_accuracy | 0.7987 [0.7389, 0.8256] | 0.7589 |
| hchi_to_rmsf | physical_0p10 | surrogate_direction_accuracy | 0.5876 [0.5775, 0.6113] | 0.5952 |
| hchi_to_rmsf | physical_0p10 | intra_direction_accuracy | 0.7918 [0.7203, 0.8229] | 0.7453 |
| hchi_to_rmsf | physical_0p10 | full_minus_surrogate_accuracy | 0.2112 [0.1415, 0.2431] | 0.1637 |
| hchi_to_rmsf | physical_0p10 | full_minus_intra_accuracy | 0.0069 [-0.0085, 0.0344] | 0.0136 |
| rmsf_to_hchi | percentile_0p25 | target_true_abs_delta_mean | 0.3331 [0.3090, 0.3435] | 0.3184 |
| rmsf_to_hchi | percentile_0p25 | full_direction_accuracy | 0.7132 [0.6526, 0.7450] | 0.6643 |
| rmsf_to_hchi | percentile_0p25 | surrogate_direction_accuracy | 0.5666 [0.5576, 0.5896] | 0.5783 |
| rmsf_to_hchi | percentile_0p25 | intra_direction_accuracy | 0.7081 [0.6460, 0.7403] | 0.6590 |
| rmsf_to_hchi | percentile_0p25 | full_minus_surrogate_accuracy | 0.1466 [0.0736, 0.1851] | 0.0860 |
| rmsf_to_hchi | percentile_0p25 | full_minus_intra_accuracy | 0.0051 [-0.0035, 0.0144] | 0.0054 |
| rmsf_to_hchi | physical_0p10 | target_true_abs_delta_mean | 0.3349 [0.3032, 0.3490] | 0.3153 |
| rmsf_to_hchi | physical_0p10 | full_direction_accuracy | 0.6975 [0.6198, 0.7342] | 0.6418 |
| rmsf_to_hchi | physical_0p10 | surrogate_direction_accuracy | 0.5445 [0.5267, 0.5885] | 0.5642 |
| rmsf_to_hchi | physical_0p10 | intra_direction_accuracy | 0.7023 [0.6265, 0.7385] | 0.6479 |
| rmsf_to_hchi | physical_0p10 | full_minus_surrogate_accuracy | 0.1530 [0.0525, 0.2025] | 0.0777 |
| rmsf_to_hchi | physical_0p10 | full_minus_intra_accuracy | -0.0048 [-0.0251, 0.0122] | -0.0061 |

所有方向一致率都先在同类型残基对内求均值，再对类型和记录依次等权。真实差为零的对仅从方向分母剔除，预测平局计0.5。代用输出指用Full的预测RMSF排列真实Hχ，或反向用预测Hχ排列真实RMSF；专属输出优于代用输出，才支持“仅靠另一读数不够”的增量判断。不能把匹配对数量当作独立样本数。

四组合分类沿用每记录每类型的25%与75%分位，低于或等于25%分位为低，高于或等于75%分位为高；分位重合、类型少于4位点或处于中段均保留为未分类。全3063记录原分类与原召回均逐项重现。此类高低为相对秩分类，不能解读为绝对物理状态；异向组合的召回必须独立报告，不能用同向组合的较高召回概括全部互补能力。

主要表格：primary_matched_summary.csv；全部来源、原标签/重复均值敏感性、失败记录及各规则结果见matched_summary_family_cluster_bootstrap.csv和weak_or_negative_records_all_rules.csv。四组合汇总见four_combinations_summary_family_cluster_bootstrap.csv。matched_per_record_AA.csv.gz、matched_per_record.csv及uniform_rep1_matched_pairs.csv.gz提供逐层源表。所有原始数据和Word均未改动。
