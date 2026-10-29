# SASA与RMSF同时相近条件下的侧链差异

质量与匹配规则见sasa_extension_protocol.json；规则在读取特征值前冻结。主结果只使用初始结构侧链完整、χ定义所需原子完整、标准残基且SASA有限的位置。每记录同类型合格位置少于4不进入匹配。SASA相近指同类型百分位差不超过0.25；真实RMSF另采用百分位差不超过0.25或绝对差不超过0.10 Å。匹配选择不依赖模型预测。

| 规则 | 指标 | 记录宏平均 [整簇95%区间] | 记录数 | 家族簇数 | 去最大簇 |
|---|---|---:|---:|---:|---:|
| sasa_q25_and_rmsf_physical_0p10 | target_true_abs_delta_mean | 0.3136 [0.2974, 0.3272] | 83 | 42 | 0.3116 |
| sasa_q25_and_rmsf_physical_0p10 | full_direction_accuracy | 0.6717 [0.5903, 0.7071] | 83 | 42 | 0.6185 |
| sasa_q25_and_rmsf_physical_0p10 | surrogate_direction_accuracy | 0.5272 [0.5084, 0.5737] | 83 | 42 | 0.5483 |
| sasa_q25_and_rmsf_physical_0p10 | intra_direction_accuracy | 0.6738 [0.5882, 0.7103] | 83 | 42 | 0.6197 |
| sasa_q25_and_rmsf_physical_0p10 | full_minus_surrogate_accuracy | 0.1445 [0.0372, 0.1927] | 83 | 42 | 0.0702 |
| sasa_q25_and_rmsf_physical_0p10 | full_minus_intra_accuracy | -0.0022 [-0.0283, 0.0250] | 83 | 42 | -0.0012 |
| sasa_q25_and_rmsf_q25 | target_true_abs_delta_mean | 0.2967 [0.2783, 0.3048] | 83 | 42 | 0.2887 |
| sasa_q25_and_rmsf_q25 | full_direction_accuracy | 0.6815 [0.6128, 0.7161] | 83 | 42 | 0.6294 |
| sasa_q25_and_rmsf_q25 | surrogate_direction_accuracy | 0.5413 [0.5280, 0.5750] | 83 | 42 | 0.5548 |
| sasa_q25_and_rmsf_q25 | intra_direction_accuracy | 0.6793 [0.6102, 0.7137] | 83 | 42 | 0.6273 |
| sasa_q25_and_rmsf_q25 | full_minus_surrogate_accuracy | 0.1402 [0.0576, 0.1823] | 83 | 42 | 0.0747 |
| sasa_q25_and_rmsf_q25 | full_minus_intra_accuracy | 0.0022 [-0.0112, 0.0151] | 83 | 42 | 0.0021 |

该分析可回答：在同一复合物记录中、同一种氨基酸、暴露程度和主链活动相近时，模型是否仍能区分侧链占据多样性。它不是去除所有结构混杂的因果证明，也不能将使用真实MD RMSF的分层伪装成部署条件。SASA百分位相近不意味着绝对面积完全相同，故另给primary_physical_sasa_balance_by_record.csv；若残留差较大，仍需结合强结构基线评价。
