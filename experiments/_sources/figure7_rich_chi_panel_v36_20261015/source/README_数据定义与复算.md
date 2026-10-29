本目录保存同一条既有 A003R 泛素 MD 轨迹的全部有效侧链位点数据。位点范围严格等于原冻结 Full 预测文件中 A 链 `valid[:,1]` 为真的 66 个位置；没有按预测高低、MD 曲线形状或误差筛选位置。保留 P19、P37、P38 三个 Pro 的原联合熵零参考。时间协议与 v27 及群体实验相同：0–100 ns，每 200 ps 一个观测，共 501 帧，对应原轨迹索引 0、2、4、……、1000。这是已有轨迹分析，没有新 MD、训练或模型推理。

66 个位置共定义 151 个侧链角：χ1 有 66 个，χ2 有 52 个，χ3 有 23 个，χ4 有 10 个。四角统一保存的 264 个槽位中，其余 113 个是化学上未定义的角。未定义角的角度、计数、概率、F*、宽度和细箱熵全部保持缺失。已定义角中的未访问箱则保存计数 0、概率 0、F* 缺失，两种情况通过 `chi_defined` 区分。不能把未定义角填为 0，不能把没有观测到的箱填成某个人工高 F* 值。

`ubiquitin_all66_chi_501frames.npz` 是完整可复算数值源。`chi_degrees` 的排列为 501 个时间点 × 66 个位置 × 4 个角，`chi_defined` 标出实际存在的角，`chi_count` 是各位置的定义角数量；`residue_index`、`chain`、`author_residue_id`、`amino_acid`、`residue` 保存原身份映射。`predicted_Hchi_nats` 和 `frozen_MD_joint_Hchi_nats` 逐值复制原冻结模型与标签，不更改定义。`counts_origin0`、`probability_origin0`、`Fstar_origin0` 按 66 个位置 × 4 个角 × 72 个角度箱排列；`bin_centers_degrees` 为 −177.5°、−172.5°、……、177.5°。`width_per_origin_degrees` 与 `fine_entropy_per_origin_nats` 末维分别对应原点 0°、1.25°、2.5°、3.75°。

`all66_angles_long_501frames.csv` 保存全部 132264 个时间点、位置和角槽组合，包括 75651 个实际角度和 56613 个缺失角度。`all66_histograms_origin0.csv` 保存全部 19008 个位置、角槽和箱组合，包括 10872 个有效角度箱记录和 8136 个缺角箱记录。此表可绘制逐位置逐角占据图或 F* 曲线。`all66_per_angle_metrics.csv` 保存 264 个角槽的逐角指标及四个原点各自结果。`all66_residue_summary.csv` 保存 66 个位置的身份、预测、原 MD 联合熵参考和只对实际定义角取平均的汇总指标。

每个角用圆周等宽 5°箱统计，原点 0°时箱为 [−180°,−175°)、……、[175°,180°)，+180°回到 −180°。概率为箱计数除以 501；经验投影 F*=−ln(p/pmax)，只在已访问箱中定义。每个角的低 F* 宽度为全部 F*≤1 已访问箱的总角度覆盖，可以包含彼此分开的片段，并非某一个连续势阱的宽度。`width_origin0_degrees` 对应原点 0°曲线；`width_four_origins_mean_degrees` 平均四种原点的宽度。`fine_entropy_origin0_nats` 是 5°占据概率的 Shannon 熵（自然对数）；`fine_entropy_four_origins_mean_nats` 平均四种原点的细箱熵。这些逐角细箱指标与模型预测的原联合 Hχ 有不同定义，不能把逐角熵相加后标成原标签，也不能将 MD F* 或 MD 占据概率标成模型输出。

原始结构与轨迹的完整 SHA256 与 v13 来源相同。MDAnalysis 只读取此目录中临时复制并核验的轨迹，临时副本关闭后删除；ILE 的 CD 原子命名与 HIS 别名只在任务拓扑副本中按原校正代码规范化，原子顺序和轨迹坐标不变。所有 66 个位置都有其完整定义角。角度使用原未修改 `_dihedral_degrees` 函数，并由 MDAnalysis 双精度函数独立核对，最大圆周角误差 2.73×10⁻⁵°。全部六个 GLN 的 9018 个角值与 v27 501 帧结果逐值完全相同；原两位点 432 个箱的计数、概率、F* 也相同。

运行 `prepare_all_ubiquitin_chi.py` 生成数据，再运行 `verify_all66_from_saved_sources.py` 独立核验。核验脚本不导入提取脚本，从保存角度另用 Counter 与 math 重新统计，检查全部 604 个逐角逐原点组合、10872 个有效角度箱、132264 条长角度记录及缺失槽。最大 F* 重算差为 8.88×10⁻¹⁶，最大宽度或细箱熵差为 1.33×10⁻¹⁵；原预测身份、三个 Pro 零参考、全部来源哈希和旧 GLN 数值均通过。冻结计划见 `analysis_plan.json`，提取及原子映射见 `extraction_audit.json`，完整独立核验见 `independent_verification.json`。
