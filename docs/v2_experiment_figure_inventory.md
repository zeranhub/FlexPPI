# FlexPPI v2 实验与图表重跑清单

盘点日期：2026-09-06。项目根目录：`E:\Project\FlexPPI`。本文只记录已查明的实验、脚本、数据依赖和改造约束，不表示实验已经重跑，不改变旧结果，也不预设新版胜过基线。文中相对路径均相对于项目根目录；新结果应写入独立的 v2 结果目录，不能覆盖旧缓存或论文图。

## 1. 评估边界和五折规则

新版是仅预测 RMSF 与侧链构象熵 Hχ 的五折模型。主评估应使用每个体系对应的留出折模型，即该体系没有参与训练的那一个模型；候选权重为 `artifacts/performance_first_v2/fold_N/epoch_040.pt`，执行前仍须核验实际 epoch、配置和 SHA-256。不能让每个样本通过全部五个模型后取平均作为 OOF 测试，因为其中四个模型训练过该样本。五折也不是同一测试集上五个随机种子的重复实验。

必须分别命名并独立报告三类结果：**历史交集复测（overlap）**只比较与旧论文队列相同的结构和残基，用于追踪版本变化；**完整 OOF**覆盖新版冻结划分的全部合格测试体系，是新版主结果；**真正的外部/零样本测试（external / zero-shot）**要求评估体系未进入对应模型训练，并审计训练重叠。ITrace 已进入新版训练来源，ITrace OOF 不得继续称为旧版的 zero-shot 或完全外部测试。第三方模型的开发域与训练重叠状态同样应单独记录；缓存有预测不等于训练无重叠。

每个样本的 full、pMHC-only、TCR-only、isolated-chain 输入都必须沿用原完整体系的同一个留出折，而不能对子结构重新分配 fold。原病例如果属于 OOF 清单就用其留出折；不在清单中的结构需要先审计与训练集的对应关系，不能自动当成外部测试。每行预测至少记录原体系 ID、fold、condition、checkpoint SHA、输入结构 SHA、残基身份和标签来源；缓存复用也验证这些标识，不能仅凭文件存在或列名齐全。

## 2. 用户五张图的确切映射

| 图组 | 原脚本或表格 | 主要输入 | 原输出 | v2 所需改造 |
|---|---|---|---|---|
| 1：DynaRepo 复合体、DynaRepo 单链、ATLAS 单链的七模型 RMSF 横评表 | `paper/tables/rmsf_results.tex` 为与截图完全匹配的手写数值表；统计生成器为 `scripts/evaluate_additional_rmsf_baselines.py` | `artifacts/evaluation/additional_baselines/results/static_summary_metrics.csv`、`static_per_system_metrics.csv`、`static_common_residue_predictions.csv`；上游为 `artifacts/evaluation/fastprotflex_comparison/` 和 `additional_baselines/raw/` | `paper/tables/rmsf_results.tex`；统计 CSV 在 `additional_baselines/results/` | 从新 OOF 逐残基预测重算，不手动挪用旧数值；分别导出 overlap 与 full OOF；补每折指标、体系等权汇总和覆盖率。`scripts/generate_benchmark_table.py` 只覆盖原 FlexPPI / FastProtFlex 两模型，不是这张七模型表的完整生成器。 |
| 2：Hχ、RMSF 预测—真值 hexbin 及边际分布 | `scripts/make_flexppi_joint_scatter_pair.py` | Hχ：`output/itrace_hchi_benchmark/raw/hchi_consensus_per_residue.csv`；RMSF：`artifacts/evaluation/additional_baselines/results/itrace_common_residue_predictions.csv` | `output/flexppi_performance_figures/fig_flexppi_hchi_rmsf_joint.{pdf,png}`，另有 `fig_flexppi_rmsf_panel_b` | 改用 OOF 行；图注标明 pooled-residue 展示与体系等权指标的不同。原两面板各自 dropna，若要作共同残基两任务比较则应另冻结两标签均有效的交集。修正旧 `\\AA` 字符串显示为 Å，并报告显示范围外点数。 |
| 3：A000K、1H16 残基曲线与方法×位置热图 | `paper/figures/make_rmsf_baseline_comparison_figures.py` | `additional_baselines/results/static_common_residue_predictions.csv`；历史选择记录在 `paper/figures/fig_rmsf_baseline_comparison_statistics.json` | `paper/figures/fig_static_residue_profiles.{pdf,png,svg}` 及 `fig_static_residue_profiles_source_data.csv` | 固定旧样本和旧窗口：`A000K` A 链 209–298，`1h16_A` 666–755。不能再调用 `select_advantage_segment` 按新版优势重选窗口。可另外提供按预先定义中位配对差选择的典型胜/负案例，但须与历史固定窗口区分。 |
| 4：按 Peptide / TCR beta / TCR alpha / MHC / beta2m 的 RMSF 与 Hχ 上下文消融 | RMSF：`scripts/make_flexppi_rmsf_context_figure.py`；Hχ：`scripts/make_flexppi_context_ablation_figure.py` | RMSF：`output/itrace_partner_ablation/per_chain_rmsf_metrics.csv` 与 `output/itrace_zero_shot_context/metrics/per_chain_metrics.csv`；Hχ：`output/itrace_hchi_multi_benchmark/metrics/per_chain_condition_hchi.csv` | `output/flexppi_performance_figures/fig_flexppi_rmsf_context_ablation.{pdf,png}`、`fig_flexppi_context_ablation_panel_a.{pdf,png}` | 使用同体系同 fold 的配对输入预测；共同有效体系上做配对统计，显示各角色样本数。旧 Hχ 图的 0.882 先验线必须移除并改为训练折先验实测值；旧结论标题不能照抄。 |
| 5：MD 与预测各自的 RMSF—Hχ 耦合 hexbin | `scripts/make_flexppi_rmsf_hchi_coupling_figure.py` | `output/itrace_residue_biology/raw/all_condition_predictions.csv`、Hχ consensus、`output/itrace_zero_shot_context/raw/common_context_predictions.csv`；按 `pdb_id,chain,residue_number,insertion_code` 合并，选择 full 且两个真值与两个预测均有效 | `output/flexppi_performance_figures/fig_flexppi_rmsf_hchi_coupling.{pdf,png}` | 严格一对一残基连接，保留共同有效遮罩；将 pooled-residue 耦合与每链/每体系耦合分开，避免长度与 AA 类型组成造成总体相关。不得将低相关直接表述为数学上的独立或正交。 |

原图 3 的 `select_advantage_segment` 遍历连续窗口，要求 MD 峰位于中央、稳健变化范围至少为两个链内标准差，再最大化 FlexPPI 相对最强基线的 z-MAE 优势。旧图和 README 已注明其属于优势案例，不是无偏代表性样本。新版纵向复测应固定上表窗口并如实呈现退步或优势消失。

## 3. 统计口径必须分开

原七模型横评在每条链内计算 Spearman、Pearson、z-MAE、top-20% AUROC/AP/P@K 等指标，然后在复合体内平均链，最后使每个体系等权。对于多链体系，这个 **chain-macro → system-macro** 口径与将整个复合体所有残基合并后计算一个相关系数的 **whole-system** 口径并不相同。新版 `scripts/aggregate_performance_first_v2_results.py` 汇总 fold 内全残基和整个体系指标，不能直接把其均值填进原七模型表。需要从逐残基 OOF 预测统一重算，并在结果字段、图注和文件名中区分两种口径。

五折均值±标准差用于展示不同划分的波动，体系/分组 bootstrap CI 用于展示测试集合上的统计不确定性，两者不是一回事。应保留每折结果、完整 OOF 体系等权结果和 paired difference；不要将五折当成大量独立样本，也不要只凭不同误差条是否重叠判断显著性。现有聚合器使用 `np.std` 默认 `ddof=0`，新导出若调整为样本标准差需明确说明且保持一致。所有模型必须采用相同合格残基遮罩；预测缺失不能为不同方法悄悄更换测试子集，应报告模型覆盖率，并区分全覆盖主榜与有限重叠补充榜。

FastProtFlex、ProDy、RMSF-net-PDB 原协议仅作为相对涨落预测，绝对 Å 误差不应与有绝对 RMSF 输出的模型直接混排。标准化应按冻结的整条链及固定遮罩计算，再展示案例窗口，不要在不同窗口独立标准化后与原表混用。复合体输入也需记录：旧协议 BackFlip/FastProtFlex 对复合体逐链预测，FlexPPI/ProDy/RMSF-net-PDB 读完整复合体；这是实际使用方式横评，不等于输入信息完全对称的机制消融。

## 4. 统一上下文预测器的兼容问题

`scripts/analyze_itrace_residue_biology.py` 的 `condition_specs` 已具备每体系八种输入：full、pmhc、tcr、五条 isolated 链，适合作为新版统一上下文数据源，避免多个实验重复推理。但其 `MODEL_METRICS = ("rmsf", "chi_entropy", "sasa")`，`run_predictions` 会在新版缺少 SASA 时直接报错，且多个下游循环与布局假定三项指标。`scripts/analyze_fixed_partner_ablation.py` 也将 `RESPONSE_METRICS` 写死为三项。需要以新版实际的两个头为准改造所有下游数据与布局，不能给 SASA 填零、借旧 SASA 或声称新版仍支持 SASA。

`scripts/evaluate_itrace_flexppi_context.py` 原先通过一个 checkpoint SHA 确保完整复合体与孤立链来自同一模型；五折应改为每体系 checkpoint 映射校验，不能简单删掉该一致性检查。旧 `analyze_itrace_residue_biology.py --reuse` 只查文件与指标列，不能区分新旧权重，必须补来源校验或完全隔离新缓存。固定坐标伙伴删除实验还要验证保留原子的坐标未变化，保留编号、插入码、链角色和界面遮罩；其含义仍是模型对输入上下文的敏感性，而非真实 apo/holo 动力学变化。

`scripts/build_itrace_hchi_benchmark.py` 与 `scripts/evaluate_itrace_hchi_multi.py` 的旧 AA-type prior 从被评估队列本身计算 AA 类型均值；新版必须从每个 fold 的训练标签计算，再用于对应测试折。旧绘图硬编码的 `0.882` 不能沿用。除了普通 Hχ 指标，建议补充同 AA 类型内排序和减去训练折 AA 均值后的残差排序，明确区分类型基线与结构上下文贡献；预测和真值同时减去同一先验后的 MAE 与原 MAE 数学上相同，不能将其作为新的独立提升证据。

## 5. 补充实验与图表支线

| 支线 | 原分析/绘图脚本 | 数据与原输出 | 重跑注意点 |
|---|---|---|---|
| ITrace 13 区域 × 7 模型、逐体系分布与典型胜负案例 | `scripts/evaluate_additional_rmsf_baselines.py` → `paper/figures/make_rmsf_baseline_comparison_figures.py` | `additional_baselines/results/itrace_*` → `paper/figures/fig_rmsf_baseline_comparison.*` | 与图 1 共享 baseline 预测；固定区域与统计规则，报告 fold 和链角色覆盖率。 |
| AlphaFlow 小样本 overlap 与帧数收敛 | `scripts/evaluate_alphaflow_ensembles.py` → 同上绘图脚本 | `additional_baselines/` 内 AlphaFlow 结果 → `paper/figures/fig_alphaflow_overlap.*` | 旧图为六体系配对与 10–250 帧收敛，不能混成 244 体系全覆盖主榜；保留独立补充结果。 |
| Hχ 总体、区域表现、MD 三重复可靠性 | `scripts/build_itrace_hchi_benchmark.py`、`scripts/make_flexppi_hchi_performance_figure.py`、`make_flexppi_hchi_panel_a.py` | `output/itrace_hchi_benchmark/` → `output/flexppi_performance_figures/fig_flexppi_hchi_performance.*` 等 | MD 重复可靠性不随模型变化，可校验后复用；AA prior 改为 train-only。 |
| TCR CDR1/2/3、framework 的 RMSF/Hχ 上下文消融 | `scripts/make_flexppi_tcr_region_context_figures.py` | all-condition 预测 + 两标签 + region → `fig_flexppi_tcr_region_rmsf_context.*`、`fig_flexppi_tcr_region_hchi_context.*` | 复用统一八条件 OOF；短区域相关系数的最少残基数/常数处理/缺失规则需冻结。 |
| CDR3β 失败模式、动态幅度、位置和 AA 偏差 | `scripts/make_flexppi_cdr3b_diagnosis_figure.py`、`make_flexppi_cdr3b_hchi_diagnosis_figure.py` | all-condition 与标签 → `fig_flexppi_cdr3b_diagnosis.*`、`fig_flexppi_cdr3b_hchi_diagnosis.*` | 与新版总体提升并列展示局部失败，不预设旧失败或旧机制解释仍成立。 |
| 固定坐标伙伴删除、接触距离响应与界面定位 | `scripts/analyze_fixed_partner_ablation.py` | `output/itrace_residue_biology/raw/matched_residue_context_responses.csv` → `output/itrace_partner_ablation/`、`paper/data/itrace_partner_ablation/`、`paper/figures/fig_itrace_partner_ablation.*` | 两头化；原距离/坐标身份分析可核验复用；按体系/相关分组进行配对统计与多重比较校正。 |
| 孤立链敏感性完整主图 | `scripts/evaluate_itrace_flexppi_context.py` → `paper/figures/make_itrace_context_ablation_figure.py` | `output/itrace_zero_shot_context/` → `paper/figures/fig_itrace_context_ablation.*` | 新输出不应沿用 zero-shot 语义；改为每体系 fold/SHA 一致性校验。 |
| 肽位置柔性、接触比例、AA 效应 | `scripts/analyze_itrace_residue_biology.py`、`scripts/summarize_tcr_peptide_residue_biology.py` → `paper/figures/make_itrace_residue_biology_figure.py`、`make_itrace_biology_figure.py` | `output/itrace_residue_biology/`、`paper/data/itrace_residue_biology/` → 对应 peptide biology 图 | 删除新版不支持的 SASA 面板；位置、AA 规律和接触解释作为相应证据等级的描述，不将敏感性自动解释为机制。 |
| 6UON 与 2PTT 静态结构案例 | `paper/figures/make_case_study_figures.py` | ITrace benchmark / conditions、DynaRepo common predictions；本地 `case_studies/render_3dmol.js` 渲染 → `fig_6uon_case_study.*`、`fig_2ptt_case_study.*` | 脚本硬编码旧 `CHECKPOINT`；需要替换为对应 OOF 预测来源。可复用结构对齐与相机；保留 2PTT His42 无 MD 映射的遮罩，不填零。 |
| 7L1D / 7RRG 静态 PyMOL 上下文案例 | `scripts/analyze_tcr_pmhc_context.py`、`scripts/analyze_tcr_binding_pmhc.py`、`scripts/create_tcr_pmhc_pymol_comparison.py` | `artifacts/evaluation/7L1D_7RRG_tcr_pmhc_flexppi/` 与 `artifacts/evaluation/tcr_pmhc_pymol/` | 默认 checkpoint 为旧 pretrain；重跑前审计病例所属 fold，保存新版着色、预测与运行元数据，不覆盖旧 PyMOL 会话。 |
| TCR3D 结构代理评测 | `scripts/compare_tcr3d_fastprotflex.py` 及项目原 TCR3D 评估入口 | `artifacts/evaluation/tcr3d_static/`、`tcr3d_full/`、`fastprotflex_comparison/tcr3d_static/` | 无 MD 真值，仅能作为结构代理或定性支线，不混入主 RMSF 真值横评。旧列名可能保留 finetune 字样，以报告实际 checkpoint 为准。 |
| 旧配对响应、Zenodo、apo/holo | `scripts/evaluate_zenodo_6638504.py` 及原配对评估入口 | `artifacts/apo_holo_validation/`、`artifacts/zenodo_6638504/evaluation/` | 新版关闭 response head，原 response-head 任务不能原样复跑；如新增两状态预测之差，应另定协议并明确与原响应头不同。 |
| 纯数据诊断 | `scripts/analyze_dynarepo_interface_dynamics.py`、训练偏差/ ThermoPCD 审计脚本 | `artifacts/dynarepo_interface_dynamics*`、`training_bias_with_zenodo_6638504/`、`thermopcd_audit/` | 不依赖模型权重，输入/标签未变时无需重复计算；ThermoPCD 缺标签候选不能悄悄纳入有监督主评估。 |
| 历史实验总览 | `paper/figures/make_completed_experiments_figure.py` | 硬编码 pretrain、finetune、pair validation、TCR3D 报告 | 按新版两头 OOF 任务重新组织总览；不能只替换局部数值后继续沿用旧三任务结论。 |

另有 `scripts/make_flexppi_metric_orthogonality_figure.py` 读取 `per_complex_rmsf_vs_hchi.csv` 和残基级数据，可作为图 5 的体系级诊断补充，但同样要区分低相关、类型组成效应与真正的独立性。`scripts/organize_all_experiment_results.py` 对旧结果进行了归档，`output/all_experiment_results/` 不是新模型重跑输入的自动权威版本；读取时以明确来源和校验为准。

## 6. 推荐执行依赖顺序与交付检查

1. 冻结 run ID、五折模型 SHA、OOF manifest、overlap 队列、结构/标签/残基遮罩、第三方模型输入模式与训练重叠审计。确认用户选择的主图/补充图范围，再启动相应计算。
2. 建立统一的逐残基两头 OOF 输出与可验证缓存。已有远程逐残基 OOF 若能验证身份与遮罩，可复用；不能只凭 fold 汇总 JSON 代替原图所需残基数据。
3. 对 ITrace 按原体系 fold 执行八条件推理，先小样本验证完整体系与保留链残基映射及坐标一致性，再跑全队列。完整输入若已等同 OOF 输入可复用，避免重复推理。
4. 校验已有 baseline 缓存的结构、残基与运行配置；对固定历史 overlap 优先复用一致预测，只为 full OOF 扩展队列和新增模型补跑缺失部分。新增预测失败与不支持的结构须列明，不暗中改变全模型共同队列。
5. 从统一残基表生成 chain-macro、whole-system、fold 明细、完整 OOF 汇总、配对差与不确定性；Hχ train-only AA prior、CDR/角色/接触距离分析全部从同一冻结表派生。
6. 重画五张主线图，再根据选定范围生成区域、失败模式、界面距离、结构案例等支线；每图同时保存源数据、统计定义、输入/权重来源、图注和视觉 QA。旧图不可覆盖。
7. 交付独立新版图表包，并列出已完成、未运行、缺数据/环境、旧任务与新版不兼容的项目。不能把待跑计划、缓存存在或绘图脚本成功执行当成实验已完成。

MD RMSF/Hχ 标签、三重复可靠性、链角色/CDR 编号、接触距离及固定坐标子结构在输入与定义未变时可校验复用；只有模型依赖的预测、对应统计和基于它们的图表需要更新。新绘图代码应通过参数指定输入/输出，不再使用旧脚本中大量硬编码的 `E:\Project\FlexPPI\output\...` 路径、旧 checkpoint 名或固定结论标题。
