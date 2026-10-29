# ITrace 测评整合说明

> 历史记录：以下说明描述早期版本，不是 2026-09-18 当前双目标正文的模型、划分与指标口径。当前稿已将 ITrace 纳入五折训练/测试，允许链家族共享，不称为独立外部测试。请以 `README.md`、`WRITING_REVISION_NOTES.md` 和主文件实际加载的 `sections/appendix_current.tex` 为准；本文保留用于追溯旧实验。

本目录中的 ITrace 内容来自 `E:\Project\FlexPPI\output\itrace_benchmark`。
正文主分析使用排除精确训练 PDB 重叠 3KXF 后的 244 个 TCR–pMHC 复合体，
每个复合体包含三条独立 MD 重复；FastProtFlex 对各链分别推理，FlexPPI 对完整
复合体一次推理。主标签是对每条链分别做刚体拟合后计算的 Cα RMSF，三条重复
在残基层面平均。指标先在链内计算，再对同一复合体的链等权平均，最后对复合体
等权平均，防止长链或长复合体支配结果。

论文使用的可审计数值保存在 `data/itrace/`、`data/itrace_context/`、
`data/itrace_biology/` 与 `data/itrace_residue_biology/`。正文的区域级生物学
主图为 `figures/fig_itrace_biology.pdf`：
它把每条链内标准化后的 MD、完整复合体预测和孤立链预测并列，显示 MHC 接触面
相对刚、TCR 接触面尤其 $\beta$ 链相对柔，以及完整模型把这种不对称压成普遍
刚化。图由 `figures/make_itrace_biology_figure.py` 重建，统计定义与 QA 记录
分别位于同名 statistics JSON 和 QA Markdown。算法基准与标签敏感性图仍为
`figures/fig_itrace_main.pdf` 和 `figures/fig_itrace_sensitivity.pdf`。
实验流程在正文中按“MD 真值生成、完整复合体预测、同模型孤立链消融”直接说明，
不再使用协议要素表。正文明确把该队列称为“精确 PDB 外部队列”：
它能排除同一 PDB 结构直接重叠，但尚未完成对训练集的序列近同源和结构拓扑
近同源隔离，因此不能据此宣称严格的新折叠泛化。

同一冻结 FlexPPI 检查点的完整复合体/结合态孤立链上下文消融来自
`E:\Project\FlexPPI\output\itrace_zero_shot_context`。该消融没有参数更新，
完整和孤立输入使用相同目标链坐标；主队列仍排除 3KXF。论文图
`figures/fig_itrace_context_ablation.pdf` 的可审计源数据保存在
`data/itrace_context/`，绘图与数据整理脚本为
`figures/make_itrace_context_ablation_figure.py`。完整减孤立只能解释为
输入上下文响应，不代表有真实 apo 对照的结合诱导 $\Delta$RMSF。

残基级扩展来自 `output/itrace_residue_biology`。每个复合体使用同一冻结检查点
计算完整体、pMHC、TCR 二聚体和五条孤立链共八种固定坐标输入，每次同时保存
RMSF、$H_\chi$ 和 SASA。`figures/fig_itrace_residue_biology.pdf` 只把 ITrace
轨迹提供的 RMSF 称为 MD 真值；$H_\chi$/SASA 明确标注为模型上下文响应。
Ala/Gly 不作 $H_\chi$ 解释，接触使用 $5\,\text{\AA}$ 重原子阈值，统计以
复合体 cluster bootstrap。全队列推理脚本为
`scripts/analyze_itrace_residue_biology.py`，同肽 7L1D/7RRG 残基表由
`scripts/summarize_tcr_peptide_residue_biology.py` 整理。
