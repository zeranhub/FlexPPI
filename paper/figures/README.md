# Figures

将论文图片放在本目录中，优先使用 PDF、PNG 或 JPG 格式。

## 扩展 RMSF 基线图

`make_rmsf_baseline_comparison_figures.py` 从
`artifacts/evaluation/additional_baselines/` 的冻结结果生成双栏论文图：

- `fig_static_residue_profiles.{pdf,png,svg}`：用统一的逐链标准化坐标展示 DynaRepo
  多链复合体与 ATLAS 单链中 MD、FlexPPI 及六种基线的高变化连续片段；每个体系同时
  给出曲线图和“方法乘残基位置”矩阵热图。候选窗口要求 MD 高柔性峰位于中央且稳健
  变化范围至少为 2 个链内标准差，再按 FlexPPI 相对最强基线的标准化 MAE 优势选择；
  因此该图明确属于优势案例可视化，队列总体结论仍由完整残基统计给出。
- `fig_rmsf_baseline_comparison.{pdf,png,svg}`：244 个 ITrace 复合体的分布、
  FlexPPI 与 BackFlip-seq 成对散点、完整的 13 区域乘 7 方法热图，以及按预先定义
  规则选择的典型胜/负案例残基曲线。
- `fig_alphaflow_overlap.{pdf,png,svg}`：独立处理的 AlphaFlow 六体系配对比较、
  10--250 帧采样收敛和 6TGK-C 绝对尺度案例。

对应的 `*_source_data.csv` 保存每一个绘图值；
`fig_rmsf_baseline_comparison_statistics.json` 记录样本定义、bootstrap 设置、配对检验、
案例选择规则和输入文件来源。在项目根目录执行：

```powershell
python paper\figures\make_rmsf_baseline_comparison_figures.py
```

正文中的插图示例：

```tex
\begin{figure}[htbp]
  \centering
  \includegraphics[width=0.85\textwidth]{model_overview.pdf}
  \caption{FlexPPI 模型整体结构}
  \label{fig:model-overview}
\end{figure}
```

## 已完成实验汇总图

`fig_completed_experiments.pdf` 是论文使用的矢量主文件，
`fig_completed_experiments.png` 是 300 dpi 预览图。图中（a）只包含 RMSF
和侧链构象熵，不包含 SASA；（b）比较 FlexPPI 与 FastProtFlex；（c）比较
三种配对响应估计；（d）展示 TCR3D 全部 1,187 条可计算链的相关系数分布。

在项目根目录执行下面的命令可由归档的 JSON/CSV 结果重新生成图、源数据表和
统计说明：

```powershell
python paper\figures\make_completed_experiments_figure.py
```

配套文件为 `fig_completed_experiments_source_data.csv`、
`fig_completed_experiments_statistics.json` 和
`fig_completed_experiments_QA.md`。

## ITrace 固定坐标伙伴消融主图

`fig_itrace_partner_ablation.{pdf,png,svg}` 使用完整五链对 pMHC-only/TCR-only
的受控比较，不拆散目标侧内部组成。六个面板依次展示结合态 RMSF 的逐复合体
Spearman、九项准确性中的七项无量纲配对优势、三类输出随伙伴距离的衰减，以及
复合体内距离--响应相关。对应分析同时逐原子验证 488 组子复合物坐标与完整拓扑
完全一致，并输出复合体 bootstrap、Wilcoxon 检验和 FDR 校正。

在项目根目录执行：

```powershell
$env:PYTHONPATH=(Resolve-Path -LiteralPath 'src').Path
python scripts\analyze_fixed_partner_ablation.py
```

完整源数据位于 `paper/data/itrace_partner_ablation/`，正文数值表为
`paper/tables/itrace_partner_ablation.tex`。该图只表示固定结合态坐标的模型
上下文敏感性，不表示 apo/holo 变化。

## ITrace 孤立链敏感性消融图

`fig_itrace_context_ablation.pdf` 为论文使用的矢量主文件，`.png` 为
300 dpi 预览，`.svg` 保留可编辑文字。四个面板依次展示逐复合体 Spearman
配对、七项指标优势、链角色异质性和孤立短肽的绝对 RMSF 尺度外推。

在项目根目录执行下列命令可重新生成图、精简源数据与统计记录：

```powershell
python paper\figures\make_itrace_context_ablation_figure.py
```

配套文件为 `fig_itrace_context_ablation_statistics.json`、
`fig_itrace_context_ablation_QA.md` 和 `paper/data/itrace_context/`。

## ITrace peptide 残基级生物学图

`fig_itrace_residue_biology.pdf` 是双栏矢量主文件，`.svg` 保留可编辑文字，
`.png` 为 300 dpi 预览。三个面板依次直接比较完整复合物 FlexPPI 与九肽 MD
真值的位置柔性、展示加入 TCR 上下文后的逐位置 RMSF/$H_\chi$/SASA 与接触
比例，以及去除长度/位置均值后的氨基酸类型效应；7L1D/7RRG 的 ALHGGWTTK
个案叠加在队列曲线上。面板（a）回答模型能否预测规律，面板（b,c）只生成待
外部验证的机制假设。

```powershell
python paper\figures\make_itrace_residue_biology_figure.py
```

配套文件为 `fig_itrace_residue_biology_statistics.json`、
`fig_itrace_residue_biology_QA.md` 和 `paper/data/itrace_residue_biology/`。

## 6UON 与 2PTT 残基级结构案例图

`make_case_study_figures.py` 生成两张双栏机制图：`fig_6uon_case_study.{pdf,png,svg}` 将 ITrace MD、完整复合物预测、6ULK/6UON 肽构象叠合、CDR3 局部接触和固定坐标伙伴上下文响应放入同一证据链；`fig_2ptt_case_study.{pdf,png,svg}` 将 DynaRepo MD、2PTT/2PTU/2PTV 叠合、跨链电荷锚点和 CC-prime 环残基配对放入同一图版。分子结构由 `case_studies/render_3dmol.js` 调用本地归档的 3Dmol.js 无界面渲染，定量图与最终合成由 Matplotlib 完成。

```powershell
python paper\figures\make_case_study_figures.py
```

配套源数据为 `fig_6uon_case_study_source_data.csv` 与 `fig_2ptt_case_study_source_data.csv`，完整编号映射、结构对齐、局部距离和统计定义保存在 `fig_case_studies_statistics.json`，质量检查记录为 `fig_case_studies_QA.md`。2PTT 的原生 2B4 His42 在结构中存在但缺少 DynaRepo 残基映射，因此数值面板将其标记为缺失而不是零柔性。

## FlexPPI 设计/筛选 pipeline 总览图

`make_flexppi_design_pipeline_figure.py` 生成四面板双栏总览图 `fig_flexppi_design_pipeline.{pdf,svg,png}`（178 mm × 95 mm，300 dpi PNG）：(a) 从头设计与突变筛选双入口泳道汇合；(b) 多 seed apo/holo 结构集合与可开发性、折叠/复合物置信度、静态 ΔΔG 三道收窄 gate；(c) 紫色高亮的 FlexPPI 动力学感知重排层，输出逐残基 RMSF/Sχ/SASA 色带、partner-induced Δ 发散色带和不确定性灰带；(d) Pareto 示意散点 + MD/湿实验验证 + 沿底部回到 (a) 的朱红实验反馈弧。图形语法与文献画法依据见 `docs/pipeline_figure_literature_review.md`，完整设计见 `docs/flexppi_design_screening_pipeline.md`。

```powershell
python paper\figures\make_flexppi_design_pipeline_figure.py
```

配套 QA 记录为 `fig_flexppi_design_pipeline_QA.md`。该图是应用/展望示意图，不包含数据统计；Pareto 散点为示意点，无真实坐标刻度。
