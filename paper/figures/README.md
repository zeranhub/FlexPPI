# Figures

将论文图片放在本目录中，优先使用 PDF、PNG 或 JPG 格式。

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

## ITrace 冻结检查点上下文消融图

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
