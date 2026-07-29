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
