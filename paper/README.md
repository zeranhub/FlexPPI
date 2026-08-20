# FlexPPI 论文目录

- `paper.tex`：论文主文件，负责组织章节。
- `preamble.tex`：宏包和统一排版设置。
- `sections/`：日常撰写的正文文件。
- `tables/`：复杂表格文件。
- `figures/`：论文图片。
- `data/itrace/`：论文中 ITrace 基准表格对应的精简原始 CSV。
- `data/itrace_biology/`：结合态界面与 CDR 链间不对称图的逐复合体源数据。
- `data/itrace_residue_biology/`：完整模型--MD 残基预测一致性、peptide 位置、
  重原子接触、$H_\chi$/SASA 上下文响应和 7L1D/7RRG 残基个案的源数据。
- `references.bib`：参考文献数据库。

ITrace 正文生物学主图、留存的算法基准图和标签敏感性图分别为
`figures/fig_itrace_biology.pdf`、`figures/fig_itrace_main.pdf` 与
`figures/fig_itrace_sensitivity.pdf`。生物学主图直接展示 pMHC/TCR 界面
动力学不对称、TCR $\alpha/\beta$ CDR 环差异，以及完整上下文在界面排序上
反而退化的现象；执行 `python paper\figures\make_itrace_biology_figure.py`
可重建 PDF、SVG、300 dpi PNG、统计 JSON、QA 报告及
`data/itrace_biology/` 的逐复合体源数据。实验流程直接写在正文，避免用协议表
重复叙述；算法基准图不再占用正文页面，总体指标与细分数值表移至附录，正文
优先解释生物学模式。

peptide 残基级图为 `figures/fig_itrace_residue_biology.pdf`，对应可编辑 SVG、
300 dpi PNG、统计 JSON 和 QA 报告均采用同名文件。该图首先直接检验完整复合物
模型能否恢复 154 个九肽的 MD 位置柔性，再展示固定坐标完整 TCR--pMHC 减
pMHC 的 RMSF/$H_\chi$/SASA 响应、TCR 重原子接触比例、同肽 7L1D/7RRG
个案，以及去除位置均值后的氨基酸类型假设。执行
`python paper\figures\make_itrace_residue_biology_figure.py` 可由
`paper/data/itrace_residue_biology/` 重建全部图件；全队列预测和接触源数据由
`scripts/analyze_itrace_residue_biology.py` 从冻结检查点生成。

固定坐标伙伴消融的正文主图为
`figures/fig_itrace_partner_ablation.pdf`；主要比较完整五链与 pMHC-only/TCR-only，
不拆散目标侧内部组成。`scripts/analyze_fixed_partner_ablation.py` 从冻结预测重建逐链
和逐复合体 RMSF 准确性、接触/距离定位、FDR 检验、488 组逐原子坐标一致性检查、
正文表格及 PDF/SVG/PNG 图件，精简源数据保存在
`data/itrace_partner_ablation/`。该比较只定义模型上下文敏感性，不表示 apo/holo
变化。

进一步拆成单链的敏感性消融图为
`figures/fig_itrace_context_ablation.pdf`，300 dpi 预览和可编辑矢量版本分别为
同名 `.png` 与 `.svg`。执行
`python paper\figures\make_itrace_context_ablation_figure.py` 可从
`output/itrace_zero_shot_context/` 重新生成图、统计 JSON 和
`data/itrace_context/` 下的全部绘图源数据。对应正文表格为
`tables/itrace_context_ablation.tex` 和
`tables/itrace_context_sensitivity.tex`。

在 VS Code 中打开任意 `sections/*.tex` 文件都可以写正文；这些文件顶部的
`% !TeX root = ../paper.tex` 会让 LaTeX Workshop 编译完整论文。

编译时优先选择 `latexmk (xelatex)`，然后使用 LaTeX Workshop 的 PDF 查看命令
预览。当前本机的 TeX Live 2026 `latexmk` 启动包装器存在 Lua 路径错误；若遇到
同样错误，可在本目录依次执行下列命令完成正文、参考文献、公式编号和交叉引用编译：

```powershell
xelatex -interaction=nonstopmode -halt-on-error paper.tex
bibtex paper
xelatex -interaction=nonstopmode -halt-on-error paper.tex
xelatex -interaction=nonstopmode -halt-on-error paper.tex
```
