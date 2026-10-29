# 最终论文六张图

本目录逐张对应最终稿 `paper_v58_20261022_中文修订_图2至5配色统一_三线表.docx` 中实际显示的图1、图2、图3、图4、图5与图S1。图1采用v56确认的框架，图2至图5采用v58，图S1沿用v39归档的相同像素版本。没有采用后来未进入这份Word的v57替代图。现有最终PNG、PDF、SVG和Draw.io均保持原文件字节；`manifest.json`记录来源、Word嵌入图片的SHA256、可编辑文件和重绘入口。

| 目录 | 内容与可编辑入口 | 数据重绘入口 |
|---|---|---|
| `01_模型框架` | `fig1_framework.drawio`；蛋白图片内嵌，文字、模块、箭头与曲线可在Draw.io编辑 | Draw.io桌面版打开源文件后导出，或统一命令选择`1` |
| `02_数值预测与高动态位点识别` | `fig2_dense_capability.svg`；`source_archive/`包含冻结CSV与NPZ | `plot.py`，读取本目录`palette.json` |
| `03_KIX与泛素局部动态验证` | `fig3_kix_ubiquitin.svg`；`source_archive/kix/`、`ubiquitin/`保存PSE、PML、真实PDB与MD代表帧 | 根目录`compose_figure.py`，读取包内表和结构面板 |
| `04_MEL5跨界面接触与双动态响应` | `fig4_mel5_compact_joint.svg`；`structure/MEL5_affinity_maturation.pse` | `plot.py`，读取本地`source_*.csv`；`structure/render_structures.py`可重渲染真实PDB |
| `05_联合读数与构象占据` | `fig5_color_consistent_joint_readouts_conformations.svg` | `plot.py`，读取`source/`中逐记录与逐角度数据 |
| `S1_标签均值与观测标准差` | `figS1_observed_label_variation.svg`；本地CSV保留标签来源分层 | `plot_label_variation.py`；最终图注与表述已并入源码 |

在项目根目录执行`python paper/figures_v58/reproduce.py --check`，可只核查最终图件的完整性。执行`python paper/figures_v58/reproduce.py 3`即可从本地数据重新合成图3；将`3`替换为`2`、`4`、`5`或`S1`可重绘对应图。结果统一写入`paper/figures_v58/reproduced/`；`--output-dir`可指定其他位置。统一入口使用临时副本运行绘图代码，保留这里已确认的最终文件。图3也支持直接运行`python paper/figures_v58/03_KIX与泛素局部动态验证/compose_figure.py`，默认输出到该图目录的`reproduced/`。常规图形重绘仅需NumPy、Pandas、Matplotlib、SciPy、Pillow与PyMuPDF；分子结构面板还可用已保存的PyMOL会话调整，重新渲染时需要PyMOL。`requirements.txt`列出常规Python依赖。

图3和图4的SVG中，结构面板是原始渲染图片，坐标轴、曲线及文字保留矢量；编辑蛋白的视角、表示或原子选择请打开对应PSE或运行结构源代码。图S1的`active_v38_image7.png`与Word内图片字节完全一致，与同目录标准命名PNG逐像素一致，两个文件的压缩编码不同，因此保留前者作为Word嵌入图的核对依据。S1的逐残基大表不参与绘图，绘图程序仅校验并读取随包的三份汇总CSV；完整数据的原始校验信息保存在`label_analysis_audit.json`。

各图只需本目录内源文件即可重新合成，不读取被移走的旧论文目录。历史实验与原始归档另外位于项目`experiments/_sources/`，其路径仅作为来源记录，不是绘图运行依赖。绘图软件与字体版本可能影响重新导出的压缩编码或文字细节；需要直接延续已接受排版时编辑现有SVG或Draw.io即可。
