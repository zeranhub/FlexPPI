# 蛋白质柔性“包络图”与结构云参考

本目录收集能把蛋白质动态压缩到一张静态三维图中的代表性视觉方案。最接近“结构外面包着一层动态包络”的是 `kozlikova_2017_fig12_variable_tube_and_flexibility_isosurface.png`：左侧以主链管径表达局部柔性，右侧以三维柔性等值面表达不同区域所覆盖的空间。所有文件均已下载并通过图像解码检查；这些图片首先用于内部调研和构图参考，公开转载前仍应按原文许可逐项核对署名和第三方素材说明。

| 优先级 | 本地图片 | 文献或官方来源 | 图在表达什么 | 对 FlexPPI 的借鉴与限制 |
|---:|---|---|---|---|
| 1 | [kozlikova_2017_fig12_variable_tube_and_flexibility_isosurface.png](kozlikova_2017_fig12_variable_tube_and_flexibility_isosurface.png) | Kozlíková et al., *Visualization of Biomolecular Structures: State of the Art Revisited*, CGF 2017, [DOI/原文](https://doi.org/10.1111/cgf.13072), Figure 12 | 同一离子通道以“可变半径骨架管”和“柔性等值面”两种方式呈现；红色柔性区占据更大空间，蓝色刚性区更集中 | 与所问“包络图”最一致，也与逐残基预测最匹配。该图是视觉综述中的示例，不是 FlexPPI 的生物学验证结果 |
| 2 | [lee_varshney_2002_crambin_fuzzy_surface.gif](lee_varshney_2002_crambin_fuzzy_surface.gif) | Lee & Varshney, *Representing Thermal Vibrations and Uncertainty in Molecular Surfaces*, VDA 2002, [项目页](https://www.cs.umd.edu/gvil/projects/probabilistic_surfaces.shtml) | 用多层半透明表面表达热振动或原子位置不确定性，左为普通表面，右为 fuzzy surface | 非常适合借鉴“置信层/概率层”的画法；原网页标注 University of Maryland 版权所有，仅作为内部参考 |
| 2 | [lee_varshney_2002_hiv_protease_fuzzy_surface.gif](lee_varshney_2002_hiv_protease_fuzzy_surface.gif) | 同上 | HIV protease 的普通表面与多层模糊表面对照 | 适合说明表面不应被当成确定边界；不能把纯逐残基标量画成有方向的真实运动 |
| 2 | [lee_varshney_2002_surface_representation.gif](lee_varshney_2002_surface_representation.gif) | 同上 | 从原子/空间填充表示过渡到平滑分子表面 | 适合解释“结构表面”和“动态包络”不是一回事 |
| 3 | [rico_2013_fig3_bfactor_and_nmr_ensemble.jpg](rico_2013_fig3_bfactor_and_nmr_ensemble.jpg) | Rico et al., *Mechanics of proteins with a focus on atomic force microscopy*, J Nanobiotechnology 2013, [原文](https://doi.org/10.1186/1477-3155-11-S1-S3), Figure 3 | B-factor 着色结构、24 个 NMR 模型叠合以及其他柔性测量并列 | 很适合支持“静态坐标看似刚性，但实验给出位置分布”的背景；B-factor 混合热运动、静态无序和晶体堆积，不能直接等同 RMSF |
| 4 | [lemmin_2017_fig1_glycan_sampling_volumes.jpg](lemmin_2017_fig1_glycan_sampling_volumes.jpg) | Lemmin et al., *Conformational Heterogeneity of the HIV Envelope Glycan Shield*, Scientific Reports 2017, [原文](https://doi.org/10.1038/s41598-017-04532-9), Figure 1 | 以三维采样体积包住代表结构，表示糖链构象所覆盖的空间 | 是很直观的“动态占据体积”实例，但对象主要是糖链柔性，不宜当成一般蛋白主链 RMSF 的直接证据 |
| 5 | [fuxreiter_fuzzy_protein_interaction_ensemble.jpg](fuxreiter_fuzzy_protein_interaction_ensemble.jpg) | Fuxreiter Lab, [Fuzzy Protein Interaction Ensembles](https://fuxreiterlab.github.io/research_fuzzy.html) | 大量构象叠合形成由无序到较有序的结构云 | 适合作为构象集合和 fuzzy binding 的概念图；偏向内在无序蛋白，网页未见明确再利用许可 |
| 6 | [janowski_2014_fig5_md_bundle_heterogeneity.jpg](janowski_2014_fig5_md_bundle_heterogeneity.jpg) | Janowski et al., *X-ray refinement significantly underestimates the level of microscopic heterogeneity in biomolecular crystals*, Nature Communications 2014, [原文](https://doi.org/10.1038/ncomms4220), Figure 5 | 电子密度、各向异性 B-factor 椭球和 100 个 MD 构象束的局部对照 | 说明同一平均结构下面仍有局部构象分布；不是整蛋白包络图，适合作为“为什么单点坐标不够”的证据 |
| 7 | [pymol_bfactor_putty_example.png](pymol_bfactor_putty_example.png) | [PyMOL Cartoon/Putty 文档](https://wiki.pymol.org/index.php/Cartoon) | 管径与颜色同时随 B-factor 改变 | 是最容易复现的绘图语法，但属于软件示例而非生物学证据；FlexPPI 可把 B-factor 字段替换成模型预测值 |

## 推荐给 FlexPPI 的图形定义

最稳妥的主图可以把目标蛋白画成灰白色或浅青色 ribbon，在其外部绘制一个半透明的预测柔性包络，同时把逐残基预测值编码为 backbone tube 的半径或颜色。包络必须在图注里给出可计算定义：如果来自 MD/NMR 构象集合，它可以是经过结构对齐后由原子占据概率产生的 50%、75% 或 90% 等值面；如果 FlexPPI 当前只有每个残基的一个标量值，那么更诚实的叫法是“predicted flexibility envelope / schematic uncertainty field”，并采用围绕 Cα 轨迹的各向同性粗细或球形扩张。单个 RMSF 数值只说明平均偏离大小，不给出运动方向，所以没有协方差、主模态或构象样本时，不应画带方向的叶瓣，也不应声称重建了真实构象集合。

## 随图保存的原文

- [kozlikova_2017_visualization_biomolecular_structures.pdf](../source_pdfs/kozlikova_2017_visualization_biomolecular_structures.pdf)
- [lee_varshney_2002_fuzzy_molecular_surfaces.pdf](../source_pdfs/lee_varshney_2002_fuzzy_molecular_surfaces.pdf)

## 下载和版权说明

Kozlíková Figure 12 从公开全文 PDF 的对应页面中提取；Rico、Lemmin 与 Janowski 图片来自出版社的 Figure 图像端点；Lee–Varshney、Fuxreiter 和 PyMOL 图片来自其项目或文档页面。下载到本地不自动获得公开转载权。论文正式作图时，最佳做法是保留这些文件作为视觉先例，然后用自己的 PDB、MD 数据或 FlexPPI 输出重画，并在方法或图注中引用视觉/计算定义的来源。
