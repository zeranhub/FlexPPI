# 蛋白质柔性与动态结构：背景开场图片索引

这组图片专门用于重写 FlexPPI 引言的最前部。正确顺序不是先讲“逐残基 RMSF 有什么用”，而是先建立：蛋白质在生理条件下占据具有不同概率的构象集合；构象之间的转换跨越多个时间和空间尺度；这种动态结构参与识别、催化和变构；静态坐标模型只描述集合的一部分。随后才能自然引出“完整构象集合很昂贵，因此需要可扩展的逐残基动态描述”。`contact_sheet.png` 是 12 张图的缩略总览。

需要区分三个概念。**结构**描述原子在某一构象中的空间组织；**柔性**描述某个区域能够访问多大范围的构象；**动态**还包含这些构象之间如何以及以多快速度转换。RMSF 只刻画位置波动幅度，是动态结构的一个边缘统计量，而不是完整能量景观、运动方向、耦合关系或转换速率。因此，FlexPPI 应定位成“从静态复合物估计一组可操作的逐残基动态统计量”，而不是“恢复完整动态结构”。

## 推荐优先级

| 优先级 | 图片 | 原文与用途 | 最适合放在哪里 | 使用边界 |
|---:|---|---|---|---|
| 1 | [henzlerwildman_kern_2007_fig1_energy_landscape.jpg](henzlerwildman_kern_2007_fig1_energy_landscape.jpg) | Henzler-Wildman & Kern, *Dynamic personalities of proteins*, Nature 2007，[原文](https://www.nature.com/articles/nature06522)。上半部分把 native basin 画成包含多个层级能垒的能量景观，下半部分连接局部运动、集体运动、时间尺度和观测方法。 | 引言第一段的首选参考图；可以学习其逻辑后重画 | Nature 版权图，只作内部调研，投稿时不应直接转载。 |
| 2 | [vandenbedem_fraser_2015_fig1_space_time_scales.jpg](vandenbedem_fraser_2015_fig1_space_time_scales.jpg) | van den Bedem & Fraser, *Integrative, Dynamic Structural Biology at Atomic Resolution—It’s About Time*, Nature Methods 2015，[原文](https://pmc.ncbi.nlm.nih.gov/articles/PMC4457290/)。同时展示蛋白运动的时间/长度尺度，以及 crystallography、NMR、MD 等方法的覆盖范围。 | 从“动态跨尺度”过渡到“任何单一技术都不能完整覆盖”的第二段 | 作者手稿中的 Nature 图，公开复用前核验版权。 |
| 3 | [fenwick_2011_fig1_structure_vs_ensemble.jpg](fenwick_2011_fig1_structure_vs_ensemble.jpg) | Fenwick et al., *Understanding biomolecular motion, recognition, and allostery by use of conformational ensembles*, European Biophysics Journal 2011，[原文](https://pmc.ncbi.nlm.nih.gov/articles/PMC3222826/)。将 X-ray 单结构、NMR 平均结构和不同时间尺度的 ensemble 并列。 | 具体说明“一个平均结构”与“构象分布”的信息差 | CC BY-NC；图片较横向、文字较小，适合参考而非直接缩放进双栏论文。 |
| 4 | [boehr_2009_fig3_ubiquitin_ensemble_recognition.jpg](boehr_2009_fig3_ubiquitin_ensemble_recognition.jpg) | Boehr, Nussinov & Wright, *The role of dynamic conformational ensembles in biomolecular recognition*, Nature Chemical Biology 2009，[原文](https://pmc.ncbi.nlm.nih.gov/articles/PMC2916928/)。游离泛素的 NMR ensemble 覆盖了与多个伙伴结合时观察到的构象，结合会稳定不同亚态。 | 把“构象集合”连接到 PPI 的最强实例之一 | Nature 版权内容/作者手稿，仅作内部参考。 |
| 5 | [nussinov_2023_fig2_energy_landscape_population_shift.jpg](nussinov_2023_fig2_energy_landscape_population_shift.jpg) | Nussinov et al., *Protein conformational ensembles in function: roles and mechanisms*, RSC Chemical Biology 2023，[原文](https://pmc.ncbi.nlm.nih.gov/articles/PMC10619138/)。从 folding funnel 放大到 native basin 的 active/inactive substates，并画出激活导致 population shift。 | 现代“sequence–ensemble–function”范式的概念图 | 许可在 PMC XML 中未明确列出；投稿复用前需向 RSC 页面核验。 |
| 6 | [henzlerwildman_kern_2007_fig7_hierarchy_catalysis.jpg](henzlerwildman_kern_2007_fig7_hierarchy_catalysis.jpg) | 同一 Nature 综述的 adenylate kinase 案例，综合 X-ray、MD、NMR 和 single-molecule FRET，说明不同方法观察同一催化循环的不同动态层次。 | 说明动态不是抽象口号，而是可被多种实验/模拟交叉观察 | Nature 版权图，只作内部参考。 |

## 其余图片

| 图片 | 主要信息 | 对 FlexPPI 的启发 |
|---|---|---|
| [boehr_2009_fig1_conformational_selection_induced_fit.jpg](boehr_2009_fig1_conformational_selection_induced_fit.jpg) | conformational selection 与 induced fit 位于同一热力学循环，并非必须二选一。 | PPI 背景中可写“结合改变构象概率分布并伴随局部优化”，避免把所有变化简单归为诱导契合。 |
| [fenwick_2011_fig6_selection_vs_induced_fit.jpg](fenwick_2011_fig6_selection_vs_induced_fit.jpg) | 用集合重叠解释 conformational selection/induced fit。 | 可作为更直观但较简化的教学图。 |
| [fenwick_2011_fig2_motion_timescales_methods.jpg](fenwick_2011_fig2_motion_timescales_methods.jpg) | 生物运动与实验/计算技术覆盖范围。 | 与 van den Bedem Figure 1 内容相近；分辨率较低，只作备选。 |
| [motlagh_2014_fig2_dynamic_continuum_allostery.jpg](motlagh_2014_fig2_dynamic_continuum_allostery.jpg) | Motlagh et al., *The ensemble nature of allostery*, Nature 2014，[原文](https://pmc.ncbi.nlm.nih.gov/articles/PMC4224315/)。变构可从刚体运动、侧链/主链动态延伸到局部解折叠和无序。 | 说明只预测主链 RMSF 会遗漏侧链和无序层面的响应，为 $H_\chi$ 等多描述符铺垫。 |
| [bahar_2010_fig1_equilibrium_motions_scales.jpg](bahar_2010_fig1_equilibrium_motions_scales.jpg) | Bahar et al., *Normal Mode Analysis of Biomolecular Structures*, Chemical Reviews 2010，[原文](https://pmc.ncbi.nlm.nih.gov/articles/PMC2836427/)。运动幅度/长度与时间尺度的层级图。 | 可用于解释不同方法为什么只能覆盖动态空间的一部分。 |
| [bahar_2010_fig2_hierarchy_energy_substates.jpg](bahar_2010_fig2_hierarchy_energy_substates.jpg) | native state 在不同分辨率下可被视为单一 minimum、多个 substates 或大量 microstates。 | 最适合解释 FlexPPI 的逐残基统计是对微观集合的压缩，而非完整重建。 |

## 建议采用的背景图叙事

如果论文只允许在引言放一张 Figure 1，不要直接转载上述综述图，而应把它们的思想压缩成 FlexPPI 的真实 running example。左侧先画一个蛋白复合物的单一坐标，明确标成“representative coordinate model”；中间画同一体系的若干淡色构象和能量/概率分布，说明真实对象是动态集合；右侧再把集合压缩成 RMSF、$H_\chi$ 和 SASA 三条逐残基 profile，并提出问题：能否仅从静态多链结构快速恢复这些可操作统计量？这样 Figure 1 的逻辑是“static model → inaccessible ensemble → residue-wise dynamic summary”，模型架构应留在后续方法图。

## 下载与核验

PMC 图片通过原文页面提供的 NCBI 图片 CDN 下载；Henzler-Wildman & Kern 两图通过 Nature 正文公开的 Springer Nature 图片端点下载。所有图片均已本地解码，且汇总图经过视觉检查。原图用于研究、比较和重画参考；是否允许公开转载取决于原论文许可，不能因为文件可下载就默认可用于投稿。
