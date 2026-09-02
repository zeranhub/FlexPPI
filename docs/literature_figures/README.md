# FlexPPI 可参考文献图片索引

本目录保存上一轮 FlexPPI 背景、引言和图片叙事调研中实际建议参考的图片。`figures/` 中是可直接浏览的原图或从官方论文 PDF 页面裁出的 Figure；`source_pdfs/` 保存能够从官方会议来源稳定取得的完整论文，便于核对图注和上下文；`contact_sheet.png` 是十张图的缩略总览。PDF 提取图只做了边界裁剪，没有重绘或修改原始内容。

> 这些图片用于论文调研和内部设计参考。若要在论文、幻灯片或公开材料中直接转载，需要根据原论文的开放许可核对署名、图注和再利用条件；更推荐学习其叙事结构后自行重绘 FlexPPI 图，而不是直接拼贴原图。

## 快速入口

- [十张图总览](contact_sheet.png)
- [图片目录](figures/)
- [官方原文 PDF](source_pdfs/)
- [“逐残基柔性有什么用”证据图总览](utility_figures/contact_sheet.png)
- [必要性/下游用途图片索引](utility_figures/README.md)
- [蛋白质柔性与动态结构背景图总览](dynamic_structure_background/contact_sheet.png)
- [动态结构开场图片索引](dynamic_structure_background/README.md)
- [柔性包络、等值面与结构云图片目录](flexibility_envelopes/)
- [柔性包络图片索引与 FlexPPI 构图建议](flexibility_envelopes/README.md)

## 图片清单与借鉴重点

| 本地文件 | 原论文与官方来源 | 主要内容 | FlexPPI 最值得借鉴 | 使用边界 |
|---|---|---|---|---|
| [alphaflow_2024_fig1.png](figures/alphaflow_2024_fig1.png) | AlphaFlow, ICML 2024，[PMLR 页面](https://proceedings.mlr.press/v235/jing24a.html) | 从先验构象到最终构象样本的单一概念流程 | 同一个蛋白保持一致配色；用少量面板讲清一个核心机制 | FlexPPI 不生成构象集合，不应借用 flow field、去噪或采样语义 |
| [backflip_flips_2025_fig1.png](figures/backflip_flips_2025_fig1.png) | BackFlip/FliPS, ICML 2025，[PMLR 页面](https://proceedings.mlr.press/v267/viliuga25a.html) | 目标柔性 profile、结构生成、柔性预测与 MD 验证闭环 | “结构 + residue profile + validation”组合非常适合未来应用图 | 目前若没有完整设计闭环，不能把设计重排序画成已验证主贡献 |
| [flexpert_2025_fig1.png](figures/flexpert_2025_fig1.png) | Flexpert, ICLR 2025，[OpenReview 页面](https://openreview.net/forum?id=L238BAx0wP) | Flexpert-Seq 嵌套在 Flexpert-3D 中，区分固定与训练模块 | 共享、嵌套和冻结/训练模块的一眼可读表达 | 适合方法图，不适合承担 FlexPPI 的生物学动机 |
| [rmsfnet_2024_fig1.png](figures/rmsfnet_2024_fig1.png) | RMSF-net, Nature Communications 2024，[期刊页面](https://www.nature.com/articles/s41467-024-49858-x) | 数据准备、MD 监督标签、网络和推理后处理 | 清楚分离“训练期才需要的 MD 标签支路”和“静态输入推理支路” | 图很密且是 cryo-EM 体素任务，不能照搬其体素和密度图细节 |
| [pesto_2023_fig1.png](figures/pesto_2023_fig1.png) | PeSTo, Nature Communications 2023，[期刊页面](https://www.nature.com/articles/s41467-023-37701-8) | 局部原子消息、多尺度几何 Transformer 与结构着色输出 | 几何架构和真实结构输出并列；从原子邻域自然汇聚到残基输出 | PeSTo 预测界面概率，不预测动态状态或 apo/holo 响应 |
| [dynaprot_2026_fig1.png](figures/dynaprot_2026_fig1.png) | DynaProt, ICLR 2026，[arXiv 页面](https://arxiv.org/abs/2509.01038) | 动态方法的表达能力–效率概念图 | 可用二维任务空间快速说明方法类别差异 | 坐标必须有可操作定义，不能主观把 FlexPPI 放到 Pareto 最优位置 |
| [atlas_2024_fig3.jpeg](figures/atlas_2024_fig3.jpeg) | ATLAS, Nucleic Acids Research 2024，[期刊页面](https://academic.oup.com/nar/article/52/D1/D384/7438909) | 结构、重复 MD、RMSF/Neq profile、接触图和结构着色 | 最接近 FlexPPI “轨迹标签—逐残基 profile—结构着色”的数据视觉语法 | 该图是数据库页面示例，不是机器学习模型架构 |
| [dynarepo_2026_graphical_abstract.jpeg](figures/dynarepo_2026_graphical_abstract.jpeg) | DynaRepo, Nucleic Acids Research 2026，[期刊页面](https://academic.oup.com/nar/article/54/D1/D393/8322558) | 复合物筛选、三次 500 ns MD、分析与数据库的完整数据链 | 适合说明复杂物数据来源、重复模拟和离线标签构建 | 不要把数据库生成流程画成 FlexPPI 推理时所需步骤 |
| [dynarepo_2026_fig2.jpeg](figures/dynarepo_2026_fig2.jpeg) | DynaRepo, Nucleic Acids Research 2026，[期刊页面](https://academic.oup.com/nar/article/54/D1/D393/8322558) | RMSD、SASA、fluctuation、聚类、电势与口袋等多种分析 | 支持 RMSF、SASA 等不同输出是互补动态维度，而非一个总分 | 信息非常密，适合数据/补充图，不适合 Figure 1 动机图 |
| [seqdance_2025_fig1.jpg](figures/seqdance_2025_fig1.jpg) | SeqDance/ESMDance，[PMC 页面](https://pmc.ncbi.nlm.nih.gov/articles/PMC11507661/) | sequence → structure ensemble → function → conservation，以及动态属性预训练 | 可借“表征信息来源”和“residue-level / pairwise outputs”的分层方式 | FlexPPI 的中心输入是多链静态结构，不应变成 sequence-only 叙事 |

## 建议的参考优先级

若只看四张，优先看 AlphaFlow Figure 1、BackFlip/FliPS Figure 1、RMSF-net Figure 1 和 PeSTo Figure 1：它们分别解决“概念故事如何简化”“结构与 profile 如何组合”“MD 标签与推理如何分开”“几何网络与残基输出如何并列”。Flexpert 和 SeqDance 更适合重画方法总览时参考；ATLAS 与 DynaRepo 更适合数据和标签构建图；DynaProt 只适合作为任务空间构图灵感。

## 对应的 FlexPPI 图序建议

1. **Figure 1 - Motivated example.** Apo/holo 或有/无伙伴的同一目标结构，旁边放 RMSF、$H_\chi$ 和 SASA profile，再比较 partner-agnostic predictor 与 FlexPPI；可借 AlphaFlow 的单一主叙事和 BackFlip 的 structure/profile/validation 组合。
2. **Figure 2 - Data and absolute-state pretraining.** MD 轨迹只在训练期生成标签，静态结构在推理期直接进入模型；可借 RMSF-net、ATLAS 和 DynaRepo 的训练/数据视觉语法。
3. **Figure 3 - Partner-conditioned response.** Apo 与 holo 两路使用 shared encoder，输出绝对状态、状态差与 uncertainty；可借 Flexpert 的共享/嵌套表达和 PeSTo 的多链几何表达。
4. **Results figures.** Overall、interface、peptide、partner 分区，并同时展示成功案例与 over-rigidification failure；不要用 DynaProt 式主观 Pareto 图代替真实测量。

## 原文 PDF

| 文件 | 来源 |
|---|---|
| [alphaflow_icml2024.pdf](source_pdfs/alphaflow_icml2024.pdf) | PMLR 官方论文 PDF |
| [backflip_flips_icml2025.pdf](source_pdfs/backflip_flips_icml2025.pdf) | PMLR 官方论文 PDF |
| [flexpert_iclr2025.pdf](source_pdfs/flexpert_iclr2025.pdf) | ICLR 官方会议 PDF |
| [dynaprot_iclr2026.pdf](source_pdfs/dynaprot_iclr2026.pdf) | 作者公开论文 PDF；图意同时用 arXiv 页面核验 |
| [tzeng_kalodimos_2012_nature_conformational_entropy.pdf](source_pdfs/tzeng_kalodimos_2012_nature_conformational_entropy.pdf) | 作者公开的 Nature 论文 PDF；用于核验“构象熵可增强/抑制结合”的 Figure 2/4 |
| [kozlikova_2017_visualization_biomolecular_structures.pdf](source_pdfs/kozlikova_2017_visualization_biomolecular_structures.pdf) | 公开全文；Figure 12 给出可变管径主链和柔性等值面 |
| [lee_varshney_2002_fuzzy_molecular_surfaces.pdf](source_pdfs/lee_varshney_2002_fuzzy_molecular_surfaces.pdf) | 作者机构公开论文；用多层透明表面表达热振动和位置不确定性 |

## 下载与核验说明

RMSF-net 和 PeSTo 使用期刊页面提供的高分辨率 Figure 1；ATLAS 与 DynaRepo 使用 Oxford Academic 页面指向的原始图片；SeqDance 使用 PMC XML 标注的原始 Figure 1 CDN 文件；AlphaFlow、BackFlip/FliPS、Flexpert 和 DynaProt 从本地保存的官方论文 PDF 对应页面裁出。所有图片均已用图像解码器打开检查，四个 PDF 也通过 PDF 元数据和目标页面渲染检查。
