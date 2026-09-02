# “逐残基柔性有什么用”图片索引

本目录不是普通的模型架构图合集，而是专门为回应“预测逐残基柔性没有实际用途”这一质疑整理的证据图。图片按一条可检验的因果链选择：静态结构或整体 RMSD 会遗漏局部动态差异；这些差异在界面、侧链转子、暴露度和变构位点上具有空间定位意义；一张逐残基动态图因此能够改变突变、口袋筛选、柔性对接、结构生成或后续 MD 的决策。`contact_sheet.png` 是全部图片的缩略总览；原图仅供调研、组图设计和内部讨论，公开转载前必须按下表核对许可证。最稳妥的论文做法是引用原文并自行重画 FlexPPI 的真实案例，而不是拼贴这些图片。

## 最值得先看的六组图

| 本地图片 | 原文 | 它能支持的必要性论点 | 最适合借鉴的位置 | 许可/边界 |
|---|---|---|---|---|
| [tzeng_kalodimos_2012_fig2_conformational_entropy_enhances_binding.png](tzeng_kalodimos_2012_fig2_conformational_entropy_enhances_binding.png)；[fig4](tzeng_kalodimos_2012_fig4_entropy_dominates_binding.png) | Tzeng & Kalodimos, *Protein activity regulation by conformational entropy*, Nature 2012，[DOI](https://doi.org/10.1038/nature11271) | 结构上处于相似活性状态并不足以解释亲和力；结合引起的内部动态/构象熵变化可以增强或抑制结合。逐残基着色把“动态发生在哪里”直接映射回结构。 | 引言第一段的反直觉例子；Figure 1 的视觉原型 | Nature 版权内容，仅作内部参考；FlexPPI 的 $H_\chi$ 是转子占据的 Shannon 熵，不能直接写成热力学 $\Delta S$ 或 $\Delta G$。 |
| [clark_2019_fig1_backbone_change.jpg](clark_2019_fig1_backbone_change.jpg)；[fig3](clark_2019_fig3_sidechain_rotamers.jpg)；[fig7](clark_2019_fig7_sasa.jpg) | Clark et al., *Inherent versus induced protein flexibility*, PLOS Computational Biology 2019，[原文](https://pmc.ncbi.nlm.nih.gov/articles/PMC6370239/)，[DOI](https://doi.org/10.1371/journal.pcbi.1006705) | 305 个蛋白家族中 apo–holo 的额外骨架位移通常小于 0.5 Å，但结合位点侧链 $\chi_1$ 可进入 apo 结构未覆盖的取向；所以“整体结构几乎不变”不等于“逐残基状态没变”。 | 背景第二段；用 Fig. 1 与 Fig. 3 成对讲最有效 | CC BY；可在规范署名后转载或改编。 |
| [cukuroglu_2012_fig1_interface_flexibility.jpg](cukuroglu_2012_fig1_interface_flexibility.jpg)；[fig2](cukuroglu_2012_fig2_profilin_actin.jpg) | Cukuroglu et al., *Roles of residues in the interface of transient protein-protein complexes before complexation*, Scientific Reports 2012，[原文](https://pmc.ncbi.nlm.nih.gov/articles/PMC3312204/)，[DOI](https://doi.org/10.1038/srep00334) | 界面不是一个统一刚性的区域：约 65% 的核心界面残基近似构象不变，外围和非刚性界面残基表现不同。逐残基定位才可能区分锚点、核心与可适应外围。 | PPI 特异的必要性段；Figure 1 的“同一界面、不同残基角色” | CC BY-NC-SA 3.0；投稿图中改编前应确认期刊的商业属性要求。 |
| [dobbins_2010_fig4_bound_unbound_mobility.jpg](dobbins_2010_fig4_bound_unbound_mobility.jpg) | Dobbins et al., *Comparing interfacial dynamics in protein-protein complexes: an elastic network approach*, BMC Structural Biology 2010，[原文](https://pmc.ncbi.nlm.nih.gov/articles/PMC2927602/)，[DOI](https://doi.org/10.1186/1472-6807-10-26) | 22 个二聚体中，不同界面类型的结合前后 mobility 关系不同：结合并不是对所有界面统一“压硬”。因此伙伴条件化响应不能由单链或一个全局柔性分数代替。 | 定义 FlexPPI 的 partner-conditioned setting | CC BY 2.0。 |
| [panjkovich_daura_2012_fig4_allosteric_pockets.jpg](panjkovich_daura_2012_fig4_allosteric_pockets.jpg)；[fig1](panjkovich_daura_2012_fig1_method.jpg) | Panjkovich & Daura, *Exploiting protein flexibility to predict the location of allosteric sites*, BMC Bioinformatics 2012，[原文](https://pmc.ncbi.nlm.nih.gov/articles/PMC3562710/)，[DOI](https://doi.org/10.1186/1471-2105-13-273) | 在 91 个非冗余变构蛋白中，70% 在变构配体结合后出现显著柔性变化；柔性变化与保守性结合后，可在候选口袋中定位变构位点，严格设置的阳性预测值为 65%。这是“动态图改变口袋选择”的直接例子。 | 引言中“预测后能做什么”的第一项下游任务 | CC BY 2.0。 |
| [seqdance_esmdance_2025_fig4_mutation_effect.jpg](seqdance_esmdance_2025_fig4_mutation_effect.jpg) | Hou et al., *Learning Biophysical Dynamics with Protein Language Models* 预印本版本，[PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC11507661/)；正式版为 PNAS 2026，[PubMed](https://pubmed.ncbi.nlm.nih.gov/41576089/) | 对野生型和突变体预测动态属性并取相对变化，可以零样本预测 412 个蛋白的稳定性突变效应；正式版 ESMDance 将动态与进化信息结合后优于相同规模 ESM2。它把“预测动态”明确连到“预测突变后果”。 | 引言中“动态表征是下游变量，而非终点” | 下载的是 CC BY-NC 4.0 预印本 Figure 4；正式版图和数值应以 PNAS 版本复核。 |

## 与 FlexPPI 输出最直接相关的图

| 本地图片 | 说明 | 对 FlexPPI 的具体启发 |
|---|---|---|
| [cole_warwicker_2002_fig1_sidechain_entropy_interface.jpg](cole_warwicker_2002_fig1_sidechain_entropy_interface.jpg)；[fig6](cole_warwicker_2002_fig6_residue_entropy_contribution.jpg) | Cole & Warwicker, *Side-chain conformational entropy at protein–protein interfaces*, Protein Science 2002，[原文](https://pmc.ncbi.nlm.nih.gov/articles/PMC2373749/)，[DOI](https://doi.org/10.1110/ps.0222702)。界面侧链熵图和单残基贡献图。 | 直接说明为什么 $H_\chi$ 需要逐残基输出：界面低熵区域可辅助识别界面，论文报告真界面在同源二聚体/大型异源亚基表面 patch 中排名第一的比例分别为 68%/66%。注意计算定义与 FlexPPI 标签并不相同，不能把结果直接移植。 |
| [flexpert_2025_fig2_design.png](flexpert_2025_fig2_design.png) | Kouba et al., *Flexpert: Learning to Predict and Control Protein Flexibility*, ICLR 2025，[OpenReview](https://openreview.net/forum?id=L238BAx0wP)。 | 展示逐位置柔性指令如何进入 inverse folding。论文的增柔实验中，Flexpert-Design 的中位 enrichment ratio 为 1.52、增柔突变比例为 0.83，基线 ProteinMPNN 为 1.07/0.61；但柔性主要由预测器评估，不能当成湿实验验证。 |
| [backflip_flips_2025_fig1_target_profile_to_backbone.png](backflip_flips_2025_fig1_target_profile_to_backbone.png) | Viliuga et al., BackFlip/FliPS, ICML 2025，[PMLR](https://proceedings.mlr.press/v267/viliuga25a.html)。 | 最清楚的“逐残基 target profile → 候选骨架 → 再预测 → MD 验证”闭环。论文报告三种长度设置中生成结构与目标 profile 的 MD 相关约为 0.70、0.78 和 0.79，说明 residue-wise flexibility 可以成为设计约束。 |
| [petrovic_kamerlin_2018_fig1_dynamics_evolvability.jpg](petrovic_kamerlin_2018_fig1_dynamics_evolvability.jpg)；[fig6](petrovic_kamerlin_2018_fig6_arylesterase_evolution.jpg)；[fig8](petrovic_kamerlin_2018_fig8_flexibility_enables_design.jpg) | Petrović & Kamerlin, *Conformational dynamics and enzyme evolution*, Interface Focus 2018，[原文](https://pmc.ncbi.nlm.nih.gov/articles/PMC6073641/)，[DOI](https://doi.org/10.1098/rsif.2018.0330)。 | Fig. 1 适合学习“主态/次态—新功能”的概念表达；Fig. 6/8 用真实进化和设计案例说明远端突变、柔性及可设计性之间的联系。它们是背景证据，不应被写成 FlexPPI 已经完成酶设计。 |
| [flexppi_6uon_running_example.png](flexppi_6uon_running_example.png) | 本项目已生成的 6UON 案例图。 | 最适合取代文献拼图成为 FlexPPI Figure 1：固定结合态坐标下的 full-minus-pMHC 只能称为“伙伴上下文敏感性”，而 ITrace 三条轨迹均值给出 MD 目标；不要把固定坐标的 partner deletion 写成真实 apo/holo 动态变化。 |

## 推荐阅读/看图顺序

先看 Clark Fig. 1+3，建立“整体看起来没变，但局部侧链状态变了”；再看 Cukuroglu Fig. 1 与 Dobbins Fig. 4，建立“界面不同残基、不同复合物的响应不一样”；随后看 Tzeng Fig. 2/4 与 Panjkovich Fig. 4，建立“局部动态会影响结合热力学并能帮助定位变构口袋”；最后看 SeqDance Fig. 4、Flexpert Fig. 2 和 BackFlip/FliPS Fig. 1，建立“逐残基动态可以进入突变评估和设计闭环”。真正用于 FlexPPI Figure 1 时，建议沿这条逻辑重画 6UON，而不是把上述图全部塞入正文。

## 下载与视觉核验

Panjkovich、Cukuroglu、Dobbins、Clark、Cole、SeqDance 和 Petrović/Kamerlin 的图片均从 PMC/NCBI 的原始图片 CDN 下载；Tzeng/Kalodimos 两图从作者公开的 Nature PDF 以 160 dpi 渲染后仅做矩形裁切；Flexpert Figure 2 从官方会议 PDF 中提取其嵌入原图；BackFlip/FliPS Figure 1 来自本目录既有的 PMLR 官方 PDF 提取图。所有文件均已实际解码并通过缩略总览检查。
