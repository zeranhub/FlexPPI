# FlexPPI：逐残基蛋白质柔性为什么值得预测——引言论证、文献与图片

调研日期：2026-08-30

> 更新说明（2026-08-31）：本文保留“逐残基预测为何有用”的下游论证；论文背景与引言开头请优先采用新版的[“从蛋白质柔性与动态结构开始”报告](protein_flexibility_dynamic_structure_background_review_2026.md)。新版先建立构象集合、能量景观与跨尺度动态，再引出逐残基动态摘要和 FlexPPI，避免一开始就陷入对预测用途的辩护。

## 结论先行

面对“预测逐残基柔性没有用”的质疑，最弱的回答是“蛋白质不是刚体，所以柔性很重要”，第二弱的回答是“模型与 MD 的 RMSF 相关很高”。前者没有说明预测结果会改变什么决定，后者只证明模型拟合了标签。更有说服力的定位是：**逐残基柔性是一种低成本、空间可定位的动态中间表征，它把单张静态结构与昂贵构象集合之间的缺口压缩成一张可操作的 residue map，用来回答应当在哪些残基允许运动、优先突变、寻找变构口袋、进行柔性对接，或追加高成本 MD。**“逐残基”是关键，因为一个全局平均数会把局部刚化和软化相互抵消；“伙伴条件化”同样关键，因为复合物形成可以在不同界面类型、界面核心/外围乃至远端位置产生不同响应；“多描述符”也是关键，因为主链 RMSF 小不代表侧链转子和溶剂暴露不变。

这一定义也限定了 FlexPPI 的承诺：它不是 MD 的替代品，不输出运动方向、残基耦合或严格的结合自由能；它是一个用于大规模初筛、假设生成和采样分配的动态图层。RMSF 表示位置波动幅度，$H_\chi$ 表示侧链转子占据的不确定性/多样性，SASA 表示溶剂暴露，三者共同给出“哪里会动、以何种局部方式变化”的紧凑描述。尤其要避免把 $H_\chi$ 直接称为热力学构象熵：它可以与构象熵相关，却不等于完整的 $\Delta S$，更不能单独推出 $\Delta G$。

## 审稿人真正会追问的三件事

第一，他会问为什么结构、pLDDT/B-factor 或一个整体 RMSD 不够。Clark 等对 305 个蛋白、2369 个 holo 结构和 1679 个 apo 结构的分析表明，apo–holo 的额外骨架位移通常小于 0.5 Å，但结合位点侧链 $\chi_1$ 会频繁进入 apo 结构没有覆盖的取向；这说明“整体骨架几乎不变”与“局部状态没有变化”不是同一件事。[Clark Figure 1](literature_figures/utility_figures/clark_2019_fig1_backbone_change.jpg) 与 [Figure 3](literature_figures/utility_figures/clark_2019_fig3_sidechain_rotamers.jpg) 应当成对看。第二，他会问为什么一定要逐残基。Cukuroglu 等指出，约 65% 的核心界面残基在结合前就近似构象不变，而核心、外围和非刚性界面残基承担不同角色；Dobbins 等又发现 22 个二聚体的不同界面类型具有不同的结合前后 mobility 关系。因此一个复合物级总分无法告诉实验者哪个残基是刚性锚点、哪个是适应性外围、哪个是远端响应点。[Cukuroglu Figure 1](literature_figures/utility_figures/cukuroglu_2012_fig1_interface_flexibility.jpg) 与 [Dobbins Figure 4](literature_figures/utility_figures/dobbins_2010_fig4_bound_unbound_mobility.jpg) 正好支撑这一层。第三，他会问得到 residue map 后能做什么。已有研究已经把局部动态用于变构口袋定位、突变效应预测和柔性条件化设计：Panjkovich & Daura 在 91 个非冗余变构蛋白中观察到 70% 的案例在配体结合后有显著柔性变化，将柔性变化与保守性结合后，在严格设置下定位变构口袋的阳性预测值为 65%；SeqDance/ESMDance 将野生型与突变体的预测动态差用于稳定性突变效应；FliPS 则把目标逐残基柔性 profile 作为生成约束，并以 MD 验证生成结构与目标 profile 的一致性。

## 推荐的六段引言逻辑

**第一段：用一个反直觉实例建立问题，而不是从“PPI 很重要”开始。** 可以先说，结构生物学与结构预测通常提供一个代表性坐标，但分子识别由构象集合完成；两个复合状态的平均骨架可以十分相似，结合亲和力和残基内部运动却显著不同。Tzeng & Kalodimos 的 CAP–DNA 体系最适合做文献引子：相似的结构/活性构象占比并不能解释不同变体的结合亲和力，而结合诱导的构象熵变化可以增强或抑制结合。随后马上落到 FlexPPI 的真实 running example——6UON：同一结合态坐标和接触图不能直接告诉读者肽段 p6/p8 是 MD 中的柔性 ridge，也不能说明删除 TCR 上下文后哪些残基的预测最敏感。这样第一段就从“蛋白会动”升级为“静态几何缺少局部、条件化的动态分辨率”。可参考 [Tzeng Figure 2](literature_figures/utility_figures/tzeng_kalodimos_2012_fig2_conformational_entropy_enhances_binding.png) 和本项目的 [6UON 图](literature_figures/utility_figures/flexppi_6uon_running_example.png)。

**第二段：解释逐残基分辨率为什么必要。** 一次结合不是把整条链均匀压硬。界面核心常由预组织、低柔性的锚点稳定，外围残基更可能承担适应性重排；不同类别复合物还可能出现结合后刚化、近似保持或更复杂的局部重分配。更重要的是，骨架位移、侧链转子与溶剂暴露不必同步。于是全局 RMSF、单一置信度或接触概率都无法回答“应当动哪一段”。这一段可引用 Cukuroglu 2012、Dobbins 2010、Cole & Warwicker 2002 和 Clark 2019，并自然引出 RMSF、$H_\chi$、SASA 的互补性。

**第三段：把预测目标写成一个会改变下游决策的中间变量。** 建议直接定义：给定静态多链复合物，目标是预测每个残基的动态状态及其对伙伴上下文的响应，使研究者能够在大量候选中选择后续的柔性 docking 区域、突变位点、变构口袋、设计约束或 MD 采样对象。Panjkovich 2012 提供变构口袋证据，SeqDance/ESMDance 提供突变效应证据，Flexpert 和 BackFlip/FliPS 提供设计闭环证据。这里应写明预测器的价值不是产生“好看的热图”，而是降低把昂贵模拟平均撒到所有体系、所有残基上的成本。

**第四段：用不超过三条限制收束现有方法。** 第一，NMA/弹性网络与快速物理近似虽然便宜，但对复杂侧链与伙伴特异响应的表达有限；第二，现有学习式柔性预测多以单链 RMSF 或伙伴无关状态为主，难以解释“同一目标链因不同伙伴而如何变化”；第三，构象生成模型信息更丰富，却计算更重，最后仍需把集合归纳为可比较的局部统计量。由此得到的缺口不是“还缺一个 RMSF 回归器”，而是**缺少静态多链几何条件下、可分辨残基和动态维度、并能报告伙伴上下文响应与不确定性的快速预测框架**。

**第五段：引入 FlexPPI 时保持模块与挑战一一对应。** 可以写成：FlexPPI 以多链结构为输入，通过界面感知的几何表征区分目标链内部环境与跨链接触；以共享编码器估计完整复合物和受控伙伴上下文下的逐残基 RMSF、$H_\chi$ 与 SASA；再将绝对状态、上下文差值和校准后的不确定性分开输出。对应的三个挑战分别是多链几何条件化、异质动态标签的联合学习，以及绝对状态与伙伴响应的可辨识/校准。若当前实现尚未包含其中一项，不应在引言中先行承诺。

**第六段：贡献必须包含用途证据，而不仅是相关系数。** 建议三项贡献分别是任务/数据、模型/预测和下游/机制验证。最关键的是第三项：至少用一个外部下游任务证明动态信息相对静态基线带来增益，或证明预测图能富集已知功能残基并节约高成本采样。仅有 RMSF Spearman、MAE 和几个结构着色案例时，“有用性”质疑仍然成立；这不是写作问题，而是证据链缺了一环。

## 可直接改写进论文的中文底稿

蛋白质–蛋白质识别通常由动态构象集合完成，而晶体、冷冻电镜或结构预测模型只提供其中一个代表性坐标。因而，近乎不变的整体骨架并不意味着结合前后具有相同的局部状态：大规模 apo/holo 结构比较显示，配体结合带来的额外骨架位移通常小于 0.5 Å，但结合位点侧链仍会占据 apo 结构未覆盖的转子状态；NMR 研究也表明，内部动态和构象熵的重新分配能够在平均结构相近时显著改变分子识别热力学。对蛋白复合物而言，这些变化还具有明显的空间异质性——预组织的界面核心、适应性外围以及远端响应残基可能表现出不同甚至相反的动态变化。因此，一个复合物级柔性分数或一张静态接触图不能回答实验与建模中真正可操作的问题：究竟应当在哪些残基允许构象采样、优先突变、寻找变构调控，或投入昂贵的分子动力学模拟。

逐残基动态预测为这一缺口提供了一种紧凑的中间表征。它不试图替代完整轨迹，而是将轨迹中与局部运动幅度、侧链转子多样性和溶剂暴露有关的信息压缩为可映射回结构的 residue-wise profile，从而支持大规模候选的筛选与排序。已有工作已经证明这种表征可以服务于下游决策：柔性变化可与序列保守性结合定位变构口袋，预测动态属性的突变前后变化可改善稳定性突变效应预测，而目标柔性 profile 可以作为蛋白序列或骨架生成的控制条件。然而，现有方法通常以单链或伙伴无关的状态为中心，难以刻画复合物中同一残基随伙伴上下文改变而产生的局部响应。

基于这一观察，我们提出 FlexPPI：从单个静态多链结构快速预测每个残基的主链波动、侧链转子熵和溶剂暴露，并显式建模伙伴上下文。FlexPPI 的目标不是用一个标量宣称某个蛋白“更柔”或“更硬”，而是生成一张可校准、可定位并可用于下游选择的动态地图。我们进一步评估该地图是否能够复现实验/MD 所揭示的界面异质性、伙伴相关响应及功能残基富集，并量化它在突变排序或模拟预算分配中相对于静态结构特征的增量价值。

## Figure 1 应如何画

Figure 1 应该是 motivated example，而不是网络架构。推荐继续使用 6UON，因为仓库中已有 ITrace 三条轨迹均值、完整五链预测、固定坐标的 partner-context ablation、接触和肽段位置映射。面板 a 只画同一复合物与静态接触，提出问题“静态图能否指出应关注哪个肽残基？”；面板 b 画 MD target 与模型逐残基 profile，突出 p6/p8 柔性 ridge，同时展示一个只看静态/伙伴无关信息时会错过的残基；面板 c 用 full-minus-pMHC 的条形图说明伙伴上下文响应，并明确标注这是 fixed-coordinate context sensitivity，不是真实 apo/holo 轨迹；面板 d 只放一个下游动作，例如“优先选择这些残基做突变或局部 MD”，不要把完整方法流程塞进来。图注第一句直接写 finding，例如：*A static contact map does not identify the peptide residues that dominate the bound-state dynamic response, whereas FlexPPI recovers a localized flexibility ridge and partner-context sensitivity at residue resolution.*

参考图片的视觉作用各不相同：Tzeng Figure 2 教你把逐残基动态变化映射到三维结构并连接热力学结果；Clark Figure 1+3 教你用一对图制造“骨架小变化—侧链大重排”的反差；Cukuroglu Figure 1 教你比较界面核心、外围与非界面残基；Dobbins Figure 4 教你画 bound-versus-unbound/不同类型的响应异质性；BackFlip/FliPS Figure 1 教你把 profile 变成一个闭环行动变量。完整本地索引见 [utility_figures/README.md](literature_figures/utility_figures/README.md)，缩略总览见 [contact_sheet.png](literature_figures/utility_figures/contact_sheet.png)。

## 哪些下游验证最能消除“没用”的质疑

优先级最高、与 FlexPPI 最贴近的是突变或界面 residue ranking：在 SKEMPI 2.0、ProteinGym 的 PPI/稳定性子集或其他具有单点突变标签的数据上，将 FlexPPI 的绝对状态和 partner-response 特征加入只含序列/静态几何的基线，比较按蛋白划分后的 Spearman、AUROC/average precision 和 calibration；必须报告静态基线、单一 RMSF、多描述符以及 partner-context ablation。第二选择是已知功能位点富集：使用 Allosteric Database 或 Catalytic Site Atlas，检验高响应/高不确定性区域是否比 SASA、接触度、B-factor 或 pLDDT 更能富集已知位点。第三选择是采样预算实验：固定总 MD 预算，让 FlexPPI 决定哪些复合物或哪些局部区域进入长模拟，比较与随机选择、置信度选择和静态几何选择在发现高波动/状态转换残基上的召回率。只要完成其中一项，论文就能把“预测标签”升级为“改变决策”；若暂时都不能完成，应把用途写成可检验的目标/潜在应用，不能写成已经证明的生物学价值。

## 核心参考文献与最适合借鉴的引入方式

| 文献 | 可引用的事实 | 引言中的角色 | 本地图片 |
|---|---|---|---|
| Frederick et al., *Conformational entropy in molecular recognition by proteins*, Nature 2007, [DOI](https://doi.org/10.1038/nature05959) | calmodulin 结合不同肽时内部动态变化与结合热力学相关 | 证明平均结构之外的内部动态参与分子识别 | 未下载；Tzeng 图更适合直接参考 |
| Tzeng & Kalodimos, *Protein activity regulation by conformational entropy*, Nature 2012, [DOI](https://doi.org/10.1038/nature11271) | 构象熵变化可增强或抑制 CAP–DNA 结合 | 最强的反直觉开场案例之一 | [Fig. 2](literature_figures/utility_figures/tzeng_kalodimos_2012_fig2_conformational_entropy_enhances_binding.png)、[Fig. 4](literature_figures/utility_figures/tzeng_kalodimos_2012_fig4_entropy_dominates_binding.png) |
| Cole & Warwicker, *Side-chain conformational entropy at protein–protein interfaces*, Protein Science 2002, [DOI](https://doi.org/10.1110/ps.0222702) | 界面通常比其余表面侧链更受限，侧链熵可帮助识别真实界面 | 为 $H_\chi$ 和 residue-wise interface map 提供直接先例 | [Fig. 1](literature_figures/utility_figures/cole_warwicker_2002_fig1_sidechain_entropy_interface.jpg)、[Fig. 6](literature_figures/utility_figures/cole_warwicker_2002_fig6_residue_entropy_contribution.jpg) |
| Dobbins et al., *Comparing interfacial dynamics in protein-protein complexes*, BMC Struct Biol 2010, [DOI](https://doi.org/10.1186/1472-6807-10-26) | 不同界面类别、核心/外围的结合前后 mobility 关系不同 | 说明 partner response 不是统一刚化 | [Fig. 4](literature_figures/utility_figures/dobbins_2010_fig4_bound_unbound_mobility.jpg) |
| Cukuroglu et al., *Roles of residues in the interface of transient protein-protein complexes before complexation*, Sci Rep 2012, [DOI](https://doi.org/10.1038/srep00334) | 约 65% 核心界面残基近似构象不变，刚性与非刚性界面残基角色不同 | 说明为什么要逐残基而非整个界面一个分数 | [Fig. 1](literature_figures/utility_figures/cukuroglu_2012_fig1_interface_flexibility.jpg)、[Fig. 2](literature_figures/utility_figures/cukuroglu_2012_fig2_profilin_actin.jpg) |
| Panjkovich & Daura, *Exploiting protein flexibility to predict the location of allosteric sites*, BMC Bioinformatics 2012, [DOI](https://doi.org/10.1186/1471-2105-13-273) | 70% 变构蛋白出现显著柔性变化；严格设置 PPV 65% | 最直接的“预测柔性→选择口袋”证据 | [Fig. 4](literature_figures/utility_figures/panjkovich_daura_2012_fig4_allosteric_pockets.jpg) |
| Campbell et al., *The role of protein dynamics in the evolution of new enzyme function*, Nat Chem Biol 2016, [DOI](https://doi.org/10.1038/nchembio.2175) | 磷酸三酯酶到芳基酯酶的定向进化伴随构象亚态和非生产性运动的重排 | 把动态连接到酶工程和功能进化 | 相关综述图 [Fig. 6](literature_figures/utility_figures/petrovic_kamerlin_2018_fig6_arylesterase_evolution.jpg) |
| Petrović & Kamerlin, *Conformational dynamics and enzyme evolution*, Interface Focus 2018, [DOI](https://doi.org/10.1098/rsif.2018.0330) | 总结动态集合、远端突变与可进化性/可设计性 | 设计应用的背景综述 | [Fig. 1](literature_figures/utility_figures/petrovic_kamerlin_2018_fig1_dynamics_evolvability.jpg)、[Fig. 8](literature_figures/utility_figures/petrovic_kamerlin_2018_fig8_flexibility_enables_design.jpg) |
| Clark et al., *Inherent versus induced protein flexibility*, PLOS Comput Biol 2019, [DOI](https://doi.org/10.1371/journal.pcbi.1006705) | 整体骨架变化小，但结合位点侧链取向可显著重排 | 支撑 RMSF+$H_\chi$+SASA 的多描述符必要性 | [Fig. 1/3/7](literature_figures/utility_figures/README.md) |
| Kouba et al., *Flexpert*, ICLR 2025, [OpenReview](https://openreview.net/forum?id=L238BAx0wP) | 逐位置柔性指令可进入 inverse folding 并改变设计分布 | 证明 flexibility profile 可以是控制变量 | [Fig. 2](literature_figures/utility_figures/flexpert_2025_fig2_design.png) |
| Viliuga et al., *BackFlip/FliPS*, ICML 2025, [PMLR](https://proceedings.mlr.press/v267/viliuga25a.html) | 目标 profile 驱动骨架生成，结果以 300 ns MD 验证 | 当前最强的直接“逐残基柔性有用”图例 | [Fig. 1](literature_figures/utility_figures/backflip_flips_2025_fig1_target_profile_to_backbone.png) |
| Hou et al., *Protein language models trained on biophysical dynamics inform mutation effects*, PNAS 2026, [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC12846831/) | 动态属性改善突变效应预测，尤其对缺少进化信息的设计/病毒蛋白 | 证明动态表征能给静态/进化表征带来增量信息 | 预印本 [Fig. 4](literature_figures/utility_figures/seqdance_esmdance_2025_fig4_mutation_effect.jpg) |

## 自动检索记录（全部结果，不等于全部纳入正文）

使用 `paper-search` 在 2000–2026 年执行两个聚焦查询。OpenAlex 因本地导入错误不可用，DBLP 与 OpenReview 未返回结果；Crossref、Semantic Scholar 与 arXiv 共返回 48 条记录。下面保留全部记录，避免只展示支持预设结论的论文。标题明显偏离“逐残基柔性用途”的记录仍列出并标为低相关。

### 查询 A：`per-residue protein flexibility functional sites binding interfaces mutation effects enzyme design molecular dynamics`

| 来源 | 年份 | 标题 | 引用快照 | 处理 |
|---|---:|---|---:|---|
| Crossref | 2021 | Review for “Structural insights on effects of mutation on interactions between 14-3-3ζ and a neurodegenerative disease-associated protein” v1 | 0 | 审稿记录，排除 |
| Crossref | 2021 | 同题 v2 | 0 | 审稿记录，排除 |
| Crossref | 2021 | Decision letter for the same manuscript v2 | 0 | 决定信，排除 |
| Crossref | 2021 | Decision letter for the same manuscript v1 | 0 | 决定信，排除 |
| Crossref | 2005 | Cancer-related mutations in BRCA1-BRCT studied by molecular dynamics | 0 | 主题邻近，非核心 |
| Crossref | 2022 | Conformational flexibility and local frustration in SARS-CoV-2 Spike variants | 3 | 突变/动态案例，次要 |
| Crossref | — | Structural insights on effects of mutation on interactions between 14-3-3ζ and a neurodegenerative disease-associated protein | 0 | 元数据不完整，次要 |
| Crossref | 2006 | Ligand-dissociation increases protein flexibility in TRAP | 0 | 结合/柔性案例，次要 |
| Semantic Scholar | 2024 | Mechanism of mutation-induced effects on catalytic function of TEV protease: a molecular dynamics study | 8 | 支持突变—动态—催化，次要 |
| Semantic Scholar | 2025 | Accurate characterization of allosteric energy landscapes in oncogenic KRAS | 3 | 预印本，支持变构应用 |
| Semantic Scholar | 2026 | Molecular modeling protein–peptide complexes in the Factor Xa–Prothrombin system | 0 | 新近、次要 |
| Semantic Scholar | 2022 | Deciphering isoniazid drug resistance in KatG using dynamics and residue networks | 19 | 动态残基网络/耐药案例 |
| Semantic Scholar | 2025 | Atomistic profiling of KRAS interactions with mutational scanning | 1 | 预印本，次要 |
| Semantic Scholar | 2014 | Compensatory role of N348I/M184V in nevirapine response | 19 | 突变/动态案例，次要 |
| Semantic Scholar | 2009 | Exploring molecular design of protein interaction sites with MD and free energy | 54 | 早期设计应用综述/方法 |
| Semantic Scholar | 2026 | Brunner syndrome mutations disrupt serotonin binding and active-site stability | 0 | 新近病例，次要 |
| arXiv | 2017 | ISLAND: sequence-based protein-binding affinity prediction | — | 与柔性不直接，排除 |
| arXiv | 2016 | Chitin–protein interface study | — | 低相关，排除 |
| arXiv | 2014 | Predicting evolutionary site variability from structure: buriedness, packing, flexibility and design | — | 支持柔性与进化位点 |
| arXiv | 2021 | Protein–protein interactions enhance thermal resilience of SpyRing enzymes: an MD study | — | 复合/稳定性案例 |
| arXiv | 2025 | PRIMRose: per-residue energy metrics for double insertions/deletions | — | 新近、邻近任务 |
| arXiv | 2024 | Generative enzyme design guided by functionally important sites and substrates | — | 功能位点设计，非柔性核心 |
| arXiv | 2025 | UniZyme cleavage-site predictor with active-site knowledge | — | 低相关，排除 |
| arXiv | 2014 | Conformational selection in protein binding and function | — | 概念背景，纳入阅读 |

### 查询 B：`residue flexibility predicts catalytic residues allosteric sites protein engineering RMSF B-factor`

| 来源 | 年份 | 标题 | 引用快照 | 处理 |
|---|---:|---|---:|---|
| Crossref | 2011 | A geometry-based generic predictor for catalytic and allosteric sites | 29 | 功能位点预测背景 |
| Crossref | 2012 | Exploiting protein flexibility to predict allosteric sites | 77 | 核心纳入 |
| Crossref | 2017 | Coupled residue–residue dynamics in protein allosteric mechanisms | 1 | 说明 RMSF 不含耦合，边界文献 |
| Crossref | 2012 | Accurate prediction of catalytic residues by side-chain orientation and contact density | 13 | 静态功能位点基线 |
| Crossref | 2026 | Figure 5 object from a peptide paper | 0 | 数据对象，排除 |
| Crossref | 2022 | Editor evaluation of coevolution-based allosteric-residue work | 0 | 评审记录，排除 |
| Crossref | — | GPCR allosteric regulation preprint | 6 | 主题邻近、元数据不完整 |
| Crossref | — | Automated machine learning for allosteric sites (meeting abstract) | 0 | 会议摘要，排除 |
| Semantic Scholar | 2024 | Conformational flexibility of a C-terminal region influences a distal active site | 3 | 远端动态—活性案例 |
| Semantic Scholar | 2024 | S–S bridge mutation enhances activity and thermostability | 4 | 工程案例，次要 |
| Semantic Scholar | 2020 | Engineering a heparin-binding pocket to enhance catalytic efficiency | 20 | 工程案例，次要 |
| Semantic Scholar | 2020 | Key dynamic residues in human cyclophilins | 0 | residue-level dynamics，邻近 |
| Semantic Scholar | 2026 | TriFit: dynamics-informed mutation fitness prediction | 0 | 最新邻近工作，需单独复核 |
| Semantic Scholar | 2026 | Theaflavin–MMP molecular dynamics study | 0 | 配体 MD 个案，排除 |
| Semantic Scholar | 2020 | Allosteric-site hits for PfProRS | 19 | 变构应用案例 |
| Semantic Scholar | 2015 | Mutational analysis of cholera toxin A1 interaction | 13 | 突变/相互作用案例，次要 |
| arXiv | 2022 | Prediction of allosteric signaling pathways and functional residues with optimized paths | — | 动态/变构位点方法 |
| arXiv | 2006 | Ligand binding, protein fluctuations and allosteric free energy | — | 理论背景 |
| arXiv | 2018 | Blind prediction of B-factor/protein flexibility | — | 预测基线背景 |
| arXiv | 2014 | Control of catalytic efficiency by a co-evolving network | — | 催化网络背景 |
| arXiv | — | AlphaFold knot-related record | — | 误召回，排除 |
| arXiv | 2018 | General characteristics of allosteric communication residues | — | residue-level allostery 背景 |
| arXiv | 2012 | Residue mobility decreased during protein evolution | — | 演化/柔性背景 |
| arXiv | 2012 | Co-evolving site pairs as a contact predictor | — | 非柔性核心，排除 |

## 搜索小结与局限

两次检索共 48 条自动记录，其中最直接支持本论证的是 Panjkovich & Daura 的变构位点工作；其余较强证据主要来自随后用 DOI、PMC、会议官网和作者公开 PDF 手工核验的 Cukuroglu、Dobbins、Clark、Tzeng、Frederick、Cole & Warwicker、Campbell、Flexpert、BackFlip/FliPS 与 SeqDance/ESMDance。自动检索的确切异常为：`[open_alex] unavailable (import failed: unsupported operand type(s) for |: 'type' and 'NoneType'); skipping this source.`；多个源出现 `Rate limited. Waiting 3 seconds...`；OpenAlex、DBLP 与 OpenReview 最终各返回 0 条。按 `agent-reach` 路径尝试通过 Jina Reader 读取两个 PMC 页面时均返回 HTTP 401，因此事实和图片改用 Europe PMC XML、PMC 原页、期刊/会议页面及公开 PDF 交叉核验。引用计数是检索时快照，只用于排序阅读优先级，不用于学术结论。
