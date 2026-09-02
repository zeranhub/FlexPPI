# FlexPPI 背景、引言与图片叙事调研（2026-08-30）

本报告涉及的十张核心参考图、四篇原文 PDF 和逐图用途索引已保存到 [`docs/literature_figures/`](literature_figures/README.md)。

## 一、结论先行：FlexPPI 应该被定义成什么

FlexPPI 最稳妥、也最容易获得审稿人认可的定位，不是又一个蛋白–蛋白相互作用位点预测器，也不是 AlphaFlow、BioEmu 一类构象集合生成器，而是一个 **partner-aware、rigid-motion-invariant、multi-task residue-dynamics predictor**：输入一个静态受体结构以及可选的伙伴链环境，输出每个残基的绝对动态状态，包括主链位置波动 RMSF、侧链转子状态多样性 $H_\chi$ 和溶剂可及表面积 SASA，并在有真实 apo/holo 配对数据时预测有符号的状态变化。中文可概括为“从静态多链结构预测伙伴条件化的逐残基动态状态”。如果实验中的 apo 与 holo 使用同一套结合态坐标、只删除伙伴链，那么该实验只能称作“伙伴上下文响应”或“固定坐标伙伴消融”，不能称为真实的“结合诱导构象响应”，因为它没有观察未结合态重新松弛后的结构与轨迹。论文的新意应落在“伙伴条件化、多描述符、逐残基、带不确定性、并显式区分绝对状态与响应”这一组合上，而不要未经系统查新声称“首次”。

一句适合放在摘要或引言末段的英文定位是：**FlexPPI is a partner-aware, rigid-motion-invariant, multi-task predictor of residue-wise dynamic states from static protein complexes, coupled with paired and controlled-context evaluations of binding-associated response.** 若当前主要实验仍是固定结合态坐标删链，则最后半句应改成更保守的 **controlled evaluations of partner-context sensitivity**。

## 二、背景应该怎样逐层收窄

建议将背景写成六段，而不是从“PPI 很重要”开始罗列大词。第一段先建立具体矛盾：蛋白质功能来自构象集合，而晶体、冷冻电镜或结构预测通常提供一个平均或代表性构象；两个状态的平均骨架甚至界面几乎不变，结合仍可能重新分配局部主链波动、侧链转子占据和溶剂暴露。这里可引用关于构象集合与分子识别的综述，以及 Frederick 等和 Tzeng & Kalodimos 关于蛋白识别中内部动态/构象熵作用的工作，但必须说明 FlexPPI 的 $H_\chi$ 是从 MD 转子占据计算的 Shannon 熵，并不等于热力学结合熵或自由能项。第二段再说测量瓶颈：NMR、重复全原子 MD 可以提供丰富动态信息，却难以在大规模复合物上常规使用；[ATLAS](https://academic.oup.com/nar/article/52/D1/D384/7438909) 和 [DynaRepo](https://academic.oup.com/nar/article/54/D1/D393/8322558) 等标准化重复 MD 数据库才让监督式动态预测成为可行问题，DynaRepo 对复合物和结合/未结合状态的覆盖尤其贴近 FlexPPI。

第三段把已有方法分成三条路线，而不是混成一张长清单：第一条是 NMA/GNM、弹性网络和 CABS-flex 等物理近似，速度快但表示能力受近似限制；第二条是从序列、结构或密度图直接回归紧凑动态描述符的模型，包括 [RMSF-net](https://www.nature.com/articles/s41467-024-49858-x)、[Flexpert](https://openreview.net/forum?id=L238BAx0wP)、[FastProtFlex](https://academic.oup.com/bioinformatics/article/42/4/btag175/8654525)、[DynaProt](https://arxiv.org/abs/2509.01038)、SeqDance 与 RocketSHP，这类方法适合蛋白组规模筛选，但多数以单链或伙伴无关的状态表征为主；第三条是 [AlphaFlow](https://proceedings.mlr.press/v235/jing24a.html)、[BioEmu](https://www.microsoft.com/en-us/research/project/biomolecules/) 一类构象集合生成器，输出更丰富，却计算更重，而且最终仍需把集合归纳成可以直接比较的逐残基统计量。第四段指出真正空缺：已有 PPI 几何网络如 [PeSTo](https://www.nature.com/articles/s41467-023-37701-8) 和 DeepInteract 擅长识别界面或残基接触，但“哪里接触”不等于“伙伴使哪里变得更硬、更软、转子更多或暴露改变”；已有动态预测又常把目标链孤立处理。因此真正欠缺的不是另一个整体 RMSF 分数，而是静态多链几何条件下、可分辨残基和动态维度的状态/响应预测及其外部复合物评估。

第五段在介绍模型前先把任务语义讲清楚。绝对状态任务是 $\hat y=f(X,C)$，其中 $X$ 是目标链静态结构，$C$ 是可选伙伴环境，$y$ 可以是 RMSF、$H_\chi$ 或 SASA；真实配对响应是 $\Delta y=y_{\mathrm{holo}}-y_{\mathrm{apo}}$，要求 apo 和 holo 各自有独立结构或轨迹；固定坐标伙伴消融则比较 $f(X,C)$ 与 $f(X,\varnothing)$，它回答的是“模型是否利用了伙伴上下文”，不是“真实未结合态会怎样松弛”。这一区分最好在引言和第一张图中都出现一次，否则审稿人很容易把 partner deletion 误读为因果性的 binding-induced dynamics。

第六段再引出 FlexPPI：共享的固定局部框架和 invariant point attention 编码多链几何，绝对状态头同时输出 RMSF、$H_\chi$ 和 SASA，配对响应头读取 apo 表示、holo 表示、二者差值与交互项，并用“直接响应约等于两种绝对状态之差”的一致性约束连接两个任务；预测方差与重复轨迹/时间窗的标签方差共同表达不确定性。结果段可以提前点出一个比“整体相关性更高”更有研究价值的发现：伙伴上下文可以改善整体排序，但模型可能把伙伴作用学成近似单向约束，从而在某些真实界面或肽段区域错误地普遍压低动态性。这一失败模式正好支持为什么要做分区、分链和响应方向的评估。

## 三、可以直接参考的引言组织方式

| 文献 | 引言的核心推进方式 | FlexPPI 可借鉴之处 | 不应照搬之处 |
|---|---|---|---|
| RMSF-net, Nature Communications 2024 | 动态重要 → 静态 PDB/密度图缺少动态 → MD 贵 → 既有密度模型数据小且结构/密度贡献不清 → 联合输入模型 | “功能需求—观测缺口—计算瓶颈—数据机会—方法”的五步收窄 | RMSF-net 是冷冻电镜体素任务，不能把密度图需求移植到 FlexPPI |
| Flexpert, ICLR 2025 | 生成式蛋白设计成功 → 当前设计偏静态 → 柔性对功能重要 → 标签昂贵且异质 → 三项方法贡献 | 从“静态设计能力强”自然过渡到“缺少可控动态属性” | FlexPPI 不是柔性条件化序列/结构设计器，不要用 designability 作为主要任务 |
| AlphaFlow, ICML 2024 | 蛋白是构象集合 → 单结构预测不能表示异质性 → 通过扰动 MSA 采样受限 → 生成模型学习集合 | 清楚地区分单结构、逐残基统计量和完整构象集合三种输出层级 | FlexPPI 不生成轨迹或自由能景观，不能暗示恢复了完整平衡分布 |
| DynaProt, ICLR 2026 | 生成集合表达力强但昂贵 → 单一 RMSF 太弱 → NMA 有方向但近似 → 明确提出“表达力–效率”问题 | 用任务空间而非逐篇罗列相关工作；解释三个输出为何比单 RMSF 更丰富 | 若没有可操作定义，不要画主观的“我们的模型位于 Pareto 前沿” |
| BackFlip/FliPS, ICML 2025 | 静态蛋白设计成熟 → 刚性限制功能 → 目标柔性 profile → 生成候选并由预测器闭环验证 | 适合未来“FlexPPI 用于动态感知重排序”的应用叙事 | 目前若没有生成/实验闭环，不应把重排序愿景当作已验证主贡献 |
| DynaRepo, NAR 2026 | 静态结构数据库成熟 → 动态仍是下一瓶颈 → 高质量重复长轨迹数据库补空缺 | 直接支撑复杂物、重复轨迹和外部评估需求 | 数据库论文不是模型精度证据 |
| PeSTo, Nature Communications 2023 | 多种相互作用都可统一为局部原子几何与界面概率 | 支撑多链几何和伙伴特异性的编码需求 | 界面分类不等同于动态状态预测 |

## 四、Figure 1 最值得采用的方案

### 4.1 图类型与叙事范式

图类型应是 **Motivated Example**，采用 **running example + failure case** 范式，而不是把当前预训练架构图直接升格为 Figure 1。最佳故事是“同一个目标蛋白在没有伙伴和存在伙伴时，平均结构可以相近，但逐残基 RMSF、侧链转子熵和 SASA 的变化具有局部性、符号性和维度差异；伙伴无关模型只能给一个通用柔性 profile，FlexPPI 则显式输出两个绝对状态及其响应，并允许检查哪些区域预测正确、哪些区域被过度刚化”。这样一张图同时回答四个审稿问题：任务是什么、为什么单静态结构不够、伙伴信息究竟改变什么、模型输出和实验验证是什么。

备选范式中，纯“现有方法 vs 我们”过于抽象，无法展示当前 ITrace/TCR 结果暴露出的方向性失败；纯性能 teaser 只显示相关性，反而重复了论文正在批评的全局指标；诱导契合/构象选择自由能漏斗虽然生物学上经典，但 FlexPPI 不预测动力学路径和自由能，因此作为中央示意会造成过度主张。

### 4.2 建议版式

建议双栏宽 178 mm、图高约 75–82 mm，三栏从左到右：

**(a) Binding redistributes local dynamics.** 左上画同一受体的 apo 静态结构，左下画 holo 复合物，目标链始终用蓝色，伙伴用橙色；两幅结构旁边放同一条序列坐标上的三条简化 profile，分别是 RMSF、$H_\chi$ 和 SASA。用一个实线框标出界面处刚化，用另一个虚线框标出远端或伙伴侧增柔，明确写出“not uniformly suppressed”。若真实数据暂时只支持一处增柔/一处减柔，就只画有数据支撑的两处，不要为了对称杜撰机制。

**(b) What existing predictors miss.** 输入仅画孤立受体或单链结构，经过一个灰色的“partner-agnostic dynamics predictor”，输出一条通用 profile；在旁边用两个短标签指出“no partner context”和“no signed apo→holo response”。这里不需要罗列模型名称，相关文献放在图注和正文。

**(c) FlexPPI task and audit.** 上下两路分别输入 apo 目标与 holo 复合物，进入标有 shared $\theta$ 的同一个几何编码器，输出 $\hat y_{apo}$、$\hat y_{holo}$、$\Delta\hat y$ 和 uncertainty。右侧放一个很小的验证块：overall / interface / peptide / partner 分区比较，并以“✓ localized response”和“! over-rigidification failure”双重编码成功与失败。真实 apo/holo 配对应使用实线箭头；同坐标删伙伴的 context ablation 使用虚线箭头，并在图例中明确二者不是同一实验。

建议第一句图注写成：**“FlexPPI separates absolute residue dynamics from partner-context response, enabling direct tests of where a static multichain model succeeds and fails.”** 后续图注再依次解释颜色、三个描述符、实线/虚线含义和右侧验证分区。

### 4.3 视觉规范

目标链用 Okabe–Ito 蓝 `#0072B2`，伙伴用橙 `#E69F00`，模型用紫 `#CC79A7`，预测正确用蓝绿色 `#009E73`，失败警示用朱红 `#D55E00`，背景模块用中性灰；颜色之外必须同时使用实线/虚线、圆点/三角形或 ✓/!，保证黑白打印和色盲可读。最终尺寸下文字不小于 8 pt，结构、箭头和 profile 均优先导出为可编辑 SVG/PDF；若使用真实分子结构，最好从 ChimeraX/PyMOL 导出矢量轮廓或在排版软件中重绘为简化丝带，而不是嵌入低分辨率截图。Figure 1 不放损失函数、IPA 内部公式、张量维度或全部特征名称，这些内容属于方法图。

## 五、最值得借鉴的已有图片

| 图片 | 图在讲什么 | 值得借的视觉语法 | 对 FlexPPI 的使用建议 |
|---|---|---|---|
| AlphaFlow Fig. 1 | 谐波先验样本经模型逐步形成构象集合 | 同一蛋白始终保持一致配色；四个面板只讲一个核心机制；从问题到输出的方向清楚 | 借它的“单一主叙事”和一致分子视觉，不要借流场/去噪符号，因为 FlexPPI 非生成模型 |
| BackFlip/FliPS Fig. 1 | 目标柔性 profile → 生成多个骨架 → 预测回柔性 → MD 验证 | 输入目标、候选、预测、验证构成闭环；结构按柔性着色 | 适合未来应用图或设计重排序图；当前模型论文 Figure 1 可借“结构 + profile + validation”三联画 |
| Flexpert Fig. 1 | Flexpert-Seq 嵌套在 Flexpert-3D 中，结构适配器补充序列模型 | 共享/冻结/可训练模块关系一眼可见 | 在 FlexPPI 方法图中用 shared $\theta$ 显示 apo/holo 共享编码器；不要把全部网络层搬入动机图 |
| RMSF-net Fig. 1 | PDB/密度输入、MD 标签构建、U-Net 与重建流程 | 清楚分离“训练时才需要的 MD 标签支路”和“推理输入” | 当前 FlexPPI 预训练图已部分采用这一点；应再加入伙伴响应分支，否则第一方法图没有体现核心新意 |
| DynaProt Fig. 1 | 动态表示方法的表达力–效率二维概念图 | 用任务地图迅速定位研究空缺 | 只有在横纵轴有明确可测定义时才画；否则容易被质疑为主观排名 |
| PeSTo Fig. 1 | 局部原子几何 → 多尺度几何网络 → 残基界面着色 | 架构与真实结构输出并列，几何信息流简洁 | 方法图可把固定局部框架、伙伴环境和 residue-wise outputs 放在同一行 |
| ATLAS protein page / DynaRepo graphical abstract | 结构、轨迹、逐残基 profile 与重复模拟 | 数据来源和最终标签之间建立直观对应 | 背景或数据图可借；ATLAS Fig. 1 本身偏数据库组成，不是最佳动机图 |
| SeqDance Fig. 1 | sequence → ensemble → function → evolution，并比较不同表征来源 | 把“静态表示”和“动态信息”放在同一条信息链 | 若画大领域背景可借，但 FlexPPI 应把 multichain structure 而非 sequence 放在输入中心 |

## 六、对现有 FlexPPI 图片的审计

当前 `paper/figures/f1_pipeline.pdf` 是一张完成度较高的预训练流程图：MD 样本、离线标签构建、几何编码器、分布头和 noise-aware NLL 的训练/推理边界清楚；问题是它只解释绝对状态预训练，没有体现伙伴条件、多链响应或 apo/holo 共享编码器，因此作为方法 Figure 2 合适，作为整篇论文最先出现的 Figure 1 则属于 **MAJOR** 级叙事缺口。`fig_flexppi_architecture.png` 技术信息准确，但模块、公式和箭头过密，更像补充材料中的完整网络图；`fig_flexppi_pretraining.png` 同样把标签、分布头和几何模块同时铺开，缩到双栏宽后可读性不足；`fig_flexppi_design_pipeline.png` 讲候选生成、结构过滤、动态重排序和验证，适合愿景或应用论文，但当前交叉箭头较多，而且若尚无完整设计闭环数据会把论点推得过远；`fig_itrace_main.png` 的结果图较清楚，可作为实验主图，但应增加 interface/peptide/partner 分区和响应方向，而不是只依赖整体相关性。

建议最终图序为：Figure 1 动机与任务定义；Figure 2 绝对状态预训练与数据标签；Figure 3 伙伴条件化共享编码器和配对响应；Figure 4 内部/外部绝对状态结果；Figure 5 ITrace 或案例级 partner ablation 成功与失败。若篇幅只允许三张主图，则把 Figure 2 和 Figure 3 合成上下两层的训练/推理总览，把当前预训练图移入补充材料。

## 七、图设计完整性与诚信门检查

设计类型为 Motivated Example，采用 running example + failure case。中心主张是“伙伴作用会重新分配而不是简单统一压低逐残基动态；FlexPPI 把绝对状态和上下文响应分开建模并允许分区审计”。图中每个面板都有正文对应：背景段对应 (a)，相关工作缺口对应 (b)，任务定义与模型概述对应 (c)。当前设计在逻辑闭合、输入输出可追踪、颜色双重编码、线型区分真实配对与消融、以及避免自由能/轨迹生成暗示方面通过；最终排版后仍需人工检查 8 pt 字号、灰度打印、投影可读性，并确保引言、方法和实验复用同一个 running example。严禁把模型预测的 $H_\chi$ 写成“结合构象熵”或“自由能贡献”，严禁把同坐标删链写成真实 apo/holo，严禁用主观坐标把 FlexPPI 放在无定义的 Pareto 最优位置。

## 八、优先阅读清单

1. Song et al., **RMSF-net: prediction of protein flexibility using 3D convolutional neural network and cryo-EM data**, Nature Communications, 2024。直接逐残基动态预测的近邻工作，也是“MD 标签只在训练时需要”图形叙事的代表。
2. Kouba et al., **Flexpert: Learning to Predict and Control Protein Flexibility**, ICLR, 2025。最接近“从静态输入预测 RMSF”的机器学习论文，适合比较数据来源、输出定义和结构适配器。
3. Jing et al., **AlphaFlow: sampling structures from a continuous distribution of conformations**, ICML, 2024。用于在引言中界定“直接描述符预测”和“构象集合生成”的边界。
4. Lewis et al., **Scalable emulation of protein equilibrium ensembles with generative deep learning / BioEmu**, Science, 2025。代表高表达力集合生成路线。
5. Prianichnikov et al., **ATLAS**, Nucleic Acids Research, 2024；Kouba et al., **DynaRepo**, Nucleic Acids Research, 2026。用于数据机会、重复轨迹、复合物与外部评估背景。
6. Viliuga et al., **BackFlip/FliPS**, ICML, 2025。代表把柔性预测器嵌入设计闭环的应用叙事。
7. Krapp et al., **PeSTo**, Nature Communications, 2023；DeepInteract, NeurIPS AI4Science, 2022。用于说明多链几何网络解决的是 partner-specific interaction，但尚不等于 partner-conditioned dynamics。
8. Frederick et al., **Conformational entropy in molecular recognition by proteins**, Nature, 2007；Tzeng & Kalodimos, **Protein activity regulation by conformational entropy**, Nature, 2012。用于说明结合可在平均结构变化很小时重分配内部动态，但须严格限定 $H_\chi$ 的统计学含义。
9. Clark et al., **Inherent versus induced protein flexibility**, PLOS Computational Biology, 2019。支持把主链波动和侧链转子变化分开，并提醒 apo–apo、holo–holo 与 apo–holo 变异不能混为一谈。
10. DynaProt, **Learning residue-level protein dynamics with multiscale Gaussians**, ICLR, 2026；FastProtFlex, Bioinformatics, 2026。代表 2026 年紧凑动态描述符预测的新进展。

## 九、数据库检索记录与限制

检索范围为 2010–2026 或 2018–2026，核心查询包括：`protein flexibility prediction from static structure residue RMSF dynamics`、`protein protein binding induced dynamics apo holo conformational change residue`、`geometric deep learning protein dynamics conformational ensemble static structure`。Semantic Scholar、arXiv、Crossref、DBLP 共返回 54 条来源记录；其中包含明显离题条目，正文结论只采用经过标题、摘要或原文核验的核心文献。数据库的引用数只是检索时快照，不用于比较模型质量。

检索异常原文如下：`[open_alex] unavailable (import failed: unsupported operand type(s) for |: 'type' and 'NoneType'); skipping this source.`；OpenReview 多次返回 HTTP 429；Exa MCP 返回 `You've hit Exa's free MCP rate limit`。因此本报告改用出版社、会议论文集、arXiv、PMC 与作者/机构官方页面交叉核验，未把搜索摘要当作最终证据。第一次本地检索还遇到 Windows GBK 的 `UnicodeEncodeError`，随后使用 UTF-8 环境重新执行并成功获得结果。

### 9.1 查询一：静态结构到逐残基柔性（18 条原始记录）

| 来源 | 年份 | 标题 | 引用快照/备注 |
|---|---:|---|---:|
| Semantic Scholar | 2024 | RMSF-net: prediction of protein flexibility using 3D convolutional neural network and cryo-EM data | 82 |
| Semantic Scholar | 2024 | Investigating protein flexibility using an atom-level graph neural network | 1 |
| Semantic Scholar | 2025 | DynaProt: Learning residue level protein dynamics with multiscale Gaussians | 3 |
| Semantic Scholar | 2025 | De novo protein backbone generation with controlled complexity using relaxed dimensionality constraints | 8；即 BackFlip/FliPS 相关记录 |
| Semantic Scholar | 2026 | Protein Language Models Beyond Static Folds: Reasoning about Dynamics, Energetics and Function | 3 |
| Semantic Scholar | 2026 | Hybrid Protein-Protein Complex Affinity Prediction with Molecular Dynamics Simulation | 0 |
| arXiv | 2022 | Alphafold prediction is an old protein with a novel structural fold | 0；离题 |
| arXiv | 2019 | A Learning Based Approach to Identify S-Glutathionylation Sites Using Protein Sequence-Derived Properties | 0；离题 |
| arXiv | 2025 | Kinetics of stress-induced unfolding of proteins: inherent mass-balance | 0；离题 |
| arXiv | 2025 | DynaProt: Learning residue level protein dynamics with multiscale Gaussians | 0；与上项重复 |
| arXiv | 2026 | Static-to-Dynamics Protein Representation Learning: A Review of Structural Embedding Paradigms | 0；预印本综述 |
| arXiv | 2023 | Introduction to Protein Structure Prediction | 0；背景材料 |
| Crossref | 未标年 | Impact of Protein Dynamics on Residue-Residue Coevolution and Contact Prediction | 5 |
| Crossref | 未标年 | Hybrid protein-ligand binding residue prediction via protein language models: a structure-free approach | 3；非 PPI 动态 |
| Crossref | 2020 | easyAmber: A comprehensive toolbox to automate the molecular dynamics simulation of proteins | 0 |
| Crossref | 未标年 | Single-Sequence, Structure Free Allosteric Residue Prediction through Masking Protein Language Models | 9 |
| Crossref | 未标年 | CopulaNet: Learning residue co-evolution directly from multiple sequence alignment for protein structure prediction | 4；离题 |
| Crossref | 2025 | Tau protein structure and dynamics | 14；综述/章节 |

### 9.2 查询二：结合、apo/holo 与残基响应（18 条原始记录）

| 来源 | 年份 | 标题 | 引用快照/备注 |
|---|---:|---|---:|
| Semantic Scholar | 2024 | Comparing AlphaFold2-predicted and experimentally determined apo/holo structures by randomized alanine scanning | 8 |
| Semantic Scholar | 2021 | Capturing the structural and dynamic consequences of ligand binding with MD | 结果集中出现的体系性研究，需按具体体系使用 |
| Semantic Scholar | 2024 | Apo/holo structural and dynamic response studies returned for specific ligand systems | 多条；非通用 PPI 模型 |
| Semantic Scholar | 2026 | Protein complex/dynamics studies returned for specific systems | 多条；不作为核心证据 |
| arXiv | 2017 | ISLAND: in-silico proteins binding affinity prediction using sequence information | 0；亲和力而非动态 |
| arXiv | 2014 | Conformational selection in protein binding and function | 0；对应 Weikl 与 Paul 综述 |
| arXiv | 2022 | Alphafold prediction is an old protein with a novel structural fold | 0；离题 |
| arXiv | 2025 | Kinetics of stress-induced unfolding of proteins: inherent mass-balance | 0；离题 |
| arXiv | 2014 | Long-range modulation of a secondary binding site in an unbound protein | 0；特定体系 |
| arXiv | 2010 | Deciphering the Molecular Mechanism of Water Boiling at Heterogeneous Interfaces | 0；离题 |
| Crossref | 2024 | Comparing AlphaFold2-predicted and experimentally determined apo/holo structures by randomized alanine scanning | 8 |
| Crossref | 2010 | An analysis of the apo/holo conformational variability of Ca-binding protein cavities | 27 |
| Crossref | 2018 | Apo and Holo Fatty Acid Binding Protein Simulations in CHARMM and AMBER Force Fields | 18 |
| Crossref | 2010 | Conformational transitions upon ligand binding: Holo-structure prediction from apo conformations | 8 |
| Crossref | 2010 | Structure, dynamics, and spectral tuning of photoactive yellow protein | 0；特定体系 |
| Crossref | 2019 | Faculty Opinions recommendation of Inherent versus induced protein flexibility | 0；原论文另行核验 |
| DBLP | — | 本查询无返回 | 0 |
| OpenAlex/OpenReview | — | 因上述异常未获得稳定结果 | 0 |

### 9.3 查询三：几何深度学习与构象集合（18 条原始记录）

| 来源 | 年份 | 标题 | 引用快照/备注 |
|---|---:|---|---:|
| Semantic Scholar | 2025 | Scalable emulation of protein equilibrium ensembles with generative deep learning | 364；BioEmu |
| Semantic Scholar | 2025 | DeepConformer: Iterative Conformer Generation for Small Molecules | 0；小分子，离题 |
| Semantic Scholar | 2025 | SpatPPI: predicting protein-protein interactions through spatial representation learning | 2；PPI 分类，非动态 |
| Semantic Scholar | 2026 | Mavchen-1: An Open Foundation Model for Molecular Dynamics | 0；新近工作 |
| Semantic Scholar | 2026 | Hyaline: Efficient Protein Dynamics Simulation via Neural Operators | 0；新近工作 |
| Semantic Scholar | 2025 | Angular Deviation Diffuser for Molecular Conformation Generation | 0；更偏分子构象 |
| arXiv | 2026 | Rigidity-Aware Geometric Pretraining for Protein Representation Learning | 0 |
| arXiv | 2018–2026 | 其余五条由关键词误召回的几何/生成条目 | 0；不纳入结论 |
| Crossref | 2018–2026 | 六条蛋白–配体/CASP 论文的审稿或 decision-letter 记录 | 不是真正独立论文，全部剔除 |
| DBLP | — | 本查询无返回 | 0 |
| OpenAlex/OpenReview | — | 因上述异常未获得稳定结果 | 0 |

### 9.4 检索汇总

三个查询共获得 54 条数据库来源记录，真正与本课题直接相关的高价值主题有四类：静态结构到逐残基动态描述符、构象集合生成、伙伴特异几何/界面建模、以及 apo/holo 或结合相关动态。标题层面最常出现的概念是 protein、structure、dynamics/flexibility、conformation、binding 与 prediction，说明文献已经从“是否预测结构”转向“怎样表示结构周围的动态分布”，但“伙伴条件化的逐残基多描述符响应”仍没有形成像界面预测或构象生成那样成熟的标准任务。引用快照中最显著的是 BioEmu、RMSF-net、钙结合位点 apo/holo 结构分析、脂肪酸结合蛋白 MD 比较与 Tau 动态综述；引用数反映传播时间和主题规模，不能视为方法优劣。

建议阅读顺序是：先读 ATLAS/DynaRepo 明确标签和数据边界，再读 RMSF-net/Flexpert 理解直接回归范式，再读 AlphaFlow/BioEmu 理解完整集合路线，随后读 PeSTo/DeepInteract 补足多链几何，最后读 DynaProt/FastProtFlex 与最新工作定位 2026 年竞争格局。对 FlexPPI 最关键的文献缺口不是“没有人研究柔性”，而是这些研究通常分别解决单链动态、界面几何或构象生成，尚未把伙伴上下文、多个可解释逐残基动态量、成对响应与受控消融统一到同一个可审计框架中。
