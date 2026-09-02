# FlexPPI 融入蛋白从头设计与突变筛选的 pipeline 方案

> 本文的文献证据面（别人怎么画 pipeline 图、六种主流画法）见 [pipeline_figure_literature_review.md](pipeline_figure_literature_review.md)；按本方案绘制的主图已生成：`paper/figures/make_flexppi_design_pipeline_figure.py` → `paper/figures/fig_flexppi_design_pipeline.{pdf,svg,png}`，QA 见 `paper/figures/fig_flexppi_design_pipeline_QA.md`。

## 核心定位

FlexPPI 最合理的位置不是候选生成端，也不是亲和力预测器，而是位于结构可信度与静态界面筛选之后、昂贵 MD 和湿实验之前的 **dynamics-aware filter / reranker（动力学感知过滤与重排序模块）**。它接收候选的静态 apo 与 holo 结构，输出逐残基 RMSF、侧链构象熵、SASA、结合诱导变化和预测不确定性。现有模型不会生成骨架或序列，也不能把“结合后更刚”直接解释为更高亲和力；因此必须与复合物结构预测、静态界面/亲和力评分和实验验证并列使用。

## 文献中的典型画法

| 工作 | 图的组织方式 | 可借鉴之处 | 不应照搬之处 |
|---|---|---|---|
| [Cao et al., Nature 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9117152/) | 两阶段、八个编号步骤；计算与实验用不同底色 | 用阶段色带和编号步骤讲清候选如何逐渐收窄 | 细节过多，不适合把 FlexPPI 网络内部再塞进去 |
| [ProteinMPNN, Science 2022](https://pmc.ncbi.nlm.nih.gov/articles/9997061/) | 输入几何图、编码器、随机顺序解码分面展示 | 用真实输入/输出名，避免“Module A”式空标签 | 它是模型架构图，不是完整设计 campaign 图 |
| [RFdiffusion, Nature 2023](https://www.nature.com/articles/s41586-023-06415-8) | 噪声结构到设计结构的连续快照，旁边展示不同条件任务 | 用重复的小结构图标表现状态转换 | 不要把 200 个去噪步展开到应用 pipeline 中 |
| [EvoPro, Nature Communications 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10187191/) | 序列变异、AF2 结构预测、打分、遗传算法构成环 | 突变筛选若有反馈，应明确画闭环 | 没有真实实验反馈时不要为了“高级感”画闭环 |
| [EVOLVEpro, 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11275896/) | 小样本实验、蛋白语言模型、回归器、下一轮候选构成主动学习环 | 将实验读数明确回流到候选选择 | FlexPPI 本身未在线更新，回流应进入排序器/校准器而不是假装更新 FlexPPI |
| [BindCraft, Nature 2025](https://www.nature.com/articles/s41586-025-09429-6) | AF2-Multimer 生成、MPNN 优化、AF2/Rosetta 过滤的简洁横向流程 | 主图只保留三至五个大阶段；把工具名放在阶段内部 | 其 target flexibility 是结构预测过程中的灵活性，不等于 FlexPPI 的动态标签 |
| [Pythia-PPI, National Science Review 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12199698/) | 局部结构图输入、模型架构、多任务/自蒸馏训练三面板 | 突变效应模型与 FlexPPI 应作为互补分支，而非互相替代 | 其输出是结合亲和力变化，不应与 FlexPPI 的 RMSF 变化共用同一颜色或单位 |
| [DynamicMPNN, 2025/ICLR 2026](https://arxiv.org/abs/2507.21938) | 多构象集合经共享编码汇聚；实线为生成，虚线为结构自洽性评估 | 明确区分主数据流与评估支路；突出多状态而非单一构象 | FlexPPI 预测的是每个静态状态的动态描述符，不是直接生成多状态序列 |
| [BinderFlow, PLOS Computational Biology 2025](https://journals.plos.org/ploscompbiol/article/figures?id=10.1371%2Fjournal.pcbi.1013747) | 编号线性流程与并行批处理架构并列 | 适合把论文主图与工程执行 DAG 分开 | 不要在论文主图中展示调度器、文件名和 GPU 作业细节 |
| [Deep learning–guided design of dynamic proteins, Science 2025](https://www.science.org/doi/full/10.1126/science.adr7094) | “设计 → 结构预测 → 动态评估 → 实验”四段式总览，动态指标作为显式 gate | 动态评估是一级独立 gate，与 FlexPPI 的槽位一致 | 其目标是构象切换/内在动态；FlexPPI 的 gate 内容是结合诱导的局部涨落变化 |
| [Hybrid AI/physics pipeline: BRD3, ChemComm 2025](https://pubs.rsc.org/cc/article/61/96/19028/883845/Hybrid-AI-physics-pipeline-for-miniprotein-binder) | AI 生成（RFdiffusion+MPNN+AF2 过滤）→ 物理/MD 层重排序 → 实验验证 | FlexPPI 与这条管线 MD 层槽位一致；图上画成“AI 动力学层 + MD 灰色复核支路” | MD 仍是更昂贵的复核，图上是“先粗筛后复核”，不是并列二选一 |
| [Dynamics-informed RFdiffusion (KRAS), JBSD 2025](https://www.tandfonline.com/doi/full/10.1080/07391102.2025.2609688) | MD 构象集合作为扩散模型生成条件 → 序列 → 过滤 | 动态信息出现的位置必须自明（生成端 vs 筛选端） | FlexPPI 在筛选端，图上不能暗示它参与生成 |
| [ProteinDJ, Protein Science 2025](https://www.ovid.com/journals/pros/fulltext/10.1002/pro.70464~proteindj-a-high-performance-and-modular-protein-design) / [Ovo, bioRxiv 2025](https://www.biorxiv.org/content/10.1101/2025.11.27.691041v1.full) | 模块化 pipeline：骨架 → 序列 → 验证 → 过滤，模块可替换 | FlexPPI 画成“可插拔 reranker 插槽”，强调模块边界 | 不要画成只有一个固定实现的黑盒 |
| [nf-binder-design, Zenodo 2025](https://zenodo.org/records/16809705) | Nextflow 工程 DAG（fold-and-dock、设计、过滤、聚类） | 工程执行图与论文主图分离 | 论文主图不画调度器、文件名和 GPU 细节 |
| [PaperDraw, 2026](https://www.sciencedirect.com/science/article/pii/S2468502X2600029X) | 自动生成论文 pipeline 图的设计原则：标准化尺寸、一致颜色语义、分层布局、统一箭头方向 | 可作为本图的自检清单（与 `docs/main_figure_design.md` 质量门一致） | 不要机械套用模板而丢掉本图的双入口语义 |

## 推荐的计算 pipeline

### 1. 定义设计任务与动态目标

输入包括目标蛋白结构、目标表位或界面、允许突变的位置、必须保留的催化/保守残基，以及一个可检验的动态目标。动态目标不能只写成“越刚越好”，而应具体到区域，例如“结合后界面核心适度降低 RMSF，但保留 hinge/CDR 环的顺应性；远端结构域不发生大范围异常重排”。如果没有生物学上合理的目标轮廓，FlexPPI 输出应作为候选描述特征和探索性排序依据，而不应成为硬淘汰门槛。

### 2. 两条候选生成入口

**从头设计入口**：目标结构与表位约束 → RFdiffusion、BindCraft 或同类方法生成候选骨架/界面 → ProteinMPNN、SolubleMPNN 或 FAMPNN 赋予序列。FlexPPI 不参与这一步。

**突变筛选入口**：野生型复合物 → 界面、热点、变构位点或允许区域定义 → 穷举单点突变、组合突变、ProteinMPNN/蛋白语言模型建议或主动学习建议。若候选很多，先做便宜的序列过滤，再构建结构。

### 3. 为每个候选建立可比较的 apo/holo 结构集合

每个候选至少需要一份游离态结构和一份复合物结构。对从头设计 binder，可分别建立“游离 binder → binder–target 复合物”和“游离 target → target–binder 复合物”两对输入，从两侧读取伙伴诱导响应。对突变体，不能只在野生型坐标上替换残基名称；应对 mutant apo 和 mutant holo 进行一致的结构预测或局部松弛。建议每个候选使用多个随机种子或多个结构模型，以免 FlexPPI 排名主要反映上游单次结构预测的偶然误差。

### 4. 分层过滤，先便宜后昂贵

第一层是序列与可开发性过滤，如异常疏水斑块、净电荷、聚集风险、化学不稳定基序和表达风险。第二层是折叠与复合物可信度，如 pLDDT、界面 PAE、复合物模型一致性、界面接触与几何互补。第三层是静态物理或监督式亲和力信号，如 Rosetta/FoldX 类能量、Pythia-PPI 或其他经过同源隔离验证的 \(\Delta\Delta G\) 模型。只有通过这些门槛的候选进入 FlexPPI；这样既节约计算，也避免给明显错误的复合物结构解释“动态机制”。

### 5. FlexPPI 动力学重排序

对每个候选输出三类信息：绝对 apo/holo 状态、带正负号的 partner-induced response、预测方差。界面可用 8 Å 重原子距离定义，也可沿用项目当前的接触定义。推荐保留以下独立特征，不要一开始压成一个不可解释的总分：

1. **界面动态匹配**：候选的逐残基 \(\Delta\mathrm{RMSF}\)、\(\Delta H_\chi\)、\(\Delta\mathrm{SASA}\) 与预先定义的目标轮廓之差。
2. **远端扰动惩罚**：非界面功能区或需要保持稳定的结构域上，预测响应绝对值是否异常增大。
3. **双侧响应**：分别评价 target 和 binder；不要只看其中一侧。
4. **不确定性**：高预测方差的候选不应因为均值好看就被排到最前。
5. **内部一致性**：直接响应头与 `holo state − apo state` 的差异可作为模型分歧指标；分歧大时送入 MD 复核，而不是直接判优。
6. **结构种子稳定性**：同一候选在多个预测结构上的动态特征均值与离散度。

若有期望动态轮廓 \(\tau_i\)，可使用一个易解释的标准化偏差：

\[
L_{\mathrm{dyn}}
=\frac{\sum_i w_i\,|\widehat{\Delta y}_i-\tau_i|/
\sqrt{\widehat{\sigma}_i^2+\epsilon}}
{\sum_i w_i}.
\]

这里 \(i\) 是残基，\(w_i\) 表示界面核心、热点、hinge 或远端保护区的重要性，\(\widehat{\sigma}_i\) 是模型不确定性。这个式子表示“偏离目标越多、且模型越自信，惩罚越大”；它不是结合自由能。

突变筛选应比较突变体相对于野生型对“结合诱导响应”的二次变化：

\[
\delta_{\mathrm{mut}}(i)
=\widehat{\Delta y}_{\mathrm{mut}}(i)
-\widehat{\Delta y}_{\mathrm{WT}}(i).
\]

这比直接比较 mutant holo 与 WT holo 更接近“该突变如何改变结合前后动态重分配”，因为两者都先扣除了各自游离态基线。

### 6. 多目标选择，不把所有信号粗暴相加

最终选择至少同时考虑折叠可信度、复合物可信度、静态界面/亲和力分数、可开发性、FlexPPI 动态特征与不确定性。首选 hard gates + Pareto front：先淘汰明显不折叠、界面不可信或不可开发的候选，再保留在多个目标上都没有被另一候选全面压倒的 Pareto 候选。若必须给单一总分，应在一套独立验证集上学习权重，并报告去掉 FlexPPI 特征后的消融结果。

### 7. 分级验证与反馈

对动态排名靠前、模型分歧大或具有新颖机制的候选做短 MD/增强采样；随后进入表达、单体稳定性、结合动力学与特异性实验。若论文声称“动态机制”，还应选少量候选用 HDX-MS、NMR、时间分辨实验或足够重复的 MD 验证逐区域变化。实验数据可以训练一个轻量候选排序/校准模型，并驱动下一轮候选生成；除非重新训练并严格验证，否则不要在图中画成“实验数据直接更新 FlexPPI”。

## 推荐的论文主图

### 图类型与范式

- 类型：Solution Overview。
- 范式：**六阶段圆形迭代 pipeline**。
- 目的：把“从头生成—第一次柔性诊断—高柔性位点改造—AlphaFold 回折—第二次柔性诊断—下一轮”画成明确闭环。
- 与现有图的关系：`fig_flexppi_architecture` 保留为模型方法图；新图单独作为应用/展望图。两张图不要合并。
- 旧的四面板线性版本 `fig_flexppi_design_pipeline` 可保留为补充材料或工程实现图；正文主图采用用户确认的圆形版本。

### 版式草图

画布建议为双栏通栏 178 mm × 124 mm。六张圆角卡片沿顺时针圆环排列，中央椭圆只保留一句总逻辑和一句设计原则：

```text
                         ① 从头设计起点
                    ↗                         ↘
    ⑥ FlexPPI 再次柔性分析                      ② 现有工具生成候选
              ↑                                      ↓
              │   FlexPPI 引导的迭代蛋白设计          │
              │   生成→诊断→定点改造→复评             │
              │   降低非期望柔性，保留功能性柔性       │
              ↑                                      ↓
    ⑤ AlphaFold 结构预测与过滤                  ③ FlexPPI 第一次柔性分析
                    ↖                         ↙
                  ④ 高柔性候选位点筛选与替换
```

圆形版本已经生成：`paper/figures/fig_flexppi_iterative_design_pipeline.drawio` 为可编辑源文件，另有 PDF、SVG 和 PNG 预览；质量记录位于同名 `.audit.md` 文件。

### 六个步骤该画什么

**① 从头设计起点**：目标结构与表位、功能约束和期望柔性轮廓。这里使用“从头设计”而不是“从头合成”，避免让读者误以为已经进入湿实验合成。

**② 现有工具生成候选**：列出 RFdiffusion、BindCraft 和 ProteinMPNN 等代表性工具，输出多个不同的初始蛋白候选。

**③ FlexPPI 第一次柔性分析**：输出 RMSF、侧链熵、SASA 和不确定性，用蓝—灰—红逐残基色带标记柔性由低到高。

**④ 高柔性候选位点筛选与替换**：先取 Top-k 高柔性候选位点，再排除催化位点、关键结合位点和功能性铰链，最后做约束性氨基酸替换。图中必须写出这道保护规则，不能表达成“最柔的残基全部替换”。

**⑤ AlphaFold 结构预测与过滤**：预测改造后蛋白及其复合物结构，并检查置信度、clash 与界面是否保持；只有结构可信的候选进入下一步。

**⑥ FlexPPI 再次柔性分析**：比较改造前后的柔性轮廓，同时检查预测不确定性与结构一致性。达到目标的候选进入备选池，未达到目标者通过朱红虚线箭头进入下一轮设计。

### 图形语法

- 顺时针主数据流使用 3 px 深灰实线和 block 箭头；⑥→① 的“下一轮设计”使用朱红虚线，与主数据流明确区分。
- 起点蓝 `#0072B2`，生成工具橙 `#E69F00`，两次 FlexPPI 紫 `#CC79A7`，位点替换朱红 `#D55E00`，AlphaFold/结构通过绿 `#009E73`，连接线与说明文字灰 `#6B7280`。颜色必须同时配合编号、标签和线型。
- 在 1400 单位宽源画布中，最小文字为 17 px；以 178 mm 宽插入论文时约为 8.2 pt。若缩窄插图，必须同步放大源字号。
- 输出 PDF/SVG；流程图内使用简化矢量蛋白轮廓，不嵌入低分辨率截图。真实结构渲染另放案例图。
- 不使用 3D 阴影、渐变底框、发光效果、装饰性 DNA 图标或超过六种语义颜色。

### 建议图注

> **FlexPPI guides an iterative protein-design cycle by identifying candidate high-flexibility residues and reassessing the redesigned structures after AlphaFold prediction.** Existing generative tools first produce de novo candidate proteins, whose residue-wise RMSF, side-chain entropy, SASA, and uncertainty are estimated by FlexPPI. High-flexibility positions are considered for constrained substitution only after catalytic, binding-critical, and functionally flexible residues are protected. AlphaFold then predicts and filters the redesigned structures before a second FlexPPI pass compares the original and redesigned flexibility profiles; candidates that do not meet the predefined structural and dynamic criteria enter the next design round.

## 验证与消融

1. 在家族、界面模板和近同源均隔离的测试集上比较 `static scores only` 与 `static scores + FlexPPI`。
2. 报告候选排名的 enrichment factor、top-k hit rate 和 wet-lab success rate，而不仅是 FlexPPI 自身的残基相关系数。
3. 对去掉不确定性、去掉 apo、只用 holo 绝对状态、只用 state subtraction、只用直接响应头分别做消融。
4. 检查被 FlexPPI 淘汰的实验阳性候选，避免把“过度刚化”规则系统性地误当成好设计。
5. 从头设计候选属于明显的分布外样本；必须按序列新颖性、折叠新颖性和结构预测置信度分层报告。

## Figure Designer 质量门

- Vector format：设计通过；导出后仍需确认没有嵌入低分辨率位图。
- Font size：设计通过；按最终 178 mm 宽度检查所有文字不低于 8 pt。
- Colour-blind safe：通过；采用 Okabe–Ito 并配合文字/线型双编码。
- Self-contained caption：通过；首句直接说明 FlexPPI 的位置和边界。
- Honest axis range：主流程图不适用；若面板 (d) 使用 Pareto 散点，必须标清轴名、单位和阈值来源。
- No chartjunk：通过。
- Integrity gate：方案层面 7 项通过；成图后的字体、矢量属性、颜色对比和实际箭头交叉仍需视觉复核。
