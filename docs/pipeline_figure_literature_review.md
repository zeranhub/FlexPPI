# FlexPPI 融入设计/筛选 pipeline 的文献检索报告：别人怎么画 pipeline 图

> 检索日期：2026-08-24。检索渠道：本会话内联 web 检索（Exa MCP 免费额度触发 429 限流，直接抓取被沙箱网络策略阻断，因此图注与流程描述来自论文图注页、摘要、正文片段和公开解读，结论均以原文为准）。检索主题：de novo binder 设计 pipeline、ML 突变筛选/定向进化 loop、动力学感知设计、pipeline 图形画法。

## 一、结论先行：主流画法可以归纳为六种模式

| 模式 | 典型代表 | 图形特征 | 适合讲什么 |
|---|---|---|---|
| ① 线性漏斗 | [BindCraft, Nature 2025](https://www.nature.com/articles/s41586-025-09429-6/figures/1)、[nf-binder-design, Zenodo 2025](https://zenodo.org/records/16809705)、[BRD3 hybrid AI/physics pipeline, ChemComm 2025](https://pubs.rsc.org/cc/article/61/96/19028/883845/Hybrid-AI-physics-pipeline-for-miniprotein-binder?silentauthchecked=true) | 从宽到窄的一串 gate，每个 gate 写过滤规则（pLDDT/iPTM/RMSD/接触/能量），可选候选数标注；实验验证在末端 | 候选如何逐级收窄、FlexPPI 插在第几级 |
| ② 双/多入口汇合 | [RFdiffusion, Nature 2023](https://www.nature.com/articles/s41586-023-06415-8)（多任务条件面板）、[ProteinDJ, Protein Science 2025](https://www.ovid.com/journals/pros/fulltext/10.1002/pro.70464~proteindj-a-high-performance-and-modular-protein-design)、[Ovo, bioRxiv 2025](https://www.biorxiv.org/content/10.1101/2025.11.27.691041v1.full) | 上下两条泳道（从头设计 / 突变或条件任务），在某一级汇合 | 一个模型同时服务两条入口 |
| ③ 闭环主动学习 | [EVOLVEpro, Science 2025](https://www.science.org/doi/10.1126/science.adr6006)、[AI 生物铸造厂综述, 2025](https://www.sciencedirect.com/science/article/pii/S0958166925001247)、[ML 酶工程综述, ACS Cent. Sci. 2023](https://pubs.acs.org/doi/10.1021/acscentsci.3c01275) | 圆环或弧形反馈箭头，轮次编号；实验读数从末端流回起点 | 实验反馈回流（FlexPPI 应回流到排序器，不回流到模型权重） |
| ④ 主数据流 + 虚线评估支路 | [DynamicMPNN, ICLR 2026 / arXiv:2507.21938](https://arxiv.org/abs/2507.21938) | 实线 = 生成主流程；虚线 = 结构自洽/评估回路；多构象输入用共享编码器汇聚 | 区分"生产流"与"校验流" |
| ⑤ 阶段色带 + 编号步骤 | [Cao et al., Nature 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9117152/)（两阶段八步骤）、多数 Nature 方案图 | 整条 pipeline 被浅色大色带切成阶段，每阶段内编号步骤，计算/实验底色不同 | 层次化叙事，避免细节淹没主线 |
| ⑥ 工程 DAG | [nf-binder-design, GitHub](https://github.com/Australian-Protein-Design-Initiative/nf-binder-design) | 调度器级流程（进程名、文件、并行分支） | 只用于技术文档/仓库 README，不用于论文主图 |

## 二、逐篇证据：每一篇的图具体怎么画

### 1. BindCraft（Nature 2025）—— Fig. 1 "De novo binder design using BindCraft"

- 画法：横向单线漏斗。输入目标结构 → 循环的骨架生成（RFdiffusion 或 AF2-Multimer 播种）+ 序列设计（ProteinMPNN）→ AF2-Multimer 结构预测 → 一串显式过滤（pLDDT、iPTM、RMSD、界面接触、能量项）→ 实验验证（SPR/ELISA）→ 命中。
- 借鉴点：工具名直接写在阶段框内；漏斗是"one-shot"而非闭环；Extended Data Fig. 1 单独放 in silico 分析细节，主图不塞超参。
- 注意：target flexibility 是其结构预测流程中的概念，不等于 FlexPPI 的 MD-RMSF 动态标签。

### 2. RFdiffusion（Nature 2023）—— Fig. 1 模型概览

- 画法：去噪扩散过程用"连续快照"表现（噪声结构逐步变成设计结构），旁边一列展示不同条件任务（monomer、binder、motif scaffolding、对称 oligomer）。
- 借鉴点：用重复的小结构图标表达"状态演化/构象转换"；用并排任务列表达"一个生成器多任务"。
- 注意：不要把 200 步去噪展开到应用 pipeline 图里。

### 3. EVOLVEpro（Science 2025）—— Fig. 1 in silico 定向进化

- 画法：闭环。初始变体 → 小样本实验读数 → 蛋白语言模型（ESM）嵌入 → 回归模型训练 → 全库打分/预测 → 选下一轮变体 → 回到实验，轮次编号（round 1、2、3…）。
- 借鉴点：实验反馈回流画成弧线回起点；PLM 与回归器是两个独立框，不合并。
- 注意：FlexPPI 未在线更新，回流应画进"候选排序/校准器"，而不是假装更新 FlexPPI。

### 4. Pythia-PPI（National Science Review 2025）

- 画法：三面板模型图——输入局部结构图（界面残基子图）→ 多任务 + 自蒸馏训练框架 → 输出 ΔΔG。
- 借鉴点：突变效应模型与 FlexPPI 是互补分支：Pythia-PPI 出结合亲和力变化，FlexPPI 出结合诱导动态变化，两者在图上应并列而非合并。
- 注意：它的 ΔΔG 与 FlexPPI 的 ΔRMSF 单位/语义不同，不要共用同一颜色或色带。

### 5. DynamicMPNN（ICLR 2026，arXiv:2507.21938）—— 多状态序列设计

- 画法：多构象集合分别进独立编码通道 → Deep Set pooling 汇聚 → 序列解码；另画一条虚线回路表示结构自洽性评估（预测结构返回检验）。
- 借鉴点：**实线 = 主数据流，虚线 = 评估/校验**，这一图形语法直接可用；突出"多状态"而不是单一构象。
- 注意：FlexPPI 输出的是每个静态状态的动态描述符，不生成多状态序列，框内文字要写明边界。

### 6. Deep learning–guided design of dynamic proteins（Science 2025, adr7094）

- 画法：Fig. 1 是"通用方法"总览——把"动态"作为设计目标：序列/结构设计 + 结构预测 + 动力学验证（MD/实验），并用动态指标（RMSF/构象切换）评估设计是否达到目标动态。
- 借鉴点：这是目前最接近"为动态而设计"的画法：**设计 → 预测 → 动态评估 → 实验**四段式，动态评估是一级显式 gate。
- 注意：该文以构象切换/内在动态为主要目标；FlexPPI 更偏"结合诱导的局部涨落变化"，可以沿用四段式但换掉 gate 内容。

### 7. Hybrid AI/physics pipeline for miniprotein binder prioritization: BRD3（ChemComm 2025）

- 画法：AI 生成（RFdiffusion + ProteinMPNN + AF2 过滤）→ **物理层（MD 模拟/自由能/动力学打分）优先排序** → 实验验证。这是"AI 生成 + 物理动力学重排"双层的标准画法。
- 借鉴点：**FlexPPI 的槽位与这条管线的 MD 层完全一致**——在静态 AI 打分之后、湿实验之前；图里把 AI 动力学层画成紫色高亮、物理 MD 层画成灰色虚线验证支路即可同时体现"替代/预筛"与"复核"。
- 注意：MD 仍是更昂贵的复核手段，图上 FlexPPI 与 MD 是"先粗筛后复核"关系，不是并列二选一。

### 8. Dynamics-informed RFdiffusion for KRAS（J. Biomol. Struct. Dyn. 2025）

- 画法：MD 构象集合作为扩散模型的输入条件（ensemble-conditioned generation）→ 序列设计 → 过滤 → 实验。
- 借鉴点：动态信息可以放在**生成端之前**（作为条件）也可以放在**筛选端**（FlexPPI 作为重排器）；图上一旦出现"动态输入"就要说明是哪一端，避免审稿人以为 FlexPPI 参与生成。

### 9. ProteinDJ（Protein Science 2025）与 Ovo（bioRxiv 2025）

- 画法：模块化 pipeline 图：设计骨架 → 设计序列 → 结构验证 → 过滤，各模块可替换；ProteinDJ 强调 fold conditioning 模式与参数搜索，Ovo 强调端到端开源生态。
- 借鉴点：模块边界用"可替换插槽"画法，FlexPPI 就是一个可替换/可插拔的 reranker 模块。

### 10. nf-binder-design（Australian Protein Design Initiative，Nextflow，Zenodo 2025）

- 画法：工程 DAG——fold-and-dock（AF2-Multimer）、binder 设计、序列设计、结构预测、过滤指标、聚类，进程名 + 文件流。
- 借鉴点：工程执行图与论文主图分离；论文主图不画调度器、文件名和 GPU 细节。

### 11. PaperDraw（2026）—— 论文 pipeline 图自动生成系统

- 内容：专门讨论"科学论文 pipeline 图"的生成与设计原则：标准化尺寸与缩放、一致的颜色语义、分层/分阶段布局、箭头方向统一、元素对齐、避免装饰。
- 借鉴点：可以作为本仓库 pipeline 图的自检清单来源（与 `docs/main_figure_design.md` 的质量门一致）。

### 12. 综述类

- [AI-driven protein design, Nat. Rev. Bioeng. 2025](https://www.nature.com/articles/s44222-025-00349-8) 与 [ML-assisted enzyme engineering, ACS Cent. Sci. 2023](https://pubs.acs.org/doi/10.1021/acscentsci.3c01275) 的 DBTL（design–build–test–learn）循环图：大圆环 + 四象限，测试/学习阶段明确标注数据回流。适合作为"闭环"画法的权威参照。

## 三、对 FlexPPI pipeline 主图的直接启示

1. **主图采用"双入口汇合 + 分层漏斗 + 单一反馈弧"**：模式②（de novo / mutation 双泳道，BindCraft/RFdiffusion 风格）+ 模式①（漏斗 gate，BindCraft 风格）+ 模式③（一条实验反馈弧，EVOLVEpro 风格）。不要画成工程 DAG（模式⑥）。
2. **FlexPPI 的槽位：静态结构/亲和力过滤之后、MD 与湿实验之前**，与 BRD3 混合管线的 MD 层一致；图上要有一个明确的"这里是动力学感知重排层"的高亮框。
3. **实线 = 主数据流，虚线 = 结构自洽/MD 复核，弧线 = 实验反馈**（DynamicMPNN 语法）；颜色用 Okabe–Ito 并做文字/线型双重编码（与 `docs/main_figure_design.md` 一致）。
4. **候选数要么写真实数字，要么不写**（honest funnel）：BindCraft 有完整计数，本仓库如果没有可靠的 10⁴→10² 数据就不标数字。
5. **工具名写进框内**（BindCraft 风格）：RFdiffusion、ProteinMPNN、AF2-Multimer、Rosetta/FoldX、Pythia-PPI、FlexPPI 各自成框，避免 "Module A" 式空标签。
6. **动态信息出现的位置要自明**：KRAS 例子把动态放生成端，FlexPPI 放筛选端；图上必须让读者一眼看出 FlexPPI 不参与生成。
7. **FlexPPI 内部网络不进这张图**：架构图（`fig_flexppi_architecture`）已单独存在；pipeline 图只画输入/输出边界（RMSF / Sχ / SASA / Δ / uncertainty）。
8. **边界声明直接写在图上**：`orthogonal dynamics signal; not ΔG`（Pythia-PPI 与 FlexPPI 并列但不共用单位/颜色）。

## 四、检索局限与核查建议

- Exa MCP 免费额度 429 限流；pwsh 沙箱阻断出站连接（curl 全部 000），因此本次以检索摘要 + 图注页 + 二手解读交叉验证，未逐张打开原图。
- 上图描述中标注了每篇的官方链接；写论文时请直接打开各篇 Fig. 1 原图核对版式细节（漏斗计数、箭头类型、色带样式）。
- 项目已有 `docs/flexppi_design_screening_pipeline.md`（设计）与 `docs/main_figure_design.md`（架构主图规范），本报告只补充"别人怎么画"的证据面，不改动既有设计结论。
