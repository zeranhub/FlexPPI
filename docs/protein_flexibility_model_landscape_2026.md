# 蛋白质柔性预测模型全景清单与 FlexPPI 预训练评测路线

> 检索与核查截止日期：2026-07-23。本文所说的“全景清单”是按可操作定义进行的高召回检索：凡是以蛋白质序列、单个静态结构、实验密度图或预测构象为输入，并直接预测逐残基柔性、输出可换算为柔性的物理量、定位刚性/柔性区域，或生成可用于计算柔性的构象集合的方法，均纳入清单。配套 CSV 共收录 82 个模型、模型族或可复现协议。由于“蛋白质柔性”没有统一真值，论文会不断出现，早期服务器也会下线，所以互联网上不存在一个能够被证明穷尽的有限集合；本文通过任务定义、别名、年份、论文和代码交叉核对来尽量减少漏项，并把不能与 FlexPPI 公平比较的方法明确隔离。

## 一、结论先行

FlexPPI 预训练柔性头的目标不是泛泛的“柔性”，而是从一份静态目标链及可选伙伴结构出发，预测分子动力学（molecular dynamics，MD）轨迹中每个残基的绝对均方根涨落（root mean square fluctuation，RMSF），单位是 Å。对残基 \(i\)，先把每一帧轨迹按同一条链做刚体对齐，消除整条链的平移和转动，再计算

\[
\mathrm{RMSF}_i =
\sqrt{\frac{1}{T}\sum_{t=1}^{T}
\left\|\mathbf r_i(t)-\overline{\mathbf r_i}\right\|^2}.
\]

这里的 \(\mathbf r_i(t)\) 可以理解为第 \(t\) 帧中该残基 Cα 原子的三维坐标，\(\overline{\mathbf r_i}\) 是它在整段轨迹中的平均位置。数值越大，表示该残基围绕平均位置晃动的幅度越大。这个定义决定了最公平的核心对手是 **FlexPred-SVR、RMSF-net 的 PDB-only 版本、Flexpert、PEGASUS、LoopFlex、BackFlip、FastProtFlex/GDV、DynaProt、RocketSHP、SeqDance/ESMDance**，以及能够从静态结构产生逐残基涨落的 **ANM/GNM、CABS-flex、ENCoM、FRI**。其中只有一部分输出绝对 Å；FastProtFlex 和 SeqDance 等主要给相对柔性，BackFlip 使用局部对齐 RMSF，RMSF-net 原版还使用 cryo-EM 密度，因此它们不能未经转换就与 FlexPPI 共用绝对 MAE 排名。

晶体 B-factor、NMR \(S^2\)、无序概率、铰链位置、刚性常数和 AlphaFold 置信度也都包含柔性信息，但不是同一个监督目标。B-factor 混合了热运动、晶格接触、分辨率、精修模型和实验误差；\(S^2\) 主要刻画 N–H 键方向在较快时间尺度上的受限程度，而且数值越大通常越刚；无序预测回答的是“在生理条件下是否缺少稳定结构”，不是“围绕平均构象晃动多少 Å”；铰链预测只定位运动枢轴；pLDDT/PAE 是结构预测置信度而不是物理轨迹。因此，这些方法适合放在代理任务或廉价基线中，不能用来宣称 FlexPPI 已超过某个 MD-RMSF 模型。

目前最值得立即实现的主基线集合是：**常数/氨基酸类型/二级结构与 SASA 基线、反向 pLDDT、ProDy-ANM、ProDy-GNM、CABS-flex、FastProtFlex、FlexPred-SVR、PEGASUS、Flexpert-Seq、Flexpert-3D、BackFlip、RocketSHP、RMSF-net-PDB**。DynaProt 的方法价值很高，但本轮没有确认到与论文完全对应、可直接复现的官方代码，应作为“代码公开后补充”的基线。现有 FlexPPI 在内部验证集和 ATLAS/DynaRepo 上的结果可以用于调试，却不能作为最终泛化结论，因为 ATLAS 和 DynaRepo 已参与训练或域适配；论文级结论必须来自冻结的、经过序列与结构双重去泄漏的外部 MD 测试集。

## 二、任务地图：哪些模型可以和 FlexPPI 比什么

| 任务层 | 典型真值或输出 | 与 FlexPPI 的关系 | 正确比较方式 |
|---|---|---|---|
| A. 逐残基绝对 MD-RMSF | Cα RMSF，Å | 同任务，最重要 | 蛋白等权 Pearson/Spearman、MAE、RMSE、CCC、校准斜率；必须统一轨迹对齐与原子定义 |
| B. 逐残基相对柔性 | 每条蛋白内部标准化或归一化曲线 | 只比较曲线形状和排序 | 每蛋白 Spearman/Pearson、top-20% AP；禁止比较 Å-MAE |
| C. 晶体或 cryo-EM B-factor | 实验坐标不确定性/位移参数 | 代理任务 | 每条链标准化后比较 Spearman；单列实验条件与分辨率 |
| D. NMR \(S^2\)、RCI、异核 NOE | 快时间尺度局部运动或化学位移代理 | 相关但不同物理量 | 单独数据集、单独指标；不要把相关系数写成 RMSF 性能 |
| E. 刚性、铰链、构象开关、无序 | 分类、分段或运动模式 | 只能验证定位能力 | AUROC、AUPRC、F1、hinge distance 等分类/定位指标 |
| F. 构象集合或轨迹生成 | 多个三维构象、概率分布或时间序列 | 可从样本间接计算 RMSF | 固定样本数后重算同定义 RMSF，并额外比较成对 RMSD、PCA/Wasserstein、接触和碰撞 |
| G. 结合/突变引起的柔性变化 | \(\Delta\)RMSF、\(\Delta S^2\)、\(\Delta\Delta S\) | 更接近 FlexPPI 后续响应头，不是本次预训练主任务 | apo/holo 或 wild-type/mutant 配对评测，保留变化正负号 |

## 三、直接预测 MD-RMSF 或可直接换算 RMSF 的学习模型

| 模型 | 年份 | 输入 → 输出 | 训练/评测依据与特点 | 对 FlexPPI 的可比性与公开状态 |
|---|---:|---|---|---|
| [PredyFlexy](https://pmc.ncbi.nlm.nih.gov/articles/PMC3394303/) | 2012；2019 更新 | 序列/PSSM → 归一化柔性轮廓、三类柔性 | 把局部序列映射到结构原型，融合晶体 B-factor 和 MD-RMSF 信息 | 相对柔性基线；有网页服务，输出不是统一绝对 Å |
| [FlexPred-SVR（Kihara/Jamroz/Kolinski）](https://pmc.ncbi.nlm.nih.gov/articles/PMC3323684/) | 2012；2017 协议 | 单个静态 PDB → 逐残基 MD 涨落，Å | 以接触数等结构特征训练支持向量回归；592 个非同源蛋白，论文报告 Pearson 约 0.669、RMSE 约 1.04 Å | **高度可比**；[实验室页面/源码](https://kiharalab.org/flexPred/)。注意不要与 2008 年同名 FlexPred 混淆 |
| [DEFMap](https://www.nature.com/articles/s42256-020-00290-y) | 2021 | cryo-EM 密度 → 逐残基 RMSF | 三维卷积网络直接从局部分辨率与密度环境学习 MD 柔性 | 需要实验密度，属于额外信息轨道；[代码](https://github.com/clinfo/DEFMap) |
| [RMSF-net](https://www.nature.com/articles/s41467-024-49858-x) | 2024 | PDB 体素 + cryo-EM 体素；也有 PDB-only 消融 → RMSF，Å | 335 个体系、30 ns MD；双 U-Net 融合原子结构和密度，论文报告平均相关约 0.746 | 原版信息量高于 FlexPPI；**PDB-only 版本可比**；[代码](https://github.com/XintSong/RMSF-net) |
| [Flexpert-Seq / Flexpert-3D](https://openreview.net/forum?id=L238BAx0wP) | 2025 | 序列；序列+结构+ANM → RMSF，Å | 基于 ATLAS 三条 100 ns 重复，按拓扑划分；论文中 Seq/3D 的逐蛋白 Pearson 约 0.76/0.83 | **核心基线**；与 FlexPPI 共享 ATLAS 风险，必须用外部测试；[代码](https://github.com/KoubaPetr/Flexpert) |
| [PEGASUS](https://pmc.ncbi.nlm.nih.gov/articles/PMC12267886/) | 2025 | 四种蛋白语言模型表示 → RMSF、二面角标准差、平均 LDDT | ATLAS 上 ECOD fold-aware 交叉验证；RMSF Pearson 约 0.75、Spearman 约 0.66、MAE 约 0.82 Å；另测 CASP15 | **最重要的 sequence-only 基线**；[代码与 Docker](https://github.com/DSIMB/PEGASUS) |
| [LoopFlex](https://openreview.net/forum?id=nGezW51Mrj) | 2025 | ESM3 多模态序列/结构/二级结构/SASA 表示 → RMSF | 以 40% 序列簇划分，论文报告 pooled RMSE 约 1.184 Å、MAE 约 0.562 Å；还用于单环序列重设计 | 任务接近，但主要是 workshop 证据；本轮未找到可直接复现官方代码 |
| [BackFlip](https://arxiv.org/abs/2508.18211) | 2025 | 主链几何 → 局部对齐 RMSF，Å | 等变网络；论文测试 Pearson 约 0.80、MAE 约 0.17 Å，并做 de novo 骨架外测 | **需重定义标签后才公平**：局部对齐会删除部分域间运动；[代码](https://github.com/graeter-group/flips) |
| [FastProtFlex / GDV](https://pmc.ncbi.nlm.nih.gov/articles/PMC13099407/) | 2026 | 全原子静态结构的 graphlet degree vector → 相对逐残基柔性 | 15 维图元特征的线性模型；ATLAS、NMR 与大型 cryo-EM/MD 外测；速度快、可解释 | **相对排序核心基线**，不是绝对 Å；[代码](https://github.com/jure-praznikar/FastProtFlex)。本仓库已经完成接口与初测 |
| [DynaProt](https://openreview.net/forum?id=uKn9PdREBA) | ICLR 2026 | 序列+结构 → 每残基 \(3\times3\) 位移协方差与残基对协方差 | RMSF 可由协方差矩阵三个对角元素之和开平方得到；在 Flexpert 测试划分上报告逐蛋白 Pearson 中位数约 0.865 | 目标比 RMSF 更丰富，理论上高度可比；本轮未确认到完整官方复现仓库 |
| [RocketSHP](https://www.biorxiv.org/content/10.1101/2025.06.12.659353v1) | 2025 预印本 | 序列或结构、温度条件 → RMSF、广义相关、结构异质性 | 同时建模局部涨落和残基耦合，并用实验异核 NOE 做外部检查 | **值得纳入**；[代码](https://github.com/flatironinstitute/RocketSHP)，需冻结具体提交与权重 |
| [SeqDance / ESMDance](https://pmc.ncbi.nlm.nih.gov/articles/PMC11507661/) | 2024/2025 预印本 | 序列语言模型 → 归一化 RMSF、SASA、二级结构/二面角分布、接触和共运动 | 从约 6.5 万个模拟体系蒸馏多种动力学属性 | 可比较每蛋白排序，不能直接比较绝对 Å；[代码](https://github.com/ShenLab/SeqDance) |
| [ENSEMBITS](https://arxiv.org/abs/2605.13789) | 2026 预印本 | 单帧结构 → 动力学 token；MLP probe → RMSF | 在 mdCATH/MISATO 上按 CATH 超家族去重，探查单帧表示是否编码动态信息 | 更像表示探针而非校准预测器；可以作为研究性基线 |
| [Gaussian Integrals + Deep Learning](https://pmc.ncbi.nlm.nih.gov/articles/PMC12485704/) | 2026 | 静态结构 → 蛋白整体平均/慢模 RMSF | 预测的是每个蛋白的汇总柔性而非逐残基曲线 | 不进入逐残基主榜，只可做蛋白级辅助任务 |
| [ProFlex](https://www.nature.com/articles/s41467-025-64103-9) | 2025 | 结构描述/语言桥接 → NMA 相关柔性量 | 主要输出全局、最小、中位和平均 RMSF及模式描述 | 不是逐残基直接对手，只作方法学邻域 |
| [用结构预测不确定性估计柔性](https://www.mlsb.io/papers_2024/Estimating_protein_flexibility_via_uncertainty_quantification_of_structure_prediction_models.pdf) | 2024 | 结构预测模型的坐标分布 → 位置协方差/柔性 | 使用扩展 AlphaFold 结构模块学习 aleatoric uncertainty，验证规模较小且以短环为主 | 概念相关，暂不适合作为全蛋白主基线 |

这里最容易造成误判的有三处。第一，BackFlip 的“局部 RMSF”先在残基附近做局部对齐，所以一个刚性结构域相对另一个结构域摆动时，它可能给出较小值，而 FlexPPI 的链级对齐会保留这种内部域运动；必须让两种模型在同一套局部标签和链级标签上各跑一次。第二，FastProtFlex 先在蛋白内部标准化后预测“哪里更柔”，因而即使相关系数很高，也不能说明它知道一条蛋白整体应当是 0.7 Å 还是 2.0 Å。第三，RMSF-net 原版获得了 cryo-EM 密度，FlexPPI 没有这份实验信息；公平主榜应使用 PDB-only 分支，而完整 RMSF-net 只能进入“附加实验信息”榜。

## 四、晶体 B-factor、NMR \(S^2\) 与早期逐残基柔性模型

这些方法构成了柔性预测的历史主线，仍适合做跨域代理评测和低成本 sanity check。它们不应与 MD-RMSF 的绝对 Å 主榜混排。

| 模型 | 年份 | 输入 → 输出 | 价值与限制 |
|---|---:|---|---|
| [Yuan SVR：Prediction of protein B-factor profiles](https://pubmed.ncbi.nlm.nih.gov/15645415/) | 2005 | 序列 → B-factor 轮廓 | 早期支持向量回归，论文报告相关约 0.53；实验 B-factor 与 MD-RMSF 不等价 |
| [PROFbval](https://pubmed.ncbi.nlm.nih.gov/16455751/) | 2006 | 序列 → 归一化 B-value、柔/刚分类 | 代表性的纯序列 B-factor 服务器；旧服务可用性需要实际测试 |
| [FlexRP](https://pmc.ncbi.nlm.nih.gov/articles/PMC1863424/) | 2007 | 序列中间隔氨基酸对 → 多结构中柔/刚残基 | 预测构象变化中参与运动的二分类，论文准确率约 79.5%；不输出振幅 |
| [Trott–Siggers–Rost–Palmer](https://pubmed.ncbi.nlm.nih.gov/18313957/) | 2008 | 三维结构描述 → NMR \(S^2\) | 16 个蛋白上的神经网络，论文 Pearson 超过 0.70；\(S^2\) 越高通常越刚，方向与 RMSF 相反 |
| [FlexPred（Kuznetsov & McDuffie）](https://pmc.ncbi.nlm.nih.gov/articles/PMC2639688/) | 2008 | 序列或结构 → 构象开关残基 | 预测两个构象间发生切换的位置；与 2012 年 FlexPred-SVR 同名但任务不同 |
| [CFP / GLP](https://pmc.ncbi.nlm.nih.gov/articles/PMC2859570/) | 2009 | 序列局部倾向 → 柔性片段 | generalized local propensity 的阈值分段方法；适合做简单可解释分类基线 |
| Pan–Shen 两阶段 SVR+RF | 2009 | 序列 → 稳健 B-factor 轮廓 | 支持向量回归后接随机森林，属于 B-factor 预测线 |
| [Hirose 等随机森林](https://pmc.ncbi.nlm.nih.gov/articles/PMC3245509/) | 2010 | 序列+预测二级结构 → NMA 内部/外部运动 | 用随机森林蒸馏正常模态运动，曾用于 PPI 位点分析；输出含义不是绝对 RMSF |
| [Sonavane 等 SVR](https://pmc.ncbi.nlm.nih.gov/articles/PMC3569600/) | 2013 | 滑动序列窗口 → 归一化 B-value | 论文平均相关约 0.47；适合重现早期序列特征上限 |
| [DynaMine](https://www.nature.com/articles/ncomms3741) | 2013 | 序列 → 主链 N–H \(S^2\)，0–1 绝对标度 | 速度快、解释清楚；预测的是快速局部主链动力学，不是 Å；[服务](https://bio2byte.be/dynamine/) |
| [RCI](https://pubmed.ncbi.nlm.nih.gov/17406296/) | 2007 | **实验 NMR 化学位移** → RMSF/\(S^2\) 代理 | 使用了实验观测，不能作为“仅凭序列/结构”公平基线；可作为实验上限参考 |
| [FLEXc](https://pmc.ncbi.nlm.nih.gov/articles/PMC5009531/) | 2016 | 序列上下文+预测结构特征 → B-factor 分类/轮廓 | 代表传统机器学习与特征工程阶段 |
| [Blind B-factor prediction / MWCG](https://pmc.ncbi.nlm.nih.gov/articles/PMC6192411/) | 2018 | 原子结构、加权彩色图特征及晶体元数据 → B-factor | 随机森林、梯度提升与卷积模型；结构与实验元数据强，信息条件需单列 |
| [MEDUSA](https://pubmed.ncbi.nlm.nih.gov/33972018/) | 2021 | 序列进化/理化特征 → 2/3/5 类 B-factor | 卷积分类器；[代码](https://github.com/DSIMB/medusa)，只能参加分类代理任务 |
| DihProFle | 2022 | 二面角与序列特征 → DNA-binding/冠状病毒蛋白柔性 | SVM 与全连接网络组合，应用域较窄；需避免把专域结果外推到一般蛋白 |
| [Pandey/Keten LSTM](https://pmc.ncbi.nlm.nih.gov/articles/PMC10499862/) | 2023 | 序列 → Cα B-factor | 在 2,442 个蛋白上评测的 LSTM；[代码](https://github.com/keten-group/B_factor_codes_repository) |
| [BiLSTM B-factor profile](https://chemrxiv.org/engage/api-gateway/chemrxiv/assets/orp/resource/item/6502b58fb6ab98a41c6bccf1/original/prediction-of-protein-b-factor-profiles-based-on-bidirectional-long-short-term-memory-network.pdf) | 2023/2024 | 序列特征 → B-factor 轮廓 | 预印本报告平均 PCC 约 0.71；证据强度低于正式发表和外部盲测 |
| [Point-cloud CNN B-factor](https://www.sciengine.com/PBB/doi/10.16476/j.pibb.2020.0383) | 2022 | 原子点云结构 → B-factor | PointNet 风格几何网络；属于结构到实验位移参数 |
| [Meta-GNN](https://arxiv.org/abs/2408.12519) | 2024 | 原子图结构 → B-factor | 在超过 4,000 个蛋白、约 1,700 万原子上测试，论文报告 Pearson 约 0.71 |
| [Persistent Sheaf Laplacian](https://pmc.ncbi.nlm.nih.gov/articles/PMC11844605/) | 2024/2025 | 多尺度拓扑和几何 → B-factor | 数学物理特征模型，论文称相对 GNM 与既有盲预测有提升；仍是 B-factor 任务 |

如果要建立一张晶体代理榜，应先在每条链内部对预测值和 B-factor 分别做 z-score，再算链等权 Spearman，并按分辨率、晶体接触、二级结构、N/C 端和缺失残基分层。绝不能把 B-factor 通过一个全局常数换算成 RMSF 后，直接与 MD 的 Å-MAE 混排；理想各向同性热运动下虽有 \(B = 8\pi^2\langle u^2\rangle/3\)，但实际 PDB 的 B-factor 并不只含热运动。

## 五、基于物理、网络与低成本采样的结构模型

| 模型族/工具 | 输入 → 输出 | 与 FlexPPI 的用途 |
|---|---|---|
| 全原子 NMA（normal mode analysis） | 势能在单个最低能构象附近的二阶近似 → 正常模态和均方涨落 | 物理基线；对初始结构和能量最小化敏感 |
| [GNM](https://www.bahargroup.org/prody/manual/reference/dynamics/index.html) | Cα 接触网络 → 各残基各向同性相对涨落 | 便宜、可解释，适合主榜相对排序 |
| [ANM](https://prody.readthedocs.io/en/latest/apps/prody/anm.html) | Cα 弹性网络 → 三维各向异性位移和模态 | 比 GNM 保留方向信息；需要统一模态数与整体缩放 |
| Tirion ENM、β-Gaussian model | 简化弹簧网络 → 模态/涨落 | 可作为 ANM/GNM 的历史变体，不必全部放进最终主表 |
| [ENCoM](https://pmc.ncbi.nlm.nih.gov/articles/PMC4489264/) | 序列敏感弹性网络 → B-factor、模态、突变柔性变化 | 对突变和侧链化学更敏感；绝对尺度仍需校准 |
| [FRI / fFRI](https://pmc.ncbi.nlm.nih.gov/articles/PMC4098052/)、mFRI、[gFRI](https://arxiv.org/abs/1609.04306) | 距离核与多尺度刚性指数 → 相对 B-factor/柔性 | 无需矩阵对角化，速度快；很适合大型复合物低成本基线 |
| [CABS-flex 2.0](https://pmc.ncbi.nlm.nih.gov/articles/PMC6031000/) | 静态 PDB → 粗粒化 Monte Carlo 构象集合与 RMSF | **核心结构物理基线**；可处理复合物，但要固定随机种子、采样步数和重复数 |
| [FlexServ](https://academic.oup.com/bioinformatics/article/25/13/1709/197147) | PDB → DMD、NMA/ANM 或布朗动力学结果 | 是多种运动方法的统一服务层，不是一个单独学习模型 |
| [ProPHet](https://bioserv.rpbs.univ-paris-diderot.fr/services/ProPHet/) | 粗粒化弹性网络+布朗动力学 → 每残基有效力常数/刚性 | 可做刚性轮廓代理；数值越大越刚，需反向后比较排序 |
| [HingeProt](https://pubmed.ncbi.nlm.nih.gov/17847101/) | GNM/ANM → 刚性结构域与铰链残基 | 适合验证 FlexPPI 高梯度区域是否靠近铰链，不参加 RMSF 回归 |
| FIRST/FRODA | 几何约束网络+低频运动 → 刚性簇和可达构象 | 适合构象可达性与刚性分析，输出不是直接 RMSF |
| ElNémo、WEBnm@、iMODS | 在线 NMA/ENM → 模态、涨落和构象动画 | 方便人工核查，但版本与服务器状态不够稳定；正式 benchmark 更应使用 ProDy 脚本化复现 |

物理模型的数值标度通常依赖温度、弹簧常数、模态截断和经验缩放。最稳妥的做法是同时报告两个轨道：原始输出只算每蛋白 Spearman；另在训练集内拟合一个全局线性校准 \(y=a x+b\)，冻结 \(a,b\) 后在测试集计算 MAE。不能为每个测试蛋白单独拟合缩放，因为那相当于提前知道了该蛋白的真实振幅。

## 六、无序、铰链、构象切换与突变响应模型

| 模型/类别 | 实际回答的问题 | 是否进入 FlexPPI 预训练主榜 |
|---|---|---|
| [IUPred2A](https://iupred2a.elte.hu/) | 序列片段能否形成足够稳定的相互作用，即内在无序概率 | 否；可作为高柔性区域的弱代理 |
| [SPOT-Disorder2](https://sparks-lab.org/server/spot-disorder2/) | 深度学习逐残基无序概率 | 否；单列分类任务 |
| [flDPnn](https://pmc.ncbi.nlm.nih.gov/articles/PMC8295265/) | 无序及其功能位点 | 否；功能性无序与轨迹涨落不同 |
| DISOPRED3、MobiDB-lite、metapredict、ESpritz、PONDR、NetSurfP-3.0 | 不同来源的无序或可及性预测 | 否；最多作为廉价输入/代理基线 |
| 2008 FlexPred、FlexRP、CFP/GLP | 构象开关、柔/刚残基或柔性片段 | 否；可在有二状态结构的数据上做定位评测 |
| HingeProt | 铰链与准刚性结构域 | 否；可测到高 \(|\nabla \mathrm{RMSF}|\) 区域的距离 |
| [DynaMut2](https://pmc.ncbi.nlm.nih.gov/articles/PMC7737773/) | 突变引起的稳定性和振动熵变化 | 否；更适合 FlexPPI 的突变/响应扩展任务 |
| ENCoM mutation、ProPHet mutation 分析 | 突变前后模式、柔性或刚性变化 | 否；放入 \(\Delta\) 柔性轨道 |

无序残基往往也很灵活，但“有高 RMSF”并不推出“无序”，例如稳定折叠蛋白的表面环、结构域铰链或末端可以在保持整体折叠的同时具有高 RMSF；反过来，实验条件、伴侣结合或有限轨迹也可能使无序区域在某段模拟中看起来不特别活跃。因此无序模型适合检查召回率，不适合给 FlexPPI 的绝对校准背书。

## 七、从生成构象集合间接得到柔性的模型

下列模型不直接输出一条 RMSF 曲线，而是生成多构象或轨迹；对它们采样后可以按与真值相同的残基映射、链级对齐和 Cα 定义重新计算 RMSF。它们计算成本通常高得多，且柔性只是生成质量的一部分，所以应单独建“ensemble track”，不能只用 RMSF 一项判定优劣。

| 模型族 | 主要机制/输出 | 评测注意事项 |
|---|---|---|
| AlphaFold2 dropout、MSA subsampling | 通过随机失活或减少 MSA 约束获得多构象 | 很便宜但样本分布不是显式动力学分布；固定 MSA 深度、随机种子数与回收次数 |
| AF-Cluster、SPEACH_AF、AFsample/AFsample2 | 聚类/扰动 MSA 或加强采样以寻找替代状态 | 适合离散状态，未必复现局部热涨落 |
| [AlphaFlow / ESMFlow](https://arxiv.org/abs/2402.04845)、[AlphaFlow-Lit](https://arxiv.org/abs/2407.12053) | 把结构预测模型转为流模型，学习 MD/实验构象分布 | 需固定每蛋白样本数；评估覆盖与精度的权衡 |
| [BioEmu](https://github.com/microsoft/bioemu) | 从序列快速采样平衡构象集合 | 可用其公开 benchmark 工具；除 RMSF 外必须测自由能盆地和接触 |
| ConfDiff、DiG、Str2Str、ESMDiff | 不同形式的结构扩散/几何生成 | 统一输入信息、采样预算、后处理和能量最小化 |
| MDGen | 生成 MD 风格轨迹或动力学序列 | 时间连续性与动力学统计要单独验证 |
| [TokenMD](https://arxiv.org/abs/2511.10056) | 通过结构 token 同义替换生成构象集合 | 逐蛋白 RMSF 可较好，但跨蛋白绝对振幅可能较弱；须同时看全局 RMSF |
| DynaProt sampler | 从预测的联合高斯位移协方差采样构象 | 适于局部连续涨落；高斯假设可能漏掉多峰跃迁 |
| [DynaMode](https://arxiv.org/abs/2607.04134) | 在频域生成有时间顺序的轨迹 | 2026-07 新预印本；速度快但仍需原子碰撞、能量和时间统计审计 |
| Mac-Diff、ConfRover、TEMPO | 2025–2026 年的宏观构象/轨迹生成方法 | 版本很新，应冻结论文版本和官方权重，并以复现状态而非论文最优值为准 |
| IDPFold、AFflecto、idpGAN、DynamICE | 内在无序蛋白或柔性区的构象集合 | 适合 IDP 专门测试，不宜直接代表折叠蛋白主榜 |

公平的 ensemble track 应给每种生成器完全相同的样本数，例如每个蛋白 250 个构象，使用同一套对齐与残基掩码，再报告逐蛋白 RMSF Pearson/Spearman/MAE、全局振幅相关、成对 RMSD 分布、PCA 子空间或 2-Wasserstein 距离、弱/瞬时接触、二面角分布、原子碰撞、推理时间和失败率。只从 1,000 个样本里挑出一个最像真值的构象，会把“采样预算大”误写成“条件控制准确”。

## 八、为 FlexPPI 建议的冻结评测方案

### 8.1 三张互不混排的榜

**主榜 A：静态结构 → 绝对 MD-RMSF。** 输入一份目标链结构；若模型支持复合物上下文，再同时给完整伙伴结构。输出必须映射到同一残基集合，以 Cα、Å、链级刚体对齐的 RMSF 为真值。主指标是蛋白等权 Spearman、蛋白等权 Pearson、MAE、RMSE、CCC（concordance correlation coefficient，一致性相关系数）、校准斜率和截距；补充最柔 20% 残基的 AUPRC 与 Precision@K。绝对输出模型参加全套指标，相对输出模型只参加相关、排序和 top-k 指标。

**代理榜 B：实验柔性。** 晶体 B-factor、cryo-EM 多构象/密度异质性、NMR \(S^2\)、异核 NOE 和 NMR ensemble 分开建表。每一种实验标签都在链内标准化，按蛋白等权汇总，并明确其不是 MD-RMSF 真值。TCR3D 应继续保留在这一榜，而不是写成独立 MD 测试。

**生成榜 C：构象集合/轨迹。** 给 AlphaFlow/BioEmu/TokenMD/DynaProt sampler/DynaMode 等相同采样预算，先从生成样本计算同定义 RMSF，再额外评价分布、多模态覆盖、接触、碰撞和速度。FlexPPI 本身不生成构象，因此它只能作为“由单结构直接预测 RMSF”的速度与曲线参照，不能与生成器争夺完整生成质量名次。

### 8.2 测试集与泄漏控制

当前配置把 ATLAS 与 DynaRepo 用于预训练或域适配，现有验证集也参与模型选择，所以这些数据上的性能只能称为 internal validation 或 in-domain diagnostic。外部测试应从 mdCATH、CASP15 公开 MD、BioEmu benchmark 的 MD 体系或自行新跑的重复 MD 中构建，随后把每个候选测试链同时与 FlexPPI 的全部训练链、所有公开基线的已知训练集做 **MMseqs2 序列搜索**和 **Foldseek/CATH/ECOD 结构拓扑搜索**。建议至少剔除序列一致性高于 30% 且覆盖度足够的近同源，并单列“新序列但旧折叠”和“新折叠”两个难度层；如果只能执行其中一种去泄漏，结论必须降格。

MD 真值本身有随机噪声。每个测试体系最好有至少三条独立初速度的重复轨迹，在相同温度、力场、水模型、离子条件、平衡流程和生产时长下运行。训练标签取重复均值，评测时同时计算“重复 A 预测重复 B”的相关与误差，把它作为经验噪声上限；如果两个独立 MD 重复之间只有 0.75 的相关，就不应把模型 0.78 与另一个模型 0.80 的差异解释成确定性进步。

### 8.3 输入公平与复合物上下文

建议把主榜再拆为四个信息条件：sequence-only、target-structure-only、target+partner structure、structure+experimental density。FlexPPI 与 Flexpert-3D、BackFlip、FastProtFlex、ANM/GNM、CABS-flex应在 target-structure-only 条件比较；FlexPPI 的伙伴感知版本进入 target+partner 条件；RMSF-net 完整版进入 density 条件。对于不支持多链的模型，应输入剥离伙伴后的同一目标链，不能把它的缺失上下文隐藏起来。每个体系同时跑 FlexPPI 的“仅目标链”和“目标链+伙伴”版本，才能量化伙伴结构到底带来多少增益。

### 8.4 必须统一的标签细节

1. 固定原子为 Cα，除非另开全原子/主链重原子轨道；不能把原子 RMSF 与残基质心 RMSF 混排。
2. 固定链级 Kabsch 刚体对齐，并额外生成 local-window 对齐标签供 BackFlip 专门比较。Kabsch 对齐就是旋转和平移整条链，让对应 Cα 的总平方距离最小，它只去掉整条链在盒子里的漂移，不删除链内部结构域的相对运动。
3. 对所有模型使用相同的缺失残基、插入码、替代构象、非标准残基和链映射规则；指标只在共同有效残基交集上计算，同时报告覆盖率。
4. termini、长环、无序段和结构缺口不能静默删除，应分别报告，因为它们既最难也最可能决定“柔性预测”是否有用。
5. 物理模型若需经验缩放，只能在训练/校准集拟合一次全局缩放并冻结；禁止按测试蛋白单独缩放。

### 8.5 主指标、统计检验与失败分析

主表应先在每个蛋白内计算 Pearson、Spearman、MAE、RMSE、CCC，再让每个蛋白等权平均；所有残基混在一起的 pooled 指标只能放补充材料，否则长蛋白和高柔末端会支配结论。使用按蛋白配对 bootstrap 给出 95% 置信区间，并用配对 Wilcoxon 或置换检验比较 FlexPPI 与每个基线；同时给效果量，而不是只给 p 值。分层至少包括 helix/sheet/loop、N/C 端、埋藏/暴露、界面/近界面/远端、单域/多域、长度、膜蛋白、是否含配体、实验结构/预测结构以及输入结构质量。还要报告每种模型的参数量、CPU/GPU、单蛋白耗时、峰值内存、安装失败和推理失败比例。

## 九、推荐的实际执行优先级

| 优先级 | 模型 | 原因 | 进入哪张榜 |
|---|---|---|---|
| P0 | 常数、氨基酸类型、二级结构、SASA、反向 pLDDT | 检查模型是否只学到暴露度、末端或置信度 | A 的全部/相对指标 |
| P0 | ProDy-GNM、ProDy-ANM | 最稳定、廉价的物理基线 | A；原始相对值+全局校准值 |
| P0 | CABS-flex | 低成本采样且可处理复合物 | A 与轻量 ensemble |
| P0 | FastProtFlex | 已在本仓库接通、速度快、ATLAS 上强 | A 的相对指标 |
| P0 | PEGASUS | 强 sequence-only、Docker 可复现 | A 全套 |
| P0 | Flexpert-Seq / 3D | 与 ATLAS-RMSF 主任务最接近 | A 全套；外部集决定结论 |
| P0 | FlexPred-SVR | 经典静态结构到绝对 Å 基线 | A 全套 |
| P1 | BackFlip | 现代几何模型且公开代码 | A；链级和局部对齐双标签 |
| P1 | RocketSHP | 同时预测涨落与耦合，可扩展 FlexPPI 分析 | A，冻结预印本版本 |
| P1 | RMSF-net-PDB | 结构体素深度网络 | A；完整模型另进 density 榜 |
| P2 | PredyFlexy、DynaMine、MEDUSA、Pandey LSTM | 历史/实验代理覆盖 | B 或 A 的相对 sanity check |
| P2 | ENCoM、FRI、ProPHet | 扩展物理与刚性解释 | A 相对榜、B 或刚性辅助分析 |
| 等待 | DynaProt | 预测对象丰富，论文结果强 | 官方可复现实现确认后进 A/C |
| 研究性 | AlphaFlow、BioEmu、TokenMD、DynaMode 等 | 成本高且任务更广 | C，不与直接回归混排 |

第一轮不需要把清单里的每个旧服务器都跑一遍。最有信息量的最小实验是：在一个真正冻结的外部 MD 测试集上，同时跑 `mean/pLDDT、GNM、ANM、CABS-flex、FastProtFlex、FlexPred-SVR、PEGASUS、Flexpert-Seq、Flexpert-3D、FlexPPI`，并严格区分绝对与相对指标。第二轮再加入 BackFlip、RocketSHP、RMSF-net-PDB 和可复现的 DynaProt。B-factor、无序、NMR 与生成模型应该在主结论稳定后作为外部证据和扩展实验加入。

## 十、FlexPPI 当前结果应怎样表述

仓库现有记录显示，FlexPPI 预训练模型在 276 个内部验证体系上的 RMSF 蛋白等权 Spearman 为 0.811、MAE 为 0.535 Å、pooled Pearson 为 0.701、最柔残基 AP 为 0.786；在当前共同 ATLAS 子集上，FlexPPI 预训练的 Spearman 约 0.8205，FastProtFlex 约 0.7984；在 DynaRepo 共同子集上，FlexPPI 静态预训练约 0.7892，FastProtFlex 约 0.7507。这些数字说明实现有有效信号，也支持继续做外部评测，但 **不能称为独立测试或 state-of-the-art**：ATLAS 全部 1,938 个体系已经进入预训练流程，DynaRepo 又参与域适配；FastProtFlex 本身也在 ATLAS 训练，而且它输出的是蛋白内相对柔性。当前最严谨的说法是“in-domain diagnostic 中表现出较强的逐残基排序能力”，等待外部冻结测试后再讨论绝对泛化。

## 十一、检索覆盖、命名冲突和可复现性说明

本次检索交叉使用了论文题名/摘要检索、官方论文页、PubMed/PMC、arXiv/OpenReview、实验室服务器页和 GitHub 仓库，关键词覆盖 `protein flexibility prediction`、`per-residue RMSF prediction`、`B-factor prediction`、`order parameter S2`、`normal mode flexibility`、`rigidity/hinge/disorder`、`conformational ensemble generation` 等组合。机器聚合检索存在明显噪声：OpenAlex 客户端因本地 Python 类型兼容错误不可用；Semantic Scholar 曾触发限流与 SSL EOF；一次包含非 ASCII 作者名的输出因 Windows GBK 编码中断；`agent-reach` 命令和 Exa MCP 后端在当前环境不可用；GitHub CLI 未登录。因此本文对核心条目均回到论文、官方仓库或实验室页面人工核查，而不是直接接受聚合搜索排序。完整机器检索审计和报错保存在历史合并快照 `allinone.md`；该文件已于 2026-09-05 移入项目外清理归档，位置与恢复方法见 [备份说明](../BACKUP_MANIFEST.md)。

还必须警惕同名：**FlexPred** 至少指 2008 年的构象开关残基分类器和 2012 年 Kihara 团队的静态结构到 MD 涨落 SVR；本文统一写成“FlexPred（Kuznetsov）”与“FlexPred-SVR（Kihara）”。“RMSF prediction”也可能指在某个配体或突变 MD 研究中用回归预测一个全蛋白汇总数，而不是逐残基曲线；这类工作只有在输出粒度与 FlexPPI 一致时才进入直接模型表。

与本文配套的可排序清单见 [`flexibility_model_inventory.csv`](flexibility_model_inventory.csv)。它包含模型、年份、类别、输入、输出、绝对标度、论文、代码、与 FlexPPI 的可比等级和复现备注，后续可以直接用作基线实现看板。已有的更聚焦 RMSF 与生成模型评测细节仍保留在 [`rmsf_prediction_and_generation_literature.md`](rmsf_prediction_and_generation_literature.md)，现有 FlexPPI 数据泄漏边界与指标定义见 [`evaluation_protocol.md`](evaluation_protocol.md)。
