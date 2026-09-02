# FlexPPI 在 Zenodo 6638504 严格 apo/holo PPI 数据集上的实验报告

本报告记录 2026-08-25 在本机完成的可复现实验。数据来自 Zenodo record 6638504，对应 Martin 与 Frezza 2022 年的 8 个经典二元蛋白复合物研究；原论文使用 Amber99SB-ILDN、TIP3P 水、0.15 M KCl、300 K 和 1 bar，在完成平衡后为 1AK4 与 1PVH 生成 500 ns 生产轨迹，其余体系生成 300 ns 生产轨迹，总计 8.4 μs。论文原文和数据记录分别见 [Frontiers 论文](https://www.frontiersin.org/journals/molecular-biosciences/articles/10.3389/fmolb.2022.970109/full) 与 [Zenodo 6638504](https://zenodo.org/records/6638504)。本次不是简单复述原论文，而是把每个复合物的游离 A、游离 B 和结合态轨迹转换成 FlexPPI 可用的严格配对标签，用冻结的 FlexPPI checkpoint 做状态预测、结合响应预测、伙伴删除消融、简单结构基线和 leave-one-complex-out 校准，再检查能否得到或验证可复核的生物学洞见。

## 1. 数据适配与质量控制

真实数据目录是 `E:\data\protein_md_datasets\zenodo_6638504_ppi_md`。8 个体系为 1AK4、1AY7、1BRS、1EMV、1GCQ、1PVH、2OOB、3SGB；每个体系都有 `A/start.gro + A/md.xtc`、`B/start.gro + B/md.xtc` 和 `complex/without_water/start.gro + md.xtc`。全部轨迹每 40 ps 保存一帧，300 ns 对应 7,501 帧，500 ns 对应 12,501 帧。适配器用游离 A、B 的精确原子数确定结合态链边界，验证结合态原子数严格等于 A+B；它保留 XTC 所需的原子顺序，把 `CYX/ASH/GLH/HID/HIE/HIP` 等力场名称规范成标准氨基酸名，但保留 `ACE/NME` 端基原子并排除它们作为模型残基。所有 16 个链侧的 apo→holo 序列身份和覆盖均为 100%，最终得到 24 个有标签状态、16 个 apo/holo 链侧配对和 1,759 个逐残基对应。

标签按 100 ns 不重叠窗口计算，因此两个 500 ns 体系的每个状态有 5 窗，其余状态各有 3 窗，总计 84 个状态窗口。RMSF 是逐链刚体拟合后的 Cα 涨落，这避免把结合伙伴的整体相对转动误当成蛋白内部柔性；χ 熵根据每个可定义侧链二面角的三态旋转异构体联合分布计算。窗口均值作为标签，窗口间方差保留为不确定性。这里与原论文有一个需要明确的口径差异：原论文报告的是逐链拟合后对主链原子 RMSF 取残基平均，本次用的是与现有 FlexPPI/DynaRepo 管线一致的 Cα RMSF，所以可以验证趋势和模型能力，但不应声称逐数值复现论文 Supplementary Figure S6。SASA 没有在全轨迹逐帧计算，因为当前本地 Shrake–Rupley 实现对 20 多万帧成本很高；因此本报告没有把模型 SASA 输出当成 MD 真值验证。

## 2. 泄漏边界与统计单位

已有 `artifacts/finetune/best.pt` 的成对微调清单没有这 16 个 apo/holo 配对，但它继承的状态预训练 checkpoint 曾见过 DynaRepo 中 1AK4、1AY7、1EMV、1GCQ、1PVH 的复合物。因此全部 8 个体系不是严格的“完全未见外部集”。本报告把这 5 个体系标成 `exact_pdb_overlap`，把 1BRS、2OOB、3SGB 标成 `exact_pdb_unseen`，核心泛化结果同时报告后 3 个体系。该审计只排除了精确 PDB ID，没有做 MMseqs2 或结构同源簇排除，所以“未见”不能等同于“无同源泄漏”。所有置信区间以 PDB 复合物为独立单位：先在 A/B 两个链侧内分别计算指标，再在同一 PDB 内等权平均，最后对 8 个复合物 bootstrap 2,000 次；3 个未见 PDB 的区间只能看成探索性范围，不能当作高功效统计检验。

## 3. 已运行实验矩阵

| 实验 | 问题 | 模型/对照 | 主要统计量 |
|---|---|---|---|
| E0 数据 QC | GRO/XTC 是否能严格配对 | 原子数、序列、帧数、窗口方差 | 24/24 状态成功，16/16 配对成功 |
| E1 绝对状态 | 静态结构能否预测 apo/holo 柔性画像 | FlexPPI、GNM、反接触数、端点先验 | 每复合物 Spearman、MAE、top-20% |
| E2 结合响应 | 能否预测 holo−apo 的逐残基变化 | 响应头、状态相减、零变化、界面距离先验 | ΔRMSF/Δχ 的 Spearman、MAE、方向、top-20% |
| E3 伙伴删除 | 信号是否真的使用对方链 | 完整复合物 vs 保持目标坐标不变并删除伙伴 | 结合态画像变化与响应变化 |
| E4 LOCO 校准 | 排序信号能否转成可用幅度 | 每次留出 1 个 PDB，在其余 7 个拟合线性校准；另测距离与混合模型 | 留出 PDB 的 MAE、Spearman |
| E5 空间生物学 | 结合效应是否局域在界面 | 0–5 Å、5–12 Å、>12 Å 重原子距离分层 | 平均变化、增/减比例、窗口不确定性 |
| E6 原论文分层 | 模型是否受界面重排影响 | 稳定界面 1GCQ/2OOB、重定向 1AK4/3SGB、其余多亚态 | 每复合物 ΔRMSF Spearman |

## 4. 模型结果：强项是绝对画像和响应排序，不是幅度

FlexPPI 对结合态绝对画像的预测较强：8 个复合物等权的 RMSF Spearman 为 0.683，95% bootstrap 区间为 [0.591, 0.761]；χ 熵为 0.870，[0.829, 0.904]。精确 PDB 未见的 3 个体系上分别为 0.686 和 0.844，说明至少按精确 ID 分层没有出现性能崩溃。不过这不能证明 FlexPPI 在 RMSF 上优于简单结构模型：结合态反接触数基线为 0.692，GNM 为 0.669；游离态 FlexPPI 为 0.685，而 GNM 达到 0.783。更准确的结论是 FlexPPI 能给出不错的绝对柔性排序，并且 χ 熵这一侧链描述符目前最有优势，但 RMSF 排序在该小数据集上与经典弹性网络/堆积度基线相当或更弱。

对结合诱导 ΔRMSF，FlexPPI 响应头的 8 复合物宏平均 Spearman 为 0.335，[0.227, 0.453]；精确 PDB 未见 3 体系为 0.337，[0.186, 0.612]。模型识别真实变化绝对值最大的 top-20% 残基时，平均 precision 为 0.484，而随机期望为 0.2，说明“哪里更可能明显改变”的排序有实质信号。伙伴删除提供了很强的机制性检查：保持目标链的结合态坐标完全不变，只删除对方链后，响应 Spearman 降到 −0.003，说明响应排序确实依赖显式伙伴上下文，不只是比较两个单链几何状态。

但是，简单界面距离先验 `−exp(−distance/5)` 的 ΔRMSF Spearman 为 0.433，[0.315, 0.546]，高于 FlexPPI；未见 3 体系上为 0.494。FlexPPI 输出本身与伙伴距离的宏平均 Spearman 为 0.438，几乎等于真实 ΔRMSF 与距离的 0.433，这说明模型学到的主要规律之一就是“越靠近伙伴越倾向刚化”。它比纯距离模型多出的非局部信息目前还不足以稳定超过这个基线。数值幅度也未校准：复合物等权 MAE 为 0.177 Å，零变化基线为 0.170 Å；仓库原生的残基池化口径同样是 0.163 Å 对 0.154 Å。LOCO 线性校准把 FlexPPI MAE 降到 0.164 Å，未见 3 体系从 0.197 降到 0.179 Å，略优于相应的零变化 0.185 Å；但校准距离模型仍略好，混合 FlexPPI+距离的全 8 体系排序为 0.441，在未见 3 体系却没有超过纯距离。因此可用结论是“响应头含有可校准的弱幅度信号”，而不是“已经能准确量化 ΔRMSF”。

Δχ 熵的结论更保守：响应头 Spearman 为 0.154，界面距离先验为 0.243；响应头 MAE 为 0.157 nats，零变化为 0.133 nats，并且 8/8 复合物都更差，配对 Wilcoxon `p=0.0078`。校准后响应头 MAE 仍为 0.140，未超过零变化。因此当前 checkpoint 的绝对 χ 熵画像很强，但结合诱导 Δχ 不应作为可靠的定量预测。一个有价值的伙伴上下文现象是：完整复合物相对删伙伴输入对结合态 χ 熵排序平均提高 0.014，8/8 体系方向一致，`p=0.0078`；RMSF 平均提高 0.023，但只有 6/8 为正，`p=0.3125`。这提示现有网络更稳定地利用伙伴环境修正侧链状态画像，而不是主链 RMSF 画像。

## 5. 得到或验证的生物学洞见

第一，结合界面的刚化是本数据集中最强、最可重复的效应。按起始结合态最小跨链重原子距离分层，0–5 Å 界面残基的平均 ΔRMSF 为 −0.280 Å，[−0.371, −0.216]，5–12 Å 近界面为 −0.094 Å，[−0.126, −0.060]，>12 Å 远端均值为 −0.001 Å，[−0.077, 0.093]；64.1% 的界面残基下降超过 0.1 Å，而上升超过 0.1 Å 的只有 2.8%。χ 熵呈同样梯度：界面平均 −0.240 nats，[−0.284, −0.198]，近界面 −0.021，远端约为零。这验证了结合限制局部主链与侧链构象空间的经典观点，也解释了为什么简单距离先验很强。

第二，远端均值接近零不代表没有变构响应，而是正负效应互相抵消。>12 Å 区域中，每个复合物等权后仍有约 11.7% 的残基 ΔRMSF 增加超过 0.1 Å，20.8% 下降超过 0.1 Å；ΔRMSF 与 Δχ 熵的残基内 Spearman 只有 0.266，[0.189, 0.347]，说明主链移动和侧链状态改变部分耦合但不是同一个量。原论文也报告结合可在远离界面的位点增加柔性，并将其解释为变构通信；本次重新计算支持“远端异质性”而不是“远端整体变柔”这一更精确表述。

第三，3SGB 给出了模型当前会漏掉的典型远端增加案例。抑制剂 B 链 N 端 2、3 位距伙伴约 19 Å，结合后 ΔRMSF 分别为 +2.575 和 +2.632 Å，相对于 apo/holo 窗口方差合成标准差的比值为 1.49 和 1.41；FlexPPI 却分别预测 −0.188 和 −0.128 Å。1BRS 的 A2 也在 16.9 Å 处增加 +0.708 Å，效应为 2.28 个窗口标准差，模型预测 −0.299 Å；A43 则在 14.6 Å 处真实下降 −0.354 Å、效应为 −2.23 个窗口标准差，模型正确预测下降。这些残基是后续做轨迹可视化、动态互相关网络或突变验证的优先候选，但因为每个状态只有一条长轨迹、窗口不是独立随机种子，仍不能把它们直接称为已证实的变构位点。

第四，模型性能似乎与界面动态复杂度相关。原论文报告 2OOB 和 1GCQ 的初始界面接触最稳定且没有分出多个界面簇，而另外 6 个体系出现多个界面亚态；1AK4 与 3SGB 还发生明显伙伴重定向。本次 FlexPPI ΔRMSF Spearman 在 1GCQ、2OOB 分别为 0.598、0.612，二者均值 0.605；在没有大重定向但有多亚态的 4 个体系均值 0.281，在 1AK4、3SGB 只有 0.160、0.186，均值 0.173。这是一个样本极小、来自事后分层的探索性洞见，但它提出了明确可证伪假设：当前静态起始结构模型更擅长稳定界面，难以预测轨迹中通过伙伴重定向或亚态切换产生的响应；未来应加入界面亚态/ensemble 表示，而不是只继续扩大单结构编码器。

## 6. 不能越过的解释边界

这 8 个体系是 rigid-body docking benchmark 的精选案例，功能类别也偏向酶/抑制剂，不能代表所有 PPI。每个状态只有一条生产轨迹，100 ns 窗口提供的是时间分块稳定性而不是独立重复；apo 与 holo 是独立模拟，holo−apo 是状态差而不是同一分子沿真实结合路径的因果时间序列。界面距离取自生产轨迹起点的静态结构，无法表达原论文观察到的动态接触重排；本次没有重新计算水桥、ASA、Gap index 或界面聚类，不能用当前结果验证水介导相互作用和结合自由能。5 个 PDB 有状态预训练身份重叠，3 个未见 PDB 又不足以形成高功效外部基准。最重要的是，任何远端候选都需要多随机种子重复、动态网络分析和最好有突变/实验数据支持后，才能升级为机制性生物学结论。

## 7. 复现命令与结果文件

```powershell
cd E:\Project\FlexPPI

# 1. 生成标准 PDB、状态清单和严格配对清单
.\.venv-md\Scripts\python -m flexppi.cli index-zenodo-ppi `
  --root "E:\data\protein_md_datasets\zenodo_6638504_ppi_md" `
  --output "artifacts\zenodo_6638504\data"

# 2. 100 ns = 2500 帧，计算 RMSF/χ 熵标签
.\.venv-md\Scripts\python -m flexppi.cli compute-md-labels `
  --manifest "artifacts\zenodo_6638504\data\zenodo_6638504_states.jsonl" `
  --output-dir "artifacts\zenodo_6638504\data\labels" `
  --window-frames 2500

# 3. 把标签路径写回 24 个状态与 16 个配对
.\.venv-md\Scripts\python -m flexppi.cli index-zenodo-ppi `
  --root "E:\data\protein_md_datasets\zenodo_6638504_ppi_md" `
  --output "artifacts\zenodo_6638504\data" `
  --label-dir "artifacts\zenodo_6638504\data\labels"

# 4. 仓库原生响应评估
.\.venv-md\Scripts\python -m flexppi.cli evaluate `
  --config "configs\evaluate_zenodo_6638504.yaml" `
  --checkpoint "artifacts\finetune\best.pt" `
  --split all --bootstrap 1000 `
  --delta-threshold rmsf=0.1 `
  --delta-threshold chi_entropy=0.1 `
  --output "artifacts\zenodo_6638504\evaluation\native_evaluation.json"

# 5. 基线、伙伴删除、空间生物学和 LOCO 校准
.\.venv-md\Scripts\python scripts\evaluate_zenodo_6638504.py `
  --pair-manifest "artifacts\zenodo_6638504\data\zenodo_6638504_pairs.jsonl" `
  --checkpoint "artifacts\finetune\best.pt" `
  --native-evaluation "artifacts\zenodo_6638504\evaluation\native_evaluation.json" `
  --output-dir "artifacts\zenodo_6638504\evaluation" `
  --bootstrap 2000
```

主要结果位于 `artifacts/zenodo_6638504/evaluation/`：`report.json` 是机器可读总报告，`summary_metrics.csv` 是按全部/重叠/未见 PDB 分层的模型指标，`residue_predictions.csv` 含 1,759 个残基的真值、预测、距离和不确定性，`biological_summary.csv` 汇总界面/近界面/远端效应，`top_remote_residues.csv` 给出远端候选及窗口标准差，`loco_calibration_*.csv` 是严格留一复合物校准结果，`source_paper_strata.csv` 是按原论文界面动态类别的探索性分析，`zenodo_6638504_audit.png/.pdf` 是四面板审计图。

## 8. 下一轮最值得做的实验

第一优先级是对 1BRS、2OOB、3SGB 做多随机种子重复或把现有长轨迹按更保守的块长重新估计不确定性，确认 A2/A43 与 3SGB B2/B3 的远端效应；第二优先级是从结合轨迹每 1–5 ns 抽取界面亚态代表结构，对同一体系输入多结构 FlexPPI，检验 ensemble 平均能否特别改善 1AK4/3SGB；第三优先级是加入动态接触持续率、界面亚态数或基于轨迹的 residue interaction network，再做严格 leave-one-PDB-out，检验模型能否超过距离先验；第四优先级才是全参数小样本微调，因为当前最明显的瓶颈是目标幅度校准和动态亚态缺失，而不是绝对状态编码完全失效。若将该数据集用于论文主结论，还需要先完成序列/结构同源簇泄漏审计，并把 8 个体系之外的 DynaRepo test53 严格配对作为更大外部验证集。
