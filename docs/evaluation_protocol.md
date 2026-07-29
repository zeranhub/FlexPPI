# FlexPPI-Bench：静态结构柔性与结合动态响应评测方案

## 1. 结论先行

当前 checkpoint 已经证明它能从静态结构中恢复相当一部分“绝对柔性排序”，尤其是侧链构象熵；但它尚未证明自己可靠掌握了“结合导致的带符号柔性变化”。最明显的证据是：微调前的绝对状态模型在 276 个内部验证系统上取得 RMSF 的蛋白等权 Spearman 0.811、侧链构象熵 0.928、SASA 0.871；结合响应微调后，这三个数分别降为 0.750、0.919 和 -0.096，其中 SASA 出现灾难性遗忘。对于 7 个内部配对验证系统，响应头把 ΔRMSF 的 MAE 从零变化基线的 0.221 Å 降到 0.209 Å，改善 5.2%，但活跃变化残基的方向准确率只有 34.3%；Δχ 熵的 MAE 0.139 反而差于零变化基线 0.133，方向准确率只有 2.7%。因此，当前最严谨的表述应是“模型较好地学习了静态结构中的绝对柔性线索，并对 ΔRMSF 提供了弱但可测的增益；尚不能声称普遍、可靠地预测结合前后所有柔性变化”。

## 2. 五支柱完整性审计

| 支柱 | 状态 | 当前内容 | 提交前必须补强的部分 |
|---|---|---|---|
| Research Gap | 部分完成 | 现有单一相关系数无法区分绝对柔性、带符号结合响应、界面效应、强变化检出与不确定性校准；“静态状态预测好”等价于“结合响应预测好”的隐含假设被本项目结果直接否定 | 需要系统检索并制作 5–8 个已有柔性/动力学 benchmark 的比较表，不能只以本项目内部评价为证据 |
| Construction Pipeline | 基本完成 | 已建立 ATLAS 三重复 RMSF、DynaRepo 窗口化 RMSF/χ 熵/SASA、test53 apo/holo 映射、TCR3D B 因子与 CDR 代理的自动流水线 | 需要冻结公开版本、输出同源簇、许可清单、失败样本和数据卡 |
| Evaluation Framework | 已实现 | 同时报残基池化与蛋白等权 MAE/RMSE/Pearson/Spearman/R²/CCC，高柔性检出，响应方向三分类，零变化与状态相减对照，界面分层和不确定性覆盖率 | 需要补充独立模型 baselines、多个随机种子与显著性检验；界面定义可从 Cα 8 Å 代理升级为全重原子 5 Å |
| Empirical Findings | 已形成 | 已得到静态能力、微调遗忘、ΔRMSF 弱增益、Δχ 失败、TCR3D 跨域排序和域编码敏感性等可行动发现 | 当前配对结论只基于 7 个模型选择验证系统，必须用未参与训练和选轮次的独立集合复核 |
| Companion Method | 已有 | FlexPPI 本身就是与 benchmark 配套的绝对状态头和响应头 | 需要通过修复 SASA 遗忘、增加方向/排序目标和独立测试，证明 benchmark 确实能驱动方法改进 |

## 3. 研究缺口、任务边界与研究问题

本 benchmark 对应“维度盲区 + 粒度失配”：只评估静态结构上的平均 RMSF 相关性，会漏掉三个实际重要的问题。第一，一个模型可能在每个状态上都给出相似的柔性轮廓，却在两个大数相减时给错结合响应；第二，残基总体 MAE 可能被大量变化很小的非界面残基主导，无法说明它是否找到了真正发生刚化或增柔的界面和远端变构残基；第三，模型给出的方差可能与实际误差不匹配，无法用于后续实验优先级排序。一个贯穿全文的 Figure 1 示例应展示同一受体的 apo/holo 两个静态结构：旧式评价只看到两个状态各自相关性较高，而 FlexPPI-Bench 同时检查 ΔRMSF 的正负号、界面与远端、响应头相对简单相减的增益以及置信区间是否覆盖真实变化。

RQ1 是“仅给一个静态结构时，模型能否在未见系统上恢复绝对 RMSF、侧链 rotamer entropy（旋转异构体熵）和 SASA，并识别最柔性的 20% 残基？”；RQ2 是“给出配对 apo/holo 静态结构时，专用响应头能否稳定超过零变化与两个状态直接相减，尤其是在界面、强变化残基和变化方向上？”；RQ3 是“这种能力能否跨到 TCR3D 的实验结构上，并且预测是否对人为的数据域编码稳定？”这里的任务边界必须写清：ATLAS/DynaRepo 的 MD 标签是主要真值；TCR3D 的晶体 B 因子和 CDR 区域只是独立代理证据，不能叫作 MD-RMSF 真值；当前也不评价生成完整构象集合、动力学时间尺度、自由能面或结合亲和力。

## 4. 数据审计与真值适用范围

| 数据源 | 本地规模与当前可用标签 | 在评测中的角色 | 关键限制 |
|---|---|---|---|
| ATLAS | 1,938 个系统，444,240 个 RMSF 残基标签；三条 100 ns 重复轨迹 | 绝对 RMSF、重复间标签方差、高柔性区域检出 | 已参与预训练；只能重建内部验证，不能当未见测试 |
| DynaRepo | 922 个系统，366,777 个 RMSF、294,358 个 χ 熵、51,459 个 SASA 标签 | 绝对状态跨数据源检验与标签不确定性 | 当前 SASA 只覆盖部分系统；数据域与 ATLAS 差异明显 |
| DynaRepo test53 | 原始 37 对，受 `max_residues=1024` 后 36 对；验证集 7 对，共 2,124 个 ΔRMSF 和 1,691 个 Δχ 有效残基 | 结合响应、响应头/状态相减/零变化对照、界面分层 | 全部 37 对已经进入微调清单；验证集还用于选择 best.pt，不是独立测试；ΔSASA 标签数为 0 |
| TCR3D | α/β 单链候选 1,206 条，排除精确训练 PDB 重叠 12 条后候选 1,194 条，成功评测 1,189 条、129,268 个残基 | 外部晶体 B 因子排序、CDR1/2/3 柔性富集、域编码敏感性 | 无 MD 轨迹；B 因子受分辨率、晶格、精修和 TLS 影响；尚未进行序列同源去污染 |
| 其他 `protein_md_datasets` | BioEmu 约 37.26 GiB 压缩训练集，mdCATH 仅一个本地样本，MISATO 只有 split 文本，GPCRmd 只有索引 | 后续独立 MD 测试候选 | 目前缺少可直接批量消费的“静态结构 + 独立轨迹真值”；不能在未解压/未下载时虚报为已评测 |

数据构建应采用“真实轨迹标签抽取 + 受控分层”的组合范式。流水线为：原始 PDB/XTC 与元数据进入许可和完整性检查；每条链单独刚体对齐后计算窗口 RMSF，χ 二面角离散成 rotamer 后计算 Shannon entropy，SASA 用 Shrake–Rupley 算法计算；apo/holo 通过序列全局比对建立残基对应；随后按同源簇而不是单个 PDB 随机拆分，冻结 train/dev/test；最后生成界面、强变化、标签噪声、长度和数据源标签。每一步必须保留输入、输出、拒绝原因和样本数量。当前代码已覆盖标签抽取、残基映射和绝大部分自动检查，但“全库同源簇拆分”仍未完成，是最终无泄漏结论的首要阻塞项。

## 5. 指标体系

绝对状态的主指标应是“每个蛋白先算 Spearman，再对蛋白等权平均”，因为它回答模型能否在一个具体蛋白内部把高柔性与低柔性残基排对次序，而且不会让长蛋白因残基多而支配结果。MAE 和 RMSE 保留绝对量纲；Pearson 衡量近似线性一致性；Spearman 衡量排序；R² 与 concordance correlation coefficient（CCC，一致性相关系数）揭示尺度和偏移；高柔性检出把每个蛋白真实值最高的 20% 作为正类，报告 top-20% recall、AUROC 和 average precision。模型方差与标签重复方差相加得到总标准差，再检查 50%/90%/95% 区间覆盖率以及“预测标准差与绝对误差”的 Spearman。

结合响应必须并列报告三种方法：专用 response head、`holo_state - apo_state` 状态相减、恒为 0 的 no-change 基线。对每种方法同时报告回归指标、绝对变化最大的 20% 残基检出，以及带 dead zone（RMSF/χ 默认 0.1，SASA 默认 2.0 Å²）的 decrease/neutral/increase 三分类；只对真实变化超过阈值的残基计算 active sign accuracy，防止接近 0 的轨迹噪声让符号指标失真。结果还要按 Cα 跨链距离不超过 8 Å 的界面代理与非界面分层。正式论文建议把界面升级为任意重原子距离不超过 5 Å，并增加 5–12 Å 邻近壳层和远端变构层。

## 6. 当前实测结果

### 6.1 绝对状态：微调前后对照

| 标签 | checkpoint | 有效系统/残基 | MAE | 残基池化 Pearson | 蛋白等权 Spearman | 高柔性 20% AP |
|---|---:|---:|---:|---:|---:|---:|
| RMSF | pretrain | 276 / 66,464 | 0.535 Å | 0.701 | 0.811 | 0.786 |
| RMSF | finetune | 276 / 66,464 | 0.592 Å | 0.679 | 0.750 | 0.743 |
| χ entropy | pretrain | 85 / 21,870 | 0.201 | 0.958 | 0.928 | 0.939 |
| χ entropy | finetune | 85 / 21,870 | 0.212 | 0.951 | 0.919 | 0.933 |
| SASA | pretrain | 11 / 4,028 | 16.80 Å² | 0.871 | 0.871 | 0.801 |
| SASA | finetune | 11 / 4,028 | 76.88 Å² | -0.168 | -0.096 | 0.234 |

**Finding 1.** 绝对状态预训练确实学到了强的局部柔性排序，侧链构象熵表现最好；RMSF 在蛋白内部的等权排序相关性明显高于跨蛋白池化相关性，说明模型更擅长判断“同一蛋白哪里更柔”，而不是把不同蛋白的绝对振幅放到统一尺度上。**Finding 2.** 结合响应微调造成多任务遗忘，SASA 几乎完全崩溃；这不是轻微数值漂移，因为高 SASA 检出已接近随机、95% 区间覆盖率从 0.927 降到 0.233。后续训练需要冻结或小学习率保护绝对状态头、加入旧任务 replay，并以每个任务独立验证指标做多目标 checkpoint 选择，而不能只用总 NLL。

### 6.2 结合前后变化

| 标签/方法 | MAE | Pearson | Spearman | 相对零变化 MAE 改善 | active sign accuracy |
|---|---:|---:|---:|---:|---:|
| ΔRMSF response head | 0.209 Å | 0.236 | 0.197 | +5.2% | 34.3% |
| ΔRMSF state subtraction | 0.275 Å | 0.293 | 0.171 | -24.7% | 37.5% |
| ΔRMSF zero change | 0.221 Å | — | — | 0 | 0%（按定义不能判断方向） |
| Δχ response head | 0.139 | 0.176 | 0.046 | -3.9% | 2.7% |
| Δχ state subtraction | 0.145 | 0.234 | 0.082 | -8.8% | 15.2% |
| Δχ zero change | 0.133 | — | — | 0 | 0% |

**Finding 3.** 响应头对 ΔRMSF 有真实但很弱的 MAE 增益，并把状态相减的 MAE 降低 24.0%；然而它对“增柔还是刚化”的判定仍不可靠，特别是 increase 类 F1 只有 0.063。模型很可能主要学会把预测向 0 收缩，从而降低平均误差，却没有学会正确的响应方向和强变化尾部。**Finding 4.** Δχ 熵尚未学成：响应头甚至不如零变化 MAE，而且非界面 Spearman 为 -0.009；界面子集 Spearman 0.407 暗示局部界面线索存在，但不足以支持全局侧链响应结论。当前 test53 没有成对 SASA 标签，因此任何“模型能预测结合导致的 SASA 变化”的说法都没有实证支撑。

### 6.3 TCR3D 外部代理

在 1,189 条成功评测的 TCR α/β 链上，预测 RMSF 与晶体 B 因子换算代理的每链等权 Spearman 为 0.389，95% bootstrap 区间为 [0.378, 0.397]。预测 RMSF 检出 CDR3 相对框架区的 pooled AUROC 为 0.612、AP 为 0.332，而随机 AP 等于正例比例 0.144；所有 CDR 合并后 AUROC 只有 0.553、AP 0.377（随机 AP 0.303），CDR1 甚至低于随机排序。χ 熵对 CDR3 近似随机，微调后的 SASA 在 CDR3 上呈反向排序。将输入域编码从 ATLAS 值 0 改成 DynaRepo 值 1 时，平均绝对预测漂移达到 RMSF 0.278 Å、χ 熵 0.287、SASA 10.73 Å²，说明新数据域的输出对一个人为条件变量较敏感。

**Finding 5.** FlexPPI 的 RMSF 头在独立实验结构上保留了中等强度的局部排序信号，并能优先标出一部分 CDR3，但这种跨域证据远弱于内部 MD 验证；域编码敏感性又表明结果会随无法客观确定的输入条件变化。正式部署应移除数据集身份这一非物理条件，或改成可测的模拟协议、温度、离子强度和时间窗，并在未知域上做显式校准。

## 7. 可复现实验命令

本地 Windows 已提供路径映射配置，下面三条命令分别重建内部绝对状态验证、内部配对验证和 TCR3D 外部代理。`--split val` 恢复训练器使用的确定性验证划分，但该验证集参与了 best checkpoint 选择，输出报告会自动写入警告；不要把它重命名为 test。

```powershell
$env:PYTHONPATH="src"
python -m flexppi evaluate --config configs/evaluate_state_local.yaml --checkpoint artifacts/pretrain/best.pt --split val --bootstrap 1000 --output artifacts/evaluation/state_validation_pretrain.json
python -m flexppi evaluate --config configs/evaluate_state_local.yaml --checkpoint artifacts/finetune/best.pt --split val --bootstrap 1000 --output artifacts/evaluation/state_validation_finetune.json
python -m flexppi evaluate --config configs/evaluate_pair_local.yaml --checkpoint artifacts/finetune/best.pt --split val --bootstrap 1000 --output artifacts/evaluation/pair_validation.json
python -m flexppi evaluate-tcr3d --checkpoint artifacts/finetune/best.pt --root "E:\Project\Data\TCR3d_database" --exclude-manifest artifacts/data/atlas_states.jsonl --exclude-manifest artifacts/data/dynarepo_states.jsonl --bootstrap 1000 --output-dir artifacts/evaluation/tcr3d_full
```

远程服务器上的清单本来就是 `/ranze/...` 路径，因此可直接用 `configs/pretrain.yaml` 作为绝对状态评测配置、`configs/finetune.yaml` 作为配对配置。服务器完成后应保存 JSON、逐残基 CSV、checkpoint 的 SHA256、代码版本和 CUDA/PyTorch 环境。若要补 ΔSASA，必须先对 test53 的 apo 与 holo 两侧都计算 SASA 标签，再重新索引配对清单；只给一侧算 SASA 仍然无法得到变化真值。

## 8. 论文逻辑骨架（Section 2–7）

Introduction 的六段逻辑应依次是：蛋白结合动态的重要性和一个 apo/holo running example；现有评价把静态状态准确性等同于动态响应的三个盲区；RQ1–RQ3；覆盖、细粒度诊断、可扩展与高质量四个设计目标；FlexPPI-Bench 的数据构建、三层证据与上述五个 findings；最后列出 benchmark、可复现流水线、系统评测和能力边界四项贡献。Section 2 定义绝对状态、结合响应和外部代理三类任务以及 G1–G4；Section 3 用 Figure 2 展示 ATLAS/DynaRepo/TCR3D 从原始文件到标签、同源拆分、分层和质量门的流水线；Section 4 介绍 companion model FlexPPI，但必须把 pretrain 与 response fine-tune 的目标冲突写清；Section 5 按 RQ1–RQ3 组织实验，并在每个子节末给出 Finding 1–5；Section 6 讨论灾难性遗忘、方向学习、域条件、独立配对数据和人类/实验验证机会；Section 7 用至少 5–8 个现有 benchmark 的比较表完成 related work，最后以结论重申“能测出能力边界”而不是“模型全面成功”。

## 9. 提交前检查与判定

Pre-submission Verdict：**NOT READY**。Critical 问题有四项：没有从训练和 checkpoint 选择中完全隔离的 MD 测试集；37 个 test53 配对已全部进入微调流程；没有 ΔSASA 成对标签；尚无已有 benchmark 的系统比较表和完整 Figure 1/Figure 2。Major 问题包括：配对验证只有 7 个系统且只跑一个训练随机种子；尚未进行 MMseqs2 序列同源和结构簇去污染；TCR3D 只提供间接代理；缺少外部模型 baselines 与人类/实验专家基线；界面当前为 Cα 8 Å 代理；SASA 任务在微调后灾难性遗忘；数据和代码许可、Croissant/数据卡、维护计划尚未整理。只有完成固定同源簇 test、重新训练且不以 test 选轮次、补齐独立 apo/holo MD 配对和 ΔSASA 后，才能将当前“内部诊断报告”升级为论文级 benchmark 结论。

