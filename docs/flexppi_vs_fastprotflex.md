# FlexPPI 与 FastProtFlex 的论文评测范围

## 1. 当前论文只回答两个定量问题

论文不再追求覆盖所有已有实验，而是集中回答两个能够由现有 MD-RMSF 标签直接支持的问题。RQ1 是：只给一个单链静态结构时，FlexPPI-static 与 FastProtFlex 谁能更准确地恢复同一蛋白内部从刚到柔的残基排序？RQ2 是：面对多链蛋白复合体时，把 FastProtFlex 分别用于每条链，与让 FlexPPI 读取完整复合体相比，哪一种方式更准确，优势主要发生在界面还是非界面区域？

本文使用 `artifacts/pretrain/best.pt` 作为 FlexPPI-static 的论文评测检查点。

## 2. 共同任务与指标

FastProtFlex 输出的是蛋白内部标准化柔性分数，FlexPPI 输出带 Å 单位的 RMSF，因此二者共同比较的是逐目标 Spearman、Pearson、蛋白内部标准化后的 z-MAE/z-RMSE，以及真实 top-20% 高 RMSF 残基的 AUROC、Average Precision 和 Precision@K。每个蛋白或复合体先独立计算指标，再做等权宏平均，避免长蛋白因为残基数多而获得更大的统计权重。误差棒使用蛋白级或复合体级 bootstrap 95% 置信区间。

## 3. 单链静态柔性预测

| 数据 | 方法 | Spearman ↑ | Pearson ↑ | z-MAE ↓ | Top-20 AUROC ↑ | Top-20 AP ↑ |
|---|---|---:|---:|---:|---:|---:|
| ATLAS，191 个蛋白 | FastProtFlex | 0.7984 | 0.7717 | 0.4870 | 0.9125 | 0.7821 |
| ATLAS，191 个蛋白 | FlexPPI-static | **0.8205** | **0.8509** | **0.3393** | **0.9230** | **0.8024** |
| DynaRepo 单链，31 个蛋白 | FastProtFlex | 0.8298 | 0.7476 | 0.5032 | 0.9134 | 0.7678 |
| DynaRepo 单链，31 个蛋白 | FlexPPI-static | **0.8517** | **0.7854** | **0.4088** | **0.9212** | **0.7835** |

ATLAS 上 FlexPPI-static 的 Spearman 平均高 0.0221，DynaRepo 单链子集上高 0.0218。两套单链数据给出一致方向，因此可以写成“FlexPPI-static 在单链静态 RMSF 排序上具有小幅、稳定的优势”。ATLAS 不是完全未见的外部测试；DynaRepo 也存在训练域不对称，所以不能把差值完全归因于网络架构。

## 4. 复合体整体预测

DynaRepo 复合体子集包含 54 个复合体、119 条蛋白链和 19,153 个共同有效残基。FastProtFlex 的正确部署方式是每条链独立预测；FlexPPI-static 则读取完整多链结构。

| 方法与输入方式 | Spearman ↑ | Pearson ↑ | z-MAE ↓ | Top-20 AUROC ↑ | Top-20 AP ↑ |
|---|---:|---:|---:|---:|---:|
| FastProtFlex，每条链独立 | 0.7599 | 0.6930 | 0.5549 | 0.8926 | 0.7374 |
| FlexPPI-static，完整复合体 | **0.7902** | **0.7506** | **0.4490** | **0.9042** | **0.7620** |

FlexPPI-static 的复合体宏平均 Spearman 高 0.0303，复合体 bootstrap 95% CI 为 [0.0066, 0.0562]。但界面分层显示，界面 Spearman 是 FlexPPI 0.7022、FastProtFlex 0.7155；非界面 Spearman 则是 FlexPPI 0.8121、FastProtFlex 0.7678。因此论文应写成“完整复合体建模改善了总体链内柔性轮廓，优势主要来自非界面区域”，而不能写成“FlexPPI 已经解决了结合界面柔性预测”。

## 5. 8PJG 定性案例

8PJG 是已经结合的 TCR–pMHC 复合体。FlexPPI 输入完整 A–E 五链结构，FastProtFlex 则将每条链单独预测。两种模型的结果被转换为各链内部的柔性百分位并写入 PDB B-factor 列，以相同的蓝–白–红色标显示。

FlexPPI 的 CDR3 偏白、FastProtFlex 的 CDR3 偏红，不应被定义为 FlexPPI 失败。CDR3 在游离状态可以具有较高内在柔性，但结合抗原肽和 MHC 后会受到界面接触限制。只有 FlexPPI 看到了这些结合伙伴，因此较低的结合态 CDR3 相对柔性与接触诱导刚化机制一致。白色只是链内中等柔性排序，不代表完全刚性；又因为缺少同一 TCR 的未结合结构或 MD 真值，这个案例只能证明两个模型对结合上下文给出了不同响应，不能证明其中一个模型准确。

建议论文表述：

> Complex-aware FlexPPI assigns lower relative flexibility to the bound CDR3 loops than chain-wise FastProtFlex, consistent with contact-induced restriction in the TCR–pMHC interface. This qualitative contrast is not treated as accuracy validation because matched unbound or MD-derived flexibility is unavailable.

## 6. 复现与论文文件

```powershell
$env:PYTHONPATH='src'
python scripts\compare_fastprotflex.py --dataset atlas --bootstrap-repetitions 2000 --output-dir artifacts\evaluation\fastprotflex_comparison
python scripts\compare_chain_context.py --bootstrap-repetitions 2000
python scripts\generate_paper_benchmark_figures.py
```

论文图表位于 `artifacts/paper/benchmark_figures`；8PJG PyMOL 文件位于 `artifacts/evaluation/tcr_pmhc_pymol/8PJG`。
