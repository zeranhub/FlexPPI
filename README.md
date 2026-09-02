# FlexPPI

FlexPPI 是一个面向蛋白–蛋白结合动态响应的研究代码库。模型把每个输入状态表示成“目标受体 + 可选伙伴环境 + 一份静态三维结构”，由同一个状态编码器预测逐残基绝对 RMSF、侧链构象熵和 SASA，再在同一受体的 apo/holo 对之间预测带正负号的变化量。这里的 RMSF（root mean square fluctuation，均方根涨落）表示一个残基在轨迹中围绕平均位置晃动的幅度；侧链构象熵表示侧链在多少种不同朝向之间切换；SASA（solvent-accessible surface area，溶剂可接触表面积）表示残基暴露给溶剂的面积。

项目参考了 [AlphaFold](https://github.com/google-deepmind/alphafold) 的 Invariant Point Attention（IPA，不变点注意力）和 [Ibex](https://github.com/prescient-design/ibex) 的 PyTorch/OpenFold 结构模块、分链 ESM 编码与批次补齐方式，但没有复制完整结构生成模型。FlexPPI 使用输入 PDB 中的 N–CA–C 原子建立固定的逐残基局部坐标架，IPA 只更新残基表示，不更新坐标；这样既保留整体旋转和平移不变性，又避免为动态标签任务引入不必要的坐标生成与扭转角重建。外部设计来源与许可证说明见 [NOTICE](NOTICE)。

## 当前实现

数据层能够直接扫描本机的 `E:\data\protein_md_datasets\ATLAS_extracted` 和 `D:\data\dynarepo_non_nucleic_full`。ATLAS 预处理器读取每个条目的三列 RMSF，保存三次重复的均值和样本方差；DynaRepo 索引器会区分 PDBbind 复合物、DynaBench、抗原、test53 复合物和 test53 单链，并先按 PDB ID 找候选，再通过 Needleman–Wunsch 全局序列比对把单链残基映射到复合物中的正确受体链。不能仅凭链名配对，因为 partner-stripped 单链和复合物中的链标识并不一定相同。

DynaRepo 轨迹标签处理是可选模块，需要 MDAnalysis。RMSF 会对复合物中的每条链单独进行刚体拟合，再计算该链内部的涨落，从而避免“两个伙伴相对转动”被误当成“残基本身更柔性”。如果本地只有一条 500 ns 轨迹，处理器会把它切成多个固定窗口，用窗口间方差近似标签不确定性；如果以后补齐多条重复轨迹，同一接口也能接受多条 `trajectory_paths`。侧链熵把每个可用的 χ 二面角分到三个旋转异构区间，然后统计联合旋转状态的 Shannon entropy；Gly 和 Ala 等没有 χ 角的残基自动保留为空标签。SASA 使用项目内的 Shrake–Rupley 实现：对每个重原子的范德华球加上 1.4 Å 水探针，再统计没有被邻近原子遮挡的球面采样点，因此不依赖 Windows 下需要额外 C++ 编译器的 MDTraj。

网络包含可切换的序列编码器、几何成对特征、固定坐标架 IPA、绝对状态头和成对响应头。默认 `sequence_backend: learned` 使用小型可学习氨基酸嵌入，方便在当前 CPU 环境运行测试；改成 `sequence_backend: huggingface_esm` 后，会按链分别调用配置中的 ESM2，再把各链表示放回复合物残基位置。状态头同时预测标签均值和模型方差，损失把模型方差与重复轨迹产生的标签方差相加，噪声很大的残基不会像高置信度残基一样主导训练。响应头读取 apo 表示、holo 表示、两者之差、逐元素乘积和绝对状态直接差值，并通过一致性项约束

\[
\widehat{\Delta y}\approx\hat y_{\mathrm{holo}}-\hat y_{\mathrm{apo}}.
\]

推理结果会同时保存响应头的 `delta_*` 和简单状态相减的 `state_subtraction_*`，因此可以直接检验响应头是否比“分别预测再相减”更好，而不是默认复杂模块一定有效。

## 安装

当前机器是 Python 3.9，并已有 CPU 版 PyTorch。核心功能可直接以可编辑方式安装：

```powershell
cd E:\Project\FlexPPI
python -m pip install -e .
```

处理 DynaRepo XTC 轨迹时安装 MD 依赖；使用 ESM2 时安装 ESM 依赖。MDAnalysis 2.7 的 Windows 二进制包需要 NumPy 1.x，而当前全局环境使用 NumPy 2.x，所以建议让轨迹处理使用隔离环境，不要直接降级全局 NumPy：

```powershell
python -m venv --system-site-packages .venv-md
.\.venv-md\Scripts\python -m pip install -e ".[md]"
python -m pip install -e ".[esm]"
```

ESM 权重会由 Hugging Face 首次下载。正式训练建议使用带 CUDA 的 PyTorch；当前环境的 `torch 2.8.0+cpu` 只适合数据检查、单元测试和小规模冒烟训练。

## 数据准备

先把 ATLAS 的现成 RMSF 转成紧凑标签和状态清单：

```powershell
flexppi index-atlas `
  --root "E:\data\protein_md_datasets\ATLAS_extracted" `
  --output "artifacts\data"
```

索引 DynaRepo 并建立 test53 配对清单：

```powershell
flexppi index-dynarepo `
  --root "D:\data\dynarepo_non_nucleic_full" `
  --output "artifacts\data"
```

第一次只建议处理一个体系和一个 100 ns 窗口，确认内存、残基数和轨迹拓扑一致。`--groups` 后面的两个英语名称分别表示 test53 单链和 test53 复合物；`--sasa` 会额外执行较耗时的 Shrake–Rupley SASA 计算，`--sasa-points` 控制每个原子球面上的采样点数，默认 96，数值越大越精确也越慢：

```powershell
flexppi compute-dynarepo-labels `
  --manifest "artifacts\data\dynarepo_states.jsonl" `
  --output-dir "artifacts\data\dynarepo_labels" `
  --groups test53_single test53_complex `
  --limit 1 `
  --window-frames 1000 `
  --max-windows 1
```

标签生成后重新执行 `index-dynarepo`，并通过 `--label-dir artifacts\data\dynarepo_labels` 把标签路径写入状态和配对清单。可以用下面的命令检查清单中有多少记录真正具备标签：

```powershell
flexppi inspect-manifest --manifest "artifacts\data\dynarepo_pairs.jsonl"
```

## 训练、评估与推理

绝对状态预训练和配对微调分别使用两个配置文件。预训练配置中的 `state_manifest` 是一个清单列表，会把 ATLAS 与已经生成标签的 DynaRepo 状态合并；尚无标签的 DynaRepo 条目会自动过滤，而不是被当成零值标签：

```powershell
flexppi train --config configs\pretrain.yaml
flexppi train --config configs\finetune.yaml
```

`configs/finetune.yaml` 默认从 `artifacts/pretrain/best.pt` 初始化。训练器以完整体系作为样本，并按 `cluster_id`、PDB ID 或体系 ID 分组划分训练/验证集，不会把同一个复合物的不同目标链随机拆到两侧。不过 ATLAS 当前本地文件没有现成同源簇；正式论文实验仍应先用 MMseqs2 或结构拓扑生成 `cluster_id`，再进行全库同源隔离，不能把当前按体系 ID 的默认划分当成最终无泄漏划分。

评估会分别报告每个标签的有效残基数、MAE、Pearson 相关和符号准确率：

```powershell
flexppi evaluate `
  --config configs\finetune.yaml `
  --checkpoint artifacts\finetune\best.pt
```

现在的评测器还会报告蛋白等权宏平均、Spearman、RMSE、R²、CCC、高柔性残基检出、不确定性覆盖率，并在结合响应任务中并列比较响应头、两个状态直接相减和零变化基线。训练器当前只有 train/validation，没有未接触 test；因此应显式使用 `--split val`，并保留报告中的模型选择偏乐观警告。本地清单来自远程服务器时可使用带 `/ranze/...` 到 Windows 路径映射的 `configs/evaluate_state_local.yaml` 与 `configs/evaluate_pair_local.yaml`：

```powershell
flexppi evaluate --config configs\evaluate_state_local.yaml --checkpoint artifacts\pretrain\best.pt --split val --output artifacts\evaluation\state_validation_pretrain.json
flexppi evaluate --config configs\evaluate_state_local.yaml --checkpoint artifacts\finetune\best.pt --split val --output artifacts\evaluation\state_validation_finetune.json
flexppi evaluate --config configs\evaluate_pair_local.yaml --checkpoint artifacts\finetune\best.pt --split val --output artifacts\evaluation\pair_validation.json
```

TCR3D 可作为独立的晶体 B 因子排序和 CDR 柔性富集代理测试；它没有 MD 轨迹，因此命令和输出都会明确标记“不是 MD-RMSF/结合响应真值”。完整评测设计、实测结果、数据泄漏边界和论文式 benchmark 骨架见 [docs/evaluation_protocol.md](docs/evaluation_protocol.md)。

只输入一个 PDB 时输出绝对状态；同时输入 apo 受体和 holo 复合物时，程序会自动做序列映射并输出成对变化：

```powershell
flexppi predict `
  --checkpoint artifacts\finetune\best.pt `
  --apo-pdb path\to\apo.pdb `
  --holo-pdb path\to\complex.pdb `
  --target-chains A `
  --output artifacts\prediction.json
```

## Zenodo 6638504 严格 apo/holo PPI 实验

本地 `E:\data\protein_md_datasets\zenodo_6638504_ppi_md` 的 8 个经典 PPI 复合物现已支持 GRO+XTC 原子顺序安全转换、100 ns 窗口标签、16 个链侧严格配对、冻结 checkpoint 评估、伙伴删除消融、GNM/接触数/距离基线以及 leave-one-complex-out 校准。完整设计、实测结果、泄漏边界、生物学解释和复现命令见 [docs/zenodo_6638504_experiments.md](docs/zenodo_6638504_experiments.md)。

```powershell
flexppi index-zenodo-ppi --output artifacts\zenodo_6638504\data
flexppi compute-md-labels `
  --manifest artifacts\zenodo_6638504\data\zenodo_6638504_states.jsonl `
  --output-dir artifacts\zenodo_6638504\data\labels `
  --window-frames 2500
```

## 代码布局与验证

`src/flexppi/data` 负责 PDB、清单、配对和 MD 标签；`src/flexppi/model` 负责序列编码、几何特征、IPA 和预测头；`losses.py` 实现噪声感知状态损失、响应损失与一致性损失；`trainer.py`、`evaluation.py` 和 `inference.py` 分别处理训练、指标和推理。运行全部测试：

```powershell
$env:PYTHONPATH="src"
pytest
```

测试覆盖 PDB/链解析、带缺口的残基对应、ATLAS 三重复均值与方差、整个状态网络在全局旋转和平移下的不变性，以及 apo/holo 前向传播、损失和反向传播。ATLAS 本地 README 标明数据采用 CC BY-NC 4.0，因此当前数据不能在未获得额外许可的情况下直接用于商业训练或商业模型发布。

本机验收结果为：ATLAS 的 1,938 个目录全部成功建立状态清单；DynaRepo 的 922 个目录全部成功索引，其中包括 406 个 PDBbind 复合物、210 个 DynaBench 复合物、215 个抗原单体、52 个 test53 单链和 39 个 test53 复合物。test53 最终得到 37 条通过“序列一致度至少 0.90、覆盖度至少 0.80”的 apo/holo 配对，另有 2 条因覆盖度只有 0.768 和 0.789 被明确拒绝。真实轨迹 A00CW/A00BR（PDB 1F6M）的闭环测试成功读取 310 残基 apo 与 834 残基 holo，映射 310 个受体残基，完成绝对状态和响应前向计算、联合损失与反向传播；A00CV 的 SASA 冒烟测试也为 295 个残基全部生成了有限标签。
