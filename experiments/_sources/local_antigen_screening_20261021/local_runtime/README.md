# 本地 CPU 推理入口

本目录只负责读取既有模型与本地结构，生成可追溯的预测。脚本不训练、不联网、不启动远程作业，也不修改输入结构或原始 NPZ。当前已验证的解释器为 `D:/Anaconda/envs/pvb_cpu/python.exe`，Torch 2.4.1，仅 CPU。代码明确从 `E:/Project/FlexPPI/src` 导入；当前 worktree 中的旧版接口缺少本次模型所需的 `information_ablation`，因此不能使用默认工作目录替代指定源码。

原计划 family 2430 权重在本地不可用且无法取得，主任务在任何本轮候选预测之前明确改用本地已有的历史 strict_full 权重；原因与旧预期SHA保存在 `checkpoint_selection_change.json`，未依据候选结果选择模型。这两套权重不是同一个模型。

当前唯一允许的权重是 历史 strict_full 2431 条训练记录、无 ITrace、epoch 40，SHA256 为 `b3eed060919df15f075329bc892f959d9355232be5feb561417edd49ea82b39f`。默认文件为上一级 `model/epoch_040.pt`。脚本同时检查该权重的 learned 序列编码、关闭条件输入、关闭配对变化预测头、4 个 IPA 模块，以及 RMSF/Hχ 两项输出。任何不同权重都会被拒绝，旧的 artifacts 权重不能补位。

## 原结构有限 smoke

```powershell
& 'D:/Anaconda/envs/pvb_cpu/python.exe' 'E:/Project/FlexPPI/paper/revision_v14_20261010/experiment_completion/local_antigen_screening_20261021/local_runtime/local_cpu_inference.py' --manifest 'E:/Project/FlexPPI/paper/revision_v14_20261010/experiment_completion/local_antigen_screening_20261021/local_runtime/smoke_manifest.json' --output-dir 'E:/Project/FlexPPI/paper/revision_v14_20261010/experiment_completion/local_antigen_screening_20261021/local_runtime/smoke'
```

该清单只含 3HG1 的原输入，用于检查全部 824 个残基、两项输出、残基身份、有限数值及有效掩码，保存实际 CPU 耗时。这里没有采用另一权重的 family NPZ 作为重现参照，也不声称复现 ff297440 模型。`smoke/run_audit.json` 保存真实耗时、各源码 SHA、模型 SHA 与输入 SHA。只有当参照确实出自当前相同权重时才可填写 `reference_prediction`，脚本才执行小于 `2e-5` 最大绝对误差的逐值重现检查。

## 候选批处理

输入 JSON 的 `candidates` 是候选列表。每个候选有唯一的英文、数字或下划线名称 `candidate_id`，以及本地绝对路径 `input_pdb`。对真实独立结构，可用 `input_mode: observed_structure`。对同长度肽序列的固定骨架扰动，提供 `chain_sequences`，例如 `{ "C": "EVAGIGILTV" }`。必须根据结构实际链名填写，脚本不硬编码某一条链为抗原。突变序列必须与该结构实际可见的该链残基数相同；不支持插入、删除或凭空补齐缺失残基。

```json
{
  "candidates": [
    {
      "candidate_id": "example_fixed_backbone",
      "input_pdb": "E:/absolute/path/template.pdb",
      "input_mode": "fixed-backbone",
      "chain_sequences": { "C": "EVAGIGILTV" }
    }
  ]
}
```

运行时只需用正式清单替换 `--manifest`，并更换 `--output-dir`。模型只加载一次，候选按顺序推理，默认每个进程使用 2 个 CPU 线程，避免与其他本地任务争抢全部线程。`--limit 1` 可先试一例；`--prepare-only` 仅检查并生成一致的模型输入，不读取权重、不进行推理。无论指定哪些目标链，模型仍接收整个复合物，输出保留全部残基并提供目标链标记。

每例输出包含 `prediction.npz`、`prediction.csv`、输入准备审计与逐例审计。预测 CSV 保留全部残基、几何有效标记及该氨基酸是否有 χ 角的化学定义。Ala/Gly 的 Hχ 数字不得用于有意义侧链熵统计；任何进一步比较都应结合 `geometry_valid`、`target` 和 `hchi_chemically_defined`，并保留完整预设区域。方差是模型输出参数，本脚本没有将其解释为已校准的置信区间。

## 固定骨架扰动的明确含义

脚本同步改写突变残基名称，并由同一 PDB 重新读取氨基酸编码；突变位点原来的侧链原子会被删除，避免让旧侧链坐标冒充突变后的真实侧链。随后逐值检查全部 N/CA/C/O 坐标、掩码、残基编号、链顺序与模板完全一致。该模型采用 learned 氨基酸嵌入，不输入 ESM 特征，因此不存在应重算却沿用 WT 的 ESM 缓存。这一层实验仅检验模型在同一骨架下对序列身份变化的响应，不能证明突变后的结构重排、稳定性、结合自由能、抗原呈递或 TCR 激活。生成的 `model_input_fixed_backbone.pdb` 也不能直接当作完整突变结构用于 MD。

## 本地结构建模能力

当前 pvb_cpu 有 BioPython 与 PyMOL，可解析结构和展示结构；没有 OpenMM、PDBFixer、PyRosetta 或 Modeller。另三个已检查的 conda 环境及项目 `.venv-md` 同样未发现这些软件包，FoldX、Rosetta、Modeller 不在 PATH。可用 ESM2-t33 缓存存在于其他项目，但目标模型不使用它。是否安装独立的任务环境以及开展真正侧链重建、受限局部松弛由主任务统一安排；本目录没有擅自安装或把序列替换称为结构重建。
