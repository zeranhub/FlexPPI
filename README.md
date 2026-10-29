# FlexPPI

## 最终论文与实验入口（2026-10-23整理）

当前论文以 [v58中文最终稿](paper/revision_v14_20261010/paper_v58_20261022_中文修订_图2至5配色统一_三线表.docx) 为唯一版本依据。文档保持原始字节与原路径。六张实际使用的图及可编辑源已整理到 [paper/figures_v58](paper/figures_v58/README.md)，实验按照论文第3章和补充材料小标题建立 [experiments入口](experiments/README.md)。目录说明见 [docs/PROJECT_LAYOUT.md](docs/PROJECT_LAYOUT.md)。

原 `experiment_completion` 中的科学实现和结果已迁至 `experiments/_sources/`；其他工作树中的最终训练、折外预测和检查点快照收拢到 `experiments/core_corrected_pipeline/`。`src/`、根 `scripts/`、`configs/`、`data/` 保留原项目功能。下面的历史技术说明包含旧论文和训练流程，不能据其判断最终论文版本。

旧稿及缓存按用户要求集中移入根目录 [待删除](待删除/说明.md)，保留原层级和移动清单。原整理清单见 [cleanup_plan_20261023.json](docs/cleanup_plan_20261023.json)。清单不是永久删除指令，实际移动状态以 [cleanup_status_20261023.json](docs/cleanup_status_20261023.json) 为准。

FlexPPI 当前主任务是从蛋白质静态结构预测每个残基的绝对 RMSF 与侧链构象熵 Hχ（代码中称为 `chi_entropy`）。RMSF 描述残基围绕平均位置的运动幅度，Hχ 描述侧链在不同朝向之间切换的多样程度。本轮采用 structure-only（只使用结构及其序列信息）设定，温度、轨迹时长和数据来源仅用于记录标签来源，不作为模型输入；SASA、apo/holo 配对、结合响应和变化量预测不属于本轮训练目标。

## 当前工作入口

训练定义以 [performance-first v2 协议](docs/performance_first_v2_training_protocol.md) 为准，数据情况见 [v2 数据划分说明](data/splits/state_rmsf_hchi_v2_performance_first/README.md)。五份当前配置是 `configs/performance_first_v2_fold0.yaml` 至 `fold4.yaml`；每次用四折训练、另一折测试，外层测试折不用于挑选训练轮次或调整参数。严格的“未见过任何同源家族”评价是独立要求，不能把现有五折结果直接当作这一结论。

目前准备流程使用 ATLAS、DynaRepo 和 ITrace 已生成的标签；ThermoPCD 仍需统一重算标签，不能把候选清单当作已完成训练数据。目录整理不代表正式训练或训练前所有验收已经完成，具体以数据审计和实际运行结果为准。

当前数据准备脚本集中在 `scripts/`：`integrate_itrace_hchi_labels.py` 整合 ITrace 标签，`prepare_thermopcd_candidates.py` 筛选 ThermoPCD 候选，`prepare_performance_first_v2.py` 构建 v2 清单与划分，`audit_pdbbind_ppi.py` 核查数据类别，`audit_performance_first_v2.py` 检查标签、划分及配置。它们不是每次训练前都要重新执行的流水线；已有冻结划分应先审计，修改划分需要单独记录版本。

## 目录与历史边界

`src/flexppi/` 是模型、数据加载和命令入口，`configs/` 存放当前五折配置与历史实验配置，`data/splits/` 存放可版本管理的数据清单，`tests/` 是回归测试，`baselines/` 是对比方法。更具体的保留范围和旧文件去向见 [项目目录与清理记录](docs/PROJECT_LAYOUT.md)。

`artifacts/` 保存本机标签、权重及审计结果，`output/` 保存基准评价和分析产物。这两处大量文件不在 Git 中，不能因为名称像“生成结果”就删除。`paper/` 保存论文、图稿及绘图素材；`docs/` 同时保留当前协议和历史研究记录。`tmp/` 还包含本机绘图辅助工具，并非可以整体删除的垃圾目录。

旧版基线的代码快照标记为 Git 标签 `benchmark-first-v1`。旧的 `pretrain.yaml`、`finetune.yaml` 和 `evaluate_*.yaml` 用于历史实验，不是当前重训的默认配置。部分旧能力仍被分析脚本和测试调用，因此保留兼容实现，但不属于当前训练目标。旧 README 与 `allinone.md` 已移入项目外清理归档，入口见 [备份说明](BACKUP_MANIFEST.md)。

## 安装与检查

在项目根目录安装开发依赖并运行测试：

```powershell
python -m pip install -e ".[dev]"
python -m pytest -p no:cacheprovider
```

不改写已冻结审计文件的复核命令：

```powershell
python scripts/audit_performance_first_v2.py --output artifacts/maintenance/v2_readiness_after_cleanup.json
```

该审计检查数据和配置，不代替训练前的前向计算、反向计算和显存验证。轨迹处理需要可选的 `md` 依赖，ESM 编码需要可选的 `esm` 依赖；现有 `.venv-md/` 保留，不在本次整理中重建。迁移电脑时需重新配置外部数据路径，详见 [备份与恢复说明](BACKUP_MANIFEST.md)。外部实现与许可证来源见 [NOTICE](NOTICE)。
