# FlexPPI

FlexPPI 从蛋白质或复合物的静态主链结构预测逐残基主链波动 RMSF 和侧链联合状态熵 Hχ。当前论文使用氨基酸类型、N–Cα–C 局部框架及主链几何；MD 用于建立监督标签和验证读数的结构含义。Hχ 是固定 χ 角分箱下的联合状态占据熵，不是总构象熵，也不直接换算为结合自由能。

当前稿件为 **v136，2026-11-01**，仓库于2026-11-02按这一版本清理。论文及附加材料集中在 [paper/current](paper/current/README.md)：[论文Word](paper/current/manuscript.docx)、[补充Word](paper/current/supplementary_material.docx)、[补充PDF](paper/current/supplementary_material.pdf)。四幅主图、五幅补充图、可编辑图源及必要源表按当前编号整理在 [figures](paper/current/figures/README.md)。

科学实现见 [experiments](experiments/README.md)，数据与环境要求见 [复现说明](docs/REPRODUCIBILITY.md)。实际论文训练使用 `experiments/core_corrected_pipeline/training_payload/` 中冻结的3,063条记录、五折清单及配置。全组成链家族隔离和ThermoPCD外部评价使用各自的独立划分。冻结清单保留原数据路径与哈希，在新机器运行前需配齐数据并配置路径映射。

```powershell
python -m pip install -e ".[dev]"
python scripts/maintenance/verify_repository.py
python -m pytest -p no:cacheprovider
python experiments/run.py training --list
```

`experiments/run.py` 默认仅列出步骤，只有显式指定 `--step` 才执行原程序。训练、MD标签重建和外部方法还需相应的轨迹、标签、检查点及软件；路径存在检查不代表完整重跑已通过。GitHub保存当前实现、冻结清单和选定源表，完整MD与模型权重另行保存。

保留标准见 [仓库范围](docs/REPOSITORY_SCOPE.md)，目录见 [PROJECT_LAYOUT](docs/PROJECT_LAYOUT.md)，恢复边界见 [备份说明](BACKUP_MANIFEST.md)。已放弃案例、旧稿、候选筛选和试绘不再列为当前结果；共享训练数据不因对应案例退出正文而删记录。许可说明见 [NOTICE](NOTICE)。
