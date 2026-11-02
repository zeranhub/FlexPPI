# 当前论文实验入口

本目录对应论文v136：正文3.1是性能与方法比较，3.2是完整结构环境的贡献，3.3是CaM–CaMKI、泛素–UBA及角度占据验证，3.4是MEL5亲和成熟。实现按工作流程登记在 [steps.json](steps.json)，不把放弃案例或中间试图列为额外结果。

```powershell
python experiments/run.py --help
python experiments/run.py performance --list
python experiments/run.py training --check
python experiments/run.py family --list
python experiments/run.py occupancy --list
python experiments/run.py mel5 --list
```

默认仅列出步骤，`--check`只核对登记路径。运行必须显式指定`--step`，原参数放在`--`之后。部分原程序没有参数解析，选择对应步骤就会开始计算，应先查看列出的用途和输入条件。

`core_corrected_pipeline/`保存冻结训练、标签、配置和评价实现；`_sources/`保留当前结果仍依赖的科学源文件。共用程序不复制为独立实验。当前图源见 [paper/current](../paper/current/figures/README.md)。完整重算所需的数据、软件及原冻结守卫要求见 [复现说明](../docs/REPRODUCIBILITY.md)。
