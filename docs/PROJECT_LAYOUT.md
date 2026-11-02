# 当前项目目录

论文v136及独立补充材料位于 [paper/current](../paper/current/README.md)，图1–4和补充图1–5位于 [figures](../paper/current/figures/README.md)。这里保存当前导出、可编辑源及必要源表。

`src/flexppi/`是共用模型和数据实现，`tests/`覆盖保留的基础行为。`experiments/run.py`和`steps.json`登记当前论文流程，`core_corrected_pipeline/training_payload/`保存实际冻结训练代码、配置及清单，`reproducibility/`保存原标签构建源码。`experiments/_sources/`保留当前结果的科学依赖；原日期或版本名用于追溯来源，不代表额外的论文案例。

`scripts/maintenance/verify_repository.py`检查发布范围、图像一致性、登记路径、入口链接和Python语法。环境要求见 [复现说明](REPRODUCIBILITY.md)，清理标准见 [仓库范围](REPOSITORY_SCOPE.md)。本机的`artifacts/`、`output/`、旧稿、完整MD和外部软件不随这份GitHub内容完整上传。
