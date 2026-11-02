# 备份与恢复

当前GitHub内容以论文v136为依据，稿件、补充材料和当前九幅图位于 [paper/current](paper/current/README.md)。[manifest.json](paper/current/manifest.json)记录所收录文件的字节数和SHA256，可运行 `python scripts/maintenance/verify_repository.py` 核对。文档、图像及保留的科学数值按原字节收录，冻结训练与分析实现继续使用原字节。

2026-11-02的清理移出过时稿件、放弃案例、候选筛选、试绘、运行状态、预览和渲染检查输出，并同步更新README和实验登记。原研究资料留在本机，逐文件清单位于本机 `tmp/github_cleanup_20261102/cleanup_plan.json`；Git历史可追溯先前已提交内容。

只克隆仓库不能恢复完整研究环境。完整MD轨迹、标签数组、模型检查点、外部方法权重、环境和部分逐样本预测另行保存，冻结清单中的来源及哈希用于识别这些输入。迁移时按 [复现说明](docs/REPRODUCIBILITY.md) 配置环境和路径，不改写冻结清单，不直接复制旧虚拟环境。恢复本机旧研究资料前应逐项核对目标，避免覆盖当前稿件或数据。
