# 备份与恢复

当前GitHub内容以论文v136为依据，稿件、补充材料和当前九幅图位于 [paper/current](paper/current/README.md)。[manifest.json](paper/current/manifest.json)记录所收录文件的字节数和SHA256，可运行 `python scripts/maintenance/verify_repository.py` 核对。文档、图像及保留的科学数值按原字节收录，冻结训练与分析实现继续使用原字节。

2026-11-02先清理GitHub，随后按相同标准实际删减本机文件：移除过时稿件、放弃案例、候选筛选、试绘、旧训练版本、运行状态、预览、渲染检查输出和废弃环境。旧研究文件已从本机删除；先前已提交的内容可从Git历史追溯，未提交的文件不能仅靠Git恢复。本机逐文件删除记录及核查结果位于 `docs/local_cleanup_20261102/`；先前GitHub清理记录位于 `tmp/github_cleanup_20261102/`。

只克隆仓库不能恢复完整研究环境。本机继续保留当前论文所需标签数组、逐样本预测、16组第40轮神经模型检查点、外部方法权重、有效运行环境和必要输入；原始MD来源及哈希按冻结清单识别。PLANET-MD部署压缩包含唯一模型权重，因此保留。图3、图4实际使用的PyMOL可编辑会话与图3输入结构已归并到 `paper/current/figures/`。迁移时按 [复现说明](docs/REPRODUCIBILITY.md) 配置环境和路径，不改写冻结清单，不直接复制旧虚拟环境。
