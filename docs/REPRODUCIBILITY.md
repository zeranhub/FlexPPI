# 当前论文的复现入口

论文v136使用冻结的3,063条ATLAS、DynaRepo和ITrace记录。实际配置位于 `experiments/core_corrected_pipeline/training_payload/configs/`，五折清单位于对应的 `splits/`。早期v1、performance-first v2和difficulty-balanced v3准备方案不作为本稿入口。

先安装 `python -m pip install -e ".[dev]"`，运行 `python scripts/maintenance/verify_repository.py` 和 `python -m pytest -p no:cacheprovider`。科学分析另需pandas、scipy、scikit-learn；图形生产脚本还使用matplotlib、Pillow、PyMuPDF、shapely、pypdf、reportlab或pdfplumber。轨迹处理与结构会话需要原分析环境中的MDAnalysis、PyMOL等；外部方法需要相应官方软件及权重。

`python experiments/run.py --help`列出当前流程。`training`对应冻结训练及种子核查，`family`对应家族隔离和外部评价，`occupancy`对应83条记录的角度占据验证，`mel5`对应MEL5结构及三重复MD方法。`--list`不执行计算，`--check`只检查登记路径；执行需明确指定`--step`，原参数放在`--`之后。

完整重算前需提供冻结清单指向的结构、标签、模型检查点及原始轨迹。固定第40轮、原体系五折、全链家族划分、外部2,430条训练记录及83条共同窗口验证不能互换。请按原程序的路径参数或`path_mappings`配置数据，不修改原清单和哈希。部分原冻结计划绑定当时的代码哈希；迁移后若触发守卫，应审核并建立新的冻结记录及独立输出，不能替换旧哈希冒充同一实验。

当前九幅图和必要源表见 [figures](../paper/current/figures/README.md)，[manifest](../paper/current/manifest.json)校验所收录文件。图源复制只适配必要路径，科学数值和绘图几何不变。完整标签、全队列逐样本预测、权重及轨迹未全部内嵌，科学重算仍需补齐原输入。本次清理验证发布内容、入口、语法和基础回归，没有重新训练或运行新MD。
