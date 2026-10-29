# 项目目录与当前论文

本项目的最终论文依据是用户指定的 v58 Word，不是历史 `paper/paper.tex` 或 `paper/paper.pdf`。Word保留原路径，整理不改正文、表格、公式或嵌入图。

| 路径 | 用途 |
| --- | --- |
| `paper/revision_v14_20261010/paper_v58_20261022_中文修订_图2至5配色统一_三线表.docx` | 唯一最终论文入口 |
| `paper/figures_v58/` | 图1–5及S1，按图号命名，含可编辑drawio/SVG、数据、代码与分子结构素材 |
| `experiments/3_*.py` | 对应论文3.1–3.5原小标题的实验入口 |
| `experiments/[A-G]_*.py` | 对应补充材料A–G的实验入口 |
| `experiments/_sources/` | 保留原内部布局的科学实现及实验结果，替代原experiment_completion路径 |
| `experiments/core_corrected_pipeline/` | 从4666工作树收拢的最终训练、折外预测、静态基线、权重及类型内分析源 |
| `experiments/_case_validation_sources/` | KIX与泛素历史验证所需的结构、真实MD帧、数值和PyMOL源 |
| `experiments/_legacy_paper_data/` | 历史论文数值源表，保留以供追溯 |
| `src/`、`scripts/`、`configs/`、`data/`、`tests/` | 项目模型、通用工具、配置、数据划分和测试 |
| `artifacts/`、`output/` | 科学数据、模型及真实结果，不按缓存整体删除 |
| `docs/` | 方法协议、调研与本次整理审计 |
| `待删除/` | 用户要求集中保存的旧版本和中间产物，按原路径分层，移动清单保存在此 |

实验入口先用 `python experiments/<对应标题>.py --list` 或 `--check` 查看步骤与依赖；仅显式指定 `--step` 才执行。最终图检查使用 `python paper/figures_v58/reproduce.py --check`。重绘写入单独输出目录，不覆盖论文当前图。

原始MD轨迹、训练环境和远程训练不因整理而自动重跑。`tmp/third_party` 含真实基线权重和结构样本，予以保留；两个PPT构建目录的node_modules是共享运行时链接，不计作可释放空间。

用户已将永久删除改为集中移动到根目录 `待删除/`。以其中 `移动清单.json` 为准，不执行永久删除。原清理计划及实际状态分别见 `cleanup_plan_20261023.json` 与 `cleanup_status_20261023.json`；路径迁移记录见 `path_migration_20261023.json`。未执行的删除不会记录为已完成。
