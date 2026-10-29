# 备份与恢复

## 2026-10-29 私有GitHub同步

同步目标为私有仓库 `zeranhub/FlexPPI`。本轮按原相对路径保存项目代码、配置、测试、研究说明、数据划分与关键分析表，以及当前中文稿 `paper/revision_v14_20261010/paper_v129_20261029_中文修订_抗原识别与抗体动态.docx`。v129已经加入6A77抗体案例；图4e的v3仍为独立预览，尚未插入Word。BCR目录中更早的“DOCX未修改”说明是历史阶段记录，应以此快照说明和v129审计为准。

当前图4保留原TCR图a–d及其v124结构输入、渲染源码和PyMOL会话，同时保存6A77的完整位点预测、MD分析表、必要NPZ证据、结构输入、图e v3源码与PNG/PDF/SVG，以及完整图4预览。6A77绘图的规范输入为 `experiments/core_corrected_pipeline/training_payload/dataset/structures/A00LZ.pdb`，图a风格参数及原会话也在本轮清单中。

逐文件范围、字节数和SHA256见 `docs/backups/2026-10-29/selected_files.csv`，选择规则和排除目录见同目录的 `selection_summary.json`，生成程序为 `scripts/maintenance/prepare_private_github_sync.py`。该备份保留源码和可核对证据，不包含完整的MD轨迹、训练数据集和模型权重，也不包含本地虚拟环境、第三方下载依赖、大批逐样本缓存或历史逐页渲染。重建分析时仍需补齐相应数据和权重，并按原配置调整机器路径。

重新生成清单：`python scripts/maintenance/prepare_private_github_sync.py --hash`。生成清单不会上传、删除或修改研究文件；提交与推送需单独执行，并在推送前确认仓库的 `isPrivate` 为 `true`。

Git 保存源代码、配置、测试、数据划分清单、文档，以及已纳入版本管理的论文与图稿。Git 不保存被忽略的 `artifacts/`、`output/`、本地虚拟环境和临时目录；尚未提交的新图稿也不受已有提交保护。因此只克隆仓库不能恢复完整研究环境，模型权重和标签必须另行备份。

## 本次清理归档

2026-09-05 的清理归档位于 `E:/Project/FlexPPI-archive/cleanup-20260905-180602-740`。其中 `receipt.json` 逐项记录原路径、归档路径、文件数、大小、操作状态和清理前 Git 提交。归档保留原有相对目录；例如 `allinone.md` 在归档根目录，旧试跑输出在归档内对应的 `artifacts/`、`output/` 或 `tmp/` 路径下。清理前的 README 与本备份说明也保存了一份。归档位于同一磁盘，是可恢复搬移，不是异地备份。

需要恢复时，先查 `receipt.json` 中的 `action`：`archive` 表示已搬移，`snapshot` 表示仅复制旧文档，`delete` 表示已删除可再生成文件。恢复归档项前确认原位置没有新文件，再复制回项目内对应相对路径；不要整批覆盖当前项目。MMseqs2 工作目录包含链接，作为中间缓存保留归档即可；若需重跑划分，优先重新生成，不依赖迁移后的旧缓存。

已删除的 Python/测试缓存和论文渲染、编译中间产物没有逐份备份，需要重新生成。其中原来被 Git 跟踪的 `paper/paper.xdv` 也可从清理前提交 `974bd51` 提取，但通常直接重新编译即可。`paper/paper.pdf`、参考文献源文件和当前图稿均保留。历史代码标签 `benchmark-first-v1` 只保存当时已提交内容，不包含本地权重及新图稿。

## 换机恢复

克隆代码后安装 `python -m pip install -e ".[dev]"`，再运行 `python -m pytest`。另行复制实际需要的 `artifacts/data/` 标签、模型权重、`output/` 分析结果及未提交的论文素材；不要直接复制旧 Python 虚拟环境作为可移植环境。

外部 ATLAS、DynaRepo、ITrace、ThermoPCD 及历史对比数据并不随仓库迁移。请检查当前 `configs/performance_first_v2_fold*.yaml` 中的 `path_mappings`，将旧机器路径映射为新机器的实际目录，再执行 [README 中的数据审计](README.md)。论文可从 `paper/paper.tex`、`paper/references.bib` 和现存图稿重新编译，所需 TeX 环境需另行安装。
