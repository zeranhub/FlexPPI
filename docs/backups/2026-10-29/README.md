# 2026-10-29 私有研究快照

本轮同步到私有仓库 `zeranhub/FlexPPI`，范围是本地源代码、配置、测试、关键说明和分析证据。逐文件本地源路径、大小和SHA256见 `selected_files.csv`；目录选择与排除规则见 `selection_summary.json`。SHA256记录的是同步前本地文件字节；Git按现有配置处理文本换行，二进制稿件、结构会话和数值归档保持原字节。

当前论文为 `paper/revision_v14_20261010/paper_v129_20261029_中文修订_抗原识别与抗体动态.docx`，已包含TCR后的6A77抗体案例。图e v3是沿用图a风格的独立预览，尚未嵌入Word；完整图4、单独e、a/e对照和源代码均在快照内。图e可复现目录为 `experiments/_sources/bcr_case_choice_20261028/figure4_extension_20261029/revision_e_v3/`，规范输入结构为 `experiments/core_corrected_pipeline/training_payload/dataset/structures/A00LZ.pdb`。图a的v124原渲染脚本、PDB和PyMOL会话也保留。

本轮没有上传完整的训练数据、模型权重与MD轨迹。已保留预测值、关键MD派生数值和审计，但在另一台机器重新执行模型推断或轨迹分析时，仍需取得相应权重与轨迹，再调整脚本中的本地路径。虚拟环境、第三方下载依赖、逐样本缓存、大量历史图像与逐页渲染不属于此快照。

`prepare_private_github_sync.py --hash` 只生成候选清单，不会自行提交、推送或删除文件。同步前后均检查GitHub仓库可见性；本次操作必须为PRIVATE。未创建计划任务，因此这是一次手动同步。
