# 运行与文件说明

先阅读“结果说明与案例筛选结论.md”。本目录保存的是实际试跑及其限制，并非已通过验证的新机制案例。screening_key_results.csv是完整模型的便于阅读摘要；全部模型、全部区域和全部逐残基输出仍在其他CSV和predictions目录，未删除不利结果。

本地运行分为三个环境，因为已有环境的科学依赖版本不同。模型前向计算使用D:\Anaconda\envs\pvb_cpu\python.exe（Python 3.10.20、PyTorch 2.4.1、NumPy 2.2.6）；编号使用mel5_absolute_affinity_20261028\affinity\venv中的Python，并加载mel5_tcr_binding_20261028\numbering\deps里的ANARCII 2.0.8；统计及最终核查使用用户Python 3.12.10，NumPy 2.4.6、pandas 3.0.5、SciPy 1.17.1。每个脚本内的R是本次结果目录，SRC是冻结模型源代码目录；它们是本机路径。搬到别处需调整这些路径、states.json的输入路径及检查点路径，不应直接解释成通用一键环境。

实际执行顺序为scripts/download_structures.py、prepare_states.py、number_antibody_regions.py、run_screen.py、analyze_screen.py、verify_and_report.py；后两个脚本分别从逐残基CSV重算全部比较、另行核对CSV/NPZ和原输入的一致性。已有MD标签的提取及重算脚本在md_reference/export_and_verify_labels.py；家族核查脚本是runtime/family_input_audit.py。完整运行记录、检查点来源及模型版本见inference_audit.json和runtime内的审计文件；原始参数与冻结代码位置在这些文件中记录。

主检查点为runtime/family_folds/runs/full/fold_0/epoch_040.pt，对照为intrachain同位置的epoch_040.pt；它们仍保存在本地完整目录，未放入结果压缩包。原MD输入PDB已按字节复制至reference/md_input_structures作为核对副本，SHA与原文件相同；原XTC轨迹和缓存标签仍在审计记录的D盘/E盘路径，压缩包未收录。重算的窗口标签、逐位MD标签和差异记录已收录，因此可以复核本报告的数值，但不能在缺少原始轨迹时重新计算全部动态真值。

压缩包内的package_manifest_sha256.json记录每个收录文件的SHA256；其自身和ZIP不参与自身校验。压缩包排除了模型权重、原始轨迹、Python缓存及环境，不改变本地完整工作目录。模型源代码在项目现有冻结目录，本包未复制一份新的源码树。
