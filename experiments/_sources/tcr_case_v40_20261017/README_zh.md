v40 TCR识别环境变化实例复现材料

最终稿位于本目录上两级，文件名paper_v40_20261017_中文修订_TCR识别环境变化实例_三线表.docx。正文新增3.5/图6；原泛素与自由能分析为3.6/图7。原稿v39的图1–5与S1及其源数据包仍适用。

figure6/：最终英文图PNG/PDF/SVG、所有13个肽位置的源表、中文图注、PyMOL会话及数据/视觉核验。plot_case.py为可直接运行的绘图脚本，使用既有plot_runtime；render_structures.py使用现有headless PyMOL。
sb27_sb47/：六条原始轨迹的逐位标签与Q7角度/占据/接触摘要、原始文件SHA256、模型输入和完整/移除TCR推理。analyze_sb27_sb47.py从E:/data/ITrace原始文件重新计算；原始XTC依ITrace许可获取。
context_inference.py：hit-gpu上的同一联合模型输入对照；运行目录/root/flexppi_completion_20261010/tcr_case_v40_20261017，环境/root/flexppi_completion_20261010/env/bin/python，原源代码/root/flexppi_revision_20261002/corrected_training_v1/src。模型检查点和原始预测已由既有完成实验备份及本包SHA记录追溯。
family_predictions/：244条ITrace当前家族留出Full/Intra预测及清单。正文实例仅用Full，未使用历史3输出模型的预测。
local_pair_screen/、change_screen/、change_screen_focus/：完整序列匹配、配对筛查及文献候选对比。biology_sources/保留所用原始文献/数据库信息；unrelated_PMC3668746_excluded.md是错误地址返回的无关论文，不用于分析，也不纳入交付压缩包。
manuscript_QA/：原稿/新稿SHA、逐数值/逐结构核验、30页页面检查记录与全文提取。可用build_v40.py重建文稿，export_v40_word.ps1与render_v40.py以隐藏Word导出后交由官方render_docx.py渲染，audit_v40.py复核数值与文档结构。
