# v63 图2入稿与文字同步

正式Word为 `E:/Project/FlexPPI/paper/revision_v14_20261010/paper_v63_20261023_中文修订_图2紧凑构图与文字同步.docx`。以v61为基稿另存，v61保留；v62不再使用。

正文第11页采用[已确认的a–c紧凑图2](../figure2_compact_layouts_20261023/fig2_compact_merged.pdf)。a为双坐标预测对应图，b为逐记录AP分布，c为三个任务的PR组合；其可编辑源、冻结数据和绘图核查均保留在[图2源目录](../figure2_compact_layouts_20261023/README.md)。原来源分层AP配套图移入补充材料C作为图S2，原图文件见[来源分层配套图](../figure2_C_selected_20261023/fig2_source_AP_companion.pdf)。

正文3.1、3.2的图号与性能描述已同步，增加图2b的真实记录级AP中位数，并说明分布重叠。图注区分了a的两类Full输出、b的四分位范围和c的PR任务；补充材料C说明三个独立连续分段线性纵轴及原始AP计分。没有新增训练、推理、显著性检验或更改原性能值。

`build_manuscript.py`通过局部OOXML编辑保留36个原生数学对象、原有表格和其他图片；`export_manuscript.ps1`用独立不可见Word实例只读导出；`render_manuscript.py`调用documents技能的渲染器与bundled Poppler生成全部28页检查图。Windows运行时为Codex bundled Python，不使用用户安装LibreOffice。核查记录保存在`document_QA/`。
