# 图3非结构子图正式文件

本文件夹包含已确认方案的10个原生Origin子图项目：b的全部25位点MD/Full/Intra三条曲线及同25位点MD/Full细节；c的66有效泛素位置一致性散点；d的全部六个Gln成对点与SASA；e的χ1–χ4四个关联图；f的两个83记录配对相关差图。绘图从v98论文问题和真实数据出发，未打开原子图作为造型参考。原始DOCX与结构图未修改。

正式分图在origin_panels_final中，每个子目录都有result.opju、result.pdf、result.png、result.tif和origin_verify_report.json。OPJU保留完整源表、绘图坐标和可编辑曲线、点、坐标轴及文字。所有项目经过保存后独立重开检查，数值完全一致，重开的PNG与首次PNG逐像素一致；原生文字对象读回Arial标题10 pt、面板字母11 pt、轴标与刻度9 pt，轴线0.6 pt。Full为#2166AC蓝色实线，Intra为#9A633D棕色虚线；散点中的Full使用同一蓝色。PDF保留矢量图形和嵌入字体，不包含整图位图；Origin打印驱动使PDF实际文字尺寸有约0.12 pt以内的量化，原生文字对象按指定字号保持可编辑，详情见final_figure_QA.json。

figure3_nonstructural_origin_collection.pdf为四页矢量合集：第1页b及b细节，第2页c/d，第3页e的四幅图，第4页f的两幅图。合集按各分图导出的物理尺寸拼接，不缩放字号；它是分图查阅集，不是已嵌入结构图的一页论文整图。figure3_origin_panels_overview.png是长图查阅索引。对单个子图继续调整时应打开相应OPJU；插入论文时按其原物理尺寸使用，避免任意缩放改变字号。

final_caption_zh.md给出对应正式图形的新版图注与正文配对差定义。f图采用横向相关差，不能沿用旧图注的绝对相关分布表述。83条记录先逐记录计算Full Hχ相关减Full RMSF相关，再记录等权汇总，区间来自来源内配对重抽记录10,000次。细箱熵差为0.189324 [0.171200,0.206497]，低F*角区总宽度差为0.161017 [0.143002,0.178873]。

source_data 中的10份CSV逐文件保持 v98 配套源数据包原值。data_audit.md/json 记录数据来源、全部数据核验、有效掩码、不同案例的留出划分与样本单位。f_paired_differences.csv 保存两个指标各83个记录内配对差，f_paired_summary.csv 保存来源分层10000次记录重抽的均值区间。以上记录、残基、角度或轨迹帧不能混作独立生物重复。

render_origin_panels.py为本次任务的Origin原生对象绘图代码，复用了EditaPlot的独立实例、串行运行、导出与读回工具；公共CLI的精确XPS样式参数不冒称支持普通散点与折线。scientific_design_proposal.json保存当时审阅的不可变提案；human_scientific_confirmation.json记录用户“确认，按此方案制作”的答复，origin_execution_approval.json绑定科学提案与实际渲染计划。环境smoke已通过。

figure3_editaplot_正式分图与源数据.zip仅收录正式分图、矢量合集、总览、源数据、代码、图注和审计。文件夹内带draft_preview的文件及origin_panels、origin_panels_v2是早期设计或排版检查记录，不属于正式交付包。重画需使用已准备的EditaPlot项目Python环境C:\Users\DELL\.agent-reach\editaplot\runtime\.editaplot-venv\Scripts\python.exe。默认调用脚本只核对数据并生成计划；正式执行要求现有确认记录和成功smoke，避免静默改变源数据。
