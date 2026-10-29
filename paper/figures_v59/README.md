# v59论文图件

本包在v58最终Word实际使用的图件上重绘图2—图5。蓝色、暖赭色和青绿色是三种彩色主调，灰色用于MD参考和背景；曲线、结构、位点标注与图例一起调整。图1和S1沿用原文件。整套预览为`figures_2_to_5_overview.png`，Word文件名为`paper_v59_20261023_中文修订_图2至5整体重绘.docx`，主项目中保存在`paper/revision_v14_20261010/`，当前工作树中保存在`paper/`。原v58不覆盖。

图2f和图3c使用完全相同的浅蓝—深蓝渐变、0—1范围与0、0.25、0.5、0.75、1刻度。图2f仍表示AP；图3c将RMSF原0—100百分位除以100并标为percentile rank。共享映射只统一视觉规则，不使两种科学指标等价。图4d表示有单位的有符号差值，因此使用零中心的蓝—浅白—暖赭渐变，分别保留Å和nats；两条色标均为−0.45至+0.45，全部值完整覆盖，灰格表示无定义或类型改变。

| 图 | 编辑与绘图入口 | 导出文件 |
|---|---|---|
| 2 数值预测与高动态识别 | `02_数值预测与高动态位点识别/plot.py` | `fig2_dense_capability.{pdf,svg,png}` |
| 3 KIX与泛素局部验证 | `03_KIX与泛素局部动态验证/compose_figure.py` | `fig3_kix_ubiquitin.{pdf,svg,png}` |
| 4 MEL5跨界面响应 | `04_MEL5跨界面接触与双动态响应/plot.py` | `fig4_mel5_compact_joint.{pdf,svg,png}` |
| 5 联合读数与构象占据 | `05_联合读数与构象占据/plot.py` | `fig5_color_consistent_joint_readouts_conformations.{pdf,svg,png}` |

所有定量图保留可编辑的矢量文字、曲线与热图。图3/4的结构层来自真实PyMOL光线追踪，包内保存可编辑PSE、真实结构和原生重渲染脚本。六个图3场景及两个图4场景的相机、原子坐标和裁切规则均保留；没有用图像像素替色。

在项目根目录执行`python paper/figures_v59/reproduce.py --check`核查图件及源文件哈希，执行`python paper/figures_v59/reproduce.py 2`重绘图2，将2替换为3、4或5重绘其他图。重绘在临时副本中运行，输出至`reproduced/`，主文件不被覆盖。`publication_style.py`为跨图共享的配色与标度源，`palette.json`保存图2分类颜色。`validate_suite.py`核对78份CSV/NPZ/PDB逐字节相同、热图范围、缺失格和输出完整性；该核对需要同项目的v58包。常规重绘依赖见`requirements.txt`，原生分子结构重渲染另需PyMOL。

`suite_validation.json`、`source_integrity.json`与各图QA记录数据核查；没有新增统计检验或改变统计结论。`document_QA/`是Word渲染及核对中间件，已检查28页概览和四张图所在的11、13、15、17页。`build_manuscript.py`使用工作区内置Python运行，仅更换Word中的四张图片及对应图注的颜色/标度用词，其他OOXML成员保持原字节。
