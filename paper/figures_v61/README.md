# v61图1输入、输出标签与分叉箭头修正

图1延续v60的轻调配色，输入框改为`Residue sequence`，静态结构下方改为`Protein or complex`；输出标签简为两行`Side-chain / entropy Hχ`。图注保留其为全部可定义χ角的联合状态熵的定义，并同步将“氨基酸类型”改为“残基序列”。

Distribution head右侧分叉竖线向左移动19个图源单位，两支末端横向箭头从24延长到38个单位，箭头头部从17缩至12，尖端与图表纵轴保留17个单位间隙。完整208个单元格中只修改7个目标单元格；配色、嵌入蛋白图、输出曲线和其余201个单元格均保持原样。

`01_模型框架/fig1_framework.drawio`为可编辑主文件，同目录保存PNG、PDF和SVG及`revision_audit.json`。`revise_fig1_v61.py`从v60源图复现该局部修改，需要Draw.io桌面版；图2—图5与S1仍沿用`paper/figures_v59/`。

当前Word为`paper_v61_20261023_中文修订_图1输入输出标注修正.docx`，主项目中保存于`paper/revision_v14_20261010/`。仅更新图1嵌图和图1图注中的三处文字，原图注的4个Word原生公式保持原样；其他OOXML文件逐字节保留。`document_QA/`保存Word渲染核对记录。旧v60保留。
