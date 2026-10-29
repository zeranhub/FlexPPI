本图将原图3中已验证的泛素材料移入补充，不新增模拟，不重新挑选代表帧，也不重复图5的66个Hχ预测剖面。a保留完整输入结构及Q41/Q62定位，b保留两处各三个实际MD保存帧，c只保留全部72个RMSF位置，d保留两处全部27个联合角度类别。原始生产脚本为 `E:\Project\FlexPPI\experiments\_case_validation_sources\compose_figure.py`，结构图包位于其 `ubiquitin` 子目录；各数值读取其 `data` 子目录的已验证源表。

Q41实际帧为2397、401、1171，对应联合类别16、9、15；Q62实际帧为924、1143、1467，对应类别3、4、13。代表帧由旧流程按前5,000个标签帧的类别占据率及圆周平均角距离选定，本轮完全沿用。完整帧出处及几何核验已复制到 `source_data/representative_frames.csv` 与 `frame_geometry_audit.csv`。联合类别覆盖合并5,000帧；Hχ标签则为五个连续窗口的均值，两种统计不混用。

为了统一新稿配色，本轮在原 `ubiquitin_2QHO.pse` 中只重定义UBA、Q41、Q62的颜色，再恢复旧相机进行PyMOL ray渲染。全部对象坐标逐一比较完全不变，局部相机与原视图相同；未人工摆动侧链或编辑结构图片。`structure_renders/recolor_render_audit.json` 保留原PSE SHA256、相机差异和新图SHA256，原素材未改。局部图直接采用旧焦点场景、完整近远裁切面，保留原PSE中显示的全部对象。结构渲染含原场景标签；新增矢量文字和坐标轴均不小于8 pt。

`compose_figS2.py` 复用原生产脚本的全位点剖面与27类累积覆盖逻辑，改为183 × 130 mm四面板布局，固定实际印刷尺寸输出PNG、PDF、SVG。量化元素保持矢量；结构为600 dpi原生PyMOL渲染。`source_data/plotted_all_72_RMSF_positions.csv` 和 `plotted_all_54_coverage_points.csv` 保存实际绘图数组，`figure_QA_and_sources.json` 记录全数值覆盖、图幅、字号、来源校验和及未重复Hχ剖面的检查。

结构重渲染使用 `D:\Anaconda\envs\pvb_cpu\python.exe recolor_existing_pymol.py`；科学绘图使用 `E:\Project\FlexPPI\experiments\_sources\local_model_benchmark_20261023\runtime_py312\Scripts\python.exe compose_figS2.py`。文档专用bundled Python未安装Matplotlib，因此该科学绘图步骤使用已验证的科学Python环境；DOCX编辑仍使用bundled Python。
