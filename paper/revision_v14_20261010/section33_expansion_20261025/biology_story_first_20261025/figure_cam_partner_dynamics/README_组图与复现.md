# CaM伙伴特异动态组图

成图为 `fig_cam_partner_specific_dynamics.png`、`.pdf` 和 `.svg`，宽183 mm、高231 mm；PNG为400 dpi，PDF保留矢量坐标、线条和文字，SVG保留可编辑文字，结构部分嵌入真实PyMOL射线渲染。图的结构层次参照论文现有图3，未修改Word。

图a展示CaM与CaMKKα/eNOS/nNOS/smMLCK肽的结构，图b放大nNOS与smMLCK的L39/L105/M144，图c/d保留四态全部87根共同N–H和44根共同甲基探针，图e/f对照五个强调残基的主链与侧链实验读数，图g保留nNOS与smMLCK之间全部44根甲基变化。全部结构位点均按NMR的成熟CaM序列编号，表中甲基轴和主链读数均来自原始BMRB记录；本图没有加入当前或旧模型预测。

源表在 `source_data/`，包括全部共同探针、误差、逐探针差值、关键位点与实验条件；结构、编号映射和原子完整性检查在 `structure_data/`。`structure_renders/cam_partner_dynamics.pse` 是可继续编辑的PyMOL工程，保存整体和局部场景；射线图、相机与引线投影也保留。`caption_zh.md`、`caption_en.md` 为对应a–g的中英文图注。

图讲的是同一CaM不同结合态中的局部运动分配。结构是沉积第一模型的定位示意，不代表运动集合；eNOS的1NIW保留原始SeMet，不能用其Met侧链化学直接替代NMR样品。nNOS和smMLCK的实验条件匹配更好，作为重点对照；结构构建体仍有末端差别。S²既不等于RMSF也不等于Hχ，不能凭此图声称预测了亲和力、酶活性或功能因果。

复现顺序为 `prepare_structure_data.py` → `build_figure_data.py` → `render_cam_structures.py` → `compose_cam_figure.py` → `audit_final_figure_data.py`。前三项读取本机既有原始PDB/NMR路径，不改原文件；若移动项目，先在准备脚本中更新原始数据路径。结构渲染使用已有PyMOL运行时 `D:/Anaconda/envs/pvb_cpu/python.exe`；数值整理与组图使用科学Python运行时 `E:/Project/FlexPPI/experiments/_sources/local_model_benchmark_20261023/runtime_py312/Scripts/python.exe`。交付包中已经包括成图、clean PDB、渲染图和源CSV，因此查看和编辑现有图不需要重新下载数据。
