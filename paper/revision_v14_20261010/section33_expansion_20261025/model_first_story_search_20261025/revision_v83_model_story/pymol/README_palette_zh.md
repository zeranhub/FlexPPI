# v83 结构图配色与复现

本文件夹只对已有的真实 PyMOL 结构场景重新配色，不改变任何原子坐标、拟合、相机、比例或裁边。CaM N 叶为中灰（#8A9299）、C 叶为浅灰（#BAC1C7）、连接区为更浅灰（#D9DDE0）；结合伙伴肽链为淡紫（#8B72A5）；重点 Met 侧链碳为青绿（#3A8D8C）、硫为金色。所有复合物使用相同的结构组分配色。这套结构配色与数值曲线的复合物状态配色分开。

whole_NHE1 / whole_IQCG 沿用 v82 的共同参考方向、同一正交相机和同一裁边框。detail_focus_M72 / detail_focus_M109 四图沿用此前局部完整链场景，两态在各自匹配对中使用同一相机和比例。局部相机裁切不等于截断肽链，场景中的完整 cartoon 选择均保留。静态结构只是坐标快照，不能代替 MD 动态或柔性预测。

复现命令：D:/Anaconda/envs/pvb_cpu/python.exe render_palette_v83.py。输入 PSE 和原始渲染 manifest 已放在 inputs 下，自包含所需结构场景。脚本使用 pymol2 真实重渲，Pillow 只按原裁边框裁图，不进行像素换色。新配色 PSE 可直接在 PyMOL 打开，已保存六个对应场景。

render_manifest.json 保留所有原始 landmarks 和 crop_box，并记录八个对象逐原子坐标完全一致的验证结果。所有对象的最大坐标变化应为 0 Å。PNG 使用透明背景输出；阅读软件若显示黑底，可覆盖在白色背景上。


全文协调依据：已实际检查论文图1、图2、图4、图5。结果主图共用蓝 #2166AC、暖赭 #B87936、青绿 #3A8D8C；对应本地 publication_style.py 和图5 compose_body.py 可核对。结构伙伴淡紫与图1方法图的 lavender 同色系；关键 Met 直接沿用全文青绿。结构组分颜色与数值状态颜色仍然互不重用。

色板源文件确切位置：E:/Project/FlexPPI/paper/revision_v14_20261010/figure2_dense_20261024/source_internal/publication_style.py（BLUE/WARM/TEAL 三常量），当前图2来源 E:/Project/FlexPPI/paper/revision_v14_20261010/figure2_relayout_20261024/relayout_fig2.py 亦使用相同三色。
