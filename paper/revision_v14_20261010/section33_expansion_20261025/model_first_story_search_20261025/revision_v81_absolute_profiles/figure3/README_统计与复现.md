# 图3：两态绝对柔性读数（v81）

NHE1用暖赭，IQCG用蓝色。b为RMSF/Hχ绝对值曲线，MD及Full上下排列，同一指标统一刻度；c为完整接触分组的两态均值，d为全部九个Met，e为两个关键位置的绝对读数。a为两态分别着色的MD与预测绝对Hχ对应，不能用该绝对对应代替对变化预测的检验；ΔHχ的MAE及排序相关仍保留正文。

使用v80完全相同的原始MD复算值与冻结Full第0折预测；没有重新训练或推理。两条记录均在该折训练中，故为结构解释案例，不是留出泛化证明。五个连续100 ns窗口来自各态同一条500 ns轨迹，细线表示最小至最大值，不是独立重复的置信区间。输入伙伴、构型与原始钙条件同时不同，不能建立单因素因果关系。Hχ缺失的Ala/Gly/Pro等位置保持断线；RMSF包含整链拟合下叶间运动贡献。

运行 `python compose_fig3.py`，需要numpy、pandas、matplotlib、Pillow。该目录内source_data与structure_assets已包含绘图所需数值及PyMOL渲染。完整可编辑PyMOL场景仍见v80源文件包。结构坐标不因这次更换纵轴而变化。

输出：450 dpi PNG、嵌入TrueType字体的PDF、可编辑文本SVG；原始物理画布183 × 245 mm。all_plotted_absolute_quantities.csv保留a–d绘制数值。figure_QA.json记录全位置数值和坐标范围检查。
