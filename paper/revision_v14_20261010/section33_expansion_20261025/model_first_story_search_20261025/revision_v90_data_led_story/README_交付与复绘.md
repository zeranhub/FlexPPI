# V90 图3数据呈现与机制推断

本版彻底撤去V89的c机制示意、回接虚线和预测活性箭头。主图仅保留真实结构背景和数据；机制解释留在正文，并明确属于依据结构、柔性与既有功能实验提出的假说。

figure3/fig3_data_led_story.pdf为183×154mm矢量主文件，SVG文字可编辑，PNG450dpi。相比V89高度178mm减少24mm。主图a保留1A06/1MXE真实结构及F298放大；b完整保留25位RMSF和21个有效Hχ的Full/Intra曲线。左下Full: expanded y scale只是将同一Full的全部25个RMSF值在0–2.3Å纵轴范围再画一次，以看清中央低值变化，不是新数据或独立证据。原双曲线全范围0–8.2Å也保留，无断轴、无裁值。

source_data/all_plotted_dynamics.csv含92个唯一读数，repeated_Full_RMSF_expanded.csv含25个重复放大读数。独立核查均与V89冻结数据精确相等。无新训练、推理或MD；Intra仍为同一结合态坐标上的链内模型预测，不代表游离肽。

figure3/compose_fig3_1MXE.py从脚本目录读取相对路径，依赖numpy/pandas/matplotlib/Pillow，可重导PDF/SVG/PNG。结构灰色、调控段紫色、功能位点青色；Full蓝圆实线、Intra暖赭方点虚线。最小7.9pt，在Word164mm宽下约7.08pt。

pymol/保存真实结构输入、渲染PNG、PSE、脚本、全部距离与投影。图a继承V89已验证资产，未改变坐标。1A06仅观测至316，1MXE不含激酶催化域；F298局部图突出较近5个邻居，完整6对象（含Y215）均在距离CSV/JSON。

narrative/3.3_定稿.md和figure3/caption_zh.md为本版最终正文和图注。build_v90_manuscript.py以V89为源，局部更新p71/p73/p75/p76与p72的图像，其他段落、6表和其他图片保持相同。Word仍29页；10–12页实际看图核查，其余26页与已完成视觉检查的V89逐像素相同。render/的包装脚本使用本机Word只读导出与内置render_docx.py生成内部QA页面。

data_audit、data_audit_only、supplement_1MXE为继承的模型和MD证据，用于复核主图及图S3。它们不属于新模型结果。核心W303–M316按结构锚点预先定义，其14位Full RMSF全部0.447–0.689Å；原文只据连续受约束和侧链位置差异结合文献[29,30]作机制推断，没有直接测得激活或抑制性回接。数据审计中的5个时间窗口来自同一条轨迹，不是独立重复。
