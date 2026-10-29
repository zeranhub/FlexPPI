# 模型支持的CaM组图

主文件为`fig_cam_model_supported_story.pdf`（183×240 mm矢量图，结构图片嵌入）与同名400 dpi PNG、可编辑SVG。原论文未改动。先看父目录`模型先行_已核实案例.md`的具体生物学问题、模型范围与全部数值。

结构来自当前模型实际读取的模拟输入PDB，不用实验结构中的另一个模型替代。`render_structures.py`只改副本链号和CaM编号，CaM4–148号逐位与成熟序列相符，去掉额外表达尾Ala149，没有重建原子或添加Ca。两整结构按N叶瓣对齐；局部图在独立副本中按对应主链对齐，同一局部对采用相同视角和裁切，因此整结构场景不会被局部叠合改写。`structure_renders/CaM_model_story.pse`保留可编辑场景。

`compose_figure.py`读取全部267行主表（145个主链位置、122个有效侧链），散点122个不下采样且全部在坐标轴范围内。四个接触集合互不重叠且数量合计122；九个Met全部保留。没有把S²混入MD比较，没有重标定模型，也没有把时间分块误画成置信区间或独立重复。全位点曲线与Met图保留模型不准确的位置。

已逐项核查结构编号、坐标保持、MD与前向NPZ身份和数值；主表对原NPZ核验及独立审查见父目录`partner_screen/final_cam_checks.json`与`regional_screen/cam_state_independent_audit.json`。原轨迹Hχ复现误差<1.2×10⁻⁷ nats、RMSF<10⁻⁶ Å。完整训练权重路径和SHA保存在`partner_screen/training_member_completion.json`。

视觉检查已在两轮PNG中完成，修正了散点坐标范围和底部注释与刻度重叠。最终文本无超出画布，PDF仅含嵌入Arial常规/粗体TrueType，SVG文字保持为可编辑文本。颜色统一，单位、Δ方向、训练范围与轨迹来源均在图内和图注中说明。图中不同探针/模型间未进行数值尺度转换。

结构复现环境：`D:/Anaconda/envs/pvb_cpu/python.exe render_structures.py`；数据与组合：`E:/Project/FlexPPI/experiments/_sources/local_model_benchmark_20261023/runtime_py312/Scripts/python.exe compose_figure.py`。使用本地Matplotlib、PyMOL，不使用图像生成工具。源表副本位于`source_data/`；整个案例包保留原始MD复算表、模型前向输出与复现代码。
