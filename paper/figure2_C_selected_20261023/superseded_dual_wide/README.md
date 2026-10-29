# 图2：选定C后的双坐标叠加

用户选择C，并明确要求两个指标画在同一块绘图区。a现采用完全重合的两个坐标变换：蓝色RMSF（Å）对应下轴、左轴，暖赭色Hχ（nats）对应上轴、右轴。两组原始带权密度以浅透明颜色叠加，趋势线置于密度上层；RMSF趋势为实线，Hχ趋势为虚线。每种颜色有单独的质量色条，采用相同原始对数质量范围；叠色不是第三种统计量。原有范围、RMSF在1 Å内线性/以外对数以及Hχ线性刻度均保留，没有换成百分位或无量纲值。

中排三个记录AP分布编号为b–d，底排三个完整PR为e–g，数据和排版保持上一轮C图件。a的颜色表示输出指标，b–g的颜色表示方法，方法图例已移到b–g上方以明确作用范围。所有记录、密度数据、AP及PR数值沿用所选C版。原C排版保存在`previous_C_layout.pdf`；上一轮误解为左右分面的排版与脚本保存在`superseded_grouped/`，已经被本次叠加图替代。

`fig2_C_selected.pdf`是矢量主文件，`fig2_C_selected.svg`可编辑，`fig2_C_selected.png`为450 dpi大图，`print_preview.png`按文稿167.64 mm宽度渲染。画布183 × 203 mm，最小字号9 pt，PDF内嵌Arial。运行`python plot.py`可从同目录`source_archive/`重建图件；数据检查见`verify.py`、`final_QA.json`及`QA_report.md`。图注和统计说明分别见`caption_zh.md`及`statistics_and_reproducibility.md`。

`panel_a_preview.png`提供a的放大预览；本次更新图件，正式v61 Word保持原样。
