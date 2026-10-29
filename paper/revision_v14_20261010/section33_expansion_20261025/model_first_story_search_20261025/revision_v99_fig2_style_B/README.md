# 图3 v99：B布局与图2统一样式

figure3/compose_fig3_v99.py是可直接运行的原生绘图脚本；同目录source_data包含完整冻结源表，publication_style.py为论文图2的真实共享配色模块。pymol为不变的真实结构资产。以安装numpy、pandas、matplotlib、Pillow的Python运行该脚本，即在figure3目录重建PNG/PDF/SVG。

本版采用PaperBanana B概念的a/b、c/d、e/f空间组织，未使用其合成分子、曲线、散点或标签。样式基准是当前v98 Word内图2的实际图像，已通过哈希匹配。figure3/figure_QA.json保存样式、字号、数值计数、裁边、嵌入字体和独立视觉审查结果；statistics_and_reproducibility.md给出统计单位与源列；figure3_caption_zh.md是包含共享图例及d空心/实心定义的图注草稿。本次不修改正文或Word。
