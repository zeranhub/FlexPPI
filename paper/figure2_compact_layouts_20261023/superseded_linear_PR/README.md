# 图2当前选定构图：a–c

当前图件为fig2_compact_merged.pdf及同名SVG/PNG，完整预览为fig2_compact_merged_preview.png。画布183 × 150 mm，较原203 mm高减少26.1%；全部可见字号至少9 pt。

a的RMSF和Hχ趋势均为实线，保留蓝色下/左轴、深赭上/右轴及无白底的图内MAE与ρ。两种密度的LogNorm范围和透明度公式相同，已在左侧合成一根双色竖向Mass标尺，只显示一个标题和一套刻度，不将两种密度质量相加。

b将三任务共12组AP分布放进共享纵轴。回收a右色条空间后，b数据区由58 mm扩至71 mm，增加22.4%；a仍为48 mm正方形。每任务内从左到右均为Full、Intra、Static、AA，与方法图例一致。

c是一个完整PR组合面板：一个总标题、一个AP行标和共用Recall说明，内部依次为RMSF、总体Hχ、同类型Hχ。每任务自己的完整曲线和真实阳性率基线保持不变，仅最左显示Precision及纵轴数字。

每次绘制逐值核对36,756个AP和1,224个PR顶点，15个源档案、3,063条记录及冻结AP均值不变，无新检验。运行python plot.py重建当前选定a–c，python verify.py核查数据、PDF/SVG、字号和编号。

上一轮5面板图件及脚本在superseded_5panel/中。旧按任务候选fig2_compact_task_panels.pdf仍保留；原7面板排版为previous_7panel_layout.pdf。当前入口只重建已选a–c，旧候选脚本保存在归档。正式v61 Word及上一轮C目录未覆盖。
