# PaperBanana Figure 3 设计输入契约

版本：2026-10-27。以已经核验的 v98 为科学内容基线，只尝试新的视觉组织。本文档是设计输入，不是新的实验结果，也不授权修改现有主图或 Word。

## 一、要解决的问题

目标是把 FlexPPI 论文的“主链柔性实例、侧链多样性实例、跨体系侧链占据验证”整合为一张紧凑、平衡、适合 Nature Communications / Genome Biology 正文的组图。用户认为当前版上方结构图仍可能显得空、整体不够有设计感，希望实际调用 PaperBanana 比较它的设计能力。请把现有图当作内容清单，不要只给现有版加边框、装饰或换色；可以重新分配整张图的面积和阅读顺序，但不能修改科学证据。

图要传达的主线是：完整复合物上下文帮助恢复 CaMKI 调控肽受约束的主链运动；泛素中同一种侧链在不同局部环境中具有不同的构象多样性，模型保留这种差异；联合 Hχ 预测与实际 MD 单角占据和跨 83 条记录的占据描述量相关。它不是结合自由能、激活机制或构象生成的直接证据。先看模型预测是否贴近参考，再看这种预测能解释哪一种局部环境/占据差异。

## 二、优先交付方式

优先让 PaperBanana 输出“完整画布版面方案/设计概念”：白底、面板边界、实际轴标签和图例的摆放、真实图块的位置和相对面积。若它不能从 CSV 精确渲染数字与曲线，就只设计版面，并把所有数据图块标为需要由原生 Matplotlib 脚本精确置换，不能把生成式曲线当作真实证据。理想流程是先得到布局参考，再在原生矢量脚本中实施，以保留数值可核查和文字可编辑。

生成式输出须单独存到本试验目录中，标为“设计试样”，不能直接覆盖 v98、投稿图或 Word。仅当所有曲线、点、坐标、误差范围及结构素材能够逐项对上源文件，才可以称为科学成图。

## 三、完整内容清单：现有 a–f 必须全部保留

| 当前面板 | 必须呈现的真实内容 | 计数与重点 | 可以调整的部分 |
|---|---|---|---|
| a | 1MXE CaM–CaMKI 与 2QHO Ubiquitin–UBA 两个真实整体丝带结构；分别标出 W303/M316 与 Q41/Q62 | 两结构放在一起；完整分子；无局部放大图 | 两结构等比放大、排位、透明边裁切、连线和文字位置 |
| b | CaMKI 全部 25 个残基的 MD、Full、Intra 绝对 RMSF 三条曲线 | 25×3=75 个纵轴值；原位点294–318；W303/M316；Full MAE=0.114859 Å、Intra MAE=2.383203 Å | 坐标框比例、线宽和图例；主线需肉眼可读 |
| c | 泛素全部 72 个序列位置，66 个有定义的侧链 Hχ；MD、Full 两条绝对曲线 | 66×2=132 个值；A/G 缺失处断开；三个 Pro 的真实零值保留；Q41/Q62；Full MAE=0.079531 nats | 同 b；不新增 RMSF 面板 |
| d | 所有 6 个 Gln 的侧链 SASA 与 MD/Full Hχ 配对散点 | 6×2=12 个值；Q41/Q62 重点标注；不能拟合暴露度单调关系 | 配对连接线、标签及与 c 的邻近关系 |
| e | 逐 χ 的 2×2 散点：同一个联合 Hχ 预测对 MD 单角 5° 细箱熵，双方已减去氨基酸类型均值 | χ1 n=66、ρ=0.50；χ2 n=52、ρ=0.45；χ3 n=23、ρ=0.52；χ4 n=10、ρ=0.13；151点×2坐标=302值；Q41菱形、Q2方形 | 整个2×2单元的位置和大小；必须保留四图、全部点及共享轴语义 |
| f | 所有 83 条记录内按氨基酸类型控制的相关性；Full Hχ 与 Full RMSF 分别对应 fine-bin entropy、low-F* width | 83×2预测量×2参考量=332 个相关系数；每组均值和既有配对 bootstrap 95%区间 | 可以改散点/区间的排列和方向；必须看得见记录级离散值及均值区间 |

主图数值总数为 853：75+132+12+302+332。不得为了“美观”删点、只选强关联记录、缩短序列、丢弃负相关或隐藏 χ4。e 中每个真实点有两个坐标，因此 302 是数值数，不能在图注里误写为302个独立样本。

## 四、可开放的版面方案空间

建议至少考虑一种与 v98 有实质区别的整体层级：把两个真实结构作为相邻且足够大的形态入口，用统一基线和接近的视觉重量连接各自的绝对序列曲线；逐χ的完整2×2作为独立整块，不能被挤成四个小贴片；6Gln环境面板和83记录总结面板可在其旁边形成较窄但不弱小的量化区。结构之间的空隙、顶部到曲线的断层、d/f信息与e的面积不平衡，都是可以重新设计的对象。不要用新增机制漫画、圆形装饰、阴影卡片或大色底填满空白。具体可超出上述建议，科学和可读性约束优先。

不必保持旧面板字母；若重排顺序，给出旧a–f到新字母的映射。不要把“当前a–f”视为固定六个相等矩形；a本身是两个相邻真实结构，e本身是一个2×2单元。最后仍须能清楚分辨每一条证据。

## 五、已确定的用户审美与出图约束

1. 结构图放在一起，整体丝带表示位置；Q41/Q62局部放大已明确删除，不要恢复。结构只能使用所给真实 PyMOL 素材，不生成、改画、补全、拉伸分子或伪造构象。
2. 主链实例只突出 RMSF；侧链实例只突出 Hχ。绝对序列纵轴从0起，不用 ΔRMSF/ΔHχ 取代。逐χ散点的偏离值有科学必要，继续保留正负值和原范围。
3. 不放表格、不加大片说明文字、不加灰色小字。背景、方法细节、机制推断写图注和正文。仅保留必要短标题、轴、图例、统计数及目标位点。
4. 与全论文现有配色协调。结构：被分析链紫色 #8B72A5、伙伴浅灰 #BAC1C7；目标位点青绿 #3A8D8C。数据：Full蓝色 #2166AC、Intra棕色 #B87936、MD近黑 #222222。Q41/Q2在逐χ图均用青绿、靠菱形/方形区分。Hχ/RMSF在83记录图可以保留同蓝色、实心圆/空心三角区分。不能让同一颜色同时含混代表复合物身份与模型方法。
5. 以白色底、黑色文字、左/下细轴、无背景网格为主；必要零线淡虚线。统一对齐、字体、线宽和面板间隙，避免粗框与阴影。
6. 当前画布179×200 mm，在Word中以164 mm宽插入后最小字约8.02 pt。新方案以双栏179–183 mm宽为目标，高度尽量接近200 mm；不可靠缩到不可读实现紧凑。最终印刷字号不少于8 pt（按164 mm宽缩放后检验），标签和图例不压住数据，结构不被切断。
7. 最终科学图由原生脚本导出可编辑文字的SVG、嵌入字体PDF和至少300 dpi PNG；PaperBanana若仅输出位图设计，不宣称可直接用于投稿或可编辑。

## 六、科学解释的硬约束

Full 与 Intra 是使用同一结合态结构输入、不同上下文范围的两个模型，不是游离态与结合态 MD；不能画“结合前→结合后”或把两条预测曲线的差称为实测结合诱导变化。1MXE只含调控短肽，没有完整 CaMKI 激酶催化域，不能用此图证明激活或自抑制解除。

Hχ是联合侧链构象状态熵；逐χ图复用同一个联合预测，没有四个独立预测头。双方做氨基酸类型均值校正是为了避免侧链类型的固有差异主导相关性。χ1–χ3为中等正相关，χ4样本少、相关弱；不能用设计把四图都表现为“强相关”。

6个Gln中的SASA与Hχ不是单调关系，不能加回归线或给出“越暴露一定越灵活”结论。Q41/Q62是局部环境对照，Q41/Q2是按同一种残基的联合Hχ预测选出的低/高位点；这两对不能混为一组。Q41/Q62都没有≤5Å的UBA直接重原子接触；不能将其差异标成直接伙伴接触效应。

low-F* width是MD角度占据推导的无量纲相对势 F*=−ln[p/pmax] 低值范围的角宽，不是结合自由能，也不是模型生成的自由能面。83条记录来自多个来源，不能改写为83个同类蛋白、83个复合物或83次独立湿实验。当前主图不新增F*曲线，既有逐角MD曲线在补充材料。

## 七、输入文件与列映射

统一基线目录：
`E:/Project/FlexPPI/paper/revision_v14_20261010/section33_expansion_20261025/model_first_story_search_20261025/revision_v98_structure_fill`

| 文件（相对于基线目录） | 用途与主要列 |
|---|---|
| `figure3/fig3_dual_case_entropy.png` | 当前完整图；SHA256=75733ef569e662e0861ccf15a499e945700033b5e4ec31d601674029b35379a1 |
| `figure3/fig3_dual_case_entropy.pdf` / `.svg` | 当前矢量主文件，保留轴/文字/线条的参考 |
| `figure3/compose_fig3_v98.py` | 精确重绘和现有面板坐标；不在试验阶段覆盖 |
| `figure3/figure_QA.json` | 计数、字号、裁边、关键MAE及科学声明审计 |
| `figure3/source_data/all_plotted_values.csv` | 当前主图全部853个纵轴/散点坐标值的可核对长表 |
| `figure3/source_data/cam_all25.csv` | b：`CaMKI_native_position`,`target_rmsf`,`pred_rmsf`,`intra_rmsf`；全部25行 |
| `figure3/source_data/ubiquitin_72_position_profile_with_original_mask.csv` | c：`residue_id`,`valid_hchi`,`MD_Hchi_nats`,`Full_Hchi_nats`；72行及原始缺失掩码 |
| `figure3/source_data/ubiquitin_all6_Gln_Hchi_SASA_plot.csv` | d：`residue`,`residue_id`,`sidechain_SASA_A2`,`MD_Hchi_nats`,`Full_Hchi_nats`；全部6行 |
| `figure3/source_data/defined_angles_long.csv` | e：`chi`,`author_residue_id`,`predicted_Hchi_nats_AA_centered`,`fine_entropy_four_origins_mean_nats_AA_centered`；151行 |
| `figure3/source_data/per_chi_correlations.csv` | e统计：筛选`predictor=predicted_joint_Hchi`,`metric=fine_entropy`；使用`chi`,`rho_AA_centered` |
| `figure3/source_data/retained_absolute_correlations_all498.csv` | f：筛选`predictor in {full_hchi,full_rmsf}`；`metric_key in {fine_entropy,lowF_width}`；每组83行；`sample_id`,`rho_type_centered` |
| `figure3/source_data/retained_panel_statistics.csv` | f区间：筛选`panel_kind=correlation`和上述预测量/目标；`mean`,`bootstrap_95_low`,`bootstrap_95_high` |
| `pymol/cam_overview.png` / `pymol/ub_overview.png` | a：原生真实结构PNG；带透明边，可等比裁透明边放大 |
| `pymol/render_manifest.json` | a：PNG/PSE/PDB关系、W303/M316/Q41/Q62屏幕定位点 |
| `pymol/cam_overview.pse` / `pymol/ub_overview.pse` | 可复现的原生PyMOL会话，仅在原生结构程序中调整相机 |
| `pymol/inputs/1MXE_CaM_native4_147_CaMKI294_318.pdb` | 真实1MXE结构、原位点编号 |
| `pymol/inputs/A003R_input_structure.pdb` | 真实2QHO对应输入结构 |

`contact_counts_retained_outside_main.csv`、`retained_prediction_groups_all498.csv`是留存的补充分析源，不属于本次必须放进主图的内容，不应无故添加面板。

## 八、通过/否决标准

通过：能指出具体改进——真实结构更充实但不挤压关键曲线；主链/侧链两案例关系清晰；e的2×2完整可读；f的83记录有足够尺寸看离散性；图例能一眼解码；留白统一而有层级；科学内容与上述853个冻结数值一致。

否决：画出新的蛋白构象；拟造或重描散点/曲线导致数值变化；丢χ4或弱/负点；把曲线改成仅Δ值；主图重新塞局部放大/表/机制漫画；为填空增加不对应数据的箭头、结合过程或活性图标；重新计算或美化ρ；用平滑/平均替代真实残基或记录；把生成式试样当作投稿成图。

## 九、可直接提供给英文模型的短输入

Design a substantially improved publication-quality multi-panel Figure 3 for a protein flexibility prediction paper. Treat the attached v98 figure as an immutable scientific content inventory, not as artwork to embellish. Aim for a compact Nature Communications / Genome Biology style, white background, aligned typography, semantic colors, meaningful hierarchy, and enough room for every real data panel. Keep the two authentic molecular ribbon views together, showing CaM–CaMKI (1MXE; W303/M316) and ubiquitin–UBA (2QHO; Q41/Q62). Do not generate molecules, conformations, mechanisms, plots, points or values. Use supplied genuine molecular renders and exact CSV data, or output only a layout concept with replaceable plot slots if exact data rendering is unavailable. Keep all 25 CaMKI residue RMSF values for MD/Full/Intra; all 72 ubiquitin sequence positions with 66 valid MD/Full Hchi pairs and original missing gaps; all six Gln environment pairs; the complete per-chi 2x2 scatter grid (n=66/52/23/10; rho=.50/.45/.52/.13); and all 83 records in each of the four correlation groups with existing means and bootstrap intervals. There are 853 frozen plotted numeric values. The same predicted joint Hchi value is compared with each MD single-angle entropy; there are no four angle-specific prediction heads. Full/Intra are models on the same bound input, not apo/bound states. Do not imply binding free energy or a proven activation mechanism. Remove no weak points and fit no SASA trend. Keep Q41/Q62 environmental and Q41/Q2 low/high prediction comparisons distinct. No tables, local structure zooms, mechanisms, gray explanatory microtext, shaded cards, decorative borders or invented experiments. Preserve purple analysed chain / gray partner / teal site highlights, blue Full / brown Intra / black MD. Give the per-chi grid substantial readable space; keep both complete sequence profiles visibly tied to the two adjacent structures. Target 179–183 mm width and about 200 mm height, with at least 8 pt text after scaling to 164 mm. Return a design trial or layout specification separately from scientific production assets. Do not overwrite v98 or the manuscript.
