# 图1审美方向复审

用户已否定v46：新图视觉不如上一版。此前独立报告确认的是无截字、无碰撞和实际尺寸可读；把这些技术验收条件当成整体美观的通过依据，是审查标准不充分。本轮不能再以“没有重叠”推导“更好看”。v47应回到用户认可的原图总体构图，改进内部组件。

原图的优势是真实蛋白、输入小图、层叠编码器、内部读出层和输出曲线共同形成充实而有层次的顶排，下方b/c/d三列各有稳定的纵向结构。v46同时减少层级、扩大空白，并把b的真实空间问题变成平面N–Cα–C三点图，虽然更简洁，却削弱了科学框架图应有的实体感和结构复杂度。v47保留原来的上下两排与下三列；顶排优先调整组件材质、几何比例、线宽、色阶和字重，不再另起一种轻量流程图。

编码器恢复前后错位的三至四层薄片，位移约1–1.2 mm即可，正面仍画IPA内部残差与FFN关系。正面留小型节点/注意力图，避免只有两个宽文字块。读出层保留LayerNorm、Linear、SiLU、Linear的内部次序，层间距与箭头统一，让这些细节成为主体轮廓内的第二层信息。Residue glyph可用节点附着短特征条；pair glyph用有明暗差别的小网格和轻错位通道面，所有通道面只作为示意，不能冒充测得注意力或新增计算。

b应以真实蛋白的局部空间为主：真实whole_chain加两处细坐标轴和i/j标识；一处连接到约20–25 mm的真实骨架球棍放大。球棍应有柔和高光、圆柱形键和前后遮挡，N/Cα/C标在外缘，坐标轴从Cα发出；不要再让一个平面直角三角形占主要视觉。可保留很淡的全局轴、一条细t_i线和变换式，但不要让全局轴穿过原子细节。长残基名可换i/j，几何必须从当前PDB和帧定义生成。

c继续用纵向的分数输入、softmax、加权汇总、返回query frame和merge/project，三个来源各配小图而不是纯文字大盒。d回到底排一列，三张相同尺度的真实MD图可轻叠或并排，连接RMSF与joint Hχ量测及噪声感知NLL；保留MD和模型两条正确来源，但不要为此横铺全宽。

## 已实际打开的可复用资产

- 真正v13原项目：`C:\Users\DELL\.codex\worktrees\4666\FlexPPI\paper\revision_20261002\figures\main_revision_v13_20261008`。其中`assets/whole_chain.pse`、`whole_chain.png`、`source_data/itrace_4g9f.pdb`、`local_frame_candidates.json`、投影receipt和原Draw.io均存在，当前蛋白与几何应从此重用。
- `E:\Project\FlexPPI\paper\figures\fig_flexppi_main_camki_dual\sidechain_0.png`已实际打开，其球棍的高光和遮挡可作为材质参照，但属于不同结构，不能直接当作当前残基帧。
- `E:\Project\FlexPPI\paper\figures\fig_flexppi_main_reconstructed_assets\frames_clean_2_target_tight.png`已实际打开，分辨率及示意几何不足，不建议当最终b的主体。
- `E:\Project\FlexPPI\docs\literature_figures\figures\pesto_2023_fig1.png`是本地已有文献图，已实际打开。可借鉴真实原子渲染呈现空间关系以及完整内部层次的方式，不复制图片或未验证几何；此处不把它冒称Genome Biology的期刊要求。

## 下一次审查必须回答

必须把原图和候选图按同一插入宽度实际对照：原图被认可的充实节奏是否保留；编码器与读出是否仍有清晰层级；b是否看起来像真实蛋白几何而非课堂直角坐标练习；真实渲染、glyph、盒子与文字是否在同一材质与线宽体系；改善是否来自组件质量而非只有删减。若任一项明显退步，审美结论应写需要修改，即使无重叠、字号合格也不能写整体通过。
