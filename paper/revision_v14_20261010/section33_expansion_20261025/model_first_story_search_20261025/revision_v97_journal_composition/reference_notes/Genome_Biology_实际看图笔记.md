# Genome Biology 实际看图笔记（2026-10-27）

按 agent-reach 的网页搜索路由使用 Exa 搜索；期刊旧域名的 Jina 读取超时，改用官方 Springer Nature Link 原文。以下明确区分“已实际打开完整图像”与“只下载”。本目录 PNG 保留期刊原图，不是本文的重绘素材。

## 最值得借鉴的两篇

**1. ENGINE：A scalable equivariant graph network framework for precise protein function prediction，Genome Biology 26:407（2025），Fig. 4。** 官方原文：https://link.springer.com/article/10.1186/s13059-025-03886-y 。原图：https://media.springernature.com/full/springer-static/image/art%3A10.1186%2Fs13059-025-03886-y/MediaObjects/13059_2025_3886_Fig4_HTML.png 。已用 view_image 实际打开整个 `engine_fig4.png`。它最贴近当前问题：左边约 36% 宽是六条等宽的全残基位置读数；右边约 62% 宽放六个蛋白的结构，采用三列、两行，各自带局部放大。所有全位置读数共用纵轴范围 0–1，数值配色直接用于结构；一次性图例放在中间的空隙，不在六个小图重复。它说明“结构集中、定量数据独立对齐”也能成立；不必一条曲线旁就重复一个庞大结构。尤其适合满足本次用户要求：两个结构集中放在同一横带，后面的两个全位置数据图用等高、等宽轴框。

但它并非每个细节都值得仿照：原图局部框内残基标签较大、占据比例高，若照搬到只有两处位置差异的泛素，会让结构比数据更抢眼；本文应按用户意见去掉 Q41/Q62 局部放大。另外它的数值是 attention activation，不是柔性，更不是因果重要性；本文不能把 RMSF/Hχ 叫 activation。只能借共同坐标和连续数值的空间映射方式。

**2. MAVE-NN: learning genotype-phenotype maps from multiplex assays of variant effect，Genome Biology 23:98（2022），Fig. 4（主参考）与 Fig. 6（不宜照抄的对照）。** 官方原文：https://link.springer.com/article/10.1186/s13059-022-02661-7 。Fig. 4 原图：https://media.springernature.com/full/springer-static/image/art%3A10.1186%2Fs13059-022-02661-7/MediaObjects/13059_2022_2661_Fig4_HTML.png 。Fig. 6 原图同路径 `13059_2022_2661_Fig6_HTML.png`。已实际打开 `mavenn_fig4.png` 与 `mavenn_fig6.png`。Fig. 4 很凝练：上、下各是一套完整案例；每行左侧约 2/3 宽放全位置效应图，右侧约 1/3 宽放对应的模型—测量关系。所有面板的文字大小、轴框高度和上下边缘完全一致，不靠大标题、粗分割线或彩色背景来制造层次。全序列位置都保留，再用少量标记解释，而不是只画挑出的热点。本文可以借它的固定网格与面积比例：中间两个宽的全位置曲线作为主证据，下部两个约等高的解释性面板；逐χ 2×2 作为一个完整单位保留，避免拆开后丢失比较关系。

Fig. 6 的强机制叙述建立在直接推断的 binding/folding energy 与实验设计上；它的状态图和 kcal/mol 参数不能移植给 FlexPPI。本文有 RMSF 与联合 Hχ，不等于结合或折叠自由能。参考它的“两种读数分工”即可，不仿造三态示意图。

## 两篇用于检验设计原则的计算模型论文

**3. SpatPPI: a geometric deep learning model for predicting protein–protein interactions involving intrinsically disordered regions，Genome Biology 26:339（2025），Fig. 4 与 Fig. 5。** 官方原文：https://link.springer.com/article/10.1186/s13059-025-03820-2 。Fig. 4 原图：https://media.springernature.com/full/springer-static/image/art%3A10.1186%2Fs13059-025-03820-2/MediaObjects/13059_2025_3820_Fig4_HTML.png ；Fig. 5 原图同路径 `13059_2025_3820_Fig5_HTML.png`。两张完整图均实际打开。Fig. 4 按三个横带组织：上行四个同尺度 embedding 视图，中行长蛋白性能曲线加样本量，下行网络例子加模型比较。这种“同一问题的多个读数放在同一行”的原则可借鉴；Fig. 5 上行消融、下行扰动与多构象，对齐也稳定。但其装饰性彩色虚线、重复大标题、多次图例、右侧饼图占位和图内字体差异，都不适合作为本文的审美模板。不能因为论文发表于 Genome Biology 就把全部版式细节当成标准。FlexPPI 不需要新增 embedding、网络或饼图来填空。

**4. High-throughput deep learning variant effect prediction with Sequence UNET，Genome Biology 24:110（2023），Fig. 4。** 官方原文：https://link.springer.com/article/10.1186/s13059-023-02948-3 。原图：https://media.springernature.com/full/springer-static/image/art%3A10.1186%2Fs13059-023-02948-3/MediaObjects/13059_2023_2948_Fig4_HTML.png 。已实际打开完整 `unet_fig4.png`。它采用主面板明显大于支撑面板的层级：上方样本规模和速度，中部物种树与一一对应的逐物种相关，下方三个解释性面板。最有用的是“同一项目顺序贯穿树和数据”，不是物种树本身；重复的文字用一次图例承担。中部有大量真实物种条目，因此允许占更多面积。本文只有 25 个 CaMKI 位点、72 个泛素位置以及 83 条记录，应该让每一类真实数据有自己的合适尺度，不能靠堆多个小图假装同样的数据规模。

## 额外反例

Cross-phyla protein annotation by structural prediction and alignment，Genome Biology 24:113（2023），Fig. 4，官方：https://link.springer.com/article/10.1186/s13059-023-02942-9 ，原图：https://media.springernature.com/full/springer-static/image/art%3A10.1186%2Fs13059-023-02942-9/MediaObjects/13059_2023_2942_Fig4_HTML.png 。已打开 `crossphyla_fig4.png`。只有两个大面板：左边用全体密度作背景，再标少量目标基因；右边是同一组基因跨物种的热图。它的优点是每一面板明确回答一个问题，不靠面板数营造信息量。当前 Fig3 若需减负，可以合并身份重复的结构为一个 a 面板，不削掉用户已认可、实际承接模型与 MD 的逐χ关联。`crossphyla_fig3.png` 已下载但本次未打开，不计入实际观图数量。

## 可直接落实到本文的版面约束

1. 两个整体结构作为同一横带，模型/伙伴身份靠直接标注，结构区面积控制在全图约 1/5–1/4；没有信息增量的 Q41/Q62 细节图删除。结构是定位器，主要证据由数据承担。
2. CaMKI RMSF 与 Ub Hχ 两幅绝对值全位置曲线应等宽、共享明确的左边缘，关键区域用浅带或轴上小标记，不画粗框。RMSF 的 CaM 接触轨道与曲线严格同一位置轴，不能另起漂移的坐标。
3. 图内只保留短对象名、坐标、图例、少量关键数字。去掉重复段落式标题、灰色解释小字、分隔色带、图内表格；方法细节和两个 OOF 协议的差别放图注。
4. 逐χ 2×2 关联作为一个完整单位，用共同横纵轴标签、相同 0 线和统一点样式；保留真实弱关联 χ4，而非只挑 χ1–χ3。旁边的全体 83 条记录统计尺度应与它平衡，不能被挤成细长条。
5. 六个 Gln 的 SASA–Hχ 图是环境解释的辅助证据，不能因为要好看强加趋势线；原始绝对值、共同的位点身份和均等字号比多样图形重要。
6. 不是要照抄哪篇的整张图，而是取 MAVE-NN 的精确对齐、ENGINE 的结构集中及空间/位置对应；明确拒绝 SpatPPI 的装饰线和重复标题。真正的新意应来自 RMSF 与 Hχ 的互补数据链，不来自没有数据支持的机制箭头。

本次实际打开：ENGINE Fig4；MAVE-NN Fig4、Fig6；SpatPPI Fig4、Fig5；Sequence UNET Fig4；Cross-phyla Fig4，共 5 篇、7 张完整组图。前四篇为主要调研对象。
