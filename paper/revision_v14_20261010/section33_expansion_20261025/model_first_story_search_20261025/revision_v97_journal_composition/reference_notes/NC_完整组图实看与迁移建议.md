# Nature Communications 模型案例组图实看记录（2026-10-27）

此次按 agent-reach 的 Exa 搜索与 Jina 网页阅读路线访问四篇官方图页，实际用 view_image 逐张打开四张完整原图。RMSF-net 与 ESMDynamic 原图此次从官方 CDN 重新获取成功；PeSTo 与 ProFlex 的 CDN 回应并非有效图片，已保留该回应并改用此前官方原始图片。四个图页本轮均成功重新获取。catalog.json 记录来源、像素尺寸、SHA256、查看状态，所有保存图都与实际查看的版本逐字节相同。下述比例是视觉估计，用于设计取舍，不冒充期刊强制规定。

## 1. RMSF-net，图4（2024）

论文：[Accurate Prediction of Protein Structural Flexibility by Deep Learning Integrating Intricate Atomic Structures and Cryo-EM Density Information](https://www.nature.com/articles/s41467-024-49858-x)。[官方图4](https://www.nature.com/articles/s41467-024-49858-x/figures/4)。本地图 rmsfnet_fig4.png。

图分成上下两个案例单元，细虚线只分隔案例。每个单元先用一条贯穿整宽的残基 RMSF 图呈现完整读数，再以三列同方向、同尺寸结构对应三种读数，结构统一使用同一种数值色标。局部图紧贴所属结构下方。a、b 只标记案例，不把每个小视图都编号；方法名称直接放在列顶。结构及局部图占多数面积，因为其颜色本身呈现真实预测数值。可借“一个案例一个完整证据单元”和空间、定量读数呼应；不能照抄三列结构，因为本文目前整体结构的颜色表示组分，连续复制不会增加信息。

## 2. PeSTo，图3（2023）

论文：[PeSTo: parameter-free geometric deep learning for accurate prediction of protein binding interfaces](https://www.nature.com/articles/s41467-023-37701-8)。[官方图3](https://www.nature.com/articles/s41467-023-37701-8/figures/3)。本地图 pesto_fig3.png。

左上是一幅性能曲线，其余大部分是结构案例。b、c、d 中同一复合物两幅结构维持同一方向、相近尺度，只改变界面预测的含义；色标只放一次。伙伴分子黄、绿与蓝白红的预测色标区分清楚。行间留白较小，标签都落在相似的左下、右下位置，结构边缘对齐形成隐形网格，几乎不用边框。可借统一尺度、共享图例、紧凑对齐；不能照抄结构占比，因为本文的静态组分结构只有定位功能，应该让数据占更大面积。

## 3. ESMDynamic，图5（2026）

论文：[ESMDynamic: Fast and accurate prediction of protein dynamic contact maps from single sequences](https://www.nature.com/articles/s41467-026-76361-2)。[官方图5](https://www.nature.com/articles/s41467-026-76361-2/figures/5)。本地图 esmdynamic_fig5.png。

两行分别讲两个案例，左列都是完整模型输出图，右列都是对应实验结构；大数据图与两个状态结构的面积大致相当。关键位置在左图直接用 B、C 或局部边界索引到右图，右侧高亮沿用同一种绿色，所以读者可以沿“数据发现—空间位置”移动视线。这里只有关键链条，没有统计图大拼盘。最值得借的是数据和结构中的同一位置使用一致标记。本文可令 W303、M316、Q41、Q62 在整体结构和绝对数值曲线中同色，且将两种结构集中成一个短行；不可借用双状态结构来暗示 Full、Intra 是游离态、结合态。

## 4. ProFlex，图6（2025，作为不照搬的例子）

论文：[ProFlex as a linguistic bridge for decoding protein dynamics in normal mode analysis](https://www.nature.com/articles/s41467-025-64103-9)。[官方图6](https://www.nature.com/articles/s41467-025-64103-9/figures/6)。本地图 proflex_fig6.png。

左侧两个分布面板共享类别顺序和横轴，右侧大规模逐位置色块图；优势来自同一批结构簇在不同读数中的对应。实际打开后，右侧颜色块极密，而与左侧箱图的视觉联系没有前几图直接，局部色块也不能单独解释本文模型。因此不建议为追求“信息丰富”照着增加大量热图或矩阵；本文每个位点只有一个或两个实数输出，应保持连续曲线和散点的可读性。

## 可直接用于下一版图3的三条规则

1. **按证据分配面积。** 用户已要求两个整体结构放在一起，应将它们集中于顶部同高、同基线的一行，作为空间定位；控制在整图高度约 20%–25%。25 位点 CaMKI 的 RMSF 与 72 位点泛素 Hχ 用绝对数值曲线占中部主要面积。不要按“一个子图一个等大盒子”平均分配。
2. **只设一套视觉语法。** 组件用紫、灰，预测方法用既有蓝、暖棕，关键残基统一青色；MD 黑虚线。相同语义只出现一个图例。所有字母统一字号和左上位置，面板标题简短，轴标题与单位保留，灰色脚注、边框、说明性小字移到图注。局部图若不能增加定量或接触证据就删除，符合用户最新要求。
3. **保留完整数据并让重点自然浮出。** 全 25 位点和全部有效泛素侧链读数不删；W303–M316 仅用很淡的带表示，同色位点标记贯穿结构与曲线。下部保留用户认可的逐 χ 2×2 散点，四小图共用轴标题、同一尺度与统一图例；83 体系验证单独占邻近的短窄面板。避免再把标题、r、n、图例各占一整行。观察到的相关程度和 χ4 的弱关联必须原样保留。

建议验收时不只看原尺寸 PNG：以 Word 实际嵌入宽度 164 mm 查看，若不能不放大读清轴字，就继续简化。期刊图有借鉴价值的是读图路径、证据层次、对齐和一致编码，不是复制它们的数据类型或增加结构图数量。

## 运行说明

agent-reach check-update 已执行：当前 v1.5.0，最新 v1.6.0。本轮未升级任何工具。该更新信息仅供维护记录，不影响本次已成功获取的官方图与图页。
