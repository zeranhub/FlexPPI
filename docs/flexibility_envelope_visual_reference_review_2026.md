# 蛋白质柔性动态“包络图”文献与 FlexPPI 构图建议

## 结论

有，而且这是分子可视化中一条相当明确的路线，只是文献不一定都使用中文直译的“包络图”。最常见的检索词是 **conformational density / conformational fuzziness / flexibility isosurface / probabilistic molecular surface / fuzzy molecular surface / occupancy volume / sausage representation / putty representation / conformational ensemble overlay**。其中最符合所问图形的是 Kozlíková 等人的综述 Figure 12：同一离子通道左边是随局部柔性改变半径的 backbone tube，右边是把柔性转成三维等值面的包络；论文图注明确说柔性红色区域占据更大空间，刚性蓝色区域更集中。这说明“用空间包络把动态压缩进一幅静态蛋白图”不是临时拼出来的视觉隐喻，而是已有分子可视化方法类别。

## 文献脉络与图形类型

第一类是 **可变管径的 sausage/putty**。每个残基仍沿代表结构的主链定位，但局部管子越粗，表示该残基的波动或不确定性越大；颜色常从蓝色刚性过渡到红色柔性。它的优点是和逐残基 RMSF 一一对应，读者能快速定位 flexible loop、terminus 或界面邻近区域；缺点是它只显示“幅度”，不显示构象往哪个方向移动。Kozlíková 等人的 Figure 12 左图和 PyMOL B-factor putty 都是这种语言。第二类是 **概率密度或占据体积等值面**。先将对齐后的 MD/NMR 构象采样到三维网格，再绘制给定占据概率或密度阈值的表面；Lee 与 Varshney用多层透明表面表现原子热振动和位置不确定性，Schmidt-Ehrenberg 等人则把动态构象采样到网格后做等值面或体渲染，Kozlíková Figure 12 右图是这一思路的直观示例。第三类是 **构象集合叠加**，也就是把 NMR 模型、MD 帧或代表构象画成 bundle/cloud；Rico Figure 3C 和 Fuxreiter Lab 的图能明显展示刚性核心彼此重合、柔性环或无序链发散。第四类是 **特定组分的采样体积**，如 Lemmin 等人把 HIV envelope glycan 的构象异质性画成围绕蛋白的三维体积。这一类最有“云团/包络”的视觉冲击力，但表达对象可能是糖链、配体或侧链，不应不加区分地移植成主链柔性证据。

## 为什么它特别适合 FlexPPI

FlexPPI 的研究对象不是给整条蛋白一个“软或硬”的标签，而是预测每个残基在特定伙伴和复合物环境下的柔性，因此三维包络能直接把任务的必要性可视化：晶体结构或预测结构给出的只是一条代表坐标，而每个残基实际允许访问的空间大小并不相同；结合伙伴遮挡、稳定或重排某些局部区域后，这个空间分布也会改变。最推荐的 motivated example 是把相同目标链在 unbound/partner A/partner B 三种条件下保持同一视角：内部放一条相同的代表 ribbon，外部放半透明包络，界面处再用轮廓线或残基标记突出。这样读者看到的不是抽象 RMSF 曲线，而是“同一平均结构周围，不同位置和不同伙伴导致可访问构象空间发生局部收缩或扩张”。这恰好回答“预测逐残基柔性有什么用”：它把一个静态复合物坐标补充成空间分辨的动态状态图，可用于界面适应、别构传播、突变效应和动态热点定位。

## 推荐的 FlexPPI Figure 1 视觉方案

建议中央使用浅灰或浅青色 representative structure，partner 使用低饱和度的另一种颜色；目标链的包络使用 20%–35% 透明度，低柔性为窄且贴近骨架的蓝紫色区域，高柔性为扩张的橙红色区域。逐残基值可以通过管径表达，也可以通过颜色表达，但不宜同时让颜色、半径、透明度和表面高度都映射同一个值，否则图看起来更“强”却难以解释。一个稳健组合是：主结构颜色表示链身份，tube 半径表示预测 RMSF；外层一个固定概率水平的半透明等值面表示构象占据空间；界面残基用细线或小球单独标记。若只有逐残基 RMSF 而没有轨迹，建议只画沿 Cα 轨迹的各向同性 envelope，并在图注写明它是从预测幅度构造的示意包络；如果有 MD 帧或生成式构象样本，则先对目标链的稳定核心做结构对齐，再在统一网格上估计原子或 Cα 占据密度，并明确等值面的概率/密度阈值。只要图注把“输入量—映射规则—物理含义”写清楚，这种图不会只是装饰。

## 三个必须避免的过度表述

第一，**B-factor 不等于纯粹的内禀动力学**。它还受到晶体堆积、模型误差和静态构象异质性的影响，因此可用来说明局部位置不确定性，但不能无条件作为溶液 RMSF 的真值。第二，**NMR 模型集合不等于按真实时间顺序采样的轨迹**。它表示与实验约束兼容的一组结构，适合说明构象异质性，但不能直接推断动力学速率。第三，**逐残基标量不包含方向信息**。如果 FlexPPI 输出的是一个数值，它可以控制包络的局部厚度，却不能支持定向箭头、椭球主轴或多叶构象云；这些需要各向异性协方差、normal mode 或实际构象集合。

## 最值得先看的本地图

1. [Kozlíková 2017 Figure 12：可变管径与柔性等值面](literature_figures/flexibility_envelopes/kozlikova_2017_fig12_variable_tube_and_flexibility_isosurface.png)：与“包络图”问题最吻合。
2. [Lee & Varshney 2002：crambin 模糊分子表面](literature_figures/flexibility_envelopes/lee_varshney_2002_crambin_fuzzy_surface.gif)：最适合学习多层透明概率表面。
3. [Rico 2013 Figure 3：B-factor 与 NMR 集合](literature_figures/flexibility_envelopes/rico_2013_fig3_bfactor_and_nmr_ensemble.jpg)：最适合支撑“蛋白不是单一静态结构”的背景。
4. [Lemmin 2017 Figure 1：构象采样体积](literature_figures/flexibility_envelopes/lemmin_2017_fig1_glycan_sampling_volumes.jpg)：最直观的动态空间占据体积。
5. [Fuxreiter：fuzzy interaction ensemble](literature_figures/flexibility_envelopes/fuxreiter_fuzzy_protein_interaction_ensemble.jpg)：最直观的结构云，但偏向无序/模糊复合物。

完整图片、来源、适用边界和版权提示见 [本地图片索引](literature_figures/flexibility_envelopes/README.md)。

## 关键参考文献

- Kozlíková, B. et al. *Visualization of Biomolecular Structures: State of the Art Revisited*. **Computer Graphics Forum** 36, 178–204 (2017). https://doi.org/10.1111/cgf.13072
- Lee, C. H. & Varshney, A. *Representing Thermal Vibrations and Uncertainty in Molecular Surfaces*. **Visualization and Data Analysis 2002** (2002). 官方项目页：https://www.cs.umd.edu/gvil/projects/probabilistic_surfaces.shtml
- Schmidt-Ehrenberg, J., Baum, D. & Hege, H.-C. *Visualizing Dynamic Molecular Conformations*. **IEEE Visualization 2002** (2002). https://doi.org/10.1109/VISUAL.2002.1183780
- Rico, F., Rigato, A., Picas, L. & Scheuring, S. *Mechanics of proteins with a focus on atomic force microscopy*. **Journal of Nanobiotechnology** 11(Suppl 1), S3 (2013). https://doi.org/10.1186/1477-3155-11-S1-S3
- Lemmin, T. et al. *Conformational Heterogeneity of the HIV Envelope Glycan Shield*. **Scientific Reports** 7 (2017). https://doi.org/10.1038/s41598-017-04532-9
- Janowski, P. A. et al. *X-ray refinement significantly underestimates the level of microscopic heterogeneity in biomolecular crystals*. **Nature Communications** 5, 3220 (2014). https://doi.org/10.1038/ncomms4220

## 检索说明

本轮同时使用分子可视化术语检索、跨数据库论文检索和出版社/作者机构页面核验。自动论文检索对这个问题噪声较高，因为 “surface”“envelope”“flexibility” 会召回很多无关几何和小分子论文，因此最终证据以论文原文、图注和官方项目页为准。Exa 免费 MCP 在本轮触发 429 限流；OpenAlex 连接器因本地类型兼容错误跳过；Semantic Scholar 多次等待限流后没有返回记录。这些错误不影响已从 Wiley、Springer Nature、University of Maryland 等原始来源核验和下载的图片。
