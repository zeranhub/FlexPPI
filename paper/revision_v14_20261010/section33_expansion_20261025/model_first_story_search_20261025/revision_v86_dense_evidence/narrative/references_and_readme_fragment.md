# 期刊图参考与可复绘说明片段

检索与核对日期：2026-10-26。按 agent-reach/Exa 定位官方期刊图页；前两篇在官方图页实际查看浏览器截图，第三篇核对官方图3说明，由主代理查看实际图。这里只借鉴数据组织，不复制原图，不把其他系统的机制结论作为CaMKI证据。

## 三篇NC实际图参考

Nguyen et al. (2025), *A domain-swapped CaMKII conformation facilitates linker-mediated allosteric regulation*, Nature Communications 16, 8461。图7官方页：https://www.nature.com/articles/s41467-025-63249-w/figures/7 ，论文：https://www.nature.com/articles/s41467-025-63249-w 。图7上方是完整复合物小插图、两个时刻的局部接触结构及相邻距离分布，下方用序列字母、电荷位置轨道与功能曲线对应。v86借鉴“具体序列位置连接结构环境与定量读数”的组织方式：c在共享残基坐标下对齐完整曲线和序列；e把每个侧链值与同位置静态遮蔽条对齐；f补足真实结构环境。CaMKII与CaMKI不同，本参考不提供本例功能或因果证据；原图中的表格数字不移植到本图。

Zhu et al. (2022), *Neural relational inference to learn long-range allosteric interactions in proteins from molecular dynamics simulations*, Nature Communications 13, 1661。图2官方页：https://www.nature.com/articles/s41467-022-29331-3/figures/2 ，论文：https://www.nature.com/articles/s41467-022-29331-3 。图2把结构、完整残基关系矩阵、结构块汇总和网络逐层呈现，矩阵边缘用细条和结构域名而非数据区内通高大框。v86的d借鉴完整矩阵加轴外叶瓣注释，将全部25×144位置保留；c以柔和核心区和细条代替深色矩形边界。v86不估计学习边、传播路径或动态接触概率，因此不采用原文网络的方向箭头与边权含义。

Song et al. (2024), *Accurate Prediction of Protein Structural Flexibility by Deep Learning Integrating Intricate Atomic Structures and Cryo-EM Density Information*, Nature Communications 15, 5538。图3官方页：https://www.nature.com/articles/s41467-024-49858-x/figures/3 ，论文：https://www.nature.com/articles/s41467-024-49858-x 。图3c以一致结构展示配合下方完整逐残基RMSF曲线，使全局空间位置与局部变化相互核查。v86保留完整25个主链及21个有效侧链位置，结构a/b、接触d、环境e共同解释曲线，不只突出两个挑选位点。原文的预测–MD对比服务准确性，本图主面板仍按用户要求只给Full/Intra绝对预测值，结合态MD核查留在S3。

## 可复制到交付README的说明

v86沿用v85冻结的CaM家族留出第1折预测和原MD标签，没有新增模型运行，也没有更换划分。主图c保留完整CaMKI294–318的25个RMSF值与21个有效侧链Hχ值的Full/Intra对照。d在同一输入结构上独立计算全部3600个残基对的最小重原子距离，以≤4.0 Å、4.0–4.5 Å、>4.5 Å三类显示；留白表示超出图示接近阈值，而非缺失数据。e保留全部21个有效侧链，按氨基酸类型排列，并把静态侧链遮蔽比例逐行对齐。遮蔽计算只比较同一靶肽坐标保留/去除CaM时的可及面积，不是游离构象集合。

f使用官方1A06自抑制结构和1MXE靶肽复合物的真实坐标。两结构共同实际有坐标的调控片段为294–316；1A06缺317–318，不补建。结构对应不是分子运动轨迹，也不是完整激酶激活态对照。图中虚线仅指向待验证的构型维持解释；预测值、结合态MD、静态接触都没有直接测量返回激酶抑制接触或激活过程。全部21侧链视图保留Lys等非单调例外；正文仍明确N308单点误差未优于Intra，避免把同类型Asn成功关系扩大为逐点占优或“埋得深就必然低熵”。

结构灰/紫/青沿用前版；数值Full蓝圆实线、Intra暖赭方虚线保持全稿一致。文献只借布局，未采用其配色或新增色彩语义。颜色定义、原始值、重排顺序及数据出处均应随主代理最终composer/source_data保留。

结构事实补充（由结构代理核实）：1A06中的M316与L45/S99在4.5 Å内接近，M316并非深埋；W303对294–316之外背景的最近重原子距离5.19 Å，因此4.5 Å内无激酶接触。1MXE中两个锚点均被CaM明显遮蔽。f可强调M316的接触对象不同，不能将两锚点都画成原先被激酶深埋再被CaM对称夺取。接触距离不涉及SASA算法差异；结构代理的额外PyMOL SASA若用于正文或图注，应独立注明算法参数，不能与v85的960球面点结果混作一套绝对值。
