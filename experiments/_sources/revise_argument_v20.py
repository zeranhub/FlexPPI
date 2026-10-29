"""Revise the manuscript argument without changing scientific results."""
from pathlib import Path
from copy import deepcopy
from collections import Counter
import hashlib, json, re, zipfile
from lxml import etree
from docx import Document
from docx.oxml.ns import qn
from docx.shared import Inches, Pt

ROOT = Path(__file__).resolve().parent
QA = ROOT / 'argument_revision_v20_20261012'
QA.mkdir(exist_ok=True)
SOURCE = ROOT.parent / 'paper_v19_20261012_中文修订_补充实验完成_三线表.docx'
OUTPUT = ROOT.parent / 'paper_v20_20261012_中文修订_联合模型_三线表.docx'
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
assert sha(SOURCE) == 'b13c4f1e1de8458609fe6ab8d0eb17ceea81b0f0992d1a5d3d0728eb056bdea6'
doc = Document(SOURCE)
paras = list(doc.paragraphs)
tables = list(doc.tables)
changes = []

def set_text(p, text):
    assert not p._p.xpath('.//m:oMath'), 'Use a targeted text substitution for native equations'
    nodes = p._p.xpath('.//w:t')
    assert nodes
    nodes[0].text = text
    nodes[0].set('{http://www.w3.org/XML/1998/namespace}space', 'preserve')
    for node in nodes[1:]:
        node.text = ''

replacements = {
2: '蛋白质的局部动态既包括主链位置的波动，也包括侧链在不同取向之间的切换；两者共同影响分子识别与功能调节，却难以由单一静态构象直接判断。我们提出FlexPPI，以氨基酸类型与静态主链结构为输入，通过共享几何编码器联合预测逐残基均方根波动（RMSF）和联合侧链状态熵（Hχ）。模型利用固定的残基局部坐标系与不变点注意力整合三维邻域，在同一框架下处理单链和多链结构。3,063条记录的五折体系留出评价中，两项目标相对静态特征对照的平均绝对误差分别降低14.8%和11.3%。ThermoPCD外部验证进一步考察两类输出在另一轨迹来源中的表现。在KIX四种组成的MD轨迹对照中，模型定位到均具有较高波动的L12连接环；在泛素–UBA复合物中，模型区分了两个主链波动相近的谷氨酰胺侧链，并与轨迹中的取向覆盖差异一致。FlexPPI通过一次结构分析同时提供主链活动区域与侧链状态多样性两种互补读数，为局部动力学模拟和功能实验筛选关注位点。',
8: '侧链动态的学习还需要区分氨基酸本身的状态范围与结构环境带来的差异。不同氨基酸具有不同数量的侧链扭转角，仅按氨基酸类型给出平均值就可能获得较高的整体相关，而无法判断同一种侧链在哪个位置受到更强约束。Bios2cor通过二面角与旋转异构状态相关分析功能位点，RING-PyMOL利用结构集合的接触变化研究残基相互作用，DynaSig-ML也展示了动态表征与下游功能预测之间的联系[12,17,18]。这些研究说明侧链及局部环境包含值得分析的信息。从静态主链同时估计位置波动与多个侧链角的联合状态多样性，需要明确的监督定义和类型内评价。我们因此将RMSF与Hχ作为互补目标，并使用同类残基比较检验模型是否学到了氨基酸类型之外的结构信息。',
10: '我们从整体预测、外部迁移和结构实例三个层面评价FlexPPI。首先，通过五折体系留出比较完整模型、仅链内模型、静态特征回归及氨基酸类型均值，分析几何表示、完整结构环境和残基类型提供的预测信息；发布模型对照、单目标训练及全组成链家族评价见补充分析。其次，在ThermoPCD外部数据上排除与训练集合关联的体系，考察两类输出在另一轨迹来源中的表现。最后，将预测对应到具体结构和轨迹：KIX连接环用于检验主链活动区域的定位，泛素–UBA中的同类侧链用于检验局部取向多样性的区分。本文建立静态主链到两类互补动态指标的联合预测框架，并通过整体误差、类型内比较和结构实例说明其在局部动力学研究中的用途。',
30: '链内模型保留相同网络架构和参数量，分别编码各条链；完整模型共同编码输入结构中的全部有效链。两者分别从头训练，用于比较链内几何与完整空间环境的预测信息。FlexPPI采用联合训练作为默认方案，RMSF与Hχ共享几何编码器，并保留各自的读出参数，使两类读数基于同一结构环境生成。读出层经归一化、线性映射和SiLU激活，分别生成指标均值参数和对数方差；均值经softplus转为非负值，对数方差限制在[−8,6]。预测方差参与残差损失，覆盖率与事后分组校准见补充材料C。',
39: '五折体系留出同时考虑生物体系和序列关系[19,20]。共享PDB编号、生物体系、完整复合物序列组合或重复标签文件的记录归为同组。非ITrace单链记录进一步按序列家族合并，多链记录按完整链家族集合合并，序列聚类阈值为30%一致性与70%双向覆盖。ITrace按完整复合物组合分组，组成链家族可出现在不同折。每次以一折为测试集，其余四折训练，每条记录由对应留出模型产生折外预测。另在固定划分上改变随机初始化，评价训练随机性的影响。',
40: '补充评价按全部组成链的序列关联构建家族划分。长度至少20个残基的链使用关闭低复杂度屏蔽的MMseqs2搜索，以至少30%一致性和70%双向覆盖定义关联；短链按精确序列连接，另合并所有长度链的精确重复及同一来源生物体系。只要任意组成链关联，相应复合物记录即通过连通关系归为同组。3,063条记录形成1,662个连通家族，按记录数确定性平衡为633、608、608、607、607条的五折，同一连通家族及生物体系保持在同一折。划分依据为序列与体系关系。32条外部记录仅用于识别关联关系。完整与链内模型在新划分上各训练五个模型，使用种子7、第40轮权重计分。全部244条ITrace位于第0折测试集合，该折训练集合来自ATLAS和DynaRepo。',
41: '外部验证使用ThermoPCD中在标签提取与模型计分前确定的32条记录。该资源收集抗体–抗原复合物在生理与发热温度范围内的MD轨迹[27,28]。纳入记录具有可追溯的初始拓扑与原始轨迹，并能够区分温度和随机种子来源。训练排除使用上述全组成链家族关系：对长度至少20个残基的链运行MMseqs2，关闭低复杂度屏蔽、敏感度为7.5，以至少30%一致性和70%双向覆盖定义序列关联；短链检查精确重复，并合并生物体系关系。外部队列固定，相关训练记录按序列及体系关系排除。',
42: '外部方案排除633条相关训练记录，保留2,430条。完整与链内模型使用全链家族第0折在该集合上训练的单模型，采用种子7、第40轮权重；静态特征与类型均值对照在同一集合拟合。外部链与保留训练集合在上述搜索及精确序列条件下无命中。32条外部记录属于一个关联家族，按记录等权汇总。原训练范围的五模型平均结果作为部署配置见补充材料B。',
44: '氨基酸类型均值对照根据训练标签估计每种残基的平均RMSF和Hχ，用于评价残基类型提供的预测信息。静态特征回归模型使用47项特征，包括中心残基及链内前后四位的类型、链长与末端位置、链组成、连续残基的主链角正余弦、N–Cα–C夹角，以及链内和其他链的接触摘要。回归器采用直方图梯度提升，固定训练120轮，学习率0.07、最多31个叶节点、最小叶节点样本数50、L2正则为1、最大分箱数63，并关闭提前停止。五折及外部方案分别使用对应训练集合拟合对照。',
46: 'PLANET-MD[31]使用公开full及sequence-only权重，对每条链独立推断后合并对应记录。full输入序列和权重对应的encoded结构特征，sequence-only输入序列。采用发布模型默认300 K与float32精度，RMSF输出从nm乘以10转换为Å。其输出目标与本文Hχ定义不同，比较采用RMSF。RMSF-net[10]使用PDB-only发布权重，将完整蛋白结构转换成4 Å分辨率、1.5 Å网格的原子序数加权高斯模拟密度，评价其相对柔性分数的排序及记录内标准化剖面误差。PLANET、RMSF-net、ANM和GNM共同覆盖3,063条记录中的961,321个RMSF有效位置。发布模型采用原权重，静态特征和类型均值对照使用本文对应训练集合拟合。',
47: '评价先在每条记录的有效残基上计算平均绝对误差（MAE）、均方根误差（RMSE）、预测减参考的平均偏差，以及Spearman和Pearson相关，再对记录等权平均。MAE评价幅值估计，RMSE对较大误差更敏感，Spearman评价高低波动残基的相对排序。另在同一来源的生物体系内先平均相关记录，再对体系等权汇总。相关仅在至少4个有效位置且预测和参考均有变化时计算。同类型残基评价先计算每条记录中相应类型的误差，再对包含该类型有效位置的记录等权平均；类型均值在同类型内部为常数，采用MAE评价。',
48: '补充五模型分析对逐位置均值等权平均，总方差为模型方差均值加模型均值之间的总体方差。固定划分的重复训练报告各模型成绩的均值与样本标准差。逐残基曲线按静态模型MAE减完整模型MAE升序排列，分别选择差值最小的记录和居中记录，并对RMSF与Hχ分别排序，展示局部预测在不同表现下的结构对应。曲线按作者编号和实际链边界排列，在无效位置断开。',
51: 'KIX模拟使用OpenMM 8.4、CHARMM36力场及其配套TIP3P水模型。在完整三链结构中补齐缺失重原子后派生各状态，以pH 7.0的构建选项加入氢原子，KIX的五个组氨酸采用HID形式，链端带电且未封端。各体系置于立方体周期性水盒，初始蛋白与盒边至少相隔1.0 nm，加入0.1 mol/L NaCl并电荷中和。长程静电使用粒子网格Ewald方法（PME），非键作用截断为1.2 nm、开关距离为1.0 nm；约束所有涉及氢原子的共价键，并使用刚性水分子，积分步长为2 fs。经两阶段L-BFGS能量最小化后，在298 K下采用Langevin恒温器，阻尼系数为5 ps⁻¹，采用恒温恒容（NVT）模拟。先运行10 ns蛋白重原子位置约束平衡，约束势为U=k d²、k=4,184 kJ mol⁻¹ nm⁻²，再解除位置约束运行100 ns生产模拟，每1 ps保存蛋白坐标。每种状态分析一条生产轨迹，以KIX全部Cα对齐生产首帧后按式(1)计算RMSF；另核查逐帧共价键几何，并比较前后半段的指标。采样诊断和结合状态差值见补充材料F。',
56: '五折体系留出评价中，FlexPPI相对静态特征回归模型的RMSF和Hχ平均绝对误差分别降低14.8%和11.3%，残基排序相关也有所提高（表1，图2）。这说明对残基间距离、方向及主链局部几何的直接学习提供了常规结构摘要之外的预测信息。两类输出共同描述主链活动幅度和侧链状态多样性，后续类型内比较与结构实例进一步分析它们的互补用途。',
59: '完整模型相对链内模型的RMSF平均误差由0.636降至0.577 Å，Hχ由0.170降至0.166 nats（表1）。图2按实际输入链数展示两种分别训练方案的差值，区分单链模型的训练差异与多链结构中的完整环境信息。全组成链家族评价中，完整模型相对链内模型的两项MAE分别降低7.5%和3.8%（补充材料H）。',
60: '发布模型在RMSF预测中呈现不同的表现。BackFlip2的总体MAE为0.537 Å、ρ为0.848；PLANET full为0.568 Å和0.828，FlexPPI为0.577 Å和0.796。PLANET full的RMSE为1.186 Å，FlexPPI为0.792 Å，表明平均幅值误差与较大误差的贡献需要分别观察。PLANET sequence-only的MAE为0.714 Å、ρ为0.771，RMSF-net PDB-only的ρ为0.658。来源分层与相对柔性基线见补充材料I。',
61: '联合训练使主链波动与侧链状态多样性在同一几何表示上读出。保持相同输入、编码器、留出划分和40轮训练预算时，联合模型的RMSF MAE为0.577 Å，单目标模型为0.582 Å；Hχ分别为0.166和0.169 nats（补充材料G）。我们采用联合模型作为默认方案，在一次推理中提供两类动态指标，便于对同一结构位置联合分析骨架活动和侧链状态。',
63: '改以生物体系等权汇总后，完整模型相对静态特征和链内对照的主要误差方向保持一致，BackFlip2仍具有较好的总体RMSF成绩。各折与体系汇总结果随源表提供，ANM/GNM对照见补充材料D。',
65: '图2：五折折外预测及输入链数分层。Full为完整结构模型，Intra为分别训练的链内模型，Static为静态特征回归，AA为氨基酸类型均值。（a,b）两项目标逐记录MAE的经验累积分布，各包含3,063条记录；RMSF横轴在近零区域线性、远端为对数。（c,d）逐记录Spearman相关分布，箱体为四分位区间，线为中位数，须线覆盖完整范围，圆点及右侧数字为记录等权均值。（e,f）按实际输入链数分为1、2、3条及以上，展示同一记录的Intra减Full MAE；正值表示Full误差较低，每点为一条记录，黑色菱形为组内记录等权均值。差值轴在零附近线性、两端为对数；1链组对应两种分别训练模型的差异。',
66: '来源分层中，FlexPPI相对静态特征模型在三个来源的两项目标上均具有较低平均误差。ATLAS主要提供单链蛋白，DynaRepo包含单链和多链结构，ITrace集中于pMHC–TCR复合物。各来源的绝对误差结合蛋白组成、动态幅度及重复轨迹或时间窗口的标签汇总方式分析。',
67: 'BackFlip2在ATLAS及整体汇总中的RMSF误差较低，FlexPPI在DynaRepo和ITrace中的幅值误差较低；DynaRepo中两者的残基排序相关接近。幅值估计与高波动区域排序呈现不同的比较结果，提供了选择柔性分析读数的依据。',
69: '固定留出划分并改变随机初始化后，完整模型相对链内模型的整体误差方向保持一致。两项目标的模型间平均差距均大于完整模型的种子间标准差，说明总体比较在这些初始化中保持一致（表S1，图S4–S6）。',
72: 'RMSF较大的记录通常具有较大的绝对误差。在ATLAS和DynaRepo中，将MAE除以参考均值后，该趋势明显减弱；ITrace仍保留较强关系（图3）。以Å表示的MAE描述幅值偏差，残基排序相关描述高波动位点的定位，两者对应不同用途。Hχ的绝对误差趋势随来源变化；参考熵接近零时，比值误差增大，因此采用原单位下的误差作为主要评价，并以相对误差辅助观察尺度效应。',
74: '侧链状态预测的收益随残基类型变化。Pro在当前分箱下99.668%的有效标签为零，类型均值已接近误差下限；Val与链内模型的误差差异小于0.001 nats。RMSF方面，完整模型相对链内、静态特征及类型均值对照在全部20种类型上均有较低平均误差。与BackFlip2的类型内比较与整体结果一致。',
76: '图3：动态幅度与预测误差的关系及同类残基内的结构信息收益。（a,b）3,063条记录的参考标签均值与Full MAE，每点为一条记录。（c,d）MAE除以相应记录的参考标签均值并乘100%，RMSF相关坐标采用对数轴；颜色与形状区分来源。（e）Full在同一氨基酸类型内的记录等权MAE，RMSF与Hχ分别使用Å和nats，以独立坐标展示；Ala和Gly的Hχ不适用。（f）在相同有效记录上相对Intra和Static的MAE降幅，定义为100%×（对照MAE−Full MAE）/对照MAE；颜色区分对照，圆点实线为RMSF，三角虚线为Hχ，百分比坐标在±10%内线性、远端为对数。类型内分数先逐记录计算，再对具有相应类型有效标签的记录等权平均。',
77: '3.3 ThermoPCD外部验证',
78: 'ThermoPCD外部验证考察两类输出在另一轨迹来源中的表现。根据方法所述的全组成链序列与体系关系，移除633条相关记录，保留2,430条训练记录。完整和链内模型使用该集合上的家族第0折模型，静态特征与类型均值对照在同集合拟合；32条外部记录用于计分（表2，图4）。',
80: '完整模型的RMSF MAE为0.474 Å，链内和静态模型均为0.471 Å；Hχ MAE分别为0.208、0.229和0.211 nats。完整模型的两项Spearman相关为0.650和0.895，均高于链内和静态对照（表2，图4）。外部结果显示，联合输出能够保留对高波动残基及侧链状态多样性的排序信息。',
81: '侧链状态熵在外部数据中的表现体现了完整结构环境的作用：相对链内模型，完整模型的Hχ误差降低约9.4%；相对静态对照降低约1.5%。RMSF幅值误差在三种模型之间相近。BackFlip2在该队列的RMSF MAE为0.749 Å、ρ为0.593，其结果与主评价中的表现不同，说明跨来源分析需要同时查看误差和排序。',
83: '图4：32条ThermoPCD记录的外部验证。（a,b）RMSF与Hχ逐记录MAE；Full、Intra、Static和AA均使用保留的2,430条训练记录建立模型或对照。灰线连接同一记录的四种方法，小圆点为逐记录误差，大菱形及方法名下方数字为记录等权均值；RMSF采用对数轴。Full与Intra为第0折、种子7、第40轮单模型。',
88: '在未结合伙伴、结合cMyb、结合MLL和双伙伴四种组成中，模型与MD的完整KIX残基排序相关为0.729–0.896。最高20%柔性位点的预测与MD集合有61%–89%重合，双伙伴状态的18个优先筛选位点有14个命中。四条100 ns轨迹中的L12均具有较高柔性，双伙伴状态下模型将其八个位点中的六个列入最高20%集合（图5c、d）。这些结果将静态结构上的预测对应到一个具体活动连接环，为转录因子识别体系的局部轨迹分析提供了关注区域。四种组成的幅值及状态差值见补充材料F。',
91: '两个实例展示了联合读数在结构分析中的用途：RMSF定位KIX的活动连接环，Hχ进一步区分骨架波动较小的泛素位点中的侧链取向多样性。这些结果可用于优先选择局部动力学模拟、突变分析和功能实验中的关注位点，并形成关于空间约束与侧链状态分布的具体假设。',
95: 'FlexPPI以共同的静态主链输入刻画两种局部动态：RMSF描述骨架位置的波动，Hχ描述侧链角度组合的多样性。共享几何编码器使两个指标基于同一残基环境生成，便于在相同结构位置上联合查看。单目标训练对照中，联合模型的两项MAE略低；本文采用联合训练作为默认方案，其用途在于同时提供互补读数。静态特征对照已包含主链角和接触摘要，几何模型仍降低两项目标的整体误差，说明残基间距离、方向与局部主链朝向提供了附加信息。',
96: '同类型比较使侧链输出的结构含义更加清楚。氨基酸类型决定可定义的侧链角数量和状态范围，类型均值已能解释部分整体排序；类型内误差比较及泛素Q41/Q62实例进一步显示，模型能够使用局部环境区分同一种侧链在不同位置的状态多样性。Q41与Q62主链波动相近，却具有不同的取向覆盖范围，说明Hχ为RMSF补充了主链幅值之外的信息。',
97: '完整结构输入在体系留出和全组成链家族评价中均降低了相对链内模型的平均误差，在外部队列中也改善了侧链状态预测。分析多链结构时，将邻近伙伴纳入同一编码器，使局部几何同时包含链内堆积和分子间接触。模型输出可按链和残基映射回结构，用于寻找高波动区域及侧链状态多样的位点。',
99: '标签由轨迹采样和处理方式共同确定。共同时间窗口分析在83条记录上统一使用100 ns、200 ps间隔、逐链对齐和第一条轨迹，完整模型相对静态对照的两项MAE收益在三个来源中均保留（补充材料E）。当前输入由残基类型和静态主链构成；进一步纳入温度、溶剂及配体状态，可将联合预测扩展到更明确的条件比较。方差输出对应残差建模，补充材料C报告其对现有标签均值的覆盖率和区间宽度。',
100: 'KIX与泛素实例分别把主链活动区域和侧链状态多样性对应到具体结构位置，为局部动力学分析提供了可操作的位点选择。后续可围绕这些位置设计局部模拟或突变实验，并比较空间约束改变后的动态读数。在抗体识别等应用中，也可将活动位点预测与结合前后的局部环区构象比较结合[16]。',
102: 'FlexPPI从静态主链结构联合估计逐残基RMSF与侧链状态熵，在共享表示中结合残基类型、主链局部方向和完整空间环境。体系留出评价显示，几何表示降低了相对静态特征对照的两项平均误差；类型内比较与泛素实例说明侧链输出能够补充主链波动之外的局部状态信息。KIX实例进一步将预测对应到活动连接环。联合模型通过一次推理提供两种互补读数，为局部动力学模拟、突变分析和功能实验筛选关注位点。',
143: '输入结构的PDB背景通过RCSB PDB Data API获取，包括条目分类、分子实体描述、EC编号和UniProt对应[29]。3,063条记录对应2,751个PDB条目，其中2,748个条目可解析，覆盖3,060条记录，其余3条保留未知。逐记录来源注释与原始API响应随复现材料提供。',
145: '图S2区分单次观测、多个已保存观测相同及观测次数未保存三种零标准差来源。单次观测的样本方差记为零占位。Hχ为零表示参与汇总的各窗口或重复内部只观测到一种联合状态。按角数归一化的熵用于尺度展示，训练和计分使用原始nats；Ala和Gly不适用。',
149: '在固定划分上改变随机初始化，完整与链内模型的整体比较方向保持一致（表S1，图S4–S6）。逐残基曲线按方法所述的误差差值排序，展示差值最小及中间位置记录的局部预测。原训练范围的ThermoPCD五模型平均采用部署配置，结果见表S2及图S3。',
151: '表S2：原训练范围的ThermoPCD五模型结果。逐位置平均五个模型的输出后，在有效残基上计分，再对记录等权汇总。完整与链内使用五模型平均，BackFlip2采用发布权重。',
154: '原训练范围的部署配置采用五模型逐位置平均，其两项目标的平均误差均低于静态特征模型。完整与链内模型的RMSF误差接近，Hχ误差以完整模型较低（表S2）。该配置保留与外部集合关联的训练记录；正文采用移除633条关联记录后的单模型方案。',
177: '事后校准使用已冻结的测试预测和序列关联分量。以种子1729与分量标识的SHA-256哈希奇偶值确定两侧，每个分量整体保留在同一侧。每折、每项目标分别拟合正的模型方差缩放系数：按记录等权最小化高斯负对数似然，其解析解为各记录有效位置上“平方误差除以原生方差”的均值，再对记录平均。方差下限为10⁻⁶。每侧记录采用另一侧拟合的系数，均值预测固定；模型方差加标签变异另作敏感性方案。覆盖率与区间宽度先在记录内计算，再按记录等权汇总，均值 ± SD为五折均值和样本标准差。',
178: '冻结五折预测的95%高斯区间，以现有标签均值为参考。模型原生方差对应的RMSF和Hχ覆盖率分别为51.76% ± 2.71%和85.19% ± 1.31%；加入标签在已有重复轨迹或时间窗之间的变异后，分别为71.85% ± 2.26%和92.46% ± 0.53%（表S4）。',
179: '事后校准在每个测试折内按序列关联分量分为两侧，一侧拟合方差缩放系数并用于另一侧，交换后汇总全部记录，均值预测保持固定。RMSF和Hχ的95%覆盖率分别为96.30% ± 0.84%和94.84% ± 0.29%；对应平均区间宽度由0.876 ± 0.085 Å增至6.469 ± 1.525 Å、由0.623 ± 0.019 nats增至0.964 ± 0.025 nats（表S4）。覆盖率与区间宽度共同描述这一校准的效果。',
180: '来源分层中，校准后Hχ的95%覆盖率在ATLAS和DynaRepo接近标称值，在ITrace为90.65% ± 1.05%；ITrace的RMSF为99.14% ± 0.65%。ITrace属于同一关联分量，其缩放系数来自另一侧的其他分量。所有区间用于现有参考标签均值的高斯诊断，保留未截断的原区间。',
182: 'RMSF宽度单位为Å，Hχ为nats。标签变异为重复或时间窗之间的样本方差；每个评价记录采用另一侧关联分量拟合的缩放系数，均值预测保持固定。',
185: '图S14：方差输出与事后分组校准。（a,b）RMSF与Hχ的标称高斯覆盖率和实际覆盖率，虚线表示相等覆盖率。（c,d）95%区间在总体及各来源的覆盖率；（e,f）对应平均区间宽度。蓝色为原生模型方差，橙色为模型方差加标签在已有重复或时间窗之间的变异，紫色为按关联分量两侧互换拟合后的方差。每条记录使用另一侧拟合的系数，均值预测固定。点与误差线为五折均值 ± 样本标准差，半透明小点为各折记录宏平均。总体为3,063条记录，RMSF和Hχ有效位置分别为961,321和824,988。覆盖率与宽度均以参考标签均值为诊断对象。',
190: 'ANM、GNM与完整模型的整体Spearman相关分别为0.715、0.663和0.796，标准化MAE为0.450、0.487和0.371，最高10%位点命中率为56.7%、57.7%和65.1%。在ATLAS中，ANM与完整模型的排序相关分别为0.785和0.788；DynaRepo为0.635和0.774，ITrace为0.453和0.938。来源分层帮助比较弹性网络与几何学习对当前MD柔性剖面的描述能力。',
193: '共同窗口分析按来源内样本标识排序，以种子20261010从ATLAS、DynaRepo和ITrace各抽取30条，共90条记录。根据原轨迹时间戳，7条DynaRepo记录只有约1 ns，按至少10 ns和时间戳单调的预设条件排除，保留30、23、30条，共83条。其203条可用轨迹均覆盖0至100 ns；采用ITrace的200 ps保存间隔作为共同观测网格，共501帧，同一时间点取首次保存帧。',
194: '主方案对每条记录使用第一条轨迹，统一窗口、采样间隔和估计方法。各链Cα逐帧Kabsch拟合到共同窗口首帧，按式(1)计算RMSF；Hχ沿用完整χ原子定义、三个角度分区、联合状态和自然对数，Ala和Gly无侧链监督，Pro采用原标签定义。完整、链内、静态、类型均值与BackFlip2预测保持固定，在原标签和共同窗口标签的相同有效位置计分，共38,599个RMSF和33,493个Hχ位置。另计算ATLAS和ITrace三条可用轨迹的共同窗口均值，DynaRepo保持一条轨迹。',
195: '表S6：83条记录中完整模型在原标签与共同窗口第一条轨迹标签上的表现。RMSF误差单位为Å，Hχ为nats；ρ为记录内Spearman相关的记录等权均值。',
196: '统一使用第一条轨迹后，三个来源的RMSF和Hχ MAE均提高，记录内排序相关均降低（表S6）；完整模型相对静态对照的两项平均误差仍较低。抽样ATLAS的RMSF在原标签上以链内误差较低，在共同窗口方案下以完整模型较低，说明时间窗口和重复数也会影响较小的模型差距。',
197: '使用三条可用轨迹均值后，ATLAS的RMSF/Hχ MAE为0.723 Å/0.172 nats，接近原标签结果；ITrace为0.158 Å/0.110 nats，DynaRepo保持一条轨迹。该分析统一估计器、窗口和采样间隔，模型预测保持冻结；逐记录标签变化、各对照配对分数和有效位置随复现材料提供。',
198: 'F KIX轨迹与状态比较',
199: '四种KIX组成各分析一条MD轨迹，每种状态含10 ns约束平衡和100 ns无约束生产。生产轨迹各含100个连续1 ns分段，每段1,000帧、间隔1 ps，共400,000个蛋白构象。全帧重算的RMSF和Hχ与原分析逐位一致，共价键几何保持连续。各状态按KIX的87个Cα对齐生产首帧，RMSF围绕各自对齐后的轨迹均值计算。模型采用正文图5相同的固定第1折留出权重及有效位置。',
200: '表S7：四种KIX组成的RMSF评价。每种状态一条100 ns轨迹；MAE在87个残基上计算，ρ为Spearman相关，命中数为最高20%位点的集合重合数。L12为614–621八个位点的算术平均；MAE与均值单位为Å。',
201: '采样分析另比较前后两个50 ns窗口。按游离、cMyb、MLL、双伙伴顺序，KIX Cα相对生产首帧的时间平均RMSD为2.843、3.902、3.232、4.070 Å，平均温度为298.88–298.90 K。L12在前后50 ns的RMSF均值依次为1.972/2.419、2.897/1.720、2.494/2.231、3.044/2.559 Å。初始界面以制备结构中距离不超过4.5 Å的KIX–伙伴重原子对定义；至少一个原始原子对仍满足阈值时，相应残基接触计为保持。cMyb、MLL单伙伴及双伙伴中的cMyb、MLL界面的平均初始残基接触保持率为46.7%、41.6%、54.7%、48.0%，该统计跟踪初始接触。',
202: '状态差值定义为表中前一状态的L12平均RMSF减后一状态，分别对MD与模型计算（表S8）。双伙伴减游离的MD差值为+1.017 Å，模型为−1.589 Å；MLL减cMyb为+0.014 Å和−1.263 Å，MD前后50 ns对应−0.402和+0.511 Å。四种组成的残基排序与状态间幅值差值呈现不同结果，正文案例采用残基排序和高波动集合定位活动区域。',
206: '在五折体系留出划分上分别训练RMSF与Hχ单目标模型，各五个，共十次训练。采用种子7、第40轮权重，与联合模型保持相同几何输入、编码器、抽样方法、优化器和训练预算。两种训练方案对应相同测试记录及有效位置，RMSF和Hχ分别覆盖961,321和824,988个位置。指标先逐记录计算，再对3,063条记录等权平均。',
207: '表S9：单目标与联合目标的记录等权结果。RMSF误差单位为Å，Hχ为nats，ρ为Spearman相关；同一目标采用相同有效位置。',
208: '联合模型的RMSF与Hχ MAE相对单目标分别降低约0.8%和1.8%。两种方案的RMSF相关接近，Hχ的MAE、RMSE和排序均以联合模型较好（表S9）。联合模型作为默认方案，以共享几何编码器生成两种互补的动态读数。',
210: '1,662个连通家族分别保留在对应折内，完整与链内模型共训练十次，每次40轮。每条记录由对应留出折模型预测，使用与体系留出相同的有效位置，并对3,063条记录等权平均。表S10列出整体及来源分层结果；ITrace的244条记录均在第0折测试。',
212: '全组成链家族评价中，完整模型相对链内模型的RMSF和Hχ MAE分别降低7.5%和3.8%，两项排序相关也较高（表S10）。来源分层结果呈现同一比较方向，进一步显示完整结构环境提供了链内几何之外的预测信息。',
214: 'PLANET full与sequence-only分别覆盖4,973条链，RMSF-net、ANM和GNM覆盖3,063条结构。七种方法在961,321个共同RMSF位置上计分，先汇总同一记录的有效链，再对记录等权平均。FlexPPI采用五折体系留出预测，PLANET与RMSF-net采用发布权重；ANM和GNM由相同静态结构计算。',
215: '表S11：共同RMSF位置上的发布模型与相对柔性基线。MAE单位为Å，ρ为Spearman相关；标准化MAE先将每条记录的参考与预测分别减去自身均值、除以自身标准差，再计算平均绝对差，为无量纲。相对分数的Å误差不适用。',
216: 'PLANET full在ATLAS、DynaRepo和ITrace的MAE/ρ分别为0.398 Å/0.876、0.704 Å/0.752和1.412 Å/0.732；sequence-only为0.499 Å/0.826、0.782 Å/0.668和2.171 Å/0.719。ITrace的预测减参考平均偏差，full为+1.381 Å、sequence-only为+2.150 Å。ATLAS占全部记录约63%，因此同时提供来源分层和逐记录结果。PLANET使用发布模型默认温度和独立链输入；RMSF-net使用PDB模拟密度和相对分数。',
}

for index, text in replacements.items():
    p = paras[index]
    changes.append({'source_paragraph':index, 'before':p.text, 'after':text})
    set_text(p, text)

# Remove repeated rhetoric and superseded external-protocol illustrations.
removed_paras = [57, 62, 68, 84, 85, 98, 139, 140, 167, 168, 169, 170, 204, 217, 218, 219]
for index in removed_paras:
    p = paras[index]
    assert not p._p.xpath('.//w:commentRangeStart | .//w:commentRangeEnd | .//w:commentReference | .//m:oMath')
    changes.append({'source_paragraph':index, 'before':p.text, 'after':None,
                    'drawing_removed':bool(p._p.xpath('.//w:drawing'))})
    p._p.getparent().remove(p._p)
tables[13]._tbl.getparent().remove(tables[13]._tbl)

# Keep the baseline table focused on its single evaluation protocol. Family
# results remain, unchanged, in the preceding dedicated table S10.
baseline = tables[12]
for row in baseline.rows:
    original_cells = list(row._tr.findall(qn('w:tc')))
    for index in [4, 3]:
        row._tr.remove(original_cells[index])
grid = baseline._tbl.tblGrid
for index in [4, 3]:
    grid.remove(grid[index])
widths = [1.80, 1.35, 1.20, 1.85]
for gridcol, width in zip(grid, widths):
    gridcol.set(qn('w:w'), str(round(width * 1440)))
for row in baseline.rows:
    for cell, width in zip(row.cells, widths):
        cell.width = Inches(width)
for cell, text in zip(baseline.rows[0].cells, ['方法', 'MAE', 'ρ', '标准化MAE']):
    set_text(cell.paragraphs[0], text)
paras[216].paragraph_format.space_before = Pt(8)

# Renumber the surviving supplementary figures; supplementary tables keep
# their existing sequence S1-S11 after removing the final historical table.
figure_map = {2:1, 3:2, 4:3, 5:4, 6:5, 7:6, 8:7, 9:8, 12:9, 13:10, 14:11}
for p in doc.paragraphs:
    if p._p.xpath('.//m:oMath'):
        continue
    original = p.text
    new = re.sub(r'(?<!表)S(\d+)', lambda m: 'S'+str(figure_map[int(m.group(1))])
                 if int(m.group(1)) in figure_map else m.group(0), original)
    if new != original:
        set_text(p, new)

doc.save(OUTPUT)
# Restore the original comments bytes rather than accepting serialization-only
# changes. All other archived parts remain available, including historical media.
with zipfile.ZipFile(SOURCE) as old, zipfile.ZipFile(OUTPUT) as new:
    old_comments = old.read('word/comments.xml')
    assert etree.tostring(etree.fromstring(old_comments), method='c14n') == etree.tostring(etree.fromstring(new.read('word/comments.xml')), method='c14n')
    saved = [(info, old_comments if info.filename == 'word/comments.xml' else new.read(info.filename)) for info in new.infolist()]
temp = OUTPUT.with_suffix('.preserved.tmp')
with zipfile.ZipFile(temp, 'w') as package:
    for info, data in saved:
        package.writestr(info, data)
temp.replace(OUTPUT)

ns = {'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main',
      'm':'http://schemas.openxmlformats.org/officeDocument/2006/math'}
with zipfile.ZipFile(SOURCE) as old, zipfile.ZipFile(OUTPUT) as new:
    assert old.read('word/comments.xml') == new.read('word/comments.xml')
    assert all(old.read(name) == new.read(name) for name in old.namelist() if name.startswith('word/media/') and not name.endswith('/'))
    a, b = [etree.fromstring(z.read('word/document.xml')) for z in (old, new)]
    anchors = lambda x: Counter((e.tag,e.get(qn('w:id'))) for e in x.iter() if e.tag in {qn('w:commentRangeStart'),qn('w:commentRangeEnd'),qn('w:commentReference')})
    equations = lambda x: [etree.tostring(e,method='c14n') for e in x.xpath('.//m:oMath',namespaces=ns)]
    assert anchors(a) == anchors(b)
    assert equations(a) == equations(b)

current = Document(OUTPUT)
assert len(current.tables) == 13
for index in range(12):
    assert [[c.text for c in r.cells] for r in current.tables[index].rows] == [[c.text for c in r.cells] for r in Document(SOURCE).tables[index].rows], index
assert [[c.text for c in r.cells] for r in current.tables[12].rows][1:] == [[r.cells[k].text for k in [0,1,2,5]] for r in Document(SOURCE).tables[12].rows][1:]

blocks=[]; p_index=0; t_index=0
for node in current._element.body:
    if node.tag == qn('w:p'):
        text=''.join(node.xpath('.//w:t/text() | .//m:t/text()'))
        blocks.append({'type':'paragraph','index':p_index,'text':text})
        p_index+=1
    elif node.tag == qn('w:tbl'):
        blocks.append({'type':'table','index':t_index,'rows':[[c.text for c in r.cells] for r in current.tables[t_index].rows]})
        t_index+=1
(QA/'v20_content.json').write_text(json.dumps(blocks,ensure_ascii=False,indent=2),encoding='utf8')
lines=[]
for block in blocks:
    if block['type']=='paragraph' and block['text']:
        lines.append('P'+str(block['index'])+': '+block['text'])
    elif block['type']=='table':
        lines.extend('T'+str(block['index'])+'R'+str(i)+': '+' | '.join(row) for i,row in enumerate(block['rows']))
(QA/'v20_content.txt').write_text('\n'.join(lines)+'\n',encoding='utf8')
(QA/'editorial_changes.json').write_text(json.dumps({'source':str(SOURCE),'output':str(OUTPUT),
    'source_sha256':sha(SOURCE),'output_sha256':sha(OUTPUT),'changes':changes,
    'removed_old_figures':['S1','S10','S11'],'removed_old_table':'S12',
    'family_table_preserved':True,'all_retained_numeric_table_cells_unchanged':True,
    'comments_and_anchors_preserved':True,'native_equations_preserved':True,
    'tables':13,'supplementary_figures':11,'figure_map':figure_map},ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps({'output':str(OUTPUT),'sha256':sha(OUTPUT),'changed_paragraphs':len(replacements),
                 'removed_paragraphs':len(removed_paras),'tables':13},ensure_ascii=False))
