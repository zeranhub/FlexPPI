# v21 可用文献与论证素材

核验日期：2026-10-12。建议核心新增四条[32–35]；若引言需要RMSF与局部弹性约束的物理背景，可再用[36]。以下句子是方法定义与研究动机，不包含尚未完成的新实验结果。

## 可直接加入文末的参考文献

[32] Gregory R. Bowman and Phillip L. Geissler. Extensive Conformational Heterogeneity within Protein Cores. *The Journal of Physical Chemistry B*, 118(24):6417–6423, 2014. doi:10.1021/jp4105823. [出版社](https://pubs.acs.org/doi/10.1021/jp4105823)，[原始全文](https://pmc.ncbi.nlm.nih.gov/articles/PMC4065209/)。

[33] Alexander D. Scouras and Valerie Daggett. The Dynameomics rotamer library: amino acid side chain conformations and dynamics from comprehensive molecular dynamics simulations in water. *Protein Science*, 20(2):341–352, 2011. doi:10.1002/pro.565. [作者原文PDF](https://depts.washington.edu/daglab/downloads/Rotamer_ProtSci_2011.pdf)，[PMC原始全文](https://pmc.ncbi.nlm.nih.gov/articles/PMC3048419/)。

[34] Rajitha Rajeshwar T., Jeremy C. Smith, and Marimuthu Krishnan. Hidden Regularity and Universal Classification of Fast Side Chain Motions in Proteins. *Journal of the American Chemical Society*, 136(24):8590–8605, 2014. doi:10.1021/ja5024783. [出版社原始记录](https://pubs.acs.org/doi/10.1021/ja5024783)，[PubMed作者与摘要](https://pubmed.ncbi.nlm.nih.gov/24844417/)。第一作者按出版商显示的完整名字保留；BibTeX中使用文字保护，避免把末尾T.误解析为另一个作者。该文在线时间为2014-06-06，正式卷期时间为2014-06-18。

[35] Robert B. Best and Gerhard Hummer. Coordinate-dependent diffusion in protein folding. *Proceedings of the National Academy of Sciences of the United States of America*, 107(3):1088–1093, 2010. doi:10.1073/pnas.0910390107. [原始全文](https://pmc.ncbi.nlm.nih.gov/articles/PMC2824289/)，[PubMed记录](https://pubmed.ncbi.nlm.nih.gov/20080558/)。该文在线发表于2009-12-28，正式卷期年份是2010，文献应写2010。

[36，可选] Ivet Bahar, Ali Rana Atilgan, and Burak Erman. Direct evaluation of thermal fluctuations in proteins using a single-parameter harmonic potential. *Folding & Design*, 2(3):173–181, 1997. doi:10.1016/S1359-0278(97)00024-2. [出版社原始记录](https://www.sciencedirect.com/science/article/pii/S1359027897000242)，[作者机构保存的原文](https://ccbb.pitt.edu/archive/PDFs/DynPro97.pdf)。

## 推荐引言段

蛋白质的动态特征来自其采样的构象分布。已有长时间模拟研究表明，即使在紧密堆积的蛋白核心，侧链也可访问多个旋转状态；状态占比与状态转换时间分别提供构象多样性和运动时间尺度的信息[32,33]。增强采样研究进一步将侧链占比、构象熵和转变过程与局部自由能面联系起来[34]。据此，RMSF与侧链离散状态熵可用于概括构象集合的不同方面：前者描述主链位置分布的宽度，后者描述侧链状态占比的多样性。我们进一步结合原始轨迹中的连续坐标分布，考察这两类读数对应的动态形态，并分析折外预测在不同分布形态中的表现。

若需要在引言中补充RMSF的物理背景，可在上段首句后加：“在结构附近的近似谐和运动中，主链涨落与局部弹性约束相连，弹性网络可由残基间的空间连接预测这种涨落[36]。”这里的谐和条件是具体模型前提，无需另加一段免责叙述。

## 方法：经验投影自由能

从每条轨迹的指定分析窗口提取连续坐标q，并以固定箱宽估计其经验概率密度p̂(q)。各箱以出现频率除以箱宽计算密度。相对经验投影自由能定义为F*(q)=−ln[p̂(q)/max p̂(q)]，最低点设为零，以无量纲kBT标度展示。该曲线显示当前采样窗口在所选坐标上的高占比区域和分布谷；未访问的箱记为缺失，不绘制为无限高峰。模型预测与这些轨迹分布特征逐记录配对评价，时间组织另通过时间序列或滞后统计描述。投影坐标的选择会改变曲线形状，运动时间尺度还涉及沿该坐标的扩散行为[35]。

实施时须在具体图注补上q的定义、箱宽、窗口长度、保存间隔、是否平滑以及分析了哪条重复。对相同等宽箱，用箱概率p̂_j代替密度得到同一条相对曲线；不同箱宽时必须用密度。p̂是有限窗口的经验占比，保留“经验投影”这一名称即可。温度未知时也能计算上述无量纲读数，无需写入假设的温度或转换为kcal/mol。若程序使用伪计数或截断纵轴，图注需写实际规则，不能把这些数值处理形成的峰高称为能垒。

## 方法：χ₁分区与圆周期

χ₁按360°圆周坐标处理。对具有常规三态旋转构象的非Pro残基，参考分区以−60°、+60°和180°附近为中心，边界为−120°、0°和+120°：g−为[−120°,0°)，g+为[0°,120°)，trans为[120°,180°)与[−180°,−120°)的并集。它对应Scouras与Daggett图2中四面体碳扭转角的常规分区[33]。原冻结标签仍采用[−180°,−60°)、[−60°,60°)、[60°,180°]区间。新增分析分别计算两种χ₁分区的边缘熵，比较其差异与排序，从而检验角度分区对读数的影响；模型训练和主计分目标保持原完整联合Hχ定义。

Pro的χ₁受环皱褶约束，原文采用围绕0°/±30°的特殊处理；Asn/Asp/His/Phe/Tyr χ₂与Gln/Glu χ₃等末端角也有不同分区。因此，非Pro χ₁敏感性结果不能写成所有χ角已经完成物理rotamer校正。χ₁边缘熵与多个χ角组合的联合Hχ应分别命名和标注。

## 可用于新图图注的共同说明

“连续坐标分布及经验投影自由能由MD轨迹计算，模型数据为对应记录的冻结折外预测。经验投影自由能为−ln(p̂/p̂max)，最低值设为零；χ₁分布采用圆周坐标。χ₁边缘熵用于分区敏感性分析，完整联合Hχ沿用原训练标签定义。”实际面板编号与样本数待新实验和图稿确定后填写，不预写结果方向。

## 论证使用规则

Bowman与Geissler支持“多个χ角联合状态的Shannon熵”这一构造，但本文原分箱边界不同，引用时不要写“与其rotamer熵定义完全一致”。Scouras与Daggett提供分箱边界、非rotamer末端角与状态占比/时间尺度分开分析的直接依据。Rajeshwar等使用ABF增强采样，分别计算熵、概率流与活化能，适合支持研究动机；本轮已核实出版社与PubMed原始摘要，未获得其付费全文，故不从该文推断具体分箱实现。Best与Hummer支持投影坐标与扩散的具体关系。RMSF/Hχ是摘要读数，经验F*来自完整轨迹分布；新图的主要问题应写为“预测读数对应何种分布形态”，避免把它直接升级成“模型预测能垒”。

研判报告中的3,068条旧OOF文件仅是当时抽查信息。v21应继续使用当前v20已冻结的3,063条修正队列及对应模型、有效掩码、窗口和标签，不因本次文献材料而切换数据版本。
