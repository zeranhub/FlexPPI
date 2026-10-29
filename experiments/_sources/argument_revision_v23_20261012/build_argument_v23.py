"""Prune duplicated analyses and organize evidence around the manuscript claims."""
from pathlib import Path
from copy import deepcopy
from collections import Counter
import hashlib,json,re,zipfile
from lxml import etree
from docx import Document
from docx.shared import Inches
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml.ns import qn
ROOT=Path(__file__).resolve().parent
BASE=ROOT.parent.parent
SOURCE=BASE/'paper_v22_20261012_中文修订_自由能分析融入正文_三线表.docx'
OUT=BASE/'paper_v23_20261012_中文修订_论证精简与补充整合_三线表.docx'
QA=ROOT/'manuscript_QA';QA.mkdir(exist_ok=True)
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
assert sha(SOURCE)=='315ca2108ab9dbdfae83774f4d9fe1556087d66f8ea09697afa9c9eb67b67bdd'
d=Document(SOURCE)
def find(start):
    pp=[p for p in d.paragraphs if p.text.startswith(start)];assert len(pp)==1,(start,len(pp));return pp[0]
def replace(p,text):
    nodes=p._p.xpath('.//w:t');assert nodes and not p._p.xpath('.//m:oMath')
    nodes[0].text=text;nodes[0].set('{http://www.w3.org/XML/1998/namespace}space','preserve')
    for n in nodes[1:]:n.text=''
def change(start,text):replace(find(start),text)
def delete(start):
    p=find(start);assert not p._p.xpath('.//w:commentRangeStart | .//w:commentRangeEnd | .//w:commentReference');p._p.getparent().remove(p._p)
def table_block(label):
    p=find(label+'：');t=p._p.getnext();assert t.tag==qn('w:tbl');return p,t
def append_table(block,new_label):
    p,t=block;replace(p,re.sub(r'^表S?\d+：',new_label+'：',p.text));d._body._body.insert(-1,p._p);d._body._body.insert(-1,t)
def para(text,style='Normal'):return d.add_paragraph(text,style)
def section(text):
    p=para(text,'Heading 1');p.paragraph_format.page_break_before=True;return p

# Capture the existing editable evidence tables before rebuilding the supplement.
keep={label:table_block(label) for label in ['表S1','表S2','表S9','表S11','表S10','表S6','表3','表S7','表S12']}
old_table_cells={label:[[c.text for c in r.cells] for r in Document(SOURCE).tables[i].rows] for label,i in [('表S1',3),('表S2',4),('表S9',11),('表S11',13),('表S10',12),('表S6',8),('表3',2),('表S7',9),('表S12',14)]}
label_fig=find('图S1：');label_pic=label_fig._p.getprevious();assert label_pic.xpath('.//w:drawing')
sampling_table_caption,sampling_table=keep['表3']
sampling_table_caption._p.getparent().remove(sampling_table_caption._p);sampling_table.getparent().remove(sampling_table)

# Remove repeated navigation and keep only the actual claims and comparisons.
change('五折体系留出评价中，FlexPPI相对','五折体系留出评价中，FlexPPI相对静态特征回归的RMSF和Hχ平均绝对误差分别降低14.8%和11.3%，残基排序相关也提高（表1，图2）。对残基间距离、方向及主链局部几何的学习提供了常规结构摘要之外的预测信息。')
change('完整模型相对链内模型的RMSF平均误差','完整模型相对链内模型的RMSF平均误差由0.636降至0.577 Å，Hχ由0.170降至0.166 nats（表1）。图2按实际输入链数展示两种训练方案的差值；全组成链家族评价中，两項MAE相对链内模型分别降低7.5%和3.8%（补充材料D）。')
change('发布模型在RMSF预测中呈现','发布RMSF模型的总体成绩以BackFlip2和PLANET full较好：MAE/ρ分别为0.537 Å/0.848和0.568 Å/0.828，FlexPPI为0.577 Å/0.796。PLANET full与FlexPPI的RMSE分别为1.186和0.792 Å。共同位置上的发布模型、弹性网络与来源分层比较见补充材料C。')
change('联合训练使主链波动与侧链状态','FlexPPI采用共享编码器的联合模型作为默认方案，一次推理同时提供两类动态读数；同输入和训练预算下的单目标及随机初始化对照见补充材料B。')
delete('改以生物体系等权汇总后')
change('来源分层中，FlexPPI相对','在三个来源中，FlexPPI的两项平均误差均低于静态特征模型。BackFlip2在ATLAS中的RMSF误差较低，FlexPPI在DynaRepo和ITrace中的幅值误差较低，DynaRepo中的残基排序相关接近。')
delete('BackFlip2在ATLAS及整体汇总中的')
delete('固定留出划分并改变随机初始化后')
change('3.2 动态幅度与氨基酸类型对评价','3.2 同类型残基的局部动态差异')
change('动态指标的数值尺度和氨基酸种类','氨基酸类型规定侧链角数与状态范围，同类型比较进一步检验模型能否利用局部结构环境区分位置差异（图3）。')
change('RMSF较大的记录通常具有','RMSF较大的记录通常有较大的绝对误差；按参考均值归一化后，这一趋势在ATLAS和DynaRepo中减弱，在ITrace中仍较明显。原单位MAE描述幅值偏差，排序相关描述活动位点的相对位置。Hχ接近零时比值误差较大，本文以原单位误差为主要评价。')
change('侧链状态预测的收益随残基类型','RMSF方面，完整模型在20种氨基酸类型中均具有低于链内、静态特征及类型均值对照的平均误差。Hχ的类型内收益随残基种类变化；Pro在原分区下99.668%的有效标签为零，Val相对链内模型的差距小于0.001 nats。')
change('全队列的配对分析显示','全队列中，RMSF与Hχ呈相关但不同步的局部变化。原值的参考与预测Spearman相关分别为0.233和0.260；减去记录内同类型均值后分别为0.391和0.517（图4a）。类型内的四分位组合进一步定位两类读数的差异。')
change('以同记录同类型的四分位区间','同记录同类型中，24,842个位点属于低RMSF／高Hχ，22,512个位点属于高RMSF／低Hχ，合计47,354个，占共同有效位置的5.740%。这些参考组合显示，单一读数不能覆盖全部局部动态形态。模型对低／低和高／高组合的记录平均召回率为54.6%和61.4%，对低／高和高／低为15.8%和16.5%（图4c，表S9）。两项目标具有不同信息，而当前模型对不同组合的恢复程度存在差异。')
change('ThermoPCD外部验证考察两类','在按2.5节排除关联训练记录的32条ThermoPCD外部记录中，比较两类读数的幅值与排序（表2，图5）。')
change('完整模型的RMSF MAE为0.474','完整模型的RMSF MAE为0.474 Å，链内和静态模型均为0.471 Å；Hχ MAE依次为0.208、0.229和0.211 nats。完整模型的两项Spearman相关为0.650和0.895，均高于链内及静态对照。其主要幅值收益来自Hχ，相对链内与静态模型的误差分别降低9.4%和1.5%；BackFlip2在该集合的RMSF MAE/ρ为0.749 Å/0.593。')
delete('侧链状态熵在外部数据中的表现')
change('为检验两类动态读数能否对应','将留出预测对应到KIX与泛素–UBA的具体结构位置，可以检验活动区域定位及同类型侧链差异。KIX在不同表面结合cMyb和MLL[23]，图6a、b显示三元复合物及连接相邻螺旋的L12环区。')
change('在未结合伙伴、结合cMyb','四种KIX组成中，模型与MD的残基排序相关为0.729–0.896，最高20%位置的集合重合为61%–89%。双伙伴状态的18个预测优先位点有14个命中，L12八个位点中六个进入最高20%集合；四条100 ns轨迹的L12均具有较高波动（图6c、d）。状态内的幅值与状态间差值见补充材料F。')
change('第二个体系是泛素与EDD','泛素–EDD UBA复合物提供同类型侧链的结构对照[24]。泛素链中，15个预测最高主链波动位点有14个命中；参考与预测同时减去训练折的氨基酸类型均值后，Hχ相关仍为0.868（图6e、g），显示类型平均值之外的位置差异。')
change('Q41与Q62均为谷氨酰胺','Q41与Q62均为谷氨酰胺，主链RMSF分别为0.469和0.615 Å，而侧链联合Hχ为1.890和2.631 nats，模型给出2.000和2.684 nats。Q41较埋藏，Q62朝向溶剂；MD代表帧及90%观测帧覆盖所需的七类和十二类联合状态显示后者更丰富的取向覆盖（图6f、h）。下节进一步通过逐角分布解释这一差异。')
delete('两个实例展示了联合读数在结构分析中的用途')
change('RMSF与Hχ概括了构象集合中的不同变化','RMSF概括对齐主链的位置散布，Hχ概括多角联合状态的占比。轨迹的相对经验投影自由能F*将常访问取向显示为低值区域，进一步展开两类预测摘要所对应的构象形态。')
change('泛素Q41与Q62的参考RMSF','泛素Q41与Q62的χ1均主要采样−60°附近，五窗口平均的原分区χ1边缘熵为0.837和0.833 nats。Q41的χ2集中在跨越±180°接缝的同一邻域，Q62的χ2主要采样−60°附近，χ3覆盖更分散；χ2边缘熵分别为0.680和0.835 nats，χ3为0.574和1.032 nats（图7）。Q62较高的联合熵因此包含后续角度及其共同组合的信息。')
change('把三个角的区间共同平移60°后','三个角的区间共同平移60°后，Q41/Q62的五窗口平均联合熵为0.528/1.399 nats，而常规χ1边缘熵为0.492/0.403 nats，两者的高低方向不同。联合类别与单角投影分别保留完整侧链的取向组合及某一旋转轴的占据；图7展开了这些不同信息。')
change('83条记录的统一窗口分析显示','在83条记录的相同采样帧上，原分区与60°平移的完整联合Hχ均值为1.484和0.962 nats（图8a）。原边界接近常见χ1峰，同一峰内的局部摆动也会分到相邻类别；连续角度与F*保持相同，这一比较说明离散熵如何读取指定区间下的占据分布。全部平移及时间采样统计见补充材料E。')
change('主链与侧链的联系还可以','以每帧联合χ状态对对齐Cα位置分组，RMSF²可以分解为状态内部均方波动与状态平均位置之间的均方波动。后者在原分区和60°平移中平均占总量的6.2%和7.8%（图8b），大部分骨架均方散布仍在同一侧链类别内部。因此，主链活动幅度与侧链类别占比能够呈现不同的局部组合。')
delete('时间组织进一步补充占据分布的形态')
fig8=find('图8：');pic8=fig8._p.getprevious();assert pic8.xpath('.//w:drawing')
for child in list(pic8):
    if child.tag!=qn('w:pPr'):pic8.remove(child)
from docx.text.paragraph import Paragraph
pp=Paragraph(pic8,d._body);pp.paragraph_format.page_break_before=False;pp.paragraph_format.keep_with_next=True;pp.alignment=WD_ALIGN_PARAGRAPH.CENTER
pp.add_run().add_picture(str(ROOT/'figures/fig8_joint_partition_CA.png'),width=Inches(6.2))
replace(fig8,(ROOT/'figures/fig8_joint_partition_CA.caption_zh.md').read_text(encoding='utf8').strip())
change('3.6.3 不同蛋白中的局部自由能形态','3.6.3 不同蛋白的局部占据与经验投影自由能')
change('为将这一解释对应到不同结构来源','按预设中位误差规则选取的2py5_A、A0022和ITrace 6p64，在完整残基剖面中呈现同类型位点的不同主链幅度与联合Hχ组合；逐角时间序列与F*进一步解释这些位置访问的取向（图9）。')
# Retain the actual angle-centre, single/two-region interpretation paragraph.
find('3.6.3 不同蛋白的局部占据').paragraph_format.page_break_before=False
fig9=find('图9：');pic9=fig9._p.getprevious()
for p in [find('按预设中位误差规则选取的'),find('逐角F*进一步区分了')]:pic9.addprevious(p._p)
Paragraph(pic9,d._body).paragraph_format.page_break_before=True
find('参考文献').paragraph_format.page_break_before=True

change('FlexPPI以共同的静态主链输入刻画','FlexPPI的主要价值是在同一个主链几何输入上估计主链波动与侧链联合状态两种摘要。静态特征、链内编码及同类型比较表明，残基间距离、方向与完整邻域提供了附加预测信息。共享编码器将两类输出对应到相同结构位置，便于共同组织局部动力学分析。')
delete('同类型比较使侧链输出的结构含义')
delete('完整结构输入在体系留出和全组成链家族评价中')
change('构象占据分析为两个预测目标','轨迹的连续分布与经验投影自由能为两项目标提供了分布解释。单个侧链旋转轴的主要占据区可以相近，而后续角度和共同组合仍不同；骨架位置散布也主要发生在同一侧链类别内部。联合熵读取指定角区间下的组合多样性，逐角F*展开高占据取向及其范围。两类摘要与连续投影共同连接了静态结构环境和轨迹访问的构象形态。')
delete('角区间比较明确了Hχ读取连续占据分布的方式')
change('KIX与泛素实例分别把主链','联合读数可以用于组织随后局部模拟、突变分析和功能实验需要关注的位置，再通过轨迹或干预实验检验堆积与取向约束。抗体识别中的局部环区也可沿这一思路比较结合前后构象[16]。进一步的分布模型可结合末端角化学对称性和模拟条件，学习更丰富的占据信息。')
change('FlexPPI从静态主链结构联合估计','FlexPPI从静态主链结构联合估计逐残基RMSF与侧链联合状态熵；五折体系留出评价中，完整几何环境降低了相对同训练集合对照的两项平均误差。类型内评价和结构案例显示，侧链读数能够补充主链幅度之外的位置差异。原始轨迹的连续角度与经验投影自由能分析进一步解释了两个目标对应的坐标散布和角度占据，为随后局部动力学分析提供两种互补的初步读数。')

# Refresh method and introductory references to the seven retained modules.
for p in d.paragraphs:
    if p.text.startswith('补充材料'):break
    if p._p.xpath('.//m:oMath'):continue
    t=p.text.replace('补充材料H','补充材料D').replace('补充材料I','补充材料C').replace('补充材料G','补充材料B')
    if t!=p.text:replace(p,t)
change('链内模型保留相同网络架构',find('链内模型保留相同网络架构').text.replace('覆盖率与事后分组校准见补充材料C。','本文主要评价均值预测，方差作为残差尺度参与训练。'))
change('轨迹分析沿用补充材料E的83条',find('轨迹分析沿用补充材料E的83条').text.replace('全部χ角的统一平移用于分区敏感性；Pro及原联合标签保持原定义。','所有原定义χ角及Pro保留在联合分区敏感性分析中，训练标签保持冻结。'))
change('补充五模型分析对逐位置均值等权平均','固定划分的随机初始化对照报告每个初始化模型的留出成绩。逐残基曲线按作者编号和实际链边界排列，无效位置断开。')
for p in d.paragraphs:
    if '原训练范围的五模型平均结果作为部署配置见补充材料B。' in p.text:
        replace(p,p.text.replace('原训练范围的五模型平均结果作为部署配置见补充材料B。',''))

# Rebuild supplement in the order of its seven questions; retain native cells.
sup=find('补充材料');e=sup._p.getnext()
while e is not None and e.tag!=qn('w:sectPr'):
    nxt=e.getnext();assert not e.xpath('.//w:commentRangeStart | .//w:commentRangeEnd | .//w:commentReference');e.getparent().remove(e);e=nxt
section('A 数据来源与标签观测').paragraph_format.page_break_before=False
para('ATLAS以重复轨迹统计建立标签，DynaRepo使用连续窗口，ITrace使用完整重复轨迹；各来源的观测数及有效Hχ位置见表S1。观测标准差描述不同重复或窗口的读数差异，单次观测的样本方差以零占位；零Hχ表示参与汇总的各窗口或重复内只出现一种联合状态（图S1）。')
append_table(keep['表S1'],'表S1')
d._body._body.insert(-1,label_pic);d._body._body.insert(-1,label_fig._p)
para('逐记录PDB来源注释及原API响应随数据提供[29]。原生模型方差的95%高斯区间对参考标签均值的RMSF和Hχ覆盖率分别为51.76% ± 2.71%和85.19% ± 1.31%；完整方差诊断保存在配套源数据，训练定义见2.4。')
section('B 训练方案与随机初始化')
para('固定原体系划分的第0折训练与留出集合，以种子7、11和23分别训练完整与链内模型，使用相同记录及40轮预算（表S2）。单目标模型每目标训练五个，与默认联合模型保持相同输入、编码器及预算，在相同961,321个RMSF和824,988个Hχ位置计分（表S3）。')
append_table(keep['表S2'],'表S2');append_table(keep['表S9'],'表S3')
replace(keep['表S2'][0],'表S2：原体系划分第0折的训练随机初始化对照。相同训练与留出记录分别以三个种子训练40轮，每行对应单独模型；MAE为Å或nats，ρ为Spearman相关。')
para('联合模型两项MAE相对单目标分别低约0.8%和1.8%，RMSF单目标的RMSE和排序略好。本文选择默认共享模型，在一次推理中提供两个读数。')
section('C 基线实现与来源比较')
para('发布模型保留原权重；PLANET以默认300 K、独立链输入推断，RMSF由nm转为Å。RMSF-net使用PDB模拟密度的相对柔性分数。ProDy 2.6.1的ANM/GNM由完整输入Cα构建，15 Å截断、弹簧常数1及最低20个非零模态，平方波动和开方为相对柔性分数[3,30]。相同961,321个有效RMSF位置上，各分数先逐记录计算再等权平均（表S4）。')
append_table(keep['表S11'],'表S4')
para('PLANET full在ATLAS、DynaRepo和ITrace中的MAE/ρ为0.398 Å/0.876、0.704 Å/0.752和1.412 Å/0.732；sequence-only为0.499 Å/0.826、0.782 Å/0.668和2.171 Å/0.719。ANM/GNM与完整模型的最高10%位点命中率为56.7%/57.7%及65.1%；完整来源分层与逐记录结果见源表。')
section('D 全组成链家族隔离评价')
para('按2.5节的全部组成链关联，将1,662个连通家族分在五折内，完整与链内模型各训练五次、40轮。使用与体系留出相同的有效位置并对3,063条记录等权汇总；ITrace的244条记录均在第0折测试（表S5）。完整模型相对链内模型的两项MAE低7.5%和3.8%。')
append_table(keep['表S10'],'表S5')
section('E 共同窗口与构象统计')
para('种子20261010在来源内各抽30条记录，按预设时长与时间戳条件排除7条约1 ns的DynaRepo轨迹，保留ATLAS/DynaRepo/ITrace 30/23/30条，共83条。其203条轨迹统一取0–100 ns、200 ps间隔的501帧，主分析每记录第一条，ATLAS/ITrace其余重复用于占据比较。共同有效位置为38,599个RMSF及33,493个Hχ，其中31,621个非Pro位置用于常规χ1分析。')
append_table(keep['表S6'],'表S6')
para('预测保持固定时，第一条共同窗口标签的误差高于原汇总标签；完整模型相对静态对照的两项平均误差仍较低。三条共同窗口轨迹的均值标签下，ATLAS的RMSF/Hχ MAE为0.723 Å/0.172 nats，ITrace为0.158 Å/0.110 nats，DynaRepo保持一条轨迹。')
append_table(keep['表3'],'表S7')
para('圆周平移将θ−s映射到[−180°,180°)，s=0保留原端点，s=60对应常规χ1中心分区；完整联合状态保留Pro。非Pro χ1原/60°的边界10°观测比例为61.3%/0.13%，保存帧类别变化为189.2/21.9次。全窗联合熵减五窗口帧数加权均值为0.135/0.170 nats，83条记录的联合状态半程JSD为0.050/0.067 nats。60条多重复记录的联合状态三对重复JSD为0.049/0.069 nats；同一批记录的非Pro χ1半程为0.021/0.040 nats、重复为0.022/0.040 nats。原角度、逐位置分解与来源统计随源数据提供。')
section('F KIX轨迹与状态比较')
para('四种KIX组成各使用一条100 ns生产轨迹，87个Cα对齐各自生产首帧后计算RMSF；MD设置见2.7，固定留出模型及相同位置的定位与幅值见表S8。')
append_table(keep['表S7'],'表S8')
para('L12六个状态差值的预测方向均与单轨迹参考相反。双伙伴减游离状态的MD与预测差值为+1.017与−1.589 Å，MLL减cMyb为+0.014与−1.263 Å。状态内排序与状态間幅值变化呈现不同结果；全部六对差值和前后半程统计保存在病例源表。')
section('G 双读数组合的完整统计')
para('正文图4采用824,988个共同有效位置。表S9提供各参考类别的完整计数、记录等权比例、误差及组合召回率；各类误差和召回率只对具有该类别的记录等权平均，未分类位置仍计入比例分母。类型内相关在ATLAS、DynaRepo和ITrace的参考/预测值分别为0.397/0.534、0.388/0.514和0.349/0.390。')
append_table(keep['表S12'],'表S9')
for p in d.paragraphs:
    if p._p.xpath('.//m:oMath'):continue
    t=p.text.replace('两項','两项').replace('状态間','状态间').replace('采样诊断和结合状态差值见补充材料F。','四种组成的定位与状态差值见补充材料F。')
    if t!=p.text:replace(p,t)

# Keep the precise original comment XML and all surviving native math.
d.save(OUT)
with zipfile.ZipFile(SOURCE) as a,zipfile.ZipFile(OUT) as b:
    comment=a.read('word/comments.xml');assert etree.tostring(etree.fromstring(comment),method='c14n')==etree.tostring(etree.fromstring(b.read('word/comments.xml')),method='c14n')
    entries=[(info,comment if info.filename=='word/comments.xml' else b.read(info.filename)) for info in b.infolist()]
tmp=OUT.with_suffix('.preserved.tmp')
with zipfile.ZipFile(tmp,'w') as z:
    for info,content in entries:z.writestr(info,content)
tmp.replace(OUT)
ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main','m':'http://schemas.openxmlformats.org/officeDocument/2006/math','a':'http://schemas.openxmlformats.org/drawingml/2006/main'}
with zipfile.ZipFile(SOURCE) as a,zipfile.ZipFile(OUT) as b:
    ax,bx=etree.fromstring(a.read('word/document.xml')),etree.fromstring(b.read('word/document.xml'))
    anchors=lambda x:Counter((e.tag,e.get(qn('w:id'))) for e in x.iter() if e.tag in {qn('w:commentRangeStart'),qn('w:commentRangeEnd'),qn('w:commentReference')})
    math=lambda x:[etree.tostring(e,method='c14n') for e in x.xpath('.//m:oMath',namespaces=ns)]
    assert anchors(ax)==anchors(bx) and math(ax)==math(bx)
    assert a.read('word/comments.xml')==b.read('word/comments.xml')
    assert len(bx.xpath('.//a:blip',namespaces=ns))==10
assert len(d.tables)==11
lines=[]
for i,p in enumerate(d.paragraphs):
    t=''.join(p._p.xpath('.//w:t/text() | .//m:t/text()'))
    if t:lines.append(f'P{i}: {t}')
for i,t in enumerate(d.tables):
    for j,row in enumerate(t.rows):lines.append(f'T{i}R{j}: '+' | '.join(c.text for c in row.cells))
(QA/'v23_content.txt').write_text('\n'.join(lines)+'\n',encoding='utf8')
audit={'document':str(OUT),'document_sha256':sha(OUT),'source_sha256':sha(SOURCE),'main_figures':9,'supplementary_figures':1,'main_tables':2,'supplementary_tables':9,'supplement_sections':list('ABCDEFG'),'native_math_and_comments_preserved':True,'removed_tables':['old_tableS3_ordinary_external','old_tableS4_posthoc_calibration','old_tableS5_duplicate_ENM','old_tableS8_all_state_deltas_replaced_with_key_values_and_source_archive'],'removed_supplement_figures':list(range(2,12)),'new_fig8_sha256':sha(ROOT/'figures/fig8_joint_partition_CA.png'),'old_table_cells':old_table_cells,'review':'pending'}
(QA/'integration_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps({k:v for k,v in audit.items() if k!='old_table_cells'},ensure_ascii=False))
