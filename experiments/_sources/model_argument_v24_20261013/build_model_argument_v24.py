"""Reorganize Chinese manuscript around evidenced model capabilities."""
from pathlib import Path
from copy import deepcopy
from collections import Counter
import hashlib,json,re,zipfile
from lxml import etree
from docx import Document
from docx.oxml.ns import qn
from docx.shared import Inches
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.text.paragraph import Paragraph
ROOT=Path(__file__).resolve().parent;BASE=ROOT.parent.parent
QA=ROOT/'manuscript_QA';QA.mkdir(exist_ok=True)
SOURCE=BASE/'paper_v23_20261012_中文修订_论证精简与补充整合_三线表.docx'
OUT=BASE/'paper_v24_20261013_中文修订_模型论证与构象占据验证_三线表.docx'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
assert sha(SOURCE)=='04e952548a9e2de9091f138d510238f93e816f5731ad801bee9af5f07fcbf4f8'
assert json.loads((QA/'predictive_occupancy_reconstruction_audit.json').read_text(encoding='utf8'))['status']=='passed'
values=json.loads((ROOT/'manuscript_values.json').read_text(encoding='utf8'))
d=Document(SOURCE);original=list(d.paragraphs)
P={f'P{i}':p for i,p in enumerate(original)}
changed=[];removed=[]
def setp(p,text):
    assert not p._p.xpath('.//m:oMath')
    nodes=p._p.xpath('.//w:t');assert nodes
    nodes[0].text=text;nodes[0].set('{http://www.w3.org/XML/1998/namespace}space','preserve')
    for n in nodes[1:]:n.text=''
def change(key,text):
    changed.append({'source':key,'before':P[key].text,'after':text});setp(P[key],text)
def remove(key):
    p=P[key];assert not p._p.xpath('.//w:commentRangeStart | .//w:commentRangeEnd | .//w:commentReference | .//m:oMath')
    removed.append({'source':key,'text':p.text});p._p.getparent().remove(p._p)
def append_after(p,text):
    new=deepcopy(P['P99']._p)
    for child in list(new):
        if child.tag!=qn('w:pPr'):new.remove(child)
    pp=Paragraph(new,d._body);pp.add_run(text)
    pp.paragraph_format.page_break_before=False;pp.paragraph_format.keep_with_next=False
    p._p.addnext(new);return pp

# Adopt model-centred Results wording, then override restructuring decisions.
proposed=json.loads((ROOT/'editorial_results/results_discussion_revisions.json').read_text(encoding='utf8'))
for key,text in proposed['replacements'].items():
    if key not in {'P179','P184','P117','P118','P120','P166','P175'}:change(key,text.replace('两項','两项'))
change('P75',P['P75'].text.replace('Hχ接近零时比值误差较大，因此后续类型内比较使用Å和nats的原单位MAE。',''))
change('P72',P['P72'].text.replace('来源分层显示不同训练资源与输入方式下的幅值和排序表现。',''))
change('P101',P['P101'].text.replace('为上述两个预测摘要提供构象解释。','为这两个位点的联合熵预测提供构象解释。'))
change('P0','FlexPPI 基于静态主链几何的残基波动与侧链联合状态熵预测')
change('P2','主链位置波动与侧链取向分布提供不同的局部动态信息，需要在同一结构上联合估计。我们提出FlexPPI，以氨基酸类型及N、Cα和C坐标为输入，通过残基局部框架与不变点注意力编码链内和链间几何，联合预测逐残基均方根波动（RMSF）与固定χ角区间下的联合状态熵（Hχ）。在3,063条记录的五折体系留出评价中，相对静态特征对照，两項平均绝对误差分别降低14.8%和11.3%；类型内比较显示，几何输入提供了残基种类之外的位置差异。在ThermoPCD外部集合中，完整模型的Hχ误差相对链内模型降低9.4%。KIX与泛素–UBA实例分别检验活动区域排序和同类侧链状态差异。进一步在83条记录中，将冻结折外预测与MD连续角度的细箱占据熵对应，侧链输出的同类型排序相关为0.510，高于主链输出的0.321；预测确定的高值组在全部记录中具有较高的真实占据熵。泛素的多角分布与经验投影自由能展开这一预测差异的具体取向形态。FlexPPI在一次推理中提供对应到相同残基位置的两个动态摘要，可用于单链及复合物的局部动力学分析。'.replace('两項','两项'))
intro=[
'蛋白质的局部构象变化影响分子识别、催化和功能调节。腺苷酸激酶的多时间尺度运动与催化过程相关，钙调蛋白的侧链运动也参与分子识别中的构象熵变化[25,26]。主链位置散布与侧链取向覆盖提供不同信息：均方根波动（RMSF）描述对齐后Cα位置围绕其均值的波动幅度；侧链联合状态熵（Hχ）依据多个χ角在指定区间内共同出现的频率，描述取向组合的多样性。长时间模拟显示，紧密堆积的侧链也可访问多个取向，联合状态占据能够表征这种局部多样性[32]。在相同残基位置估计这两类读数，有助于同时定位主链活动区域与侧链状态差异。',
'已有方法从不同输入估计蛋白质动态。弹性网络方法由残基间连接描述集体涨落，ProDy提供相应模型与分析工具[3]；PEGASUS从序列估计MD派生柔性指标，Flexpert结合序列与结构预测RMSF并用于柔性设计，BackFlip系列利用几何表示学习结构波动与运动耦合[4–6,9]。RMSF-net结合原子结构和冷冻电镜密度信息，AlphaFlow生成构象集合，PLANET-MD从序列及可选结构表示预测RMSF、动态相关和结构异质性[10,11,31]。本文聚焦从静态主链几何同时估计主链波动与侧链多角联合状态占据，使两种摘要能够对应到同一输入结构。',
'联合预测需要区分氨基酸类型先验与局部几何环境。不同残基具有不同的侧链角数和理论状态范围，类型平均值能够解释部分整体排序；模型还需要判断同一种氨基酸在不同位置的状态差异。复合物中的局部环境同时包括链内堆积和邻近伙伴接触，输入表示也需要兼容独立链与完整多链结构。几何学习在蛋白界面预测中展示了三维邻域的预测价值[14,15]。本文以残基距离、相对方向和局部主链朝向描述空间环境，并通过类型内及完整／链内模型对照检验这些信息的作用。',
'FlexPPI以N、Cα和C构建每个残基的局部坐标框架，通过不变点注意力共同聚合残基、残基对和三维点信息[7,13]。同一编码器处理输入中的全部有效链，两个读出头分别生成RMSF与Hχ的均值和残差尺度。局部坐标使邻域方向与残基主链朝向共同参与表示，完整结构输入使链内和链间接触在同一几何框架中参与计算。模型采用联合训练，在一次推理中给出对应到相同残基位置的两个动态摘要。',
'我们整合ATLAS、DynaRepo与ITrace的轨迹标签，在五折体系留出评价中比较完整模型、链内模型、静态特征与类型均值对照[1,2,21]。整体及类型内成绩检验几何输入和局部环境的预测信息，公开权重比较确定其RMSF性能位置；ThermoPCD外部队列考察另一轨迹来源中的幅值与排序。KIX和泛素–UBA实例分别对照活动区域定位与同类型侧链差异。进一步固定83条记录的折外输出，检验侧链预测能否识别实际连续角度占据的差异，并用泛素的角度密度与经验投影自由能解释预测对应的构象形态。各项评价共同检验两个结构读数的准确性及其在局部动力学分析中的含义。']
for i,text in enumerate(intro,5):change('P'+str(i),text)
for key in ['P10','P11']:remove(key)

# Methods define actual inputs and measurements positively, with no repeated defenses.
change('P14','FlexPPI以单链或多链静态结构为输入，使用氨基酸类型、链归属、链内顺序及N、Cα和C坐标，估计每个残基的RMSF与联合侧链状态熵Hχ。共享几何编码器整合局部主链朝向、残基距离和相对方向，两个读出头分别生成目标均值及用于异方差回归的残差尺度。MD轨迹与侧链原子用于建立监督标签，推理使用静态主链几何。网络输入由上述几何与类型特征构成，标签观测设计见2.2；图1概括输入、编码与读出流程。')
for node in P['P25']._p.xpath('.//w:t'):
    node.text=(node.text or '').replace('该指标描述给定角度分区与采样时长内的取向多样性，不等同于完整热力学构象熵，也不直接表示侧链转换速率。','模型估计固定角区间和观测窗口下的联合占据多样性。')
for node in P['P26']._p.xpath('.//w:t'):
    node.text=(node.text or '').replace('，而非标签均值的标准误','').replace('表示无法估计观测间变异，不能据此认为该标签没有不确定性。连续窗口也不视为独立模拟重复。','连续窗口与独立重复按各自的观测设计汇总。')
change('P31','链内对照保留相同架构和参数量，分别编码各条链；完整模型共同编码全部有效输入链。两个模型分别从头训练，比较编码范围对预测成绩的影响。默认联合模型共享几何编码器，并为RMSF和Hχ保留各自的均值与对数方差读出参数。读出层经归一化、线性映射和SiLU激活，均值经softplus变为非负值，对数方差限制为[−8,6]；方差头提供式(3)中的残差尺度。')
change('P37','式(4)对记录和具有监督的目标等权聚合；有效掩码确定各项目标参与损失的位置。')
change('P38',P['P38'].text.replace('所有模型从头训练40轮，统一使用第40轮权重，不依据测试成绩选择训练轮次，也不另设验证集。','全部模型从头训练40轮，统一以第40轮权重计分。').replace('主分析使用固定随机种子，并在固定划分上进行不同种子的重复训练，以观察训练随机性；完整配置随代码提供。','主分析使用固定随机种子；第0折随机初始化对照以种子7、11和23分别训练相同训练与留出集合，完整配置随代码提供。'))
change('P40',P['P40'].text.replace('[19,20]',''))
change('P46',P['P46'].text.replace('该模型保留原研究的训练权重，其训练数据与FlexPPI不同，因此发布模型比较与使用本研究相同训练标签拟合的对照分别解释。','发布模型保留原训练权重；静态特征及类型均值使用本文各对应训练集合拟合。')) if not P['P46']._p.xpath('.//m:oMath') else None
for node in P['P46']._p.xpath('.//w:t'):
    node.text=(node.text or '').replace('该模型保留原研究的训练权重，其训练数据与FlexPPI不同，因此发布模型比较与使用本研究相同训练标签拟合的对照分别解释。','发布模型保留原训练权重；静态特征及类型均值使用本文各对应训练集合拟合。')
change('P51',P['P51'].text.replace('选择依据为相关训练记录A00AQ位于该折测试集合，而非本次MD的预测成绩。','相关记录A00AQ位于该折测试集合。'))
change('P54','泛素侧链溶剂可及面积由完整复合物重原子计算，探针半径为1.4 Å。侧链取向叠加使用实际MD代表帧，经局部主链对齐后展示；逐残基预测、有效位置、联合状态频率及代表帧来源随配套数据提供。')
change('P55','2.8 联合预测与构象占据的验证')
change('P56','将冻结预测与同一残基的轨迹占据对应，可以检验读数所识别的动态差异。侧链标签由多个角的联合类别占比计算；将类别c的占比记为pc，最大占比记为pmax，其相对经验投影自由能以最常见类别为零点。式(2)可以改写为')
change('P58','式(5)将模型的侧链监督目标与联合类别的相对占据联系起来，熵取决于各类别的占比和最常见类别的集中程度。逐角投影进一步以5°箱中的经验概率密度计算F*(q)，高占据角度对应低F*，低值范围展示该旋转轴的采样形态[32,33,35]。泛素实例将留出预测对应到这些真实角度分布，分析同类型位点预测差异的构象来源。')
change('P60','连续占据评价使用83条冻结记录，ATLAS、DynaRepo和ITrace分别为30、23和30条，共有203条可用轨迹。全部轨迹统一取0–100 ns、200 ps间隔的501个实际保存帧；主分析使用每记录第一条轨迹，重复评价对该记录全部可用重复的同一占据指标取均值。提取全部可定义χ角并保留Pro、化学名称校正和逐残基映射，Full与Intra使用对应留出折的原始预测。共同有效位置要求原RMSF与Hχ标签、轨迹角度及Cα均有效，预测值均有限；获得33,493个位点。')
change('P61','泛素Q41与Q62沿用前5,000帧和五个标签窗口。各χ角的经验密度由5°箱出现频率除以箱宽得到，F*(q)=−ln[p̂(q)/p̂max]，最高密度位置为零，采用无量纲热能标度；零计数箱断线。密度与投影曲线合并全部帧，窗口曲线分别计算。该分析展开图6中同一留出预测位点的实际取向范围。')
change('P62','群体评价的轨迹指标为平均细箱角度占据熵：对每个可定义χ角计算5°细箱Shannon熵，箱格原点取0°、1.25°、2.5°和3.75°，先对四个原点平均，再对该残基的全部角平均，单位为nats/角。此指标细化每个旋转轴的占据范围。记录内同一种氨基酸至少具有四个共同有效位置时参与评价，保留33,142个位点；分别从该轨迹指标及Full Hχ、Intra Hχ、Full RMSF折外预测减去记录内同类型均值，在相同位置上计算Spearman相关，再对记录等权平均。')
append_after(P['P62'],'位置筛选仅使用各预测读数：同记录同类型内按预测排序，最低和最高两组各取位点数的一半并向下取整，并列值以输入残基顺序确定；奇数类型的中间位置仅不参加高低组比较。计算两组真实细箱熵之差，先对类型等权平均，再对记录等权平均。相关及组差均按来源分层重抽记录10,000次，种子20261013，报告均值的95%百分位区间；方法差值对同一记录配对。1至4个χ角的分层和可用重复均值评价采用相同指标、位置及聚合规则。')

# Model -> local prediction -> actual conformations, followed by a cohort test.
for key in ['P94','P95']:P['P101']._p.addprevious(P[key]._p)
for key in ['P102','P106','P107','P110','P111','P112']:remove(key)
pic9=P['P114']._p.getprevious();assert pic9.xpath('.//w:drawing');pic9.getparent().remove(pic9);remove('P114')
change('P105','3.6.2 折外侧链预测对角度占据差异的群体验证')
core='在83条记录的33,142个同类型有效位置上，完整模型Hχ预测与真实平均细箱角度占据熵的记录平均相关为0.510，链内Hχ为0.493，完整模型RMSF为0.321（图8a）。完整侧链输出相对主链输出的配对相关增量为0.189，95%重采样区间为0.172–0.206；相对链内侧链输出的增量为0.017，区间为0.006–0.027。比较控制了每条记录内的氨基酸类型，并使用完全相同的位置和连续角度指标，支持侧链输出识别主链幅度之外的实际占据差异。'
new1=append_after(P['P105'],core)
new2=append_after(new1,'仅按完整模型Hχ预测确定高、低半组时，真实平均细箱占据熵之差为0.293 nats/角，95%区间为0.277–0.307，83条记录均为正值（图8b）。链内侧链预测与完整模型主链预测确定的组差分别为0.281和0.198 nats/角。1至4个χ角分层中的完整侧链相关为0.529、0.503、0.528和0.591；使用全部可用重复的指标均值后，侧链／主链相关为0.541／0.342（表S7）。因此，模型从静态几何给出的侧链高低值能够用于筛选同类型中角度占据更分散的位置；泛素的自由能投影展示了这一群体联系在具体旋转轴上的形态。')
fig8=P['P109'];pic8=fig8._p.getprevious();assert pic8.xpath('.//w:drawing')
for child in list(pic8):
    if child.tag!=qn('w:pPr'):pic8.remove(child)
pp=Paragraph(pic8,d._body);pp.paragraph_format.page_break_before=False;pp.paragraph_format.keep_with_next=True;pp.alignment=WD_ALIGN_PARAGRAPH.CENTER
pp.add_run().add_picture(str(ROOT/'figures/fig8_model_angular_occupancy.png'),width=Inches(6.2))
setp(fig8,(ROOT/'figures/fig8_model_angular_occupancy.caption_zh.md').read_text(encoding='utf8').strip())
change('P117','侧链输出的物理意义由预测对照和轨迹占据共同建立。泛素留出预测首先识别类型均值之外的位置差异，并在Q41与Q62上保留联合熵的高低顺序；联合状态覆盖和后续χ角的自由能投影进一步展开该顺序对应的取向范围。83条记录中，侧链预测与细箱占据的同类型相关高于主链预测，预测高低组的真实占据差异也在全部记录中同向。这使侧链输出的价值从对固定标签的回归，延伸到同类型位点实际角度占据差异的筛选。')
change('P118','两类输出对应不同的筛选问题：RMSF排序定位主链活动区域，Hχ区分侧链联合状态多样性。KIX中的共同高波动环区和泛素中的同类型位点分别给出具体验证；全队列组合召回则量化两种信息同时用于定位时的表现，低／高与高／低组合的召回率为15.8%和16.5%。据此，联合读数可用于优先确定局部模拟和功能实验需要关注的位置，再围绕目标运动类型安排后续分析。')
change('P120','FlexPPI从静态主链几何联合估计逐残基RMSF与侧链联合状态熵。五折体系留出和类型内评价表明，几何输入提高了相对同训练集合对照的两项平均估计精度；KIX与泛素分别展示活动区域定位和同类型侧链多样性区分。冻结折外预测对连续角度占据的群体验证进一步支持侧链输出筛选实际取向差异的能力，泛素自由能投影解释这种差异的具体形态。两类输出在同一结构位置提供互补的局部动力学读数。')

# Supplement retains model comparisons; replace pure-label table with new validation.
change('P166','逐记录PDB来源注释及原API响应随数据提供[29]。标签观测方差和模型残差尺度按2.4节用于训练，逐位诊断随配套源数据保存。')
change('P179','E 共同窗口与预测驱动的构象占据验证')
change('P180',P['P180'].text.replace('ATLAS/ITrace其余重复用于占据比较。','ATLAS/ITrace其余重复用于可用重复指标均值评价。').replace('，其中31,621个非Pro位置用于常规χ1分析。','，全部可定义角及Pro保留在本节预测验证中。'))
oldtab=d.tables[8];table_data=__import__('pandas').read_csv(ROOT/'tableS7_predictive_occupancy.csv')
header=['角数分层','记录数','有效位置','Full Hχ ρ','Intra Hχ ρ','RMSF ρ','高减低']
change('P183','表S7：冻结预测与细箱角度占据的同类型评价。各角数层内沿用同一记录、同一氨基酸至少四个位点的控制；ρ为记录平均的同类型去均值Spearman相关，最后一列为仅按完整模型Hχ预测确定高、低半组的真实细箱熵之差，单位nats/角。各列使用相同有效位置，先逐记录计分再等权平均；1至4个χ角层的类型组成各不相同。')
bodyrows=[]
for _,row in table_data.iterrows():bodyrows.append(['全部' if row.scope=='all' else row.scope[1:]+'角',str(int(row.n_records)),f'{int(row.positions):,}',*[f'{row[k]:.3f}' for k in ['Full_Hchi_rho','Intra_Hchi_rho','Full_RMSF_rho','Full_Hchi_high_low']]])
# Create each row/cell from existing native table formatting; seven columns fit page.
tbl=oldtab._tbl
row_template=deepcopy(tbl.xpath('./w:tr')[1]);header_template=deepcopy(tbl.xpath('./w:tr')[0])
for e in list(tbl):
    if e.tag in {qn('w:tr'),qn('w:tblGrid')}:tbl.remove(e)
grid=etree.Element(qn('w:tblGrid'))
for width in [1300,880,1250,1270,1270,1170,1170]:etree.SubElement(grid,qn('w:gridCol')).set(qn('w:w'),str(width))
tbl.insert(1,grid)
for i,texts in enumerate([header]+bodyrows):
    row=deepcopy(header_template if i==0 else row_template)
    oldcells=row.xpath('./w:tc')
    for child in list(row):
        if child.tag==qn('w:tc'):row.remove(child)
    for j,text in enumerate(texts):
        cell=deepcopy(oldcells[min(j,len(oldcells)-1)])
        width=cell.xpath('./w:tcPr/w:tcW')[0];width.set(qn('w:w'),str([1300,880,1250,1270,1270,1170,1170][j]))
        borders=cell.xpath('./w:tcPr/w:tcBorders')[0]
        for e in borders:e.set(qn('w:val'),'single' if etree.QName(e).localname in ({'top','bottom'} if i==0 else {'bottom'} if i==len(bodyrows) else set()) else 'nil')
        nodes=cell.xpath('.//w:t');nodes[0].text=text
        for e in nodes[1:]:e.text=''
        for sz in cell.xpath('.//w:sz | .//w:szCs'):sz.set(qn('w:val'),'20')
        row.append(cell)
    tbl.append(row)
change('P184','共同窗口中的33,493个位点有351个因记录内同类型位置少于四个而不参加类型控制，保留33,142个位点和1,296个记录内类型组。高低组比较对各类型等权汇总，奇数类型的中间位置按2.8节处理。使用全部可用重复指标均值时，Full Hχ、Intra Hχ和Full RMSF的平均相关分别为0.541、0.527和0.342，按Full Hχ确定的真实高低差为0.290 nats/角。完整来源、角数、逐记录及预测分组统计随源数据提供。')

# Remove references without a current argument and consistently number in order of use.
refparas={int(re.match(r'^\[(\d+)\]',p.text).group(1)):p for p in d.paragraphs if re.match(r'^\[(\d+)\]',p.text)}
pattern=re.compile(r'\[([\d,，\s–\-]+)\]')
def numbers(s):
    ids=[]
    for item in re.split('[,，]',s):
        bounds=re.split('[–-]',item.strip())
        ids.extend(range(int(bounds[0]),int(bounds[1])+1) if len(bounds)==2 else [int(bounds[0])])
    return ids
order=[]
for p in d.paragraphs:
    if re.match(r'^\[\d+\]',p.text):continue
    for match in pattern.finditer(p.text):
        for n in numbers(match.group(1)):
            if n in refparas and n not in order:order.append(n)
renumber={old:i+1 for i,old in enumerate(order)}
def mapped(match):return '['+','.join(str(renumber[n]) for n in numbers(match.group(1)))+']'
for p in d.paragraphs:
    if re.match(r'^\[\d+\]',p.text):continue
    for node in p._p.xpath('.//w:t'):
        node.text=pattern.sub(mapped,node.text or '')
for old,p in refparas.items():
    if old not in renumber:
        assert not p._p.xpath('.//w:commentRangeStart | .//w:commentReference');p._p.getparent().remove(p._p)
    else:
        nodes=p._p.xpath('.//w:t');assert nodes[0].text.startswith('['+str(old)+']')
        nodes[0].text=re.sub(r'^\[\d+\]',f'[{renumber[old]}]',nodes[0].text)
reference_heading=P['P124']._p
for old in order:
    reference_heading.addnext(refparas[old]._p);reference_heading=refparas[old]._p
# Clean empty figure-only spacer paragraphs left from removed content.
for p in list(d.paragraphs):
    if not p.text.strip() and not p._p.xpath('.//w:drawing | .//m:oMath | .//w:commentRangeStart | .//w:commentReference'):
        p._p.getparent().remove(p._p)

d.save(OUT)
# Retain original comment XML exactly, independent of python-docx serialization.
with zipfile.ZipFile(SOURCE) as a,zipfile.ZipFile(OUT) as b:
    entries=[(info,a.read(info.filename) if info.filename=='word/comments.xml' else b.read(info.filename)) for info in b.infolist()]
tmp=OUT.with_suffix('.preserved.tmp')
with zipfile.ZipFile(tmp,'w') as z:
    for info,data in entries:z.writestr(info,data)
tmp.replace(OUT)
ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main','m':'http://schemas.openxmlformats.org/officeDocument/2006/math','a':'http://schemas.openxmlformats.org/drawingml/2006/main'}
with zipfile.ZipFile(SOURCE) as a,zipfile.ZipFile(OUT) as b:
    ax,bx=etree.fromstring(a.read('word/document.xml')),etree.fromstring(b.read('word/document.xml'))
    anchor=lambda x:Counter((e.tag,e.get(qn('w:id'))) for e in x.iter() if e.tag in {qn('w:commentRangeStart'),qn('w:commentRangeEnd'),qn('w:commentReference')})
    math=lambda x:[etree.tostring(e,method='c14n') for e in x.xpath('.//m:oMath',namespaces=ns)]
    assert anchor(ax)==anchor(bx) and math(ax)==math(bx)
    assert a.read('word/comments.xml')==b.read('word/comments.xml')
    assert len(bx.xpath('.//a:blip',namespaces=ns))==9
assert len(d.tables)==11
audit={'document':str(OUT),'document_sha256':sha(OUT),'source_sha256':sha(SOURCE),'main_figures':8,'supplementary_figures':1,'main_tables':2,'supplementary_tables':9,'references':len(order),'removed_old_references':sorted(set(refparas)-set(order)),'reference_mapping':renumber,'source_math_and_comments_preserved':True,'changes':changed,'removed_paragraphs':removed,'new_figure8_sha256':sha(ROOT/'figures/fig8_model_angular_occupancy.png'),'new_tableS7':bodyrows,'review':'pending'}
(QA/'integration_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf8')
lines=[]
for i,p in enumerate(d.paragraphs):
    t=''.join(p._p.xpath('.//w:t/text() | .//m:t/text()'))
    if t:lines.append(f'P{i}: {t}')
for i,t in enumerate(d.tables):
    for j,row in enumerate(t.rows):lines.append(f'T{i}R{j}: '+' | '.join(c.text for c in row.cells))
(QA/'v24_content.txt').write_text('\n'.join(lines)+'\n',encoding='utf8')
print(json.dumps({k:v for k,v in audit.items() if k not in {'changes','removed_paragraphs','new_tableS7'}},ensure_ascii=False))
