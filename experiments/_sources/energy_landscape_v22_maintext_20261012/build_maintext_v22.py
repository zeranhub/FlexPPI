"""Integrate the energy/distribution argument into the manuscript main text."""
from pathlib import Path
from copy import deepcopy
from collections import Counter
import hashlib,json,zipfile
from lxml import etree
from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn

ROOT=Path(__file__).resolve().parent
PREV=ROOT.parent/'energy_landscape_v21_20261012'
BASE=ROOT.parent.parent
SOURCE=BASE/'paper_v21_20261012_中文修订_构象分布与自由能分析_三线表.docx'
OUTPUT=BASE/'paper_v22_20261012_中文修订_自由能分析融入正文_三线表.docx'
QA=ROOT/'manuscript_QA';QA.mkdir(exist_ok=True)
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
assert sha(SOURCE)=='0773be657cbc977000f4326ae0eabbdea5c782ecc6b7308ccf9c80978f014ac8'
doc=Document(SOURCE)
def find(start):
    p=[p for p in doc.paragraphs if p.text.startswith(start)]
    assert len(p)==1,(start,len(p))
    return p[0]
def replace(p,text):
    nodes=p._p.xpath('.//w:t');assert nodes and not p._p.xpath('.//m:oMath')
    nodes[0].text=text;nodes[0].set('{http://www.w3.org/XML/1998/namespace}space','preserve')
    for node in nodes[1:]:node.text=''
def before(anchor,text,style='Normal'):
    p=doc.add_paragraph(text,style);anchor._p.addprevious(p._p);return p
def remove(p):p._p.getparent().remove(p._p)
def picture_for(caption):
    e=caption._p.getprevious();assert e.xpath('.//w:drawing');return e
def main_picture(anchor,picture,caption,text):
    picture.find(qn('w:pPr'))
    anchor._p.addprevious(picture)
    caption.paragraph_format.keep_with_next=False
    replace(caption,text)
    anchor._p.addprevious(caption._p)
def M(name,*children):
    e=OxmlElement('m:'+name)
    for c in children:e.append(c)
    return e
def R(text):return M('r',MT(text))
def MT(text):
    e=OxmlElement('m:t');e.text=text;return e
def sub(base,index):return M('sSub',M('e',R(base)),M('sub',R(index)))
def fstar_c():return M('sSubSup',M('e',R('F')),M('sub',R('c')),M('sup',R('*')))
def equation(anchor):
    p=before(anchor,'');p.alignment=WD_ALIGN_PARAGRAPH.CENTER
    math=M('oMath')
    math.append(fstar_c());math.append(R(' = −ln'))
    math.append(M('f',M('num',sub('p','c')),M('den',sub('p','max'))))
    math.append(R(',     '));math.append(sub('H','χ'));math.append(R(' = '))
    nary=M('nary');props=M('naryPr');char=M('chr');char.set(qn('m:val'),'∑');props.append(char)
    lim=M('limLoc');lim.set(qn('m:val'),'subSup');props.append(lim)
    hide=M('supHide');hide.set(qn('m:val'),'1');props.append(hide);nary.append(props)
    nary.append(M('sub',R('c: '),sub('p','c'),R(' > 0')));nary.append(M('sup'))
    nary.append(M('e',sub('p','c'),fstar_c()));math.append(nary)
    math.append(R(' − ln '));math.append(sub('p','max'));math.append(R('.     (5)'))
    p._p.append(math)
    p.paragraph_format.keep_with_next=True
    return p

# Keep every scientific image part and native table; promote their existing nodes.
fig7=find('图7：');pic7=picture_for(fig7)
fig8=find('图S12：');pic8=picture_for(fig8)
fig9=find('图S13：');pic9=picture_for(fig9)
table3=find('表S13：');table3_node=table3._p.getnext();assert table3_node.tag==qn('w:tbl')
heading=find('3.6 从柔性读数到连续构象分布')
discussion=find('4 讨论')
node=heading._p.getnext()
while node is not discussion._p:
    next_node=node.getnext()
    if node is not pic7 and node is not fig7._p:node.getparent().remove(node)
    node=next_node
replace(heading,'3.6 双动态读数对应的构象占据与投影自由能')
before(discussion,'RMSF与Hχ概括了构象集合中的不同变化：RMSF读取对齐后的主链位置散布，Hχ读取多个侧链角共同出现的状态及其占比。将角度占据率转换为相对经验投影自由能F*，高占据取向表现为较低的F*区域，曲线的低值范围和多个低值区域进一步显示侧链访问的角度形态。下面把这些曲线与同一位置的参考标签及折外预测共同分析，解释主链幅度相近时侧链状态多样性仍可不同的原因。')
before(discussion,'3.6.1 泛素侧链的多角度分布','Heading 3')
before(discussion,'泛素Q41与Q62的参考RMSF为0.469和0.615 Å，模型分别给出0.482和0.574 Å；两者均处于局部较小的主链波动范围。五个连续窗口的平均联合Hχ分别为1.890和2.631 nats，预测为2.000和2.684 nats。连续密度与F*把这一区别展开到各条侧链旋转轴：两个位点的χ1均主要采样−60°附近，原分区χ1边缘熵为0.837和0.833 nats；Q41的χ2主要集中于跨越±180°接缝的同一邻域，Q62的χ2主要采样−60°附近，χ3覆盖更分散。χ2边缘熵分别为0.680和0.835 nats，χ3为0.574和1.032 nats（图7）。Q62较高的联合读数因此来自多个角度及其共同组合的覆盖，单独χ1曲线没有呈现同样的高低差异。')
before(discussion,'把三个角的区间共同平移60°后，Q41/Q62的五窗口平均联合熵为0.528/1.399 nats，Q62仍较高；常规χ1分区的边缘熵却为0.492/0.403 nats，以Q41较高。联合状态和单角投影保留了不同信息：前者汇总一个完整侧链的取向组合，后者展开某一旋转轴的占据。图7的合并帧F*及各窗口曲线使这种多角度差异在具体残基上可直接观察。')
main_picture(discussion,pic7,fig7,fig7.text)

before(discussion,'3.6.2 角度占据与主链位置波动','Heading 3')
before(discussion,'83条记录的统一窗口分析显示，同一构象分布可以在不同角区间下形成不同的离散摘要。保持全部采样帧固定，完整联合Hχ的记录均值从原分区的1.484 nats变为60°平移后的0.962 nats；非Pro χ1边缘熵由0.717变为0.374 nats（表3，图8a、b）。原边界附近的观测比例为61.3%，常规χ1中心分区为0.13%；相邻保存帧的类别变化由189.2次变为21.9次（图8c）。原边界位于常见χ1峰附近时，围绕同一低F*区域的局部摆动也会分散到相邻类别，中心分区则把这些观测归入同一类。连续密度保持相同，熵值的变化具体显示了区间边界如何读取该分布。')
replace(table3,table3.text.replace('表S13：','表3：'))
discussion._p.addprevious(table3._p);discussion._p.addprevious(table3_node)
before(discussion,'主链与侧链的联系还可以通过同一批帧中的位置波动分解来检验。按侧链联合状态对链内对齐后的Cα坐标分组，RMSF的平方等于每个状态内部的均方波动，加上不同状态平均位置之间的均方波动；两项均按实际占据率加权。状态平均位置之间的贡献在原分区和60°平移后，逐位置计算比例再对记录等权平均，分别占总均方波动的6.2%和7.8%。这说明在当前状态分区中，大部分主链幅度来自同一侧链类别内部的位置散布，而Hχ随类别占比和多角度组合改变。两类读数因此能够呈现不同的局部高低组合，侧链类别之间的骨架位移也可从同一轨迹中单独读出。')
before(discussion,'时间组织进一步补充占据分布的形态。全窗联合熵比五个连续窗口熵的帧数加权平均高0.135和0.170 nats，分别对应原分区和60°平移，说明合并窗口还纳入了窗口之间的占据差异。在同一批60条多重复记录中，非Pro χ1前后半程占据JSD的记录均值为0.021和0.040 nats，全部重复对为0.022和0.040 nats（图8d）。这些比较将角度覆盖范围与采样时间内的占据变化共同对应到动态读数。')
caption8='图8：角区间和时间采样对侧链占据读数的影响。（a）83条记录在0°、30°、60°和90°整体平移下的完整联合Hχ，保留原定义的全部χ角及Pro。（b,c）非Pro χ1在原分区与以−60°、+60°和180°为中心的分区下，其边缘熵和相邻保存帧类别变化次数；±180°接缝两侧的trans合为同一类。（d）同一批60条ATLAS/ITrace记录的非Pro χ1前后半程及三条重复轨迹全部三对占据JSD，单位为nats。灰点和连线覆盖全部对应记录，彩色符号为来源均值，黑菱形为记录等权均值。统一窗口为0–100 ns、200 ps间隔的501帧，半程为0–49.8 ns和50–100 ns；（a–c）使用每记录第一条轨迹，提取与统计方法见2.8。'
main_picture(discussion,pic8,fig8,caption8)

case_heading=before(discussion,'3.6.3 不同蛋白中的局部自由能形态','Heading 3')
case_heading.paragraph_format.page_break_before=True
before(discussion,'为将这一解释对应到不同结构来源，我们按原折外两项误差的来源内平均秩最接近中位数的规则，选取ATLAS 2py5_A、DynaRepo A0022与ITrace 6p64，并在各记录内选取同类型但两读数不同步的位点对。2py5_A的Y254/Y521参考RMSF为0.630/1.850 Å、Hχ为1.573/0.358 nats；A0022的I78/I21为0.343/0.464 Å、1.694/0.794 nats；ITrace的E链T149/D链T25为0.724/1.735 Å、1.046/0.514 nats。三组均表现为主链幅度较低的位置具有较高的联合侧链熵，完整剖面同时给出这些位置的折外预测与参考值（图9）。')
before(discussion,'逐角F*进一步区分了这些联合读数所对应的形态。2py5_A的Y254和Y521在第一条轨迹中分别集中于χ1的负向和正向取向区，5°箱的最高占据位置为−57.5°和67.5°；Y254较高的多角度联合熵与χ1单个主要占据区并存。A0022的I78主要位于负向区，并有4.8%的保存帧访问正向区，I21则集中于负向区。ITrace的E链T149有约47.9%帧处于正向区、51.9%处于trans区，F*显示两个分离的高占据区域；D链T25集中于trans区，接缝两端对应同一角度邻域。这些曲线把不同位点的主链幅度、侧链多角度状态和具体旋转轴的取向约束放到一起，说明联合读数如何帮助选定位点，再由局部分布解释其构象形态。')
caption9='图9：不同来源蛋白的动态剖面与局部经验投影自由能。（a–d）ATLAS 2py5_A，（e–h）DynaRepo A0022，（i–l）ITrace 6p64。各来源在83条记录子集内，按两项原折外MAE的来源内百分位秩均值最接近0.5选取代表记录，三条记录均为0.5。每行前两面板展示完整输入残基顺序上的原MD标签与折外预测，灰色竖线为链边界；后两面板为选定同类型位点的χ1全部501个保存帧及5°箱宽F*。位点按参考Hχ秩减RMSF秩最相离选取，粗线为第一条轨迹，淡线为其他已有重复，零计数箱断开。原标签保留其观测汇总口径，角度窗口为0–100 ns、200 ps间隔；位点编号对应模型输入PDB，筛选方法见2.8。'
main_picture(discussion,pic9,fig9,caption9)
# Pair the section heading with its full scientific plate; put the connected
# interpretation on the following page rather than leaving a few orphan lines.
case_heading._p.addnext(pic9);pic9.addnext(fig9._p)
break_node=pic9.find(qn('w:pPr')).find(qn('w:pageBreakBefore'))
if break_node is not None:break_node.getparent().remove(break_node)
for p in doc.paragraphs:
    if p.text.startswith('泛素Q41与Q62的参考RMSF'):
        replace(p,p.text.replace('原分区χ1边缘熵','五窗口平均的原分区χ1边缘熵').replace('χ2边缘熵分别为','五窗口平均的χ2边缘熵分别为'))
    if p.text.startswith('把三个角的区间共同平移'):
        replace(p,p.text.replace('常规χ1分区的边缘熵','常规χ1分区的五窗口平均边缘熵'))
    if '大部分主链幅度来自' in p.text:
        replace(p,p.text.replace('大部分主链幅度来自','大部分主链均方波动来自'))

# Explain the probability connection at the point where the metric is defined.
replace(find('2.8 构象分布与双读数分析'),'2.8 构象占据与经验投影自由能分析')
pair=find('两类读数的配对分析使用全部')
before(pair,'RMSF与Hχ从轨迹采样的构象集合中概括不同信息。RMSF的平方是对齐后Cα位置相对均值的平方偏离平均，Hχ由多个侧链角联合类别的占据比例计算。将同一窗口内某个联合类别c的占比记为pc，最大占比记为pmax，其相对经验投影自由能F*c以最常见类别为零点；式(2)可以改写为')
equation(pair)
before(pair,'式(5)把联合熵与同一组类别的占据自由能联系起来：熵同时取决于各类别的相对占据和最常见类别的集中程度。逐角曲线进一步使用5°等宽箱中的经验概率密度，按相同的负对数关系绘出F*(q)，高占据角度对应低F*，低值范围展开该旋转轴上的采样形态[32,33,35]。联合类别与逐角投影分别保留整个侧链的角度组合及单条旋转轴的分布。模型输出RMSF与Hχ，轨迹提供占据概率及对应的F*曲线。')
oldf=find('连续角度分布以5°等宽箱')
replace(oldf,'连续角度的经验密度由5°箱出现频率除以箱宽得到，F*(q)=−ln[p̂(q)/p̂max]，最低点为零，采用无量纲热能标度。所有角度箱保留原横坐标，零计数箱断线。前后半程与全部重复对的占据通过Jensen–Shannon散度（JSD）比较，另计算全窗熵与五个连续窗口熵的帧数加权平均。边界观测以角度到任一边界的圆周距离不超过10°定义，相邻保存帧变化按实际采样间隔计数。为连接主链幅度与侧链占据，按每帧联合状态对链内对齐Cα坐标分组，计算各状态内部均方波动和状态平均位置之间的均方波动，按状态占比加权后相加为RMSF²。状态间部分占比先逐位置计算，再对记录内有效位置及全部记录等权平均。')
representative_method=find('连续分布案例沿用泛素Q41/Q62')
replace(representative_method,representative_method.text+' 非Pro χ1的负向、正向及trans取向区沿用上述常规中心分区；主要角度占据区与5°箱的经验峰位置从第一条轨迹计算。')

# The main narrative now explicitly includes distribution and F* interpretation.
abstract=find('蛋白质的局部动态既包括')
replace(abstract,abstract.text.replace('83条记录的203条轨迹分析进一步将离散读数对应到连续角度分布。','83条记录的203条轨迹分析结合连续角度分布、经验投影自由能及主链波动分解，解释两类读数对应的构象占据。'))
intro=find('我们从整体预测、外部迁移和结构实例')
replace(intro,intro.text.replace('本文建立静态主链到两类互补动态指标的联合预测框架，并通过整体误差、类型内比较和结构实例说明其在局部动力学研究中的用途。','进一步通过连续侧链分布、经验投影自由能与状态条件下的主链波动分解，将预测读数对应到具体构象占据。本文建立静态主链到两类互补动态指标的联合预测框架，并从整体误差、类型内比较及构象分布解释其在局部动力学研究中的用途。'))
coupling=find('全队列的配对分析显示')
replace(coupling,coupling.text.replace('同类型控制后仍存在较大未共同变化的部分，说明主链波动幅度与侧链状态多样性提供了不同的残基信息。','同类型控制后的参考排序相关为0.391，两类读数有关联而并不一致；下述四分位组合进一步定位这一差异在具体残基上的表现。'))
disc=find('离散熵与连续分布是同一构象集合')
replace(disc,'构象占据分析为两个预测目标提供了共同的分布解释。RMSF读取骨架位置的散布，Hχ读取整个侧链多个角度联合出现的比例，经验投影自由能则展开特定角度的高占据取向及其范围。泛素中χ1边缘熵接近而完整联合Hχ不同，三个来源的实例也显示单个旋转轴可以呈现一个或多个主要占据区。状态条件分解进一步表明，在这些侧链类别中，主链均方波动主要来自类别内部的位置散布。这些结果将两类读数的互补性从预测误差推进到具体构象形态，使结构筛选能够形成关于局部堆积、取向约束和主链活动的可检验解释。')
before(find('KIX与泛素实例分别把主链'), '角区间比较明确了Hχ读取连续占据分布的方式。常见峰附近的摆动在原分区下也可以分到不同类别，中心分区则将同一峰附近的观测更集中地归为一类；移动边界改变离散熵，连续角度与经验F*保持相同。因此，联合状态读数适合与逐角分布共同分析：前者概括组合多样性，后者显示这种多样性在具体旋转轴上的来源。进一步的分布模型可以结合末端角的化学对称性与模拟条件，直接学习更丰富的占据信息。')
conclusion=find('FlexPPI从静态主链结构联合估计')
replace(conclusion,conclusion.text.replace('全库配对及连续角度分析进一步显示，两个读数概括了不同的局部构象信息，并将侧链状态熵与轨迹中的占据分布相联系。','全库配对、经验投影自由能及状态条件波动分解进一步解释了两个读数与构象占据的关系，将预测对应到局部高占据取向、角度组合与骨架位置散布。'))

# Retain detailed protocols in the supplement; no energy result requires it.
replace(find('K 角区间与采样分布'),'K 角度分布与采样分析的复现细节')
details=find('轨迹分析完整覆盖83条记录和203条轨迹')
replace(details,'正文2.8及3.6报告构象分布分析。第一条共同窗口轨迹覆盖38,599个RMSF和33,493个Hχ有效位置，原统一窗口标签逐位复现；常规χ1保留31,621个非Pro位置，1,872个Pro位置在完整联合状态中保留。圆周平移将θ−s映射到[−180°,180°)，s=0保留原端点处理，s=60对应常规χ1中心分区。完整83条记录、203条轨迹的逐位置统计及四种平移的来源分层源表随复现材料提供；角区间汇总见正文表3。')
sampling=find('完整联合状态的全窗熵高于')
replace(sampling,'完整联合状态的半程JSD在原分区和60°平移下为0.050和0.067 nats；60条ATLAS/ITrace记录的全部三对重复平均为0.049和0.069 nats。正文图8d另展示同一批60条记录的非Pro χ1边缘占据JSD。状态条件分解使用RMSF与Hχ共同有效的33,493个位置，各位置的状态内与状态间均方波动之和逐位等于总均方波动；保存帧数加权、均值汇总口径及数值核验均记录在源表和审计中。')
remove(find('L 连续角度分布代表案例'))
remove(find('三个来源的记录按原折外两项目标误差'))

doc.save(OUTPUT)
with zipfile.ZipFile(SOURCE) as a,zipfile.ZipFile(OUTPUT) as b:
    comment=a.read('word/comments.xml')
    assert etree.tostring(etree.fromstring(comment),method='c14n')==etree.tostring(etree.fromstring(b.read('word/comments.xml')),method='c14n')
    entries=[(info,comment if info.filename=='word/comments.xml' else b.read(info.filename)) for info in b.infolist()]
temp=OUTPUT.with_suffix('.preserved.tmp')
with zipfile.ZipFile(temp,'w') as package:
    for info,content in entries:package.writestr(info,content)
temp.replace(OUTPUT)
with zipfile.ZipFile(SOURCE) as a,zipfile.ZipFile(OUTPUT) as b:
    assert a.read('word/comments.xml')==b.read('word/comments.xml')
    assert all(a.read(n)==b.read(n) for n in a.namelist() if n.startswith('word/media/'))
    ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main','m':'http://schemas.openxmlformats.org/officeDocument/2006/math'}
    ax,bx=etree.fromstring(a.read('word/document.xml')),etree.fromstring(b.read('word/document.xml'))
    anc=lambda e:Counter((x.tag,x.get(qn('w:id'))) for x in e.iter() if x.tag in {qn('w:commentRangeStart'),qn('w:commentRangeEnd'),qn('w:commentReference')})
    eq=lambda e:[etree.tostring(x,method='c14n') for x in e.xpath('.//m:oMath',namespaces=ns)]
    assert anc(ax)==anc(bx) and eq(ax)==eq(bx)[:len(eq(ax))]
    assert len(eq(bx))==len(eq(ax))+1
old=Document(SOURCE)
cells=lambda t:tuple(tuple(c.text for c in r.cells) for r in t.rows)
assert len(doc.tables)==15 and Counter(cells(t) for t in doc.tables)==Counter(cells(t) for t in old.tables)
audit={'document':str(OUTPUT),'document_sha256':sha(OUTPUT),'source':str(SOURCE),'source_sha256':sha(SOURCE),'main_figures':9,'supplementary_figures':11,'main_tables':3,'supplementary_tables':12,'all_table_cell_values_preserved':True,'scientific_media_parts_preserved':True,'comments_and_anchors_preserved':True,'original_native_math_preserved':True,'added_native_math':'Equation5 exact entropy/relative empirical free-energy identity','promotions':{'figureS12':'figure8','figureS13':'figure9','tableS13':'table3'},'no_new_MD_or_training':True,'final_review':'pending'}
(QA/'integration_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf8')
out=[]
for i,p in enumerate(doc.paragraphs):
    t=''.join(p._p.xpath('.//w:t/text() | .//m:t/text()'))
    if t:out.append(f'P{i}: {t}')
for i,t in enumerate(doc.tables):
    for j,row in enumerate(t.rows):out.append(f'T{i}R{j}: '+' | '.join(c.text for c in row.cells))
(QA/'v22_content.txt').write_text('\n'.join(out)+'\n',encoding='utf8')
print(json.dumps(audit,ensure_ascii=False))
