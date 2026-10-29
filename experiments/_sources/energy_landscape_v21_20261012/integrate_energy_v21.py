"""Integrate completed conformational-distribution experiments into v20."""
from pathlib import Path
from copy import deepcopy
from collections import Counter
import hashlib,json,re,zipfile
from lxml import etree
import pandas as pd
from docx import Document
from docx.shared import Inches,Pt
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn

ROOT=Path(__file__).resolve().parent
BASE=ROOT.parent.parent
SOURCE=BASE/'paper_v20_20261012_中文修订_联合模型_三线表.docx'
OUTPUT=BASE/'paper_v21_20261012_中文修订_构象分布与自由能分析_三线表.docx'
QA=ROOT/'manuscript_QA'
QA.mkdir(exist_ok=True)
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
assert sha(SOURCE)=='4b4597b42fdc8a1300ab8f4155e98a4960872a2368fc887397d133af7d264d79'
doc=Document(SOURCE)
sources={}
def csv(path):
    p=ROOT/path;sources[str(p)]=sha(p)
    return pd.read_csv(p)
def read(path):
    p=ROOT/path;sources[str(p)]=sha(p)
    return p.read_text(encoding='utf8').strip()
def find(start):
    p=[p for p in doc.paragraphs if p.text.startswith(start)]
    assert len(p)==1,(start,len(p))
    return p[0]
def replace(p,text):
    assert not p._p.xpath('.//m:oMath')
    nodes=p._p.xpath('.//w:t');assert nodes
    nodes[0].text=text;nodes[0].set('{http://www.w3.org/XML/1998/namespace}space','preserve')
    for node in nodes[1:]:node.text=''
def before(anchor,text,style='Normal'):
    p=doc.add_paragraph(text,style);anchor._p.addprevious(p._p)
    return p
def after(anchor,text,style='Normal'):
    p=doc.add_paragraph(text,style);anchor._p.addnext(p._p)
    return p
def figure_before(anchor,name,caption,width=6.3,newpage=False):
    path=ROOT/'figures'/name;sources[str(path)]=sha(path)
    p=doc.add_paragraph();p.alignment=WD_ALIGN_PARAGRAPH.CENTER
    p.paragraph_format.keep_with_next=True
    p.paragraph_format.page_break_before=newpage
    p.add_run().add_picture(str(path),width=Inches(width))
    anchor._p.addprevious(p._p)
    p=before(anchor,caption,'Caption');p.paragraph_format.keep_with_next=False
    return p
def heading(text):
    p=doc.add_paragraph(text,'Heading 1');p.paragraph_format.page_break_before=True
    return p
def table(caption,columns,rows,widths):
    p=doc.add_paragraph(caption,'Caption');p.paragraph_format.keep_with_next=True
    t=doc.add_table(rows=1,cols=len(columns));t.autofit=False
    for col,w in zip(t.columns,widths):col.width=Inches(w)
    for cell,text in zip(t.rows[0].cells,columns):cell.text=str(text)
    for row in rows:
        for cell,text in zip(t.add_row().cells,row):cell.text=str(text)
    for i,row in enumerate(t.rows):
        row._tr.get_or_add_trPr().append(OxmlElement('w:cantSplit'))
        if i==0:row._tr.get_or_add_trPr().append(OxmlElement('w:tblHeader'))
        for cell,width in zip(row.cells,widths):
            cell.width=Inches(width)
            for p in cell.paragraphs:
                p.alignment=WD_ALIGN_PARAGRAPH.CENTER
                p.paragraph_format.space_before=Pt(4);p.paragraph_format.space_after=Pt(4)
                p.paragraph_format.first_line_indent=Pt(0)
                p.paragraph_format.left_indent=Pt(0);p.paragraph_format.right_indent=Pt(0)
                for r in p.runs:
                    r.font.name='Times New Roman';r.font.size=Pt(10);r.font.bold=i==0
                    r._element.get_or_add_rPr().rFonts.set(qn('w:eastAsia'),'宋体')
    return t

# Renumber only existing main figures; supplementary labels are independent.
for p in doc.paragraphs:
    if p._p.xpath('.//m:oMath'):continue
    text=re.sub(r'图\s*([45])(?=\D|$)',lambda m:'图'+{'4':'5','5':'6'}[m.group(1)],p.text)
    # Supplementary tables follow their first appearance in the final document.
    text=re.sub(r'表\s*S([123])(?=\D|$)',lambda m:'表S'+{'1':'2','2':'3','3':'1'}[m.group(1)],text)
    text=text.replace('3.3 ThermoPCD外部验证','3.4 ThermoPCD外部验证').replace('3.4 将联合动态预测对应到具体结构位置','3.5 将联合动态预测对应到具体结构位置')
    if text!=p.text:replace(p,text)

abstract=find('蛋白质的局部动态既包括')
replace(abstract,abstract.text.replace('ThermoPCD外部验证进一步考察两类输出在另一轨迹来源中的表现。',
    '全体折外预测的配对分析显示，同类型残基中存在主链波动与侧链状态多样性不同步的局部组合；83条记录的203条轨迹分析进一步将离散读数对应到连续角度分布。ThermoPCD外部验证考察两类输出在另一轨迹来源中的表现。'))
p=find('侧链动态的学习还需要区分')
after(p,'蛋白质动态来自其采样的构象分布。长时间模拟表明，紧密堆积的蛋白核心也可访问多个侧链旋转状态，状态占比与转换时间分别描述构象多样性和运动时间尺度[32,33]。增强采样研究进一步联合分析侧链占据、熵、概率流及局部自由能面[34]。RMSF概括主链位置分布的宽度，Hχ概括多个侧链角共同出现的状态及其占比。本文在联合预测之外，结合连续角度分布和经验投影自由能，分析两类读数对应的构象形态。')

method=find('3 结果')
before(method,'2.8 构象分布与双读数分析','Heading 2')
before(method,'两类读数的配对分析使用全部3,063条完整模型折外预测，RMSF和Hχ共同有效位置为824,988个。先在每条记录内计算两个指标的Spearman相关，再对记录等权平均；同类型控制版本分别减去该记录内相同氨基酸的指标均值。局部组合按同一记录、同一种氨基酸中各指标的第25与第75百分位定义低和高，并列值全部保留。至少具有四个有效位置且两个分位不同的类型参与四类组合分类，其余位置保留为未分类。参考与预测分别使用各自的分位，召回率为恢复相同组合的交集位置数除以参考组位置数。')
before(method,'轨迹分析沿用补充材料E的83条冻结记录，ATLAS、DynaRepo和ITrace分别为30、23和30条，读取全部203条可用轨迹。统一取0–100 ns、200 ps间隔的501个实际保存帧，主分析每记录使用第一条轨迹，其余重复用于占据一致性比较。提取完整χ角并保留化学名称校正与残基映射，在原等宽联合分区之外，另将圆周分区整体平移30°、60°和90°。非Pro残基的χ1单独使用以−60°、+60°和180°为中心的常规三区间，边界为−120°、0°和+120°，接缝两侧归为同一trans状态[33]。全部χ角的统一平移用于分区敏感性；Pro及原联合标签保持原定义。')
before(method,'连续角度分布以5°等宽箱计算经验概率p̂j，相对经验投影自由能定义为F*j=−ln(p̂j/max p̂j)，最低点为零，以无量纲热能标度展示。零计数箱保留为缺失，曲线沿实际访问箱绘出。前后半程及重复轨迹的占据分布通过Jensen–Shannon散度（JSD）比较；将完整轨迹分为五个连续窗口，另比较全窗熵与按帧数加权的窗口熵。相邻保存帧的类别变化按实际采样间隔计数。F*由轨迹的所选坐标计算，RMSF与Hχ由模型预测；投影坐标与扩散行为共同影响对运动过程的描述[35]。')
before(method,'连续分布案例沿用泛素Q41/Q62的前5,000帧及五个标签窗口。另在83条子集内，每个来源选取原折外两项MAE的来源内百分位秩均值最接近0.5的记录，获得2py5_A、A0022及ITrace 6p64。在每个记录的同一种氨基酸内，选取参考Hχ秩减RMSF秩最相离的一对非Pro位点，用于展示两类运动的不同组合；该位点选择使用参考读数，完整残基曲线与三条记录的原始选择结果同时保留。')

external=find('3.4 ThermoPCD外部验证')
before(external,'3.3 主链与侧链读数的互补组合','Heading 2')
before(external,'全队列的配对分析显示，两种动态指标呈相关而不同的局部变化。原值的参考与预测Spearman相关分别为0.233和0.260；在每条记录内减去同类型均值后，分别为0.391和0.517（图4a）。同类型控制后仍存在较大未共同变化的部分，说明主链波动幅度与侧链状态多样性提供了不同的残基信息。')
before(external,'以同记录同类型的四分位区间考察局部组合，24,842个位点属于低RMSF／高Hχ，22,512个位点属于高RMSF／低Hχ，合计47,354个位置。两类组合的记录等权比例为6.105%，对应全部共同有效位置的5.740%。原始及预测类别占比见图4b，中间及不能形成分位区间的位置保留为未分类。模型对低／低和高／高组合的记录平均召回率为54.6%和61.4%，对低／高和高／低组合为15.8%和16.5%（图4c，表S12）。因此，联合查看可以发现单一读数遗漏的组合，同时组合恢复比单项排序提出了更具体的位点筛选要求。')
figure_before(external,'fig4_joint_complementarity.png',read('figures/fig4_joint_complementarity_caption_zh.md'))

discussion=find('4 讨论')
before(discussion,'3.6 从柔性读数到连续构象分布','Heading 2')
before(discussion,'将泛素侧链读数对应到连续轨迹，可以进一步说明联合熵的来源。Q41与Q62的χ1分布均主要集中在−60°附近，原分区的χ1边缘熵分别为0.837和0.833 nats，基本相近；χ2边缘熵为0.680和0.835 nats，χ3为0.574和1.032 nats。Q62较高的完整联合Hχ因此对应多个角度的组合覆盖，而单独χ1投影没有呈现同样的高低差异（图7）。概率密度与F*显示，Q41的χ2主要围绕圆周接缝180°采样，Q62则主要围绕−60°采样；Q62的χ3覆盖也更分散。这些分布给出两位点在局部空间约束下不同的取向形态。')
before(discussion,'分区比较将连续角度与离散状态的关系具体化。在83条记录的统一100 ns窗口中，完整联合Hχ的记录均值在原分区为1.484 nats，在整体平移60°后为0.962 nats；30°和90°结果分别为1.021和1.026 nats（表S13，图S12）。非Pro的χ1边缘熵由0.717变为0.374 nats，原边界10°邻域覆盖61.3%的角度观测，常规χ1中心分区为0.13%；相邻保存帧的平均类别变化由189.2次变为21.9次。原分区经常把围绕同一角度峰的摆动分到相邻箱中，因此Hχ同时描述了原区间下的局部摆动和更广范围的取向覆盖。')
before(discussion,'泛素的比较保留了这一结构差异：Q41与Q62在60°联合分区平移后的五窗口平均熵为0.528和1.399 nats，Q62仍具有较高的多角度状态多样性，而常规χ1边缘熵为0.492和0.403 nats。χ1、各角边缘分布与完整联合状态分别保留不同信息。进一步按预设中位误差规则展示的三个来源案例，将完整残基读数、具体同类型位点的角度时间序列与F*放在同一行比较（图S13），使主链活动、侧链取向覆盖和采样时间组织能够在具体蛋白中联合观察。')
figure_before(discussion,'fig7_ubiquitin_angular_landscapes.png',read('figures/fig7_ubiquitin_angular_landscapes_caption_zh.md'),width=6.1,newpage=True)

p=find('同类型比较使侧链输出')
replace(p,p.text+' 全库组合分析和连续角度投影进一步将这种互补性量化：相同类型的主链与侧链变化并不完全同步，多χ角的联合状态也可以区分χ1分布相近的位点。联合读数适合用作第一步结构筛选，选定位点的连续分布再由轨迹具体分析。')
p=find('标签由轨迹采样和处理方式共同确定')
replace(p,'离散熵与连续分布是同一构象集合的不同摘要。分区敏感性分析显示，边界穿过常规χ1峰附近时，局部摆动也会改变箱占据和保存帧类别；常规中心分区将同一峰附近的观测更集中地归入一个区间。本文的联合目标对应明确定义的角区间，连续密度及F*补充了区间之间和圆周接缝附近的形态。后续可针对不同末端χ角的化学对称性设计专门状态，并将更明确的条件输入与侧链分布目标结合。')
p=find('FlexPPI从静态主链结构联合估计')
replace(p,p.text.replace('KIX实例进一步将预测对应到活动连接环。','KIX实例将预测对应到活动连接环；全库配对及连续角度分析进一步显示，两个读数概括了不同的局部构象信息，并将侧链状态熵与轨迹中的占据分布相联系。'))

# Four new primary references inserted before the supplement.
supp=find('补充材料')
references=[
'[32] Gregory R. Bowman and Phillip L. Geissler. Extensive Conformational Heterogeneity within Protein Cores. The Journal of Physical Chemistry B, 118(24):6417–6423, 2014. doi:10.1021/jp4105823.',
'[33] Alexander D. Scouras and Valerie Daggett. The Dynameomics rotamer library: amino acid side chain conformations and dynamics from comprehensive molecular dynamics simulations in water. Protein Science, 20(2):341–352, 2011. doi:10.1002/pro.565.',
'[34] Rajitha Rajeshwar T., Jeremy C. Smith, and Marimuthu Krishnan. Hidden Regularity and Universal Classification of Fast Side Chain Motions in Proteins. Journal of the American Chemical Society, 136(24):8590–8605, 2014. doi:10.1021/ja5024783.',
'[35] Robert B. Best and Gerhard Hummer. Coordinate-dependent diffusion in protein folding. Proceedings of the National Academy of Sciences of the United States of America, 107(3):1088–1093, 2010. doi:10.1073/pnas.0910390107.']
sample=find('[31] Samuel')
for text in references:
    p=before(supp,text,sample.style)
    if sample._p.pPr is not None:
        if p._p.pPr is not None:p._p.remove(p._p.pPr)
        p._p.insert(0,deepcopy(sample._p.pPr))
    for r in p.runs:
        if sample.runs[0]._r.rPr is not None:r._r.insert(0,deepcopy(sample.runs[0]._r.rPr))

heading('J 全队列双读数组合')
doc.add_paragraph('全部3,063条记录纳入配对分析，RMSF与Hχ共同有效位置为824,988个。减去同记录同类型均值后，参考与预测的Spearman相关分别为0.391和0.517；按ATLAS、DynaRepo和ITrace分别为0.397/0.534、0.388/0.514和0.349/0.390。四类组合使用每记录内同类型的25%及75%分位，比例以该记录全部共同有效位置为分母；各类误差和召回率对具有该类位置的记录等权汇总。未分类位置包括中间区间、同类型不足四个位点及分位重合，保留在原始计数与比例中。')
combo=csv('oof_analysis/reference_class_record_macro.csv')
combo=combo[combo.dataset=='all'].set_index('reference_class')
names=[('low_RMSF_low_Hchi','低RMSF／低Hχ'),('low_RMSF_high_Hchi','低RMSF／高Hχ'),
       ('high_RMSF_low_Hchi','高RMSF／低Hχ'),('high_RMSF_high_Hchi','高RMSF／高Hχ'),('unclassified','未分类')]
rows=[]
for key,name in names:
    r=combo.loc[key]
    rows.append([name,f'{int(r.common_residues):,}',f'{100*r.reference_fraction_record_macro:.3f}',
                 f'{r.rmsf_mae_record_macro:.3f}',f'{r.hchi_mae_record_macro:.3f}',
                 f'{100*r.rank_class_recall_record_macro:.1f}' if key!='unclassified' else '无'])
table('表S12：同类型相对双读数组合的参考位置数、记录等权比例及预测表现。比例和召回率单位为%；RMSF MAE为Å，Hχ MAE为nats。各类误差和召回率使用存在该类的记录；未分类为分位定义之外的位置。',
      ['参考组合','位置数','比例','RMSF MAE','Hχ MAE','召回率'],rows,[1.5,1.0,.8,1.0,1.0,1.0])

heading('K 角区间与采样分布')
doc.add_paragraph('轨迹分析完整覆盖83条记录和203条轨迹。各记录第一条轨迹的共同有效位置为38,599个RMSF及33,493个Hχ位置，原统一窗口标签逐位复现。常规χ1分析保留31,621个非Pro位置，1,872个Pro位置在完整联合分区中保留。圆周平移先将θ−s映射到[−180°,180°)，再按120°区间编码；s=0保留原端点处理，s=60对应常规χ1中心分区。模型和主训练标签保持冻结，平移分区用于同一轨迹的读数比较。')
partition=csv('trajectory_analysis/manuscript_summary_table.csv')
partition=partition[partition.dataset=='all'].set_index('shift_degrees')
chi=csv('trajectory_analysis/nonpro_chi1_source_summary.csv')
chi=chi[(chi.dataset=='all')&(chi.replica_mode=='replica1')].set_index('shift_degrees')
rows=[]
for phase in [0,30,60,90]:
    r=partition.loc[phase];q=chi.loc[phase]
    rows.append([f'{phase}°',f'{r.joint_H_mean_nats:.3f} ± {r.joint_H_record_sd_nats:.3f}',
                 f'{q.h_chi1_mean:.3f}',f'{100*q.chi1_boundary_10deg_fraction_mean:.2f}',
                 f'{q.apparent_chi1_changes_mean:.1f}',f'{r.half_joint_JSD_mean_nats:.3f}'])
table('表S13：角区间平移的记录等权统计。主分析每记录第一条100 ns轨迹、501帧；联合Hχ保留全部原定义角和Pro。χ1熵、边界10°观测比例和保存帧变化数均采用非Pro位置；半程JSD采用完整联合状态，熵与JSD单位为nats。均值 ± SD为83条记录的均值与样本标准差。',
      ['平移','联合Hχ','非Pro χ1熵','边界比例%','变化次数','半程JSD'],rows,[.7,1.4,1.0,1.1,1.0,1.1])
doc.add_paragraph('完整联合状态的全窗熵高于五个约20 ns窗口熵的帧数加权平均，原分区与60°平移的记录平均差值为0.135和0.170 nats，半程JSD为0.050和0.067 nats。60条ATLAS/ITrace记录具有三条重复，全体三对重复的联合状态JSD为0.049和0.069 nats。以联合状态为条件，将链内对齐Cα的总均方波动分解为状态内波动和状态均值之间的波动，后者平均占6.2%和7.8%；两部分相加逐位等于总波动。分段、重复及状态条件的比较从不同方面描述了同一轨迹窗口中的占据和位置变化。')
# Figure pages deliberately keep a large scientific plate paired with its caption.
p=doc.add_paragraph();p.paragraph_format.page_break_before=True;p.paragraph_format.keep_with_next=True
p.alignment=WD_ALIGN_PARAGRAPH.CENTER
p.add_run().add_picture(str(ROOT/'figures/figS12_partition_sampling.png'),width=Inches(6.0))
sources[str(ROOT/'figures/figS12_partition_sampling.png')]=sha(ROOT/'figures/figS12_partition_sampling.png')
doc.add_paragraph(read('figures/figS12_partition_sampling.caption_zh.md').replace('图 S12','图S12'),'Caption')

heading('L 连续角度分布代表案例')
doc.add_paragraph('三个来源的记录按原折外两项目标误差的平均来源内秩接近中位数选取，2py5_A、A0022和ITrace 6p64的均值秩均为0.5。参考与预测的完整逐残基曲线和局部角度分布见图S13。选取的同类型位点对分别为Y254/Y521、A链I78/I21、E链T149/D链T25；这些编号对应模型输入PDB。例如2py5_A的Y254/Y521参考RMSF为0.630/1.850 Å、Hχ为1.573/0.358 nats，A0022的I78/I21为0.343/0.464 Å、1.694/0.794 nats，ITrace所选位点为0.724/1.735 Å、1.046/0.514 nats。它们在同类型内表现为两读数不同步的组合，原始轨迹中的χ1时间序列及F*用于具体显示这些位点的连续分布。')
p=doc.add_paragraph();p.paragraph_format.page_break_before=True;p.paragraph_format.keep_with_next=True
p.alignment=WD_ALIGN_PARAGRAPH.CENTER
p.add_run().add_picture(str(ROOT/'figures/figS13_representative_distribution_cases.png'),width=Inches(6.1))
sources[str(ROOT/'figures/figS13_representative_distribution_cases.png')]=sha(ROOT/'figures/figS13_representative_distribution_cases.png')
doc.add_paragraph(read('figures/figS13_representative_distribution_cases_caption_zh.md'),'Caption')

# Apply the user's three-line-table requirement to the two added native tables.
for t in doc.tables:
    for e in list(t._tbl.tblPr):
        if e.tag in {qn('w:tblBorders'),qn('w:shd')}:t._tbl.tblPr.remove(e)
    for i,row in enumerate(t.rows):
        for cell in row.cells:
            cp=cell._tc.get_or_add_tcPr()
            for e in list(cp):
                if e.tag in {qn('w:tcBorders'),qn('w:shd')}:cp.remove(e)
            borders=OxmlElement('w:tcBorders')
            for side in ['top','left','bottom','right','insideH','insideV']:
                outer=(i==0 and side=='top') or (i==len(t.rows)-1 and side=='bottom')
                head=i==0 and side=='bottom'
                e=OxmlElement('w:'+side);e.set(qn('w:val'),'single' if outer or head else 'nil')
                if outer or head:e.set(qn('w:sz'),str(12 if outer else 6));e.set(qn('w:color'),'000000');e.set(qn('w:space'),'0')
                borders.append(e)
            cp.append(borders)
doc.save(OUTPUT)
with zipfile.ZipFile(SOURCE) as old,zipfile.ZipFile(OUTPUT) as new:
    comment=old.read('word/comments.xml')
    assert etree.tostring(etree.fromstring(comment),method='c14n')==etree.tostring(etree.fromstring(new.read('word/comments.xml')),method='c14n')
    entries=[(info,comment if info.filename=='word/comments.xml' else new.read(info.filename)) for info in new.infolist()]
temp=OUTPUT.with_suffix('.preserved.tmp')
with zipfile.ZipFile(temp,'w') as package:
    for info,content in entries:package.writestr(info,content)
temp.replace(OUTPUT)
ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main','m':'http://schemas.openxmlformats.org/officeDocument/2006/math'}
with zipfile.ZipFile(SOURCE) as old,zipfile.ZipFile(OUTPUT) as new:
    assert old.read('word/comments.xml')==new.read('word/comments.xml')
    assert all(old.read(n)==new.read(n) for n in old.namelist() if n.startswith('word/media/') and not n.endswith('/'))
    a,b=[etree.fromstring(z.read('word/document.xml')) for z in [old,new]]
    anchors=lambda e:Counter((x.tag,x.get(qn('w:id'))) for x in e.iter() if x.tag in {qn('w:commentRangeStart'),qn('w:commentRangeEnd'),qn('w:commentReference')})
    equations=lambda e:[etree.tostring(x,method='c14n') for x in e.xpath('.//m:oMath',namespaces=ns)]
    assert anchors(a)==anchors(b) and equations(a)==equations(b)
assert len(doc.tables)==15
old=Document(SOURCE)
assert all([[c.text for c in r.cells] for r in doc.tables[i].rows]==[[c.text for c in r.cells] for r in old.tables[i].rows] for i in range(13))
audit={'document':str(OUTPUT),'sha256':sha(OUTPUT),'source_document_sha256':sha(SOURCE),'tables':15,
       'new_figures':['4','7','S12','S13'],'references':35,'scientific_sources_sha256':sources,
       'original_tables_preserved':True,'comments_and_anchors_preserved':True,'native_equations_preserved':True,
       'actual_existing_new_trajectory_experiments':203,'new_MD_production':0,'new_model_training':0,
       'render_and_content_review':'pending'}
(QA/'integration_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf8')
lines=[]
for i,p in enumerate(doc.paragraphs):
    text=''.join(p._p.xpath('.//w:t/text() | .//m:t/text()'))
    if text:lines.append(f'P{i}: {text}')
for i,t in enumerate(doc.tables):
    for j,row in enumerate(t.rows):lines.append(f'T{i}R{j}: '+' | '.join(c.text for c in row.cells))
(QA/'v21_content.txt').write_text('\n'.join(lines)+'\n',encoding='utf8')
print(json.dumps({'document':str(OUTPUT),'sha256':sha(OUTPUT),'tables':15},ensure_ascii=False))
