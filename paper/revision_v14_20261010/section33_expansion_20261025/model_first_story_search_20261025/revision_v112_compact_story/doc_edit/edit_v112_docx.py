"""Targeted OOXML update from v98 to v112, preserving unrelated ZIP parts.
Run with the bundled Codex Python. No Word/LibreOffice mutation is performed.
"""
from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED
from copy import deepcopy
from lxml import etree
from PIL import Image
import json, hashlib, re

BASE=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010')
HERE=Path(__file__).resolve().parent
SOURCE=BASE/'paper_v98_20261027_中文修订_结构区域放大.docx'
TARGET=BASE/'paper_v112_20261027_中文修订_双尺度柔性与角度占据整合.docx'
FIGURE=BASE/'figure3_v112_20261027_紧凑布局与等残基宽度.png'
NS={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main',
    'a':'http://schemas.openxmlformats.org/drawingml/2006/main',
    'wp':'http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing',
    'r':'http://schemas.openxmlformats.org/officeDocument/2006/relationships'}
W='{'+NS['w']+'}'
def text(p):return ''.join(p.xpath('.//w:t/text()',namespaces=NS))
def sha(data):return hashlib.sha256(data).hexdigest()
assert FIGURE.is_file(),FIGURE
source_sha=sha(SOURCE.read_bytes())
with ZipFile(SOURCE) as z:
    items=z.infolist(); payload={i.filename:z.read(i.filename) for i in items}
root=etree.fromstring(payload['word/document.xml'])
original=deepcopy(root)
paras=root.xpath('//w:body/w:p',namespaces=NS)
assert text(paras[68]).startswith('3.3 ')
assert text(paras[77]).startswith('3.4 ')
assert text(paras[81]).startswith('3.5 TCR')
assert len(paras[72].xpath('.//w:drawing',namespaces=NS))==1
changes=[]
def replace(index,new):
    p=paras[index]
    old=text(p)
    assert not p.xpath('.//w:fldChar|.//w:fldSimple|.//w:drawing|.//m:oMath',namespaces={**NS,'m':'http://schemas.openxmlformats.org/officeDocument/2006/math'})
    # Retain paragraph properties and the original text run's typography.
    runs=p.findall(W+'r')
    rpr=next((deepcopy(r.find(W+'rPr')) for r in runs if r.find(W+'rPr') is not None),None)
    for child in list(p):
        if child.tag!=W+'pPr':p.remove(child)
    r=etree.SubElement(p,W+'r')
    if rpr is not None:r.append(rpr)
    t=etree.SubElement(r,W+'t');t.text=new
    changes.append({'source_paragraph_index':index,'old':old,'new':new})

replace(67,'这些层面的贡献将在后续分析中落实到具体结构位置。3.3节结合整体RMSF–Hχ分布、CaM–CaMKI与泛素–UBA两个实例及MD角度占据，依次检验双尺度柔性预测的准确性、结构环境对应及侧链构象含义；3.4节以MEL5亲和成熟和肽替换为例，将两类读数进一步联系到跨界面接触网络。')
replace(68,'3.3 主链与侧链柔性预测及角度占据的结构解释')
replace(69,'主链RMSF描述对齐后残基位置的波动，侧链联合Hχ概括固定χ角区间下取向状态的占据多样性。为说明两类读数的关系，首先汇总3,063条记录中824,988个共同有效位点：MD参考中RMSF与Hχ的Spearman相关为0.209，完整模型预测中为0.262（图3b）。在相近RMSF范围内仍可观察到较宽的Hχ分布，预测也保留了这一总体特征。这种分布对应说明主链运动与侧链状态多样性有一定关联，但不是一一对应。总体分布受氨基酸类型及侧链角数影响，因此进一步在两个具体结构实例中核对完整有效位置上的预测，再检验连续识别核心与同类型侧链局部环境的对应。')
replace(70,'CaM–CaMKI复合物提供连续识别片段的主链约束实例。CaMKI的调控序列与自抑制区域重叠，CaM结合需要重新组织该片段与激酶的相互作用；1MXE中CaM两个叶瓣包绕294–318号靶肽，W303和M316分别结合C端及N端叶瓣，形成1–14识别模式[29,30]（图3a）。全组成链家族划分第1折的留出预测在靶肽全部25个位置上的RMSF MAE为0.115 Å，Spearman相关为0.803；21个有效侧链位置的Hχ MAE为0.254 nats（图3c）。图中同时保留两项绝对读数及其MD参考，使预测精度与区域差异在相同残基位置上对应。')
replace(71,'结构定义的W303–M316核心具有连续较低的主链波动，MD平均RMSF为0.553 Å，模型为0.513 Å；其余11个位点分别为1.089和0.919 Å。按4.5 Å重原子距离阈值，W303与M316分别邻接9个和6个CaM残基，模型恢复的约束跨越了两个锚点之间的整个片段，而不只局限于锚点本身。下方对齐的侧链溶剂可及面积（SASA）补充展示各位置的静态暴露环境（图3c）；该面积不等同动力学柔性，不能将曲线与柱高的局部对应直接解释为因果。连续低RMSF核心与CaM双叶瓣稳定调控片段的结构背景一致，可为调控构型的维持提供动态线索，但1MXE不含激酶催化域，不能直接检验片段回接触或酶激活。完整与链内模型的比较及全部参考核查保留在图S3；二者使用相同结合态结构分别训练，其差异不是游离与结合状态之间的实验变化。')
replace(73,'图3：主链与侧链柔性预测的结构对应及角度占据解释。（a）CaM–CaMKI（1MXE）与泛素–UBA（2QHO）输入结构；紫色为分析链，灰色为伙伴，青绿色标记W303/M316与Q41/Q62。（b）同824,988个共同有效位点的MD与完整模型RMSF–联合Hχ分布及边缘频率；MD为灰色，预测为蓝色，透明叠加，无等密度轮廓。主视野为0–10 Å，MD与预测分别有4,234和906个位点位于右侧范围外；图示相关按全部位点计算，分别为0.209和0.262。（c、d）CaMKI294–318与泛素1–72的RMSF（蓝色，左轴Å）和Hχ（赭色，右轴nats）；预测为实线，MD为虚线；青绿色位置标记对应W303/M316与Q41/Q62。两例纵轴范围一致，残基横向间距相同。下方为共用残基横轴的起始复合物侧链SASA（青灰柱，Å²），不是MD时间均值或结合前后变化；Gly无侧链所产生的零值不表示深埋。RMSF/Hχ有效位置分别为25/21和72/66，MAE分别为0.115 Å/0.254 nats和0.040 Å/0.080 nats；无定义位置断线，合法零值保留。（e）同一联合Hχ预测与MD逐χ角5°细箱熵的关联，两轴分别去除氨基酸类型均值；Q41为菱形，Q62为方形，其余位置为圆点。χ1–χ4的全部有效点数为66、52、23和10，相关分别为0.50、0.45、0.52和0.13；四图不代表四个独立预测输出。模型不生成构象或直接预测结合自由能。')
replace(74,'泛素–UBA复合物2QHO用于检验侧链环境差异是否能在低主链波动的蛋白骨架中被区分[23]。使用原体系划分第0折的留出预测，模型在泛素全部72个主链位置上的RMSF MAE为0.040 Å，在全部66个有效侧链位置上的Hχ MAE为0.080 nats，Hχ的Spearman相关为0.979（图3d）。原始掩码中无χ角的位置留空，主链参考与预测同时保留。为避免把侧链种类和χ角数的差异误认为环境识别，进一步比较全部六个谷氨酰胺，其Hχ MAE为0.105 nats，类型内排序相关为0.943。相同化学类型仍具有不同的Hχ，表明本例的预测差异超出了残基类型平均值所能提供的信息。')
replace(75,'Q41与Q62给出了这一环境差异的具体结构对照（图3a、d）。Q41位于紧密堆积的内部环境，起始结构侧链SASA为0.0 Å²，MD和预测Hχ分别为1.890和2.000 nats；Q62更暴露于溶剂，SASA为100.67 Å²，对应Hχ为2.631和2.684 nats。模型恢复的两位置差为0.684 nats，接近MD参考的0.741 nats。这一结果与局部堆积和溶剂暴露所伴随的侧链取向限制一致，说明模型能够区分同类型侧链在不同环境中的状态多样性。全部六个谷氨酰胺并不随SASA单调排序，暴露面积因而不是Hχ的单因素解释；Q41与Q62均无5 Å以内的UBA直接接触，二者反映的是泛素自身的局部环境，不能据此认定存在位点上的直接伙伴效应。联合状态占据与全部位置的补充核查见图S2。')
replace(76,'为将标量Hχ进一步联系到实际角度占据，我们把同一个联合Hχ预测与MD逐χ角的5°细箱熵对应，并分别去除两轴的氨基酸类型均值（图3e）。χ1–χ3呈中等正相关，ρ分别为0.50、0.45和0.52，包含66、52和23个有效位置；χ4仅有10个位置，相关较弱（ρ=0.13）。图中继续标出Q41与Q62，使结构环境对照与逐角占据对应到相同残基。四幅图共享同一个联合Hχ预测，不能解释为模型分别预测了四个χ角；逐角细箱熵也不同于用于监督的粗箱多角联合Hχ。总体关联支持将模型输出解释为侧链占据多样性的摘要，但不意味着每个χ角、每个位置的占据细节都被恢复。')
replace(78,'这一角度占据解释进一步在83条留出记录中进行核查，包含ATLAS、DynaRepo和ITrace的30、23和30条记录，各记录统一取第一条轨迹0–100 ns、200 ps间隔的501帧（补充材料E）。控制记录内氨基酸类型后，完整模型Hχ与逐角细箱熵的平均Spearman相关为0.510，高于其RMSF的0.321；与经验投影F*=−ln[p̂(q)/p̂max]下F*≤1角区总宽度的相关分别为0.400和0.239（图S4b）[3,24,25]。仅按Hχ预测将同类型位置分为低、高两半后，实际细箱熵均值由2.498增至2.791 nats/角，角区宽度由35.41增至43.17°/角（图S4c）。这些结果支持侧链输出定位角度占据更丰富的位置；F*角区反映同一轨迹中的相对占据范围，不是整体蛋白构象熵或结合自由能，也不是独立实验验证。逐角MD占据曲线与方法敏感性保留在图S4及补充材料E。')
for index in (77,79,80):
    p=paras[index];changes.append({'source_paragraph_index':index,'old':text(p),'new':None});p.getparent().remove(p)
replace(81,'3.4 TCR–pMHC识别中的界面柔性与亲和力差异')
# Retain the main conclusion while making the statistical and mechanistic limits explicit.
replace(91,'侧链输出进一步与角度占据多样性建立统计对应。Hχ概括固定χ角区间下多角联合状态的占据多样性，同类型控制使比较建立在相同侧链自由度上。泛素逐角细箱熵及83条记录的配对相关与预测分组，支持侧链输出定位实际取向范围更广、角度占据更多样的位置。经验投影F*=−ln(p̂/p̂max)以所分析轨迹中占据最高的角区为零，描述该角度坐标下的相对占据偏好[3,24,25]；F*≤1角区覆盖与细箱熵分别刻画主占据范围和整体占据分布。这些MD导出指标不构成独立实验验证，模型仅预测RMSF与Hχ标量。上述对应可为后续局部模拟、接触分析和功能实验确定关注的残基位置。')
replace(158,text(paras[158]).replace('图3a、b中的1MXE预测','图3a、c中的1MXE预测'))
replace(167,text(paras[167]).replace('逐角比较在泛素六个谷氨酰胺中按留出模型Hχ预测选取最低和最高两个位置；','图3e高亮结构环境对照的Q41与Q62；图S4d另在泛素六个谷氨酰胺中按留出模型Hχ预测选取最低的Q41和最高的Q2作逐角占据核查，两种选点目的不同；'))
replace(168,text(paras[168]).replace('泛素侧链的起始结构溶剂可及面积采用','图3c、d两例侧链的起始结构溶剂可及面积统一采用').replace('以完整复合物的935个重原子为环境，','以各自完整复合物为环境（泛素–UBA含935个重原子），').replace('作为预测后的结构环境解释量。','作为预测后的静态结构环境解释量。图3c的CaMKI全部25位和图3d的泛素全部72位使用相同方法，保留全部97个位置而不受Hχ有效掩码限制；Gly无侧链重原子时侧链SASA为零，不表示深埋。图3c的SASA为保留CaM的完整复合物面积，不与补充链移除对照的遮蔽比例混用。'))
replace(171,text(paras[171]).replace('完整66个Hχ位点预测见图3c','完整66个Hχ位点预测见图3d').replace('主文的Q41/Q62用于环境对照','主文图3a、d、e的Q41/Q62用于环境及角度占据对照'))
replace(179,text(paras[179]).replace('先在每条记录内计算两个指标的Spearman相关，再对记录等权平均；','图3b将全部共同有效位点直接合并计算Spearman相关，MD为0.208615、预测为0.261905；该逐位点汇总会赋予较长记录更多权重。主视野为RMSF 0–10 Å与Hχ 0–4.5 nats，全部非零密度格保留，范围外位置不从统计删除，完整0–35 Å图随源数据提供。图S4a另先在每条记录内计算两个指标的Spearman相关，再对记录等权汇总；'))

# Replace exactly the figure3 image part. Preserve the relation ID and on-page width.
figp=paras[72]
embed=figp.xpath('.//a:blip/@r:embed',namespaces=NS)
assert embed==['rId13']
rels=etree.fromstring(payload['word/_rels/document.xml.rels'])
rel=[r for r in rels if r.get('Id')=='rId13'][0]
image_part='word/'+rel.get('Target')
assert image_part=='word/media/image5.png'
old_image_sha=sha(payload[image_part])
new_image=FIGURE.read_bytes()
with Image.open(FIGURE) as image:iw,ih=image.size;image.verify()
extent=figp.xpath('.//wp:extent',namespaces=NS)[0]
cx=int(extent.get('cx'));old_cy=int(extent.get('cy'));cy=round(cx*ih/iw)
extent.set('cy',str(cy))
for node in figp.xpath('.//a:xfrm/a:ext',namespaces=NS):node.set('cx',str(cx));node.set('cy',str(cy))
for node in figp.xpath('.//wp:docPr',namespaces=NS):
    node.set('name','Figure 3 dual-scale flexibility and angular occupancy')
    node.set('descr','Two structures, pooled RMSF-Hchi distribution, aligned CaMKI and ubiquitin profiles with SASA, and per-chi occupancy associations; highlights Q41 and Q62.')
payload[image_part]=new_image
payload['word/document.xml']=etree.tostring(root,xml_declaration=True,encoding='UTF-8',standalone=True)
with ZipFile(TARGET,'w',compression=ZIP_DEFLATED,compresslevel=6) as z:
    for item in items:z.writestr(item,payload[item.filename])

# Deterministic semantic/ZIP audits.
with ZipFile(SOURCE) as a, ZipFile(TARGET) as b:
    assert a.namelist()==b.namelist()
    changed_parts=[name for name in a.namelist() if a.read(name)!=b.read(name)]
    assert set(changed_parts)=={'word/document.xml',image_part},changed_parts
    assert b.testzip() is None
    r2=etree.fromstring(b.read('word/document.xml'))
    assert sha(b.read(image_part))==sha(new_image)
source_fields=[etree.tostring(x) for x in original.xpath('//w:fldChar|//w:fldSimple',namespaces=NS)]
target_fields=[etree.tostring(x) for x in r2.xpath('//w:fldChar|//w:fldSimple',namespaces=NS)]
assert source_fields==target_fields
fulltext='\n'.join(r2.xpath('//w:t/text()',namespaces=NS))
for obsolete in ['3.5节','3.5 TCR','3.4 侧链预测','图3f','图3b中的1MXE','主图集中呈现主链结果']:
    assert obsolete not in fulltext,obsolete
assert '3.4 TCR–pMHC识别中的界面柔性与亲和力差异' in fulltext
assert fulltext.count('3.3 主链与侧链柔性预测及角度占据的结构解释')==1
assert 'Q41与最高的Q2' in fulltext # supplemental extremes remain truthful
# Original bibliography and downstream TCR paragraphs are exactly retained.
original_refs=[text(p) for p in original.xpath('//w:body/w:p',namespaces=NS) if re.match(r'^\[\d+\]',text(p))]
new_refs=[text(p) for p in r2.xpath('//w:body/w:p',namespaces=NS) if re.match(r'^\[\d+\]',text(p))]
assert original_refs==new_refs
for i in range(82,88):assert etree.tostring(paras[i])==etree.tostring(original.xpath('//w:body/w:p',namespaces=NS)[i])
assert sha(SOURCE.read_bytes())==source_sha
report={'source':str(SOURCE),'target':str(TARGET),'source_sha256':source_sha,'target_sha256':sha(TARGET.read_bytes()),'changed_zip_parts':changed_parts,'figure':str(FIGURE),'figure_sha256':sha(new_image),'old_figure_sha256':old_image_sha,'image_pixels':[iw,ih],'image_emu':[cx,cy],'old_image_emu':[cx,old_cy],'figure_width_mm':cx/36000,'figure_height_mm':cy/36000,'preserved_field_nodes':len(source_fields),'bibliography_entries_preserved':len(new_refs),'paragraph_changes':changes,'validation':'XML/content/package passed; render review pending'}
(HERE/'document_change_report.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
(HERE/'rewritten_3.3_and_figure_caption.md').write_text('\n\n'.join([changes[j]['new'] for j in range(len(changes)) if changes[j]['source_paragraph_index'] in (68,69,70,71,73,74,75,76,78)]),encoding='utf-8')
(HERE/'README_改动与核验.md').write_text('从v98另存为v112；合并原3.3和3.4，后续TCR改为3.4，更新图3五面板图注并嵌入v112真实数据图。只改变document.xml和原图3的image5.png；其他ZIP部件、三个字段节点、全部文献及TCR正文保持原样。图3沿用164 mm文内宽度并按新图比例更新高度。补充图S4保留真实Q41/Q2极值核查，主图Q41/Q62使用真实坐标。使用原段落属性与字体，原文件未改。\n\n本目录的edit_v112_docx.py可从源文件完整复现。渲染和逐页视觉检查完成后更新render_QA.json。',encoding='utf-8')
print(json.dumps({k:v for k,v in report.items() if k!='paragraph_changes'},ensure_ascii=False,indent=2))

