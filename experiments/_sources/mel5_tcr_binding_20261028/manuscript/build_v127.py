"""Revise 3.4 around TCR-specific binding predictions; preserve existing OOXML parts."""
from pathlib import Path
from copy import deepcopy
from zipfile import ZipFile,ZIP_DEFLATED
from hashlib import sha256
from lxml import etree
from docx import Document
from docx.shared import Cm
import csv,json
ROOT=Path(__file__).resolve().parent
EXP=ROOT.parent
PAPER=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010')
SOURCE=PAPER/'paper_v126_20261028_中文修订_绝对亲和力对照与肽突变验证.docx'
TARGET=PAPER/'paper_v127_20261028_中文修订_TCR结合模型对照.docx'
NS={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main','r':'http://schemas.openxmlformats.org/officeDocument/2006/relationships','wp':'http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing','a':'http://schemas.openxmlformats.org/drawingml/2006/main'}
def q(n):a,b=n.split(':');return '{'+NS[a]+'}'+b
with ZipFile(SOURCE) as z:
 infos=z.infolist();parts={x.filename:z.read(x.filename) for x in infos}
doc=etree.fromstring(parts['word/document.xml']);body=doc.find('w:body',NS);ps=body.findall('w:p',NS)
def txt(p):return ''.join(p.xpath('.//w:t/text()',namespaces=NS))
def find(prefix):return next(p for p in ps if txt(p).startswith(prefix))
main_template=find('肽单点替换提供')
body_template=find('绝对亲和力分别输入')
heading_template=find('F 绝对亲和力')
caption_template=find('图S7：')
ref_template=find('[32]')
def paragraph(value,template=body_template,page=False):
 p=etree.Element(q('w:p'))
 pp=deepcopy(template.find('w:pPr',NS))
 if pp is None:pp=etree.Element(q('w:pPr'))
 for x in list(pp):
  if x.tag in [q('w:sectPr'),q('w:pageBreakBefore')]:pp.remove(x)
 if page:etree.SubElement(pp,q('w:pageBreakBefore'))
 p.append(pp);r=etree.SubElement(p,q('w:r'))
 rp=template.find('w:r/w:rPr',NS)
 if rp is not None:r.append(deepcopy(rp))
 etree.SubElement(r,q('w:t')).text=value
 return p
with (EXP/'combined_binding_scores_wide.csv').open(encoding='utf-8-sig',newline='') as f:
 rows=list(csv.DictReader(f))
byid={r['case_id']:r for r in rows}
fmt=lambda key,model:f"{float(byid[key][model]):.3f}"
net='NetTCR2.2 known-overlap-excluded';im='ImRex downsample'
main=(
 '肽单点替换提供了亲和力显著不同而均可检测结合的对照。G4A与I7A使成熟受体的SPR解离常数由0.6 nM分别增至36与31 μM，但四态均有可测结合[28]。我们对各状态分别运行ERGO-II、NetTCR-2.2和ImRex，比较结合评分而不要求输出KD[33–35]。'
 f'NetTCR-2.2对成熟ELA、G4A、I7A的评分为{fmt("4jff",net)}、{fmt("4jfd",net)}、{fmt("4jfe",net)}，ImRex默认版为{fmt("4jff",im)}、{fmt("4jfd",im)}、{fmt("4jfe",im)}，均给出肽替换后下降的方向；ERGO-II两套权重则均未给出该方向。'
 '结合评分可以恢复部分替换方向，其响应依模型而异；这些识别倾向评分尚未描述具体位置的结合态动态。完整分数、输入、训练重叠与实验不结合对照见表S8、图S8和补充材料G；通用绝对亲和力对照见补充材料F。'
)
changes=[]
body.replace(main_template,paragraph(main,main_template));changes.append({'before':txt(main_template),'after':main})
disc=find('MEL5家族留出实例')
old=txt(disc)
new=old.replace('亲和力评分、残基动态与接触占据分别描述效应量、响应位置和轨迹中的局部组织，三者的结合为机制研究提供可检验线索；','TCR结合评分能够恢复部分肽替换的方向，残基动态与接触占据进一步描述响应位置和轨迹中的局部组织，为机制研究提供可检验线索；')
assert old!=new
body.replace(disc,paragraph(new,disc));changes.append({'before':old,'after':new})
refs=[
 '[33] I. Springer, N. Tickotsky, and Y. Louzoun. Contribution of T Cell Receptor Alpha and Beta CDR3, MHC Typing, V and J Genes to Peptide Binding Prediction. Frontiers in Immunology, 12:664514, 2021. doi:10.3389/fimmu.2021.664514.',
 '[34] M. F. Jensen and M. Nielsen. Enhancing TCR specificity predictions by combined pan- and peptide-specific training, loss-scaling, and sequence similarity integration. eLife, 12:RP93934, 2024. doi:10.7554/eLife.93934.3.',
 '[35] P. Moris, J. De Pauw, A. Postovskaya, S. Gielis, N. De Neuter, W. Bittremieux, B. Ogunjimi, K. Laukens, and P. Meysman. Current challenges for unseen-epitope TCR interaction prediction and a new perspective derived from image classification. Briefings in Bioinformatics, 22(4):bbaa318, 2021. doi:10.1093/bib/bbaa318.',
 '[36] A. Greenshields-Watson et al. ANARCII enables alignment-free antigen receptor numbering using a generalised language model. Communications Biology, 9:1085, 2026. doi:10.1038/s42003-026-10186-z.'
]
for ref in refs:find('补充材料').addprevious(paragraph(ref,ref_template))
supp=[]
def add(value,template=body_template,page=False):
 p=paragraph(value,template,page);body.insert(len(body)-1,p);supp.append(value);return p
add('G TCR结合模型对照',heading_template,True)
add('本对照采用三种TCR结合预测方法的公开预训练权重，逐条件独立输出结合评分，不使用亲和力变化专用模型，不将分数换算为KD[33–35]。四个核心状态的SPR结果均为可检测结合，G4A与I7A不能因较弱而标作不结合。为检验二分类，另纳入原论文表2中其余五种成熟受体肽替换，及正文明确报告不结合的ILAKFLHWL、SLLMWITQC、YLEPGPVTA三个SPR对照，共12行、8个可测结合和4个未检测结合[28]。其中I5A（ELAGAGILTV）为表2中的单点替换NB；三个无关肽对照原文注明数据未展示。面板及标签在预测前固定，不自设KD阈值。')
add('序列取原始RCSB晶体SEQRES，经ANARCII 2.0.8的TCR accuracy模式编号，并与已解析序列核对[36]。ERGO-II与ImRex使用含IMGT C104/F118的CDR3，NetTCR使用六个CDR且CDR3为105–117，不含锚点。三种成熟复合物六CDR完全相同；其余扫描与负对照沿用成熟受体CDR。V基因为TCR3d核验的TRAV12-2、TRBV30；J未可靠注释，ERGO-II按作者允许方式使用UNK，VDJdb模型中的无等位基因V名称也映射UNK。ERGO-II不直接输入工程化CDR1/2，ImRex仅输入CDR3β和肽；因此其比较不等于完整受体结构比较。')
add('表S8：四核心态的TCR结合评分。SPR KD以nM为单位，仅为实验背景，各模型无KD拟合目标。ERGO-II列为McPAS、VDJdb官方权重；NetTCR列为排除已确认重叠的4模型及全20模型均值；ImRex列为官方默认downsample及full数据版。分数仅在同一配置内比较。',caption_template)
methods=['ERGO-II McPAS','ERGO-II VDJdb',net,'NetTCR2.2 all20',im,'ImRex full']
headers=['PDB','SPR KD\n(nM)','ERGO-II\nMcPAS','ERGO-II\nVDJdb','NetTCR\n4模型','NetTCR\n20模型','ImRex\n默认','ImRex\nfull']
table_rows=[headers]
for key in ['3hg1','4jff','4jfd','4jfe']:
 r=byid[key];vals=[]
 for method in methods:
  v=float(r[method]);vals.append(f'{v:.3f}' if v>=.001 else f'{v:.2e}')
 table_rows.append([key.upper(),f"{float(r['experimental_KD_nM']):g}"]+vals)
tbl=etree.Element(q('w:tbl'));tp=etree.SubElement(tbl,q('w:tblPr'))
widths=[850,1100]+[1250]*6
etree.SubElement(tp,q('w:tblW'),{q('w:w'):str(sum(widths)),q('w:type'):'dxa'})
etree.SubElement(tp,q('w:tblLayout'),{q('w:type'):'fixed'})
borders=etree.SubElement(tp,q('w:tblBorders'))
for side in ['top','left','bottom','right','insideH','insideV']:etree.SubElement(borders,q('w:'+side),{q('w:val'):'single',q('w:sz'):'4',q('w:color'):'D9D9D9'})
grid=etree.SubElement(tbl,q('w:tblGrid'))
for width in widths:etree.SubElement(grid,q('w:gridCol'),{q('w:w'):str(width)})
for ri,values in enumerate(table_rows):
 tr=etree.SubElement(tbl,q('w:tr'))
 for width,value in zip(widths,values):
  tc=etree.SubElement(tr,q('w:tc'));tcpr=etree.SubElement(tc,q('w:tcPr'))
  etree.SubElement(tcpr,q('w:tcW'),{q('w:w'):str(width),q('w:type'):'dxa'});etree.SubElement(tcpr,q('w:vAlign'),{q('w:val'):'center'})
  if ri==0:etree.SubElement(tcpr,q('w:shd'),{q('w:fill'):'E8EEF4'})
  p=etree.SubElement(tc,q('w:p'));pp=etree.SubElement(p,q('w:pPr'))
  etree.SubElement(pp,q('w:ind'),{q('w:firstLine'):'0',q('w:firstLineChars'):'0',q('w:left'):'0',q('w:right'):'0'})
  etree.SubElement(pp,q('w:jc'),{q('w:val'):'center'});etree.SubElement(pp,q('w:snapToGrid'),{q('w:val'):'0'})
  etree.SubElement(pp,q('w:spacing'),{q('w:before'):'75',q('w:after'):'75',q('w:line'):'240',q('w:lineRule'):'auto'})
  for li,line in enumerate(value.split('\n')):
   rr=etree.SubElement(p,q('w:r'));rp=etree.SubElement(rr,q('w:rPr'))
   etree.SubElement(rp,q('w:rFonts'),{q('w:ascii'):'Arial',q('w:hAnsi'):'Arial',q('w:eastAsia'):'宋体'})
   etree.SubElement(rp,q('w:sz'),{q('w:val'):'17'})
   if ri==0:etree.SubElement(rp,q('w:b'))
   if li:etree.SubElement(rr,q('w:br'))
   etree.SubElement(rr,q('w:t')).text=line
body.insert(len(body)-1,tbl)
add('NetTCR统一使用pan模式，不将ELA的肽特异模型与突变肽的pan模型混用，也不加TCRbase缩放。预测前发现11条精确CDR3配对相关训练记录均位于分区0，主结果固定为t.0.v.1–4四个官方模型均值，其训练及验证不使用分区0；全20模型为敏感性对照。该处理仅排除已确认重叠，不保证所有近缘家族隔离。ImRex按官方README默认选择downsample，并保留full版；原生MEL5–ELA在其训练集中，成熟CDR3β未见。ERGO-II的VDJdb训练中存在原生β–ELA记录。三方法均非本研究FlexPPI的同一留出训练体系，不能将其当作统一家族隔离的公平性能基准。')

# Add final figure with a fresh relationship; all existing figures remain binary-identical.
rels=etree.fromstring(parts['word/_rels/document.xml.rels'])
rnum=max(int(x.get('Id')[3:]) for x in rels if x.get('Id','').startswith('rId') and x.get('Id')[3:].isdigit())+1
rid=f'rId{rnum}'
figure=EXP/'figures/figS8_TCR_binding_scores.png'
temp=Document();temp.add_picture(str(figure),width=Cm(17.0));p=etree.fromstring(etree.tostring(temp.paragraphs[-1]._p))
pp=etree.Element(q('w:pPr'));etree.SubElement(pp,q('w:jc'),{q('w:val'):'center'});etree.SubElement(pp,q('w:pageBreakBefore'));etree.SubElement(pp,q('w:keepNext'));etree.SubElement(pp,q('w:spacing'),{q('w:before'):'0',q('w:after'):'40',q('w:line'):'240',q('w:lineRule'):'auto'});p.insert(0,pp)
for blip in p.xpath('.//a:blip',namespaces=NS):blip.set(q('r:embed'),rid)
idnum=max(int(x.get('id')) for x in doc.xpath('//wp:docPr',namespaces=NS))+1
media='word/media/v127_figS8_TCR_binding_scores.png'
for x in p.xpath('.//wp:docPr',namespaces=NS):x.set('id',str(idnum));x.set('name','v127_figS8_TCR_binding_scores.png')
etree.SubElement(rels,'{http://schemas.openxmlformats.org/package/2006/relationships}Relationship',{'Id':rid,'Type':NS['r']+'/image','Target':'media/v127_figS8_TCR_binding_scores.png'})
body.insert(len(body)-1,p)
add('图S8：TCR专用结合模型对完整12行文献面板的评分。上四行为原四态，中五行为补充的成熟受体Ala替换，下三行为SPR不结合对照。加号表示原实验检测到结合，NB表示原实验未检测到结合；NB不表示严格为零的亲和力。保留三种模型全部六种配置，共72个分数，未按表现删选配置。NetTCR主四模型按已确认训练重叠在预测前固定，全20模型仅作敏感性对照。深色仅代表该配置较高结合分数，不表示跨模型经过共同概率校准；图中不以0.5切分，也不作总体性能显著性检验。',caption_template)
add('识别倾向、结合效应量和结合态动态分别解释。ImRex默认版对四核心态均给出高于0.5的分数，并恢复G4A与I7A下降方向，但I5A（实验NB）和YLEPGPVTA负对照的评分分别为0.631与0.611，在固定0.5参考切分下判为结合。full版对I5A为0.232，对I7A却高于ELA。NetTCR主结果中I7A阳性与I5A阴性的评分为0.078与0.076；ERGO-II两权重的G4A/I7A分数均高于ELA。作者推理接口主要输出连续分数，0.5只作未调参且未经本SPR实验校准的参考，不能把阈值结果直接作为方法普遍失效的证据。')
add('全部分数、输入及训练覆盖在配套CSV、逐模型审计和运行日志中提供。NetTCR对作者原CV表11行的复算最大差为1.5×10⁻⁸，与未修改predict.py交叉核验最大差为2.5×10⁻⁸；ImRex使用未修改官方CLI及特征流程；ERGO-II保留官方网络和权重，严格加载，仅适配旧Lightning推理容器及CPU放置。结合评分可以帮助识别部分替换方向；图S6–S7则进一步检验哪些残基与接触参与结合态动态响应，两类证据不互相替代。')
parts['word/document.xml']=etree.tostring(doc,xml_declaration=True,encoding='UTF-8',standalone=True)
parts['word/_rels/document.xml.rels']=etree.tostring(rels,xml_declaration=True,encoding='UTF-8',standalone=True)
with ZipFile(TARGET,'w',compression=ZIP_DEFLATED,compresslevel=6) as z:
 for item in infos:z.writestr(item,parts[item.filename])
 z.writestr(media,figure.read_bytes())
with ZipFile(SOURCE) as a,ZipFile(TARGET) as b:
 assert b.testzip() is None
 original_media=[n for n in a.namelist() if n.startswith('word/media/')]
 assert all(a.read(n)==b.read(n) for n in original_media)
 changed=[n for n in a.namelist() if a.read(n)!=b.read(n)]
 assert set(changed)=={'word/document.xml','word/_rels/document.xml.rels'}
 olddoc=etree.fromstring(a.read('word/document.xml'))
 for query in ['//w:tbl','//w:fldSimple','//w:instrText','//w:commentRangeStart','//w:commentRangeEnd','//w:sectPr']:
  oldnodes=olddoc.xpath(query,namespaces=NS);newnodes=doc.xpath(query,namespaces=NS)
  assert all(etree.tostring(x)==etree.tostring(y) for x,y in zip(oldnodes,newnodes)),query
report=dict(source=str(SOURCE),target=str(TARGET),source_sha256=sha256(SOURCE.read_bytes()).hexdigest(),target_sha256=sha256(TARGET.read_bytes()).hexdigest(),changed_original_parts=changed,original_media_preserved=len(original_media),original_tables_fields_comments_sections_preserved=True,new_media=media,changed_paragraphs=changes,references_added=refs,supplement=supp,render_validation='pending')
(ROOT/'revision_changes.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
(ROOT/'3.4_TCR对照修订正文.md').write_text(main+'\n\n'+new+'\n',encoding='utf-8')
print(json.dumps(dict(target=str(TARGET),main_paragraph_characters=len(main),media_preserved=len(original_media)),ensure_ascii=False))
