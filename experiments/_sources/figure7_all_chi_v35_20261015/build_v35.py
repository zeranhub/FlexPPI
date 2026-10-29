from pathlib import Path
import hashlib, json, zipfile
from collections import Counter
from lxml import etree
from docx import Document
from docx.oxml.ns import qn
from docx.shared import Inches
from docx.text.paragraph import Paragraph
ROOT=Path(__file__).resolve().parent
BASE=ROOT.parent.parent
QA=ROOT/'manuscript_QA'; QA.mkdir(exist_ok=True)
SOURCE=BASE/'paper_v34_20261015_中文修订_图7六谷氨酰胺顺序对照_三线表.docx'
OUT=BASE/'paper_v35_20261015_中文修订_图7全残基逐χ占据关联_三线表.docx'
STEM='fig7_all_chi_associations_and_MD_curves'
FIG=ROOT/'figures'/(STEM+'.png')
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
assert sha(SOURCE)=='da9a276b62d94fb5e7a77668a55d92100a8c9896d0a6cb586483c3670d34a0df'
d=Document(SOURCE); changes=[]
def find(prefix):
    matches=[p for p in d.paragraphs if p.text.startswith(prefix)]
    assert len(matches)==1,(prefix,len(matches))
    return matches[0]
def setp(p,text):
    before=p.text
    assert not p._p.xpath('.//m:oMath')
    nodes=p._p.xpath('.//w:t'); assert nodes
    nodes[0].text=text
    for n in nodes[1:]:n.text=''
    assert before!=text
    changes.append({'before':before,'after':text})
p=find('侧链实例使用泛素–EDD UBA复合物')
old=p.text[p.text.index('图7按残基编号'):]
new='图7a左按残基编号展开泛素链全部66个Hχ有效位置的冻结预测与主评价参考，仅标出Q2与Q41；右侧分别汇总χ1–χ4的全部有效位点。每个χ角层内，将联合Hχ预测与MD单角细箱占据熵及低自由能角区总宽度配对，各自减去相同氨基酸类型的均值后计算Spearman相关。各χ层分别包含66、52、23和10个位点。泛素六个谷氨酰胺中联合Hχ预测最低Q41与最高Q2另展示逐角MD曲线。逐角指标与曲线统一取0–100 ns、200 ps间隔501帧，指标按2.8节对四个箱格原点取平均。'
setp(p,p.text.replace(old,new).replace('以位点9 Å内局部主链对齐。','以该位点9 Å内的局部主链对齐。'))
p=find('全部六个谷氨酰胺按残基编号排列')
setp(p,'在泛素各χ角层内控制氨基酸类型后，联合Hχ预测与χ1–χ3的MD细箱占据熵均呈正相关（ρ=0.447–0.525），与低自由能角区总宽度的关联在χ3最强（ρ=0.640）；χ4的两项相关分别为0.127和−0.183（图7a右）。逐角比较将模型读数对应到具体侧链角的占据差异。在泛素六个谷氨酰胺中，模型预测最低Q41与最高Q2的Hχ分别为2.000和2.990 nats，原MD联合熵参考为1.890和3.009 nats。图7b–d进一步展开这两个预测端点的实际逐角自由能投影。两个位置的χ1主占据区均在−60°附近；Q2的χ2同时覆盖反式及−60°邻域，χ3的低值角区也更广。在共同501帧窗口和四个箱格原点平均下，Q41与Q2的逐角低自由能总宽度分别为27.50和64.17°/角。全序列对照量化预测精度，逐角关联及端点曲线进一步将侧链读数对应到实际占据差异。')
p=find('侧链输出的物理意义由预测驱动的MD检验建立')
old='泛素留出预测识别类型均值之外的位置差异；'
new='泛素留出预测识别类型均值之外的位置差异，并在χ1–χ3上与MD单角占据保持正向联系；'
assert old in p.text
setp(p,p.text.replace(old,new))
caption=find('图7：')
setp(caption,(ROOT/'figures'/(STEM+'.caption_zh.md')).read_text(encoding='utf8').strip())
pic=caption._p.getprevious(); assert pic.xpath('.//w:drawing')
for child in list(pic):
    if child.tag!=qn('w:pPr'):pic.remove(child)
picp=Paragraph(pic,d._body);picp.paragraph_format.page_break_before=False
picp.add_run().add_picture(str(FIG),width=Inches(6.2))
d.save(OUT)
with zipfile.ZipFile(SOURCE) as a,zipfile.ZipFile(OUT) as b:
    entries=[(info,a.read(info.filename) if info.filename=='word/comments.xml' else b.read(info.filename)) for info in b.infolist()]
tmp=OUT.with_suffix('.preserved.tmp')
with zipfile.ZipFile(tmp,'w') as z:
    for info,data in entries:z.writestr(info,data)
tmp.replace(OUT)
ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main','m':'http://schemas.openxmlformats.org/officeDocument/2006/math','a':'http://schemas.openxmlformats.org/drawingml/2006/main'}
with zipfile.ZipFile(SOURCE) as a,zipfile.ZipFile(OUT) as b:
    ax,bx=etree.fromstring(a.read('word/document.xml')),etree.fromstring(b.read('word/document.xml'))
    anchors=lambda x:Counter((e.tag,e.get(qn('w:id'))) for e in x.iter() if e.tag in {qn('w:commentRangeStart'),qn('w:commentRangeEnd'),qn('w:commentReference')})
    native=lambda x,query:[etree.tostring(e,method='c14n') for e in x.xpath(query,namespaces=ns)]
    assert anchors(ax)==anchors(bx) and a.read('word/comments.xml')==b.read('word/comments.xml')
    assert native(ax,'.//m:oMath')==native(bx,'.//m:oMath') and len(native(bx,'.//m:oMath'))==36
    assert native(ax,'.//w:tbl')==native(bx,'.//w:tbl') and len(native(bx,'.//w:tbl'))==11
    assert len(bx.xpath('.//a:blip',namespaces=ns))==8
assert len(changes)==4
outdoc=Document(OUT)
assert len(outdoc.paragraphs)==len(d.paragraphs)
(QA/'v35_content.txt').write_text('\n'.join(f'P{i}: '+p.text for i,p in enumerate(outdoc.paragraphs) if p.text),encoding='utf8')
audit={'document':str(OUT),'document_sha256':sha(OUT),'source_document':str(SOURCE),'source_sha256':sha(SOURCE),'figure7_sha256':sha(FIG),'changed_text_paragraphs':len(changes),'changes':changes,'all_tables_math_comments_preserved':True,'main_figures':7,'supplementary_figures':1,'tables':11,'review':'pending'}
(QA/'integration_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps({k:v for k,v in audit.items() if k!='changes'},ensure_ascii=False))
