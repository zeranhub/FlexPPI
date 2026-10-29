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
SOURCE=BASE/'paper_v36_20261015_中文修订_图7逐角散点与三子图布局_三线表.docx'
OUT=BASE/'paper_v37_20261015_中文修订_图7删去汇总表与布局精简_三线表.docx'
STEM='fig7_clean_chi_scatter_and_MD_curves'
FIG=ROOT/'figures'/(STEM+'.png')
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
assert sha(SOURCE)=='84961145b46d467f4ba262cc3c6a5ce88c4d5f4e99282c2e349b496ceca07706'
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
old='，并汇总预测与细箱熵、低自由能角区总宽度的相关'
assert old in p.text
setp(p,p.text.replace(old,''))
p=find('在泛素各χ角层内控制氨基酸类型后')
old='，与低自由能角区总宽度的关联在χ3最强（ρ=0.640）；χ4的两项相关分别为0.127和−0.183（图7b）。'
new='，χ4的对应相关为0.127（图7b）。'
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
assert len(changes)==3
outdoc=Document(OUT)
assert len(outdoc.paragraphs)==len(d.paragraphs)
(QA/'v37_content.txt').write_text('\n'.join(f'P{i}: '+p.text for i,p in enumerate(outdoc.paragraphs) if p.text),encoding='utf8')
audit={'document':str(OUT),'document_sha256':sha(OUT),'source_document':str(SOURCE),'source_sha256':sha(SOURCE),'figure7_sha256':sha(FIG),'changed_text_paragraphs':len(changes),'changes':changes,'all_tables_math_comments_preserved':True,'main_figures':7,'supplementary_figures':1,'tables':11,'review':'pending'}
(QA/'integration_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps({k:v for k,v in audit.items() if k!='changes'},ensure_ascii=False))
