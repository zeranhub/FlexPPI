from pathlib import Path
from collections import Counter
import hashlib,json,zipfile
from lxml import etree
from docx import Document
from docx.oxml.ns import qn
from docx.shared import Inches
from docx.text.paragraph import Paragraph
ROOT=Path(__file__).resolve().parent;BASE=ROOT.parent.parent
QA=ROOT/'manuscript_QA';QA.mkdir(exist_ok=True)
SOURCE=BASE/'paper_v27_20261014_中文修订_模型选择与自由能曲线验证_三线表.docx'
OUT=BASE/'paper_v28_20261014_中文修订_图7三条自由能曲线_三线表.docx'
STEM='fig7_model_selected_three_MD_free_energy_curves'
FIG=ROOT/'figures'/(STEM+'.png')
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
assert sha(SOURCE)=='2e4c566833dbf7f828f11c2a449248a7025fbe6011fc1da531c1528b816d4316'
d=Document(SOURCE);changes=[]
def find(prefix):return next(p for p in d.paragraphs if p.text.startswith(prefix))
def setp(p,text):
    assert not p._p.xpath('.//m:oMath')
    changes.append({'before':p.text,'after':text})
    nodes=p._p.xpath('.//w:t');nodes[0].text=text
    for n in nodes[1:]:n.text=''
p=find('低相对自由能角区总宽度按每个可定义χ角计算。')
setp(p,p.text.split('；范围曲线在E=')[0]+'。')
p=find('位置筛选仅使用各预测读数：')
setp(p,p.text.replace('范围曲线沿用主统计的同一组配对抽样索引，报告各阈值处的区间。',''))
p=find('按同一预测排序依据在各记录同类型内形成高、低半组，')
setp(p,p.text.replace('按同一预测排序依据在各记录同类型内形成高、低半组，测量实际MD范围随自由能阈值放宽的覆盖曲线（图7d）。E=1时，','按冻结侧链预测在各记录同类型内形成高、低半组，比较两组真实低自由能角区总宽度。F*≤1时，').replace('MD曲线及宽度用于评价这些位置','实际MD宽度用于评价这些位置'))
p=find('侧链输出的物理意义由预测驱动的MD检验建立。')
setp(p,p.text.replace('具体曲线说明模型选择了怎样的位置，群体范围曲线检验相同排序依据在更多蛋白上的表现，支持优先筛选具有较广实际取向范围的同类型位点。','具体曲线与群体预测分组结果共同支持利用侧链输出优先筛选实际取向范围较广的同类型位点。'))
caption=find('图7：');setp(caption,(ROOT/'figures'/(STEM+'.caption_zh.md')).read_text(encoding='utf8').strip())
pic=caption._p.getprevious();assert pic.xpath('.//w:drawing')
for child in list(pic):
    if child.tag!=qn('w:pPr'):pic.remove(child)
Paragraph(pic,d._body).add_run().add_picture(str(FIG),width=Inches(6.2))
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
audit={'document':str(OUT),'document_sha256':sha(OUT),'source_document':str(SOURCE),'source_sha256':sha(SOURCE),'figure7_sha256':sha(FIG),'changed_text_paragraphs':len(changes),'changes':changes,'all_tables_math_comments_preserved':True,'main_figures':7,'supplementary_figures':1,'tables':11,'review':'pending'}
(QA/'integration_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf8')
dst=Document(OUT)
assert '图7d' not in '\n'.join(p.text for p in dst.paragraphs)
(QA/'v28_content.txt').write_text('\n'.join(f'P{i}: '+p.text for i,p in enumerate(dst.paragraphs) if p.text),encoding='utf8')
oldroot=ROOT.parent/'free_energy_curves_v27_20261014'
for source,target in [('export_v27_word.ps1','export_v28_word.ps1'),('render_v27.py','render_v28.py')]:
    text=(oldroot/source).read_text(encoding='utf8').replace(SOURCE.name,OUT.name).replace('v27_word','v28_word').replace('render_v27','render_v28').replace('v27_page_map','v28_page_map')
    (ROOT/target).write_text(text,encoding='utf8')
print(json.dumps({k:v for k,v in audit.items() if k!='changes'},ensure_ascii=False))
