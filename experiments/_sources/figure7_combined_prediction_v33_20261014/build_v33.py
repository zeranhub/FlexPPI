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
SOURCE=BASE/'paper_v32_20261014_中文修订_图7横轴标注整理_三线表.docx'
OUT=BASE/'paper_v33_20261014_中文修订_图7预测合并与自由能曲线_三线表.docx'
STEM='fig7_combined_predictions_and_MD_curves';FIG=ROOT/'figures'/(STEM+'.png')
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
assert sha(SOURCE)=='301daed0997726951d5b6a958532606faffd5a74480b29170ceaabca39d95d35'
d=Document(SOURCE);changes=[]
def find(prefix):return next(p for p in d.paragraphs if p.text.startswith(prefix))
def setp(p,text):
    assert not p._p.xpath('.//m:oMath')
    changes.append({'before':p.text,'after':text})
    nodes=p._p.xpath('.//w:t');nodes[0].text=text
    for n in nodes[1:]:n.text=''
p=find('侧链实例使用泛素–EDD UBA复合物')
setp(p,p.text.replace('六个有效谷氨酰胺位点按冻结预测从低到高排列，展开各自χ1–χ3的完整角度占据；最低Q41与最高Q2另展示逐角MD曲线。逐角占据、宽度与曲线统一取0–100 ns、200 ps间隔501帧。','六个有效谷氨酰胺位点在图7a右侧按冻结预测从低到高排列，对照联合熵参考、预测与带符号误差；最低Q41与最高Q2另展示逐角MD曲线。逐角宽度与曲线统一取0–100 ns、200 ps间隔501帧。'))
p=find('全部六个谷氨酰胺按冻结Hχ预测从低到高排列')
setp(p,'全部六个谷氨酰胺按冻结Hχ预测从低到高排列，对照联合熵预测、MD参考及逐位误差（图7a右）。六个位点的预测平均绝对误差为0.105 nats，其中Q40高估0.368 nats。模型预测最低Q41与最高Q2的Hχ分别为2.000和2.990 nats，原MD联合熵参考为1.890和3.009 nats。图7b–d进一步展开这两个预测端点的实际逐角自由能投影。两个位置的χ1主占据区均在−60°附近；Q2的χ2同时覆盖反式及−60°邻域，χ3的低值角区也更广。在共同501帧窗口和四个箱格原点平均下，Q41与Q2的逐角低自由能总宽度分别为27.50和64.17°/角。完整的预测剖面与同类型逐位对照量化模型的准确性，端点曲线进一步展示预测高低对应的实际取向差异。')
caption=find('图7：');setp(caption,(ROOT/'figures'/(STEM+'.caption_zh.md')).read_text(encoding='utf8').strip())
pic=caption._p.getprevious();assert pic.xpath('.//w:drawing')
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
audit={'document':str(OUT),'document_sha256':sha(OUT),'source_document':str(SOURCE),'source_sha256':sha(SOURCE),'figure7_sha256':sha(FIG),'changed_text_paragraphs':len(changes),'changes':changes,'all_tables_math_comments_preserved':True,'main_figures':7,'supplementary_figures':1,'tables':11,'review':'pending'}
(QA/'integration_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf8')
dst=Document(OUT)
(QA/'v33_content.txt').write_text('\n'.join(f'P{i}: '+p.text for i,p in enumerate(dst.paragraphs) if p.text),encoding='utf8')
oldroot=ROOT.parent/'figure7_clean_axis_v32_20261014'
for source,target in [('export_v32_word.ps1','export_v33_word.ps1'),('render_v32.py','render_v33.py')]:
    text=(oldroot/source).read_text(encoding='utf8').replace(SOURCE.name,OUT.name).replace('v32_word','v33_word').replace('render_v32','render_v33').replace('v32_page_map','v33_page_map')
    (ROOT/target).write_text(text,encoding='utf8')
print(json.dumps({k:v for k,v in audit.items() if k!='changes'},ensure_ascii=False))
