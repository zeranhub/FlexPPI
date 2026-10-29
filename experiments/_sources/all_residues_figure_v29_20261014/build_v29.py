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
SOURCE=BASE/'paper_v28_20261014_中文修订_图7三条自由能曲线_三线表.docx'
OUT=BASE/'paper_v29_20261014_中文修订_全残基与同类型预测及自由能曲线_三线表.docx'
STEM='fig7_all_residues_Gln_and_MD_free_energy';FIG=ROOT/'figures'/(STEM+'.png')
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
assert sha(SOURCE)=='e86276be39fb476680c0f915d2f155c8be0030599789e035d50a1b83d1dd97ac'
d=Document(SOURCE);changes=[]
def find(prefix):return next(p for p in d.paragraphs if p.text.startswith(prefix))
def setp(p,text):
    assert not p._p.xpath('.//m:oMath')
    changes.append({'before':p.text,'after':text})
    nodes=p._p.xpath('.//w:t');nodes[0].text=text
    for n in nodes[1:]:n.text=''
p=find('侧链实例使用泛素–EDD UBA复合物')
setp(p,p.text.replace('图7仅按冻结Hχ预测在泛素六个有效谷氨酰胺位点中排序，选择最低Q41与最高Q2；逐角曲线统一取0–100 ns、200 ps间隔的501帧，原五窗口联合熵按主标签协议核对。','图7首先对照泛素链全部66个Hχ有效残基及全部六个有效谷氨酰胺位点的冻结预测与主评价参考，再按六个位点的Hχ预测选择最低Q41与最高Q2，展示逐角MD曲线。逐角曲线统一取0–100 ns、200 ps间隔的501帧。'))
p=find('泛素–EDD UBA复合物提供同类型侧链的留出验证')
setp(p,'泛素–EDD UBA复合物提供同类型侧链的留出验证[23]。完整模型在泛素链的15个预测最高主链波动位点中命中14个（图6e、g）。在具有侧链联合熵参考值的全部66个残基上，Hχ预测的平均绝对误差为0.080 nats，参考与预测的排序相关为0.979；分别减去训练折的氨基酸类型均值后，相关仍为0.868（图7a）。模型由静态主链输入同时定位活动区域，并区分残基类型平均值之外的侧链状态差异。')
p=find('进一步以模型预测确定曲线验证位点：')
setp(p,'同类型评价纳入泛素全部六个有效谷氨酰胺位点，联合Hχ预测的平均绝对误差为0.105 nats，逐位误差为0.019–0.368 nats，其中Q40高估0.368 nats（图7b）。进一步按这六个位点的冻结预测确定曲线验证位置：Hχ预测最低为Q41（2.000 nats），最高为Q2（2.990 nats），原MD联合熵参考分别为1.890和3.009 nats。图7c–e展示这两个预测端点的实际逐角自由能投影。两个位置的χ1主占据区均在−60°附近；Q2的χ2同时覆盖反式及−60°邻域，χ3的低值角区也更广。在共同501帧窗口和四个箱格原点平均下，Q41与Q2的逐角低自由能总宽度分别为27.50和64.17°/角。全位点评价量化模型预测的准确性，同类型对照展示各位置的预测差异，所选端点的真实曲线进一步检验这些差异对应的取向范围。')
caption=find('图7：');setp(caption,(ROOT/'figures'/(STEM+'.caption_zh.md')).read_text(encoding='utf8').strip())
pic=caption._p.getprevious();assert pic.xpath('.//w:drawing')
for child in list(pic):
    if child.tag!=qn('w:pPr'):pic.remove(child)
picp=Paragraph(pic,d._body)
picp.paragraph_format.page_break_before=False
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
assert '图7a–c' not in '\n'.join(p.text for p in dst.paragraphs)
(QA/'v29_content.txt').write_text('\n'.join(f'P{i}: '+p.text for i,p in enumerate(dst.paragraphs) if p.text),encoding='utf8')
oldroot=ROOT.parent/'curve_only_v28_20261014'
for source,target in [('export_v28_word.ps1','export_v29_word.ps1'),('render_v28.py','render_v29.py')]:
    text=(oldroot/source).read_text(encoding='utf8').replace(SOURCE.name,OUT.name).replace('v28_word','v29_word').replace('render_v28','render_v29').replace('v28_page_map','v29_page_map')
    (ROOT/target).write_text(text,encoding='utf8')
print(json.dumps({k:v for k,v in audit.items() if k!='changes'},ensure_ascii=False))
