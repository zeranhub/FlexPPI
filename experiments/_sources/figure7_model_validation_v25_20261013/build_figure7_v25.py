from pathlib import Path
from collections import Counter
from copy import deepcopy
import hashlib, json, zipfile, re
from lxml import etree
from docx import Document
from docx.oxml.ns import qn
from docx.shared import Inches
from docx.text.paragraph import Paragraph
from docx.enum.text import WD_ALIGN_PARAGRAPH

ROOT=Path(__file__).resolve().parent
BASE=ROOT.parent.parent
QA=ROOT/'manuscript_QA';QA.mkdir(exist_ok=True)
SOURCE=BASE/'paper_v24_20261013_中文修订_模型论证与构象占据验证_三线表.docx'
OUT=BASE/'paper_v25_20261013_中文修订_图7预测对照明确化_三线表.docx'
FIG=ROOT/'figures/fig7_ubiquitin_prediction_and_MD_Fstar.png'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
assert sha(SOURCE)=='de8e19c4b941dbf0217c2b539e731e5932858c64f1fb81f7c07153e1b45a242e'
d=Document(SOURCE)
changes=[]
def setp(p,text):
    assert not p._p.xpath('.//m:oMath')
    changes.append({'before':p.text,'after':text})
    nodes=p._p.xpath('.//w:t');assert nodes
    nodes[0].text=text;nodes[0].set('{http://www.w3.org/XML/1998/namespace}space','preserve')
    for n in nodes[1:]:n.text=''

site=next(p for p in d.paragraphs if p.text.startswith('Q41与Q62均为具有三个侧链角'))
setp(site,'Q41与Q62均为具有三个侧链角的谷氨酰胺。冻结留出模型的联合Hχ预测分别为2.000和2.684 nats，MD参考为1.890和2.631 nats，绝对偏差分别为0.110和0.053 nats（图7a）。模型保留Q62高于Q41的顺序，并接近两个位置的参考幅值；Q62减Q41的参考差值为0.741 nats，预测为0.684 nats。Q41较埋藏，Q62朝向溶剂，MD中覆盖90%观测帧分别需要七类和十二类联合状态（图6f、h）。这一同类型对照将模型识别的侧链状态差异对应到实际采样的取向组合。')
angles=next(p for p in d.paragraphs if p.text.startswith('逐角分布进一步解释这一预测差异'))
setp(angles,'MD的逐角经验投影自由能进一步展开这一差异的分布来源（图7b–d）。两个位点的χ1均主要访问−60°附近；Q41的χ2集中于跨越±180°接缝的同一反式取向邻域，Q62的χ2主要访问−60°附近，其χ3取向覆盖也更分散。相近的χ1主要占据区因此可以伴随不同的后续角度及联合组合多样性。联合Hχ的预测与参考对照检验模型对这两个位点的数值及高低差异的估计，逐角MD投影则解释这些读数对应的实际取向形态。')
caption=next(p for p in d.paragraphs if p.text.startswith('图7：'))
new_caption=(ROOT/'figures/fig7_ubiquitin_prediction_and_MD_Fstar.caption_zh.md').read_text(encoding='utf8').strip().removesuffix('\\n')
setp(caption,new_caption)
pic=caption._p.getprevious();assert pic.xpath('.//w:drawing')
old_image_id=pic.xpath('.//a:blip')[0].get(qn('r:embed'))
old_image_hash=hashlib.sha256(d.part.related_parts[old_image_id].blob).hexdigest()
for child in list(pic):
    if child.tag!=qn('w:pPr'):pic.remove(child)
pp=Paragraph(pic,d._body);pp.alignment=WD_ALIGN_PARAGRAPH.CENTER
pp.paragraph_format.keep_with_next=True
pp.add_run().add_picture(str(FIG),width=Inches(6.2))
d.save(OUT)
# Keep all original comments byte for byte.
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
    tables=lambda x:[etree.tostring(e,method='c14n') for e in x.xpath('.//w:tbl',namespaces=ns)]
    assert anchor(ax)==anchor(bx) and math(ax)==math(bx)
    assert a.read('word/comments.xml')==b.read('word/comments.xml')
    assert tables(ax)==tables(bx)
    math_count=len(math(bx))
    assert len(bx.xpath('.//a:blip',namespaces=ns))==9
src=Document(SOURCE);dst=Document(OUT)
media=lambda doc:[hashlib.sha256(doc.part.related_parts[p.xpath('.//a:blip')[0].get(qn('r:embed'))].blob).hexdigest() for p in doc._element.body.xpath('./w:p') if p.xpath('.//a:blip')]
ma,mb=media(src),media(dst)
assert len(ma)==len(mb)==9 and sum(x!=y for x,y in zip(ma,mb))==1
assert sha(FIG) in mb
assert len(src.paragraphs)==len(dst.paragraphs)
changed_indices=[i for i,(a,b) in enumerate(zip(src.paragraphs,dst.paragraphs)) if a.text!=b.text]
assert len(changed_indices)==3
assert len(dst.tables)==11
refs=[p.text for p in dst.paragraphs if re.match(r'^\[\d+\]',p.text)]
assert len(refs)==27 and refs==[p.text for p in src.paragraphs if re.match(r'^\[\d+\]',p.text)]
audit={'document':str(OUT),'document_sha256':sha(OUT),'source_document':str(SOURCE),'source_sha256':sha(SOURCE),'changed_paragraph_indices':changed_indices,'changes':changes,'main_figures':8,'supplementary_figures':1,'main_tables':2,'supplementary_tables':9,'references':27,'all_11_native_three_line_tables_byte_equivalent':True,'other_8_figures_byte_equivalent':True,'math_elements':math_count,'math_and_comment_anchors_preserved':True,'comments_xml_byte_equivalent':True,'old_figure7_sha256':old_image_hash,'new_figure7_sha256':sha(FIG),'layout_review':'pending'}
(QA/'integration_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf8')
lines=[]
for i,p in enumerate(dst.paragraphs):
    t=''.join(p._p.xpath('.//w:t/text() | .//m:t/text()'))
    if t:lines.append(f'P{i}: {t}')
for i,t in enumerate(dst.tables):
    for j,row in enumerate(t.rows):lines.append(f'T{i}R{j}: '+' | '.join(c.text for c in row.cells))
(QA/'v25_content.txt').write_text('\n'.join(lines)+'\n',encoding='utf8')
# Reuse the canonical render route already validated for this document.
old=ROOT.parent/'model_argument_v24_20261013'
for source,target in [('export_v24_word.ps1','export_v25_word.ps1'),('render_v24.py','render_v25.py')]:
    text=(old/source).read_text(encoding='utf8').replace(SOURCE.name,OUT.name).replace('v24_word','v25_word').replace('render_v24','render_v25').replace('v24_page_map','v25_page_map')
    (ROOT/target).write_text(text,encoding='utf8')
print(json.dumps({k:v for k,v in audit.items() if k!='changes'},ensure_ascii=False))
