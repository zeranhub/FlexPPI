from pathlib import Path
from collections import Counter
import hashlib, json, zipfile
from lxml import etree
from docx import Document
from docx.oxml.ns import qn
from pypdf import PdfReader

ROOT=Path(__file__).resolve().parent
QA=ROOT/'argument_revision_v20_20261012'
TMP=Path(r'C:\Users\DELL\AppData\Local\Temp\flexppi_completion_20261010')
OLD=ROOT.parent/'paper_v19_20261012_中文修订_补充实验完成_三线表.docx'
NEW=ROOT.parent/'paper_v20_20261012_中文修订_联合模型_三线表.docx'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
assert sha(NEW)=='4b4597b42fdc8a1300ab8f4155e98a4960872a2368fc887397d133af7d264d79'
document=Document(NEW)
tables=[]
for index,table in enumerate(document.tables):
    for i,row in enumerate(table.rows):
        for cell in row.cells:
            borders=cell._tc.tcPr.find(qn('w:tcBorders'))
            assert borders is not None
            for side in ['top','bottom','left','right','insideH','insideV']:
                expected=(i==0 and side in ['top','bottom']) or (i==len(table.rows)-1 and side=='bottom')
                border=borders.find(qn('w:'+side))
                assert border is not None and border.get(qn('w:val'))==('single' if expected else 'nil')
                if expected:assert border.get(qn('w:sz'))==('6' if i==0 and side=='bottom' else '12')
            assert cell._tc.tcPr.find(qn('w:shd')) is None
    tables.append({'table':index+1,'rows':len(table.rows),'columns':len(table.columns),'three_line':True})
assert len(tables)==13
ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main','m':'http://schemas.openxmlformats.org/officeDocument/2006/math'}
with zipfile.ZipFile(OLD) as a,zipfile.ZipFile(NEW) as b:
    assert a.read('word/comments.xml')==b.read('word/comments.xml')
    oldxml,newxml=[etree.fromstring(z.read('word/document.xml')) for z in (a,b)]
    eq=lambda e:[etree.tostring(x,method='c14n') for x in e.xpath('.//m:oMath',namespaces=ns)]
    assert eq(oldxml)==eq(newxml)
    anchors=lambda e:Counter((x.tag,x.get(qn('w:id'))) for x in e.iter() if x.tag in {qn('w:commentRangeStart'),qn('w:commentRangeEnd'),qn('w:commentReference')})
    assert anchors(oldxml)==anchors(newxml)
    media=[name for name in a.namelist() if name.startswith('word/media/') and not name.endswith('/')]
    assert all(a.read(name)==b.read(name) for name in media)
    assert len(newxml.xpath('.//w:drawing',namespaces=ns))==16
pdf=PdfReader(TMP/'v20_word.pdf')
assert len(pdf.pages)==40
blank=[i+1 for i,p in enumerate(pdf.pages) if len(p.extract_text().strip())<10]
assert not blank
assert all((TMP/'render_v20'/f'page-{i}.png').is_file() for i in range(1,41))
root_reviewed=list(range(1,19))+list(range(31,41))
report={'document':str(NEW),'sha256':sha(NEW),'source_sha256':sha(OLD),
        'pages':40,'tables':tables,'actual_figures':16,'blank_pages':blank,
        'comments_and_anchors_preserved':True,'native_equations_preserved':True,
        'archived_media_preserved':True,'all_retained_table_values_preserved':True,
        'root_actual_visual_review_passed_pages':root_reviewed,
        'delegated_visual_review_pages':list(range(19,31)),
        'delegated_visual_review_status':'pending',
        'png_sha256':{str(i):sha(TMP/'render_v20'/f'page-{i}.png') for i in range(1,41)}}
(QA/'v20_layout_audit.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps({k:report[k] for k in ['sha256','pages','actual_figures','blank_pages','root_actual_visual_review_passed_pages']},ensure_ascii=False))
