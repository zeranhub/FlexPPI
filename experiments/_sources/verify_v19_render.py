from pathlib import Path
from collections import Counter
from PIL import Image
from pypdf import PdfReader
from docx import Document
from docx.oxml.ns import qn
from lxml import etree
import hashlib, json, zipfile

ROOT = Path(__file__).resolve().parent
TMP = Path(r'C:\Users\DELL\AppData\Local\Temp\flexppi_completion_20261010')
old_doc = ROOT.parent / 'paper_v18_20261010_中文修订_复用KIX单轨迹_三线表.docx'
new_doc = ROOT.parent / 'paper_v19_20261012_中文修订_补充实验完成_三线表.docx'

def body_hash(path):
    im = Image.open(path).convert('RGB')
    return hashlib.sha256(im.crop((0, 0, im.width, im.height - 120)).tobytes()).hexdigest()

old_count = len(PdfReader(TMP / 'v18_word.pdf').pages)
old_pages = {body_hash(TMP / 'render_v18' / f'page-{n}.png'): n for n in range(1, old_count + 1)}
pdf = PdfReader(TMP / 'v19_word.pdf')
pages = [{'page': n, 'previous_identical_content_page': old_pages.get(body_hash(TMP / 'render_v19' / f'page-{n}.png'))} for n in range(1, len(pdf.pages) + 1)]
blank_pages = [i+1 for i, p in enumerate(pdf.pages) if len(p.extract_text().strip()) < 10]
assert not blank_pages, blank_pages

tables = []
for index, table in enumerate(Document(new_doc).tables):
    for i, row in enumerate(table.rows):
        for cell in row.cells:
            borders = cell._tc.tcPr.find(qn('w:tcBorders'))
            assert borders is not None
            for side in ['top', 'bottom', 'left', 'right', 'insideH', 'insideV']:
                expected = (i == 0 and side in ['top', 'bottom']) or (i == len(table.rows)-1 and side == 'bottom')
                border = borders.find(qn('w:' + side))
                assert border is not None and border.get(qn('w:val')) == ('single' if expected else 'nil')
                if expected:
                    assert border.get(qn('w:sz')) == ('6' if i == 0 and side == 'bottom' else '12')
            assert cell._tc.tcPr.find(qn('w:shd')) is None
    tables.append({'table': index+1, 'rows': len(table.rows), 'columns': len(table.columns), 'three_line': True})
assert len(tables) == 14

with zipfile.ZipFile(old_doc) as old, zipfile.ZipFile(new_doc) as new:
    assert old.read('word/comments.xml') == new.read('word/comments.xml')
    assert all(old.read(n) == new.read(n) for n in old.namelist() if n.startswith('word/media/') and not n.endswith('/'))
    def anchors(z):
        xml = etree.fromstring(z.read('word/document.xml'))
        return Counter((e.tag, e.get(qn('w:id'))) for e in xml.iter() if e.tag in {qn('w:commentRangeStart'), qn('w:commentRangeEnd'), qn('w:commentReference')})
    assert anchors(old) == anchors(new)
    def equations(z):
        xml = etree.fromstring(z.read('word/document.xml'))
        return [etree.tostring(e, method='c14n') for e in xml.iter('{http://schemas.openxmlformats.org/officeDocument/2006/math}oMath')]
    assert equations(old) == equations(new)

audit = {'document': str(new_doc), 'sha256': hashlib.sha256(new_doc.read_bytes()).hexdigest(),
         'pages': len(pdf.pages), 'tables': tables, 'blank_pages': blank_pages,
         'comments_and_anchors_preserved': True, 'original_images_preserved': True, 'native_equations_preserved':True,
         'pixel_comparison_to_previously_inspected_v18': pages,
         'changed_pages_requiring_visual_review': [p['page'] for p in pages if p['previous_identical_content_page'] is None],
         'visual_review': 'pending'}
(ROOT / 'completed_results_20261012/v19_layout_audit.json').write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding='utf8')
print(json.dumps({k: audit[k] for k in ['pages', 'sha256', 'blank_pages', 'changed_pages_requiring_visual_review']}, ensure_ascii=False))
