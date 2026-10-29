from pathlib import Path
from zipfile import ZipFile
import json
from PIL import Image, ImageChops
from lxml import etree

root = Path(__file__).resolve().parent
source = root.parent / 'paper_v14_20261010_中文修订.docx'
output = root.parent / 'paper_v15_20261010_中文修订_三线表.docx'
ns = {'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main', 'm':'http://schemas.openxmlformats.org/officeDocument/2006/math'}
with ZipFile(source) as before, ZipFile(output) as after:
    a = etree.fromstring(before.read('word/document.xml'))
    b = etree.fromstring(after.read('word/document.xml'))
    assert a.xpath('//w:t/text()', namespaces=ns) == b.xpath('//w:t/text()', namespaces=ns)
    assert [etree.tostring(x) for x in a.xpath('//m:oMathPara', namespaces=ns)] == [etree.tostring(x) for x in b.xpath('//m:oMathPara', namespaces=ns)]
    assert before.read('word/comments.xml') == after.read('word/comments.xml')
    media = [x for x in before.namelist() if x.startswith('word/media/')]
    assert all(before.read(x) == after.read(x) for x in media)
    assert len(b.xpath('//w:tbl', namespaces=ns)) == 5

old = Path(r'C:\Users\DELL\AppData\Local\Temp\flexppi_revision_20261010\render')
new = Path(r'C:\Users\DELL\AppData\Local\Temp\flexppi_completion_20261010\render')
changed = []
for page in range(1,34):
    x=Image.open(old/f'page-{page}.png').convert('RGB')
    y=Image.open(new/f'page-{page}.png').convert('RGB')
    assert x.size == y.size
    if ImageChops.difference(x,y).getbbox():
        changed.append(page)
audit = {'pages':33,'visually_changed_pages':changed,'text_preserved':True,'display_math_preserved':True,'comments_preserved':True,'all_images_preserved':True,'tables':5}
(root/'three_line_tables_verification.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(audit,ensure_ascii=False))
