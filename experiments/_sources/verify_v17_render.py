from pathlib import Path
from PIL import Image
from pypdf import PdfReader
from docx import Document
from docx.oxml.ns import qn
import hashlib,json,zipfile
tmp=Path(r'C:\Users\DELL\AppData\Local\Temp\flexppi_completion_20261010')
root=Path(__file__).resolve().parent
document=root.parent/'paper_v17_20261010_中文修订_基线与标签敏感性_三线表.docx'
def key(path):
    image=Image.open(path).convert('RGB')
    return hashlib.sha256(image.crop((0,0,image.width,image.height-120)).tobytes()).hexdigest()
old_count=len(PdfReader(tmp/'v16_word.pdf').pages)
old={key(tmp/'render_v16'/f'page-{n}.png'):n for n in range(1,old_count+1)}
count=len(PdfReader(tmp/'v17_word.pdf').pages)
pages=[dict(page=n,previous_identical_content_page=old.get(key(tmp/'render_v17'/f'page-{n}.png'))) for n in range(1,count+1)]
tables=[]
for index,table in enumerate(Document(document).tables):
    for row_index,row in enumerate(table.rows):
        for cell in row.cells:
            borders=cell._tc.tcPr.find(qn('w:tcBorders'))
            assert borders is not None
            for side in ['top','bottom','left','right','insideH','insideV']:
                expected=(row_index==0 and side=='top') or (row_index==len(table.rows)-1 and side=='bottom') or (row_index==0 and side=='bottom')
                border=borders.find(qn('w:'+side))
                assert border is not None and border.get(qn('w:val'))==('single' if expected else 'nil')
                if expected:assert border.get(qn('w:sz'))==('6' if row_index==0 and side=='bottom' else '12')
            assert cell._tc.tcPr.find(qn('w:shd')) is None
    tables.append(dict(index=index+1,rows=len(table.rows),columns=len(table.columns),true_three_line=True))
with zipfile.ZipFile(root.parent/'paper_v16_20261010_中文修订_已完成分析_三线表.docx') as oldzip,zipfile.ZipFile(document) as newzip:
    comments_preserved=oldzip.read('word/comments.xml')==newzip.read('word/comments.xml')
    old_media_preserved=all(oldzip.read(n)==newzip.read(n) for n in oldzip.namelist() if n.startswith('word/media/'))
assert comments_preserved and old_media_preserved and len(tables)==8
audit=dict(document=str(document),sha256=hashlib.sha256(document.read_bytes()).hexdigest(),pages=count,
    tables=tables,comments_preserved=comments_preserved,old_media_preserved=old_media_preserved,
    pixel_comparison_to_previously_inspected_v16=pages,changed_pages_requiring_visual_review=[p['page'] for p in pages if p['previous_identical_content_page'] is None])
(tmp/'v17_layout_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps(audit,ensure_ascii=False))
