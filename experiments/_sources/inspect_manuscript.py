from pathlib import Path
from docx import Document
root=Path(__file__).resolve().parent
doc=Document(root.parent/'paper_v15_20261010_中文修订_三线表.docx')
text='\n'.join(f'{i:04d}\t{p.style.name}\t{p.text}' for i,p in enumerate(doc.paragraphs))
(root/'manuscript_paragraphs.txt').write_text(text,encoding='utf-8')
for i,p in enumerate(doc.paragraphs):
    if any(word in p.text for word in ['不确定','功能','敏感性','局限','补充方法','图 S13','图S13']):
        print(i,p.style.name,p.text)
