from pathlib import Path
import json,hashlib
import importlib.util,sys,os
from pypdf import PdfReader
import pypdfium2 as pdfium
from PIL import Image,ImageOps,ImageDraw

ROOT=Path(__file__).resolve().parent
QA=ROOT/'manuscript_QA'
doc=ROOT.parent.parent/'paper_v50_20261022_中文修订_MEL5残基接触网络与图4精简_三线表.docx'
pdf=QA/'v50_word.pdf'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
# Packaged renderer is called with the already exported read-only Word PDF.
# Windows bundle has no LibreOffice. No system LibreOffice is used.
script=Path(r'C:\Users\DELL\.codex\plugins\cache\openai-primary-runtime\documents\26.1020.11127\skills\documents\render_docx.py')
spec=importlib.util.spec_from_file_location('packaged_docx_renderer',script)
renderer=importlib.util.module_from_spec(spec);spec.loader.exec_module(renderer)
renderer.convert_to_pdf=lambda *args,**kwargs:(str(pdf),'Owned invisible Word COM read-only PDF export')
os.environ['PATH']=r'C:\Users\DELL\.cache\codex-runtimes\codex-primary-runtime\dependencies\native\poppler\Library\bin'+os.pathsep+os.environ.get('PATH','')
sys.argv=[str(script),str(doc),'--output_dir',str(QA/'pages'),'--dpi','120']
renderer.main()
reader=PdfReader(pdf)
page_map=[]
for index,page in enumerate(reader.pages,1):
    text=page.extract_text()
    page_map.append(dict(page=index,chars=len(text),start=text[:160],end=text[-160:],MEL5_case=any(token in text for token in ['3.4 MEL5','Q72','图4：','Q72–β'])) )
(QA/'v50_page_map.json').write_text(json.dumps(page_map,ensure_ascii=False,indent=2),encoding='utf8')
imgs=sorted((QA/'pages').glob('page-*.png'),key=lambda p:int(p.stem.split('-')[-1]))
assert len(imgs)==len(reader.pages)
for begin in range(0,len(imgs),8):
    sheet=Image.new('RGB',(4*310,2*460),'#DBDBDB');draw=ImageDraw.Draw(sheet)
    for offset,path in enumerate(imgs[begin:begin+8]):
        im=Image.open(path).convert('RGB');im.thumbnail((300,424))
        x=(offset%4)*310+(310-im.width)//2;y=(offset//4)*460+24
        sheet.paste(im,(x,y));draw.text(((offset%4)*310+10,(offset//4)*460+6),f'Page {begin+offset+1}',fill='black')
    sheet.save(QA/f'contact_sheet_{begin+1}_{min(begin+8,len(imgs))}.png')
(QA/'render_source_binding.json').write_text(json.dumps(dict(document=str(doc),document_sha256=sha(doc),pdf=str(pdf),pdf_sha256=sha(pdf),page_count=len(reader.pages),renderer='Packaged render_docx.py using owned Word read-only PDF export',manuscript_saved_by_renderer=False),ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps(dict(page_count=len(reader.pages),case_pages=[p['page'] for p in page_map if p['MEL5_case']],contact_sheets=[p.name for p in QA.glob('contact_sheet*.png')]),ensure_ascii=False))
