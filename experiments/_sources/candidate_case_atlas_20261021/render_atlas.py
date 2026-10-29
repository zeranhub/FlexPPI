from pathlib import Path
import importlib.util,os,sys,json,hashlib,re
from pypdf import PdfReader
ROOT=Path(__file__).resolve().parent
TMP=Path(r'C:\Users\DELL\AppData\Local\Temp\flexppi_case_atlas_20261021')
os.environ['PATH']=r'C:\Users\DELL\.cache\codex-runtimes\codex-primary-runtime\dependencies\native\poppler\Library\bin'+os.pathsep+os.environ['PATH']
source=r'C:\Users\DELL\.codex\plugins\cache\openai-primary-runtime\documents\26.1020.11127\skills\documents\render_docx.py'
spec=importlib.util.spec_from_file_location('docx_renderer',source)
mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
mod.convert_to_pdf=lambda *a,**kw:(str(TMP/'atlas_word.pdf'),'Owned invisible Word read-only PDF export')
doc=ROOT.parent.parent/'FlexPPI_32个潜在案例图册_背景证据与局限_20261021.docx'
iteration=1
while (TMP/f'render_{iteration}').exists():iteration+=1
render_directory=TMP/f'render_{iteration}'
sys.argv=['render_atlas.py',str(doc),'--output_dir',str(render_directory),'--dpi','150']
mod.main()
reader=PdfReader(TMP/'atlas_word.pdf')
pages=[]
for i,p in enumerate(reader.pages,1):
    text=p.extract_text()
    pages.append({'page':i,'chars':len(text),'start':text[:240],'end':text[-240:],
                  'case_mentions':sorted(set(re.findall(r'C\d{2}',text))),
                  'text':text})
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
binding={'status':'rendered_pending_visual_QA','document_sha256':sha(doc),'PDF_SHA256':sha(TMP/'atlas_word.pdf'),
         'pages':len(pages),'renderer':source,'PDF_method':'Owned invisible Word COM, read-only, no save',
         'render_directory':str(render_directory),
         'PNG_SHA256':{str(i):sha(render_directory/f'page-{i}.png') for i in range(1,len(pages)+1)}}
(ROOT/'document_QA/page_map.json').write_text(json.dumps(pages,ensure_ascii=False,indent=2),encoding='utf8')
(ROOT/'document_QA/render_binding.json').write_text(json.dumps(binding,ensure_ascii=False,indent=2),encoding='utf8')
for p in pages:print(p['page'],p['chars'],p['start'][:90].replace('\n',' '))
