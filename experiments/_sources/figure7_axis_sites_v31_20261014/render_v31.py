from pathlib import Path
import importlib.util,os,sys,json
from pypdf import PdfReader
ROOT=Path(__file__).resolve().parent
tmp=Path(r'C:\Users\DELL\AppData\Local\Temp\flexppi_completion_20261010')
os.environ['PATH']=r'C:\Users\DELL\.cache\codex-runtimes\codex-primary-runtime\dependencies\native\poppler\Library\bin'+os.pathsep+os.environ['PATH']
spec=importlib.util.spec_from_file_location('docx_renderer',r'C:\Users\DELL\.codex\plugins\cache\openai-primary-runtime\documents\26.1011.11113\skills\documents\render_docx.py')
mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
mod.convert_to_pdf=lambda *a,**kw:(str(tmp/'v31_word.pdf'),'Owned invisible Word COM render')
sys.argv=['render_v31.py',str(ROOT.parent.parent/'paper_v31_20261014_中文修订_图7简化与六谷氨酰胺轴标注_三线表.docx'),'--output_dir',str(tmp/'render_v31'),'--dpi','150']
mod.main()
pages=[{'page':i+1,'chars':len(p.extract_text()),'start':p.extract_text()[:140],'end':p.extract_text()[-140:]} for i,p in enumerate(PdfReader(tmp/'v31_word.pdf').pages)]
(ROOT/'manuscript_QA/v31_page_map.json').write_text(json.dumps(pages,ensure_ascii=False,indent=2),encoding='utf8')
for p in pages:print(p['page'],p['chars'],p['start'][:70].replace('\n',' '))
