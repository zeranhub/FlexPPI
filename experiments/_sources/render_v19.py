from pathlib import Path
import importlib.util, os, sys, json
from pypdf import PdfReader
tmp = Path(r'C:\Users\DELL\AppData\Local\Temp\flexppi_completion_20261010')
os.environ['PATH'] = r'C:\Users\DELL\.cache\codex-runtimes\codex-primary-runtime\dependencies\native\poppler\Library\bin' + os.pathsep + os.environ['PATH']
spec = importlib.util.spec_from_file_location('docx_renderer', r'C:\Users\DELL\.codex\plugins\cache\openai-primary-runtime\documents\26.1011.11113\skills\documents\render_docx.py')
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)
mod.convert_to_pdf = lambda *a, **kw: (str(tmp / 'v19_word.pdf'), 'Rendered with owned invisible Word COM application')
sys.argv = ['render_v19.py', str(Path(__file__).resolve().parent.parent / 'paper_v19_20261012_中文修订_补充实验完成_三线表.docx'), '--output_dir', str(tmp / 'render_v19'), '--dpi', '150']
mod.main()
pages = [{'page': i+1, 'chars': len(p.extract_text()), 'start': p.extract_text()[:140], 'end': p.extract_text()[-140:]} for i, p in enumerate(PdfReader(tmp / 'v19_word.pdf').pages)]
(tmp / 'v19_page_map.json').write_text(json.dumps(pages, ensure_ascii=False, indent=2), encoding='utf8')
for p in pages:
    print(p['page'], p['chars'], p['start'][:65].replace('\n', ' '))
