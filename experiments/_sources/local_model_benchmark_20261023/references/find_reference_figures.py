from pathlib import Path
from pypdf import PdfReader
ROOT = Path(__file__).resolve().parent
for name in ['bioemu_accepted.pdf', 'scgpt.pdf']:
    p = ROOT / name
    if not p.exists():
        continue
    reader = PdfReader(p)
    print(name, 'pages=', len(reader.pages), flush=True)
    for i, page in enumerate(reader.pages[:32]):
        t = page.extract_text() or ''
        if any(s.lower() in t.lower() for s in ['samples functionally distinct', 'results of multi-batch', 'cell type-annotation results', 'emulates equilibrium distributions']):
            print('FOUND', name, i+1, flush=True)
            (ROOT / (p.stem + '_page_' + str(i+1) + '.txt')).write_text(t, encoding='utf-8')
