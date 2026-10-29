"""Run the bundled DOCX renderer against an owned invisible Word PDF export."""
from pathlib import Path
import argparse
import importlib.util
import json
import os
import re
import sys
from pypdf import PdfReader
from PIL import Image, ImageDraw

ROOT = Path(__file__).resolve().parent
RENDERER = Path('C:/Users/DELL/.codex/plugins/cache/openai-primary-runtime/documents/26.1022.11436/skills/documents/render_docx.py')
POPLER = Path('C:/Users/DELL/.cache/codex-runtimes/codex-primary-runtime/dependencies/native/poppler/Library/bin')

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('document', type=Path)
    parser.add_argument('pdf', type=Path)
    args = parser.parse_args()
    document = args.document.resolve()
    pdf = args.pdf.resolve()
    output = ROOT / 'document_QA' / 'pages'
    output.mkdir(exist_ok=True, parents=True)
    os.environ['PATH'] = str(POPLER) + os.pathsep + os.environ['PATH']
    spec = importlib.util.spec_from_file_location('docx_renderer', RENDERER)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.convert_to_pdf = lambda *a, **kw: (str(pdf), 'Owned invisible Word read-only PDF export')
    sys.argv = ['render_docx.py', str(document), '--output_dir', str(output), '--dpi', '130']
    module.main()
    doc = PdfReader(pdf).pages
    page_map = []
    for i, page in enumerate(doc):
        text = page.extract_text()
        page_map.append({'page': i + 1, 'text': text, 'images': len(page.images),
                         'figures': [n for n in ['1', '2', '3', '4', '5', 'S1', 'S2']
                                     if re.search(r'图\s*' + n + r'\s*[：:]', text)]})
    (ROOT / 'document_QA' / 'page_map.json').write_text(json.dumps(page_map, ensure_ascii=False, indent=2), encoding='utf-8')
    for start in range(0, len(doc), 8):
        sheet = Image.new('RGB', (1680, 1230), 'white')
        draw = ImageDraw.Draw(sheet)
        for index in range(start, min(start + 8, len(doc))):
            page = Image.open(output / ('page-' + str(index + 1) + '.png')).convert('RGB')
            page.thumbnail((400, 570))
            slot = index - start
            x, y = (slot % 4) * 420 + 10, (slot // 4) * 615 + 30
            sheet.paste(page, (x, y))
            draw.text((x, y - 20), 'Page ' + str(index + 1), fill='black')
        sheet.save(ROOT / 'document_QA' / ('contact-' + str(start + 1) + '.png'))
    print('Rendered', len(doc), 'pages; figures:', [(p['page'], p['figures']) for p in page_map if p['figures']])

if __name__ == '__main__':
    main()
