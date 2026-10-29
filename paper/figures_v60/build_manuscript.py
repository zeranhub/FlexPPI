"""Replace only Figure 1's PNG part in a byte-preserving copy of v59."""
from pathlib import Path
from zipfile import ZipFile
import hashlib
import json
from PIL import Image

ROOT = Path(__file__).resolve().parent
BASE = Path('E:/Project/FlexPPI/paper/revision_v14_20261010/paper_v59_20261023_中文修订_图2至5整体重绘.docx')
OUTPUT = ROOT.parent / 'paper_v60_20261023_中文修订_图1配色轻调.docx'
FIGURE = ROOT / '01_模型框架' / 'fig1_framework.png'

def main():
    assert Image.open(FIGURE).size == (3802,2288)
    data = FIGURE.read_bytes()
    with ZipFile(BASE) as source:
        old = {name:source.read(name) for name in source.namelist()}
        assert 'word/media/image1.png' in old
        with ZipFile(OUTPUT,'w') as target:
            for info in source.infolist():
                target.writestr(info,data if info.filename=='word/media/image1.png' else old[info.filename])
    with ZipFile(OUTPUT) as target:
        changed=[name for name,content in old.items() if target.read(name)!=content]
        assert changed == ['word/media/image1.png'], changed
    audit={'status':'source_preservation_passed_render_pending','baseline':str(BASE),
        'output':str(OUTPUT.resolve()),'changed_zip_members':changed,
        'all_text_equations_styles_and_other_figures_byte_identical':True,
        'document_xml_byte_identical':True,
        'baseline_sha256':hashlib.sha256(BASE.read_bytes()).hexdigest(),
        'document_sha256':hashlib.sha256(OUTPUT.read_bytes()).hexdigest(),
        'figure1_sha256':hashlib.sha256(data).hexdigest()}
    (ROOT/'document_QA').mkdir(exist_ok=True)
    (ROOT/'document_QA/document_edit_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
    print(OUTPUT)

if __name__=='__main__':
    main()
