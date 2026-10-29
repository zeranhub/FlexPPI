"""Preserve manuscript content and native comments while applying three-line tables."""
from pathlib import Path
import hashlib
import json
from docx import Document
from docx.oxml import OxmlElement
from docx.oxml.ns import qn

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT.parent / 'paper_v14_20261010_中文修订.docx'
OUTPUT = ROOT.parent / 'paper_v15_20261010_中文修订_三线表.docx'


def border(parent, tag, value, size=0):
    element = OxmlElement('w:' + tag)
    element.set(qn('w:val'), value)
    if size:
        element.set(qn('w:sz'), str(size))
        element.set(qn('w:color'), '000000')
        element.set(qn('w:space'), '0')
    parent.append(element)


document = Document(SOURCE)
counts = []
for index, table in enumerate(document.tables, 1):
    properties = table._tbl.tblPr
    for item in list(properties):
        if item.tag in {qn('w:tblBorders'), qn('w:shd')}:
            properties.remove(item)
    borders = OxmlElement('w:tblBorders')
    for side in ['top', 'left', 'bottom', 'right', 'insideH', 'insideV']:
        border(borders, side, 'single' if side in {'top', 'bottom'} else 'nil',
               12 if side in {'top', 'bottom'} else 0)
    properties.append(borders)
    for row_number, row in enumerate(table.rows):
        for cell in row.cells:
            cell_properties = cell._tc.get_or_add_tcPr()
            for item in list(cell_properties):
                if item.tag in {qn('w:tcBorders'), qn('w:shd')}:
                    cell_properties.remove(item)
            cell_borders = OxmlElement('w:tcBorders')
            for side in ['top', 'left', 'bottom', 'right', 'insideH', 'insideV']:
                outer = (row_number == 0 and side == 'top') or (row_number == len(table.rows) - 1 and side == 'bottom')
                header = row_number == 0 and side == 'bottom'
                border(cell_borders, side, 'single' if outer or header else 'nil',
                       12 if outer else (6 if header else 0))
            cell_properties.append(cell_borders)
    counts.append({'table': index, 'rows': len(table.rows), 'columns': len(table.columns)})
document.save(OUTPUT)
audit = {'input': str(SOURCE), 'output': str(OUTPUT), 'tables': counts,
         'line_width_pt': {'top': 1.5, 'header_bottom': 0.75, 'bottom': 1.5},
         'vertical_lines': False, 'internal_body_lines': False,
         'input_sha256': hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
         'output_sha256': hashlib.sha256(OUTPUT.read_bytes()).hexdigest()}
(ROOT / 'three_line_tables_audit.json').write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding='utf-8')
print(json.dumps(audit, ensure_ascii=False))
