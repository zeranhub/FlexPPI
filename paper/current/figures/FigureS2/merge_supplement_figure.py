"""Vector-only composition of frozen source heatmap and internal PR panels.

This is a cosmetic merge: every heatmap value, curve vertex, axis mapping,
legend entry, and colour is carried over from the original PDF/SVG assets.
Only panel letters and empty exterior margins change. No data are recalculated.
"""
from pathlib import Path
from io import BytesIO
from copy import deepcopy
import hashlib
import json
import xml.etree.ElementTree as ET

from pypdf import PdfReader, PdfWriter, Transformation
from pypdf.generic import ContentStream, NameObject, FloatObject, TextStringObject
from reportlab.pdfgen import canvas
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
import pypdfium2 as pdfium
import pdfplumber

OUT = Path(__file__).resolve().parent
HEAT = OUT / 'source_panels/source_AP'
PRS = OUT / 'source_panels/internal_PR'
PT = 72 / 25.4
W = 183 * PT
HEAT_H = 122.5  # Original content bottom is 119.93 pt; empty margin trimmed.
PR_H = 237.5    # Original content bottom is 233.84 pt; empty margin trimmed.
GAP = 2 * PT
H = HEAT_H + GAP + PR_H
X_LEFT = 0.5 * PT
X_RIGHT = 94.5 * PT
SV = 'http://www.w3.org/2000/svg'
XL = 'http://www.w3.org/1999/xlink'
ET.register_namespace('', SV)
ET.register_namespace('xlink', XL)

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def relabel_pdf_pr(page):
    stream = ContentStream(page.get_contents(), page.pdf)
    replacements = [('a', 'b', X_LEFT), ('b', 'c', X_RIGHT)]
    changed = []
    for old, new, x in replacements:
        old_x = 39.6850394 if old == 'a' else 306.141732
        matches = [i for i, (args, op) in enumerate(stream.operations)
                   if op == b'TJ' and len(args) == 1 and list(args[0]) == [old]
                   and stream.operations[i - 4][1] == b'cm'
                   and abs(float(stream.operations[i - 4][0][4]) - old_x) < 0.001]
        assert len(matches) == 1, (old, matches)
        idx = matches[0]
        args, _ = stream.operations[idx]
        # Existing embedded bold font is subsetted for its source a/b glyphs.
        # Remove these exact text objects; add all three labels with Arial below.
        args[0][0] = TextStringObject('')
        # Matplotlib's panel text uses its own q / cm / BT text block.
        cm_args, cm_op = stream.operations[idx - 4]
        assert cm_op == b'cm' and abs(float(cm_args[4]) - (39.6850394 if old == 'a' else 306.141732)) < 0.001
        cm_args[4] = FloatObject(x)
        changed.append({'old': old, 'new': new, 'x_pt': x})
    page[NameObject('/Contents')] = stream
    return changed

def write_pdf():
    heat = PdfReader(HEAT.with_suffix('.pdf')).pages[0]
    pr = PdfReader(PRS.with_suffix('.pdf')).pages[0]
    pr_labels = relabel_pdf_pr(pr)
    heat_original_h = float(heat.mediabox.height)
    pr_original_h = float(pr.mediabox.height)
    heat.cropbox.lower_left = (0, heat_original_h - HEAT_H)
    pr.cropbox.lower_left = (0, pr_original_h - PR_H)
    writer = PdfWriter()
    page = writer.add_blank_page(width=W, height=H)
    page.merge_transformed_page(heat, Transformation().translate(0, H - heat_original_h))
    page.merge_transformed_page(pr, Transformation().translate(0, -(pr_original_h - PR_H)))
    # The new heatmap letter matches the existing Arial bold panel typography.
    pdfmetrics.registerFont(TTFont('FigureArialBold', r'C:\Windows\Fonts\arialbd.ttf'))
    buf = BytesIO()
    c = canvas.Canvas(buf, pagesize=(W, H))
    c.setFont('FigureArialBold', 11)
    c.drawString(X_LEFT, H - 10.0, 'a')
    c.drawString(X_LEFT, H - HEAT_H - GAP - 14.837093, 'b')
    c.drawString(X_RIGHT, H - HEAT_H - GAP - 14.837093, 'c')
    c.save()
    buf.seek(0)
    page.merge_page(PdfReader(buf).pages[0])
    writer.add_metadata({'/Title': 'Supplementary Figure S1: source-stratified AP and internal PR',
                         '/Subject': 'Frozen native vector panels; layout and panel letters only'})
    path = OUT / 'figS1_merged.pdf'
    writer.write(path)
    return pr_labels

def prefix_ids(root, prefix):
    ids = {el.attrib['id']: prefix + el.attrib['id'] for el in root.iter() if 'id' in el.attrib}
    for el in root.iter():
        for key, val in list(el.attrib.items()):
            if key == 'id':
                el.set(key, ids[val])
            else:
                for old, new in ids.items():
                    val = val.replace('url(#' + old + ')', 'url(#' + new + ')')
                    if val == '#' + old:
                        val = '#' + new
                el.set(key, val)

def write_svg():
    root = ET.Element(f'{{{SV}}}svg', {'width': '183mm', 'height': f'{H / PT:.8f}mm',
                                     'viewBox': f'0 0 {W:.9f} {H:.9f}', 'version': '1.1'})
    defs = ET.SubElement(root, f'{{{SV}}}defs')
    for ident, height in [('clip_heat_outer', HEAT_H), ('clip_pr_outer', PR_H)]:
        clip = ET.SubElement(defs, f'{{{SV}}}clipPath', {'id': ident})
        ET.SubElement(clip, f'{{{SV}}}rect', {'width': str(W), 'height': str(height)})
    for source, prefix, y, clip_id in [(HEAT, 'heat_', 0, 'clip_heat_outer'),
                                        (PRS, 'pr_', HEAT_H + GAP, 'clip_pr_outer')]:
        original = ET.parse(source.with_suffix('.svg')).getroot()
        if source == PRS:
            texts = list(original.iter(f'{{{SV}}}text'))
            labels = [el for el in texts if el.text in ('a', 'b')]
            assert len(labels) == 2
            # Process the source a before source b so the two targets stay distinct.
            for element, new, x in zip(labels, ['b', 'c'], [X_LEFT, X_RIGHT]):
                element.text = new
                element.set('x', str(x))
                element.attrib.pop('transform', None)  # Existing transform is rotate(-0).
        prefix_ids(original, prefix)
        translated = ET.SubElement(root, f'{{{SV}}}g', {'transform': f'translate(0 {y})'})
        clipped = ET.SubElement(translated, f'{{{SV}}}g', {'clip-path': f'url(#{clip_id})'})
        for child in original:
            if child.tag != f'{{{SV}}}metadata':
                clipped.append(deepcopy(child))
    text = ET.SubElement(root, f'{{{SV}}}text', {'x': str(X_LEFT), 'y': '10',
                 'style': "font: 700 11px 'Arial', 'Helvetica', 'Liberation Sans', sans-serif; text-anchor: start"})
    text.text = 'a'
    ET.ElementTree(root).write(OUT / 'figS1_merged.svg', encoding='utf-8', xml_declaration=True)

def validate_and_preview(labels):
    target = OUT / 'figS1_merged.pdf'
    pdf = pdfium.PdfDocument(target)
    image = pdf[0].render(scale=300/72).to_pil()
    image.save(OUT / 'figS1_merged.png', dpi=(300, 300))
    thumb = image.copy()
    thumb.thumbnail((1500, 1500))
    thumb.save(OUT / 'figS1_merged_preview.png')
    expected = []
    for source in [HEAT, PRS]:
        with pdfplumber.open(source.with_suffix('.pdf')) as d:
            expected += [w['text'] for w in d.pages[0].extract_words() if w['text'] not in ('a', 'b')]
    with pdfplumber.open(target) as d:
        actual = [w['text'] for w in d.pages[0].extract_words() if w['text'] not in ('a', 'b', 'c')]
    from collections import Counter
    assert Counter(actual) == Counter(expected), Counter(expected) - Counter(actual)
    path_operators = {b'm', b'l', b'c', b'v', b'y', b'h'}
    expected_paths = Counter()
    for source in [HEAT, PRS]:
        page = PdfReader(source.with_suffix('.pdf')).pages[0]
        expected_paths.update(op.decode() for _, op in ContentStream(page.get_contents(), page.pdf).operations if op in path_operators)
    page = PdfReader(target).pages[0]
    actual_paths = Counter(op.decode() for _, op in ContentStream(page.get_contents(), page.pdf).operations if op in path_operators)
    assert expected_paths == actual_paths, (expected_paths, actual_paths)
    svg = ET.parse(OUT / 'figS1_merged.svg').getroot()
    ids = [el.attrib['id'] for el in svg.iter() if 'id' in el.attrib]
    assert len(ids) == len(set(ids))
    assert Counter(el.text for el in svg.iter(f'{{{SV}}}text') if el.text in ['a','b','c']) == Counter(['a','b','c'])
    report = {'width_mm': 183.0, 'height_mm': H / PT, 'dpi': 300, 'png_pixels': image.size,
              'composition': 'a: source-stratified AP; b: internal RMSF PR; c: internal Hchi PR',
              'data_changes': 'none; no recalculation, no curve edits, no axis edits',
              'text_validation': 'all original non-panel-letter PDF words preserved with identical multiplicity',
              'vector_path_validation': dict(actual_paths),
              'vector_path_validation_status': 'all source path operator counts preserved exactly',
              'svg_validation': 'all IDs unique; exactly one a, b, and c panel label',
              'margin_trimming': 'only empty lower exterior whitespace; original curves and all text inside clips',
              'labels': labels,
              'sources': [{'path': str(p), 'sha256': sha(p)} for stem in [HEAT, PRS] for p in [stem.with_suffix('.pdf'), stem.with_suffix('.svg')]],
              'outputs': [{'path': str(p), 'sha256': sha(p)} for p in [target, OUT / 'figS1_merged.svg', OUT / 'figS1_merged.png']]}
    (OUT / 'figS1_merged_audit.json').write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps(report, ensure_ascii=False, indent=2))

if __name__ == '__main__':
    OUT.mkdir(parents=True, exist_ok=True)
    labels = write_pdf()
    write_svg()
    validate_and_preview(labels)
