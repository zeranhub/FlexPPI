"""Reproducible targeted v79 -> v80 DOCX patch using the bundled Python runtime.

Only text selected in the reviewed JSON plan, Figures 3/5, appended references,
and supplementary CaM methods/Figure S2 are changed. The source remains intact.
OOXML targeting preserves all other parts, equations, tables, styles and section
3.4 byte-for-byte as XML elements. The root agent owns the exactly-once artifact
operation marker and the final bundled render / full-page visual QA.
"""
from __future__ import annotations
import argparse
import copy
import hashlib
import io
import json
from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED

from lxml import etree
from PIL import Image

HERE = Path(__file__).resolve().parent
NS = {
    'w': 'http://schemas.openxmlformats.org/wordprocessingml/2006/main',
    'r': 'http://schemas.openxmlformats.org/officeDocument/2006/relationships',
    'a': 'http://schemas.openxmlformats.org/drawingml/2006/main',
    'wp': 'http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing',
    'pic': 'http://schemas.openxmlformats.org/drawingml/2006/picture',
    'm': 'http://schemas.openxmlformats.org/officeDocument/2006/math',
    'pr': 'http://schemas.openxmlformats.org/package/2006/relationships',
    'ct': 'http://schemas.openxmlformats.org/package/2006/content-types',
}

def q(prefix, local):
    return '{%s}%s' % (NS[prefix], local)

def sha(data):
    return hashlib.sha256(data).hexdigest()

def xml_signature(element):
    return sha(etree.tostring(element, method='c14n'))

def ptext(paragraph):
    return ''.join(paragraph.xpath('.//w:t/text()', namespaces=NS))

def replace_text(paragraph, text):
    """Keep paragraph format and copy its first nonempty run's character style."""
    assert not paragraph.xpath('.//m:oMath|.//m:oMathPara|.//w:drawing', namespaces=NS), 'Unsafe rich paragraph replacement'
    rpr = paragraph.find('.//' + q('w', 'rPr'))
    rpr = copy.deepcopy(rpr) if rpr is not None else None
    for child in list(paragraph):
        if child.tag != q('w', 'pPr'):
            paragraph.remove(child)
    run = etree.SubElement(paragraph, q('w', 'r'))
    if rpr is not None:
        run.append(rpr)
    node = etree.SubElement(run, q('w', 't'))
    node.text = text
    if text[:1].isspace() or text[-1:].isspace():
        node.set('{http://www.w3.org/XML/1998/namespace}space', 'preserve')

def new_text_paragraph(template, text):
    paragraph = copy.deepcopy(template)
    replace_text(paragraph, text)
    return paragraph

def set_image_geometry(paragraph, path, alt_text, width_mm=None):
    with Image.open(path) as image:
        px_w, px_h = image.size
    extents = paragraph.xpath('.//wp:extent', namespaces=NS)
    assert len(extents) == 1
    cx = round(width_mm * 36000) if width_mm else int(extents[0].get('cx'))
    cy = round(cx * px_h / px_w)
    for node in paragraph.xpath('.//wp:extent|.//a:xfrm/a:ext', namespaces=NS):
        node.set('cx', str(cx)); node.set('cy', str(cy))
    for node in paragraph.xpath('.//wp:effectExtent', namespaces=NS):
        for field in ('l', 't', 'r', 'b'):
            node.set(field, '0')
    for node in paragraph.xpath('.//wp:docPr', namespaces=NS):
        node.set('descr', alt_text)
    return {'path': str(path), 'pixels': [px_w, px_h], 'display_emu': [cx, cy], 'sha256': sha(path.read_bytes())}


def paragraph_setting(paragraph, name, **attributes):
    ppr = paragraph.find(q('w', 'pPr'))
    if ppr is None:
        ppr = etree.Element(q('w', 'pPr')); paragraph.insert(0, ppr)
    node = ppr.find(q('w', name))
    if node is None:
        node = etree.SubElement(ppr, q('w', name))
    for key, value in attributes.items():
        node.set(q('w', key), str(value))


def keep_figure_with_caption(image_paragraph, caption):
    # Fit the integrated figure and its complete caption on the same A4 page.
    paragraph_setting(image_paragraph, 'keepNext', val=1)
    paragraph_setting(image_paragraph, 'keepLines', val=1)
    paragraph_setting(image_paragraph, 'spacing', before=0, after=0, line=240, lineRule='auto')
    paragraph_setting(caption, 'keepNext', val=0)
    paragraph_setting(caption, 'keepLines', val=1)
    paragraph_setting(caption, 'spacing', before=40, after=80, line=224, lineRule='exact')
    paragraph_setting(caption, 'snapToGrid', val=0)
    for run in caption.findall(q('w', 'r')):
        rpr = run.find(q('w', 'rPr'))
        if rpr is None:
            rpr = etree.Element(q('w', 'rPr')); run.insert(0, rpr)
        for name in ('sz', 'szCs'):
            node = rpr.find(q('w', name))
            if node is None: node = etree.SubElement(rpr, q('w', name))
            node.set(q('w', 'val'), '18')

def read_caption(path):
    text = Path(path).read_text(encoding='utf-8').strip()
    assert text and not text.startswith('```')
    return text.replace('\n\n', ' ').replace('\n', ' ')

def source_audit(plan):
    source = Path(plan['source_docx'])
    with ZipFile(source) as package:
        root = etree.fromstring(package.read('word/document.xml'))
        paragraphs = root.xpath('//w:body/w:p', namespaces=NS)
        original_inventory = json.loads((HERE / 'v79_paragraph_inventory.json').read_text(encoding='utf-8'))
        assert len(paragraphs) == len(original_inventory) == 176
        for index, item in enumerate(original_inventory):
            assert item['text'] == ptext(paragraphs[index]), ('Source paragraph mismatch', index)
        checked = {
            'source': str(source), 'source_sha256': sha(source.read_bytes()),
            'paragraphs': len(paragraphs), 'tables': len(root.xpath('//w:tbl', namespaces=NS)),
            'equations': len(root.xpath('//m:oMath', namespaces=NS)),
            'original_references_preserved': list(range(1, 29)),
            'section_34_indices_unchanged': list(range(76, 83)),
            'replaced_indices': sorted(map(int, plan['replacements'])),
            'docx_created': False,
        }
    (HERE / 'v83_source_audit.json').write_text(json.dumps(checked, ensure_ascii=False, indent=2), encoding='utf-8')
    return checked

def build(args):
    plan_path = Path(args.plan)
    plan = json.loads(plan_path.read_text(encoding='utf-8'))
    audit = source_audit(plan)
    if args.audit_only:
        print(json.dumps(audit, ensure_ascii=False, indent=2))
        return
    assert args.operation_mark_completed, 'Root must confirm the exactly-once artifact edit marker before authoring.'
    for name in ('fig3', 'fig5', 'fig_s2', 'fig_s3', 'caption3', 'caption5', 'caption_s3'):
        value = getattr(args, name)
        assert value and Path(value).is_file(), f'Missing finalized asset: {name}'
    source = Path(plan['source_docx'])
    output = Path(args.output or plan['output_docx'])
    assert output.resolve() != source.resolve()
    with ZipFile(source) as package:
        original_parts = {item.filename: package.read(item.filename) for item in package.infolist()}
        document = etree.fromstring(original_parts['word/document.xml'])
        rels = etree.fromstring(original_parts['word/_rels/document.xml.rels'])
        paragraphs = document.xpath('//w:body/w:p', namespaces=NS)
        body = document.find(q('w', 'body'))
        protected_34 = [xml_signature(p) for p in paragraphs[76:83]]
        protected_tables = [xml_signature(p) for p in document.xpath('//w:tbl', namespaces=NS)]
        protected_equations = [xml_signature(p) for p in document.xpath('//m:oMath', namespaces=NS)]
        protected_references = [xml_signature(p) for p in paragraphs[100:128]]
        changed = []
        for key, text in plan['replacements'].items():
            index = int(key)
            changed.append({'source_index': index, 'before': ptext(paragraphs[index]), 'after': text})
            replace_text(paragraphs[index], text)
        replace_text(paragraphs[72], read_caption(args.caption3))
        replace_text(paragraphs[89], read_caption(args.caption5))
        modified_parts = {}
        figure_report = {}
        for name, index, path, descr in (
            ('Figure 3', 71, Path(args.fig3), 'CaM matched rigid-reference orientation with right-hand local views and absolute two-state profiles'),
            ('Figure 5', 88, Path(args.fig5), 'Joint readout population validation and ubiquitin held-out prediction with angular occupancy'),
        ):
            paragraph = paragraphs[index]
            embed = paragraph.xpath('.//a:blip/@r:embed', namespaces=NS)
            assert len(embed) == 1
            relationship = rels.xpath('./pr:Relationship[@Id=$rid]', namespaces=NS, rid=embed[0])
            assert len(relationship) == 1
            target = 'word/' + relationship[0].get('Target')
            assert target.endswith('.png')
            modified_parts[target] = path.read_bytes()
            figure_report[name] = set_image_geometry(paragraph, path, descr, width_mm=164)
            keep_figure_with_caption(paragraph, paragraphs[index + 1])
        # Keep the full Figure 3 and caption on their own page, then resume the case narrative.
        for case_index in (73, 74):
            paragraph_setting(paragraphs[case_index], 'keepLines', val=1)
        for text in plan.get('insert_before_158', []):
            paragraphs[158].addprevious(new_text_paragraph(paragraphs[164], text))
        references = plan.get('references_to_append', [])
        assert len(references) == 2 and references[0].startswith('[29]') and references[1].startswith('[30]')
        for text in references:
            paragraphs[128].addprevious(new_text_paragraph(paragraphs[127], text))
        s2_path = Path(args.fig_s2)
        s2_paragraph = copy.deepcopy(paragraphs[154])
        s2_embed = s2_paragraph.xpath('.//a:blip', namespaces=NS)
        assert len(s2_embed) == 1
        s2_embed[0].set(q('r', 'embed'), 'rIdSupplementV80Ub')
        max_id = max(int(x) for x in document.xpath('//wp:docPr/@id', namespaces=NS))
        for node in s2_paragraph.xpath('.//wp:docPr|.//pic:cNvPr', namespaces=NS):
            node.set('id', str(max_id + 1)); node.set('name', 'Supplementary Figure S2 ubiquitin')
        s2_rel = etree.SubElement(rels, q('pr', 'Relationship'))
        s2_rel.set('Id', 'rIdSupplementV80Ub')
        s2_rel.set('Type', 'http://schemas.openxmlformats.org/officeDocument/2006/relationships/image')
        s2_rel.set('Target', 'media/figureS2_ubiquitin_v80.png')
        modified_parts['word/media/figureS2_ubiquitin_v80.png'] = s2_path.read_bytes()
        figure_report['Figure S2'] = set_image_geometry(s2_paragraph, s2_path, 'Supplementary ubiquitin Q41 Q62 MD representatives, RMSF and 27-state occupancy')
        s2_caption_text = read_caption(args.caption_s2) if args.caption_s2 else plan['supplement_s2_caption']
        s2_caption = new_text_paragraph(paragraphs[155], s2_caption_text)
        paragraphs[167].addprevious(s2_paragraph)
        paragraphs[167].addprevious(s2_caption)
        # Preserve the reference comparisons in a separate supplemental figure.
        s3_path = Path(args.fig_s3)
        s3_paragraph = copy.deepcopy(s2_paragraph)
        s3_paragraph.xpath('.//a:blip', namespaces=NS)[0].set(q('r', 'embed'), 'rIdSupplementV83CaM')
        for node in s3_paragraph.xpath('.//wp:docPr|.//pic:cNvPr', namespaces=NS):
            node.set('id', str(max_id + 2)); node.set('name', 'Supplementary Figure S3 CaM reference validation')
        s3_rel = etree.SubElement(rels, q('pr', 'Relationship'))
        s3_rel.set('Id', 'rIdSupplementV83CaM')
        s3_rel.set('Type', 'http://schemas.openxmlformats.org/officeDocument/2006/relationships/image')
        s3_rel.set('Target', 'media/figureS3_CaM_v83.png')
        modified_parts['word/media/figureS3_CaM_v83.png'] = s3_path.read_bytes()
        figure_report['Figure S3'] = set_image_geometry(s3_paragraph, s3_path, 'Supplementary CaM MD references and frozen-model prediction validation', width_mm=164)
        s3_caption = new_text_paragraph(paragraphs[155], read_caption(args.caption_s3))
        keep_figure_with_caption(s3_paragraph, s3_caption)
        paragraphs[167].addprevious(s3_paragraph)
        paragraphs[167].addprevious(s3_caption)
        assert protected_34 == [xml_signature(p) for p in paragraphs[76:83]]
        assert protected_tables == [xml_signature(p) for p in document.xpath('//w:tbl', namespaces=NS)]
        assert protected_equations == [xml_signature(p) for p in document.xpath('//m:oMath', namespaces=NS)]
        assert protected_references == [xml_signature(p) for p in paragraphs[100:128]]
        # Avoid a three-line overflow page at the end of the discussion.
        for index in (84, 85, 86, 87, 91, 92, 93, 95, 97, 98):
            paragraph_setting(paragraphs[index], 'spacing', after=80, line=300, lineRule='auto')
        for index in (94, 96):
            paragraph_setting(paragraphs[index], 'spacing', before=160, after=80)
        modified_parts['word/document.xml'] = etree.tostring(document, xml_declaration=True, encoding='UTF-8', standalone=True)
        modified_parts['word/_rels/document.xml.rels'] = etree.tostring(rels, xml_declaration=True, encoding='UTF-8', standalone=True)
        output.parent.mkdir(parents=True, exist_ok=True)
        with ZipFile(output, 'w', ZIP_DEFLATED) as package_out:
            for name, data in original_parts.items():
                package_out.writestr(name, modified_parts.pop(name, data))
            for name, data in modified_parts.items():
                package_out.writestr(name, data)
    with ZipFile(output) as checked:
        assert checked.testzip() is None
        out_document = etree.fromstring(checked.read('word/document.xml'))
        all_text = '\n'.join(out_document.xpath('//w:body/w:p//w:t/text()', namespaces=NS))
        main_text = all_text.split('补充材料')[0] if False else '\n'.join(ptext(p) for p in out_document.xpath('//w:body/w:p', namespaces=NS)[:130])
        paragraphs_out = out_document.xpath('//w:body/w:p', namespaces=NS)
        main_end = next(i for i,p in enumerate(paragraphs_out) if ptext(p) == '参考文献')
        main_text = '\n'.join(ptext(p) for p in paragraphs_out[:main_end])
        assert 'KIX' not in main_text and 'Q62' not in main_text.replace('Q41/Q62的MD代表帧', '')
        assert '图3e、g' not in main_text and '图3f、h' not in main_text
        assert '图5d' in main_text and '图S2' in main_text
        residual = [{'index':i,'text':ptext(p)} for i,p in enumerate(paragraphs_out) if 'KIX' in ptext(p)]
        final_audit = {
            **audit, 'docx_created': True, 'output': str(output),
            'output_sha256': sha(output.read_bytes()), 'plan_sha256': sha(plan_path.read_bytes()),
            'figures': figure_report, 'changed_paragraphs': changed,
            'section_34_preserved': True, 'equations_preserved': True,
            'tables_preserved': True, 'references_1_to_28_preserved': True,
            'KIX_remaining_only_reference_and_supplement': residual,
            'render_visual_QA_pending': True,
        }
    (HERE / 'v83_build_audit.json').write_text(json.dumps(final_audit, ensure_ascii=False, indent=2), encoding='utf-8')
    (HERE / 'v83_final_text.txt').write_text(all_text, encoding='utf-8')
    print(json.dumps({'output': str(output), 'sha256': final_audit['output_sha256'], 'render_visual_QA_pending':True}, ensure_ascii=False, indent=2))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--plan', default=str(HERE / 'v83_text_plan.json'))
    parser.add_argument('--audit-only', action='store_true')
    parser.add_argument('--operation-mark-completed', action='store_true')
    parser.add_argument('--fig3'); parser.add_argument('--fig5'); parser.add_argument('--fig-s2')
    parser.add_argument('--caption3'); parser.add_argument('--caption5'); parser.add_argument('--caption-s2')
    parser.add_argument('--fig-s3'); parser.add_argument('--caption-s3')
    parser.add_argument('--output')
    build(parser.parse_args())
