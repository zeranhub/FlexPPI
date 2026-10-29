"""Revise the information hierarchy in 3.2 without changing frozen evidence."""
from pathlib import Path
from copy import deepcopy
import hashlib
import json
import zipfile
from lxml import etree as E

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT.parent / 'paper_v72_20261024_中文修订_结果五节重构.docx'
TARGET = ROOT.parent / 'paper_v73_20261024_中文修订_AA至Full分层贡献.docx'
NS = {'w': 'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
W = '{' + NS['w'] + '}'

def text(el):
    return ''.join(el.xpath('.//w:t/text()', namespaces=NS))

def main():
    assert hashlib.sha256(SOURCE.read_bytes()).hexdigest() == 'c0a380369098ad572381bda4f9219fddaa0c5304699c6039916ab00688018b94'
    assert not TARGET.exists()
    section = json.loads((ROOT / 'section32_v73.json').read_text(encoding='utf-8'))
    with zipfile.ZipFile(SOURCE) as package:
        entries = package.infolist()
        parts = {item.filename: package.read(item.filename) for item in entries}
    root = E.fromstring(parts['word/document.xml'])
    body = root.find('w:body', NS)
    original = list(body)
    start = next(i for i, el in enumerate(original) if text(el).startswith('3.2 '))
    end = next(i for i, el in enumerate(original) if text(el).startswith('3.3 '))
    heading = original[start]
    templates = original[start+1:end]
    assert len(templates) == 4 and all(el.tag == W+'p' for el in templates)
    assert not any(el.xpath('.//w:drawing|.//w:object', namespaces=NS) for el in templates)
    before = [E.tostring(el) for el in original[:start]]
    after = [E.tostring(el) for el in original[end:]]

    def paragraph(template, value):
        el = deepcopy(template)
        for attr in list(el.attrib):
            if attr.endswith('paraId') or attr.endswith('textId'):
                del el.attrib[attr]
        rpr = el.find('w:r/w:rPr', NS)
        for child in list(el):
            if child.tag != W+'pPr':
                el.remove(child)
        run = E.SubElement(el, W+'r')
        if rpr is not None:
            run.append(deepcopy(rpr))
        node = E.SubElement(run, W+'t')
        node.set('{http://www.w3.org/XML/1998/namespace}space', 'preserve')
        node.text = value
        return el

    replacement = [paragraph(heading, '3.2 ' + section['title'])]
    replacement += [paragraph(templates[0], p) for p in section['paragraphs']]
    for el in original[start:end]:
        body.remove(el)
    for offset, el in enumerate(replacement):
        body.insert(start+offset, el)
    assert [E.tostring(el) for el in list(body)[:start]] == before
    assert [E.tostring(el) for el in list(body)[start+len(replacement):]] == after
    new_xml = E.tostring(root, xml_declaration=True, encoding='UTF-8', standalone=True)
    with zipfile.ZipFile(TARGET, 'w', zipfile.ZIP_DEFLATED) as package:
        for entry in entries:
            package.writestr(entry, new_xml if entry.filename == 'word/document.xml' else parts[entry.filename])
    with zipfile.ZipFile(TARGET) as package:
        assert package.testzip() is None
        changed = [name for name, content in parts.items() if package.read(name) != content]
        assert changed == ['word/document.xml']
    report = {
        'source': str(SOURCE), 'target': str(TARGET),
        'target_sha256': hashlib.sha256(TARGET.read_bytes()).hexdigest(),
        'changed_package_parts': changed, 'only_section_3_2_changed': True,
        'all_other_body_elements_exact': True, 'all_media_tables_math_and_captions_exact': True,
        'title': section['title'], 'paragraph_count': len(section['paragraphs']),
        'characters': sum(map(len, section['paragraphs'])),
        'concept_checks': ['Static includes partner contact summaries', 'Intra uses learned intrachain geometry', 'Full jointly encodes interchain geometry', 'same visible chain range for single-chain input', 'case links to sections 3.3, 3.4 and 3.5'],
    }
    (ROOT/'v73_patch_audit.json').write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps(report, ensure_ascii=False))

if __name__ == '__main__':
    main()
