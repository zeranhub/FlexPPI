"""Replace only Figure 2 artwork after tightening the a/b title spacing."""
from pathlib import Path
import hashlib
import json
import zipfile
from lxml import etree as ET
from PIL import Image

ROOT = Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010')
OUT = ROOT / 'figure2_title_spacing_20261024'
SOURCE = ROOT / 'paper_v67_20261024_中文修订_图2组合子图与图例重排.docx'
TARGET = ROOT / 'paper_v68_20261024_中文修订_图2标题间距优化.docx'
NS = {'w': 'http://schemas.openxmlformats.org/wordprocessingml/2006/main',
      'wp': 'http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing',
      'a': 'http://schemas.openxmlformats.org/drawingml/2006/main'}


def main():
    assert not TARGET.exists(), 'Preserve existing versions.'
    picture = OUT / 'fig2_title_spacing.png'
    assert picture.is_file()
    with zipfile.ZipFile(SOURCE) as package:
        entries = package.infolist()
        parts = {item.filename: package.read(item.filename) for item in entries}
    root = ET.fromstring(parts['word/document.xml'])
    text_before = root.xpath('//w:t/text()', namespaces=NS)
    drawing = root.xpath('//w:p[.//wp:docPr[@id="2"]]', namespaces=NS)
    assert len(drawing) == 1
    extents = drawing[0].xpath('.//wp:extent | .//a:xfrm/a:ext', namespaces=NS)
    original_extents = [[int(el.get('cx')), int(el.get('cy'))] for el in extents]
    cx = original_extents[0][0]
    with Image.open(picture) as image:
        image_size = list(image.size)
        cy = round(cx * image.height / image.width)
    for element in extents:
        element.set('cx', str(cx))
        element.set('cy', str(cy))
    assert root.xpath('//w:t/text()', namespaces=NS) == text_before
    changed = {'word/media/image2.png': picture.read_bytes()}
    if any(pair != [cx, cy] for pair in original_extents):
        changed['word/document.xml'] = ET.tostring(root, xml_declaration=True,
                                                  encoding='UTF-8', standalone=True)
    with zipfile.ZipFile(TARGET, 'w', zipfile.ZIP_DEFLATED) as package:
        for item in entries:
            package.writestr(item, changed.get(item.filename, parts[item.filename]))
    with zipfile.ZipFile(TARGET) as package:
        assert package.testzip() is None
        actual = [name for name, value in parts.items() if package.read(name) != value]
        assert set(actual) == set(changed)
        final_root = ET.fromstring(package.read('word/document.xml'))
        assert final_root.xpath('//w:t/text()', namespaces=NS) == text_before
    report = {'source': str(SOURCE), 'target': str(TARGET),
              'source_sha256': hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
              'target_sha256': hashlib.sha256(TARGET.read_bytes()).hexdigest(),
              'changed_package_parts': actual, 'all_other_parts_byte_identical': True,
              'all_document_text_unchanged': True, 'figure_pixels': image_size,
              'old_extents_emu': original_extents, 'new_extent_emu': [cx, cy]}
    (OUT / 'docx_patch_audit.json').write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
