"""Replace figures 2–5 and their colour/scale wording in a copy of v58.

Uses direct OOXML edits to preserve all other text, equations, tables and styles.
Run with the bundled workspace Python. No scientific data are recalculated.
"""
from copy import deepcopy
from pathlib import Path
import argparse
import hashlib
import json
from zipfile import ZipFile
from lxml import etree as ET
from PIL import Image

ROOT = Path(__file__).resolve().parent
SOURCE = Path("E:/Project/FlexPPI/paper/revision_v14_20261010/paper_v58_20261022_中文修订_图2至5配色统一_三线表.docx")
NS = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main",
      "wp": "http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing",
      "a": "http://schemas.openxmlformats.org/drawingml/2006/main",
      "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships"}
MATH_NS = "http://schemas.openxmlformats.org/officeDocument/2006/math"

def sha(data):
    return hashlib.sha256(data).hexdigest()

def caption(number, original):
    if number == "2":
        return original + "蓝色热图使用0—1统一色标，深色为高AP。"
    if number == "3":
        return original.replace("四种组成下完整KIX的残基排序", "四种组成下完整KIX的RMSF百分位秩（0—1，深蓝为高秩）")
    if number == "4":
        return original.replace("两项分别标示单位和色标", "两项分别标示单位，色标均为−0.45至+0.45并以零为中心").replace("蓝、红表示降低、升高", "蓝、暖赭表示降低、升高")
    if number == "5":
        return original.replace("Q2/41为紫/橙色", "Q2为青绿方形、Q41为暖赭菱形")
    raise ValueError(number)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT.parent / "paper_v59_20261023_中文修订_图2至5整体重绘.docx")
    args = parser.parse_args()
    manifest = json.loads((ROOT / "manifest.json").read_text(encoding="utf-8"))
    with ZipFile(SOURCE) as source:
        original = {name: source.read(name) for name in source.namelist()}
        document = ET.fromstring(original["word/document.xml"])
        relationships = ET.fromstring(original["word/_rels/document.xml.rels"])
        targets = {r.attrib["Id"]: "word/" + r.attrib["Target"] for r in relationships}
        changed = {}
        audit = {"source_document": str(SOURCE), "source_sha256": sha(SOURCE.read_bytes()), "figures": {}}
        for figure in manifest["figures"]:
            number = figure["figure"]
            if number not in ["2", "3", "4", "5"]:
                continue
            image = ROOT / figure["folder"] / (figure["stem"] + ".png")
            media = figure["word_media"]
            changed[media] = image.read_bytes()
            updated = 0
            for paragraph in document.findall(".//w:p", NS):
                text = "".join(paragraph.xpath(".//w:t/text()", namespaces=NS))
                if not text.startswith("图" + number + "："):
                    continue
                assert not paragraph.findall(".//{%s}oMath" % MATH_NS), "caption contains native math"
                replacement = caption(number, text)
                first_run = paragraph.find("w:r", NS)
                props = first_run.find("w:rPr", NS) if first_run is not None else None
                run = ET.Element("{%s}r" % NS["w"])
                if props is not None:
                    run.append(deepcopy(props))
                node = ET.SubElement(run, "{%s}t" % NS["w"])
                node.text = replacement
                for child in list(paragraph):
                    if child.tag != "{%s}pPr" % NS["w"]:
                        paragraph.remove(child)
                paragraph.append(run)
                updated += 1
            assert updated == 1, (number, updated)
            # Preserve the original figure width and update height only if the
            # deliberate drawing canvas changed; never stretch a new image.
            with Image.open(image) as im:
                ratio = im.height / im.width
            shapes = 0
            for drawing in document.findall(".//w:drawing", NS):
                blip = drawing.find(".//a:blip", NS)
                if blip is None or targets[blip.attrib["{%s}embed" % NS["r"]]] != media:
                    continue
                extent = drawing.find(".//wp:extent", NS)
                width = int(extent.attrib["cx"])
                height = round(width * ratio)
                extent.attrib["cy"] = str(height)
                for shape in drawing.findall(".//a:xfrm/a:ext", NS):
                    shape.attrib["cy"] = str(height)
                shapes += 1
            assert shapes == 1
            audit["figures"][number] = {"media": media, "sha256": sha(changed[media]),
                "source": str(image), "caption_updated": True, "aspect_ratio_preserved": True}
        changed["word/document.xml"] = ET.tostring(document, xml_declaration=True, encoding="UTF-8", standalone=True)
        args.output.parent.mkdir(exist_ok=True, parents=True)
        with ZipFile(args.output, "w") as destination:
            for info in source.infolist():
                destination.writestr(info, changed.get(info.filename, original[info.filename]))
    with ZipFile(args.output) as result:
        differences = [name for name, data in original.items() if result.read(name) != data]
        assert set(differences) == set(changed), differences
        for name, data in original.items():
            if name not in changed:
                assert result.read(name) == data
    audit.update({"output_document": str(args.output.resolve()), "output_sha256": sha(args.output.read_bytes()),
                  "modified_zip_members": differences, "all_other_members_byte_identical": True})
    (ROOT / "document_QA" / "document_edit_audit.json").write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8")
    print(args.output.resolve())

if __name__ == "__main__":
    main()
