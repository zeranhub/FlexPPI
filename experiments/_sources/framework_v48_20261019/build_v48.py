"""Surgical Figure 1 replacement in the v45 manuscript.

Only three things change: the Figure 1 image part, the Figure 1 caption text
(native R_i / t_i equations preserved in place) and the image extent. Every
other package part is copied byte for byte.
"""
from pathlib import Path
from copy import deepcopy
import hashlib
import json
import re
import zipfile

from lxml import etree
from PIL import Image

ROOT = Path(__file__).resolve().parent
BASE = ROOT.parent.parent
QA = ROOT / "manuscript_QA"
QA.mkdir(exist_ok=True)

SOURCE = BASE / "paper_v45_20261018_中文修订_紧凑多层模型评价图_三线表.docx"
OUT = BASE / "paper_v48_20261019_中文修订_模型框架图重绘_三线表.docx"
FIG = ROOT / "figure1" / "fig1_framework.png"
CAPTION = ROOT / "figure1" / "caption_zh.md"
SOURCE_SHA = "0e069a66c99e2196d50452cf52080eb703fc45a832dfef4e17f1d325b0cdfb63"
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
assert sha(SOURCE) == SOURCE_SHA

NS = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main",
      "m": "http://schemas.openxmlformats.org/officeDocument/2006/math",
      "a": "http://schemas.openxmlformats.org/drawingml/2006/main",
      "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
      "wp": "http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing"}

with zipfile.ZipFile(SOURCE) as z:
    infos = z.infolist()
    parts = {p.filename: z.read(p.filename) for p in infos}

root = etree.fromstring(parts["word/document.xml"])
body = root.find("w:body", NS)
txt = lambda n: "".join(n.xpath(".//w:t/text()", namespaces=NS))

caps = [p for p in body if p.tag == "{" + NS["w"] + "}p" and txt(p).startswith("图1：")]
assert len(caps) == 1
cap = caps[0]
fig = cap.getprevious()
figidx, capidx = list(body).index(fig), list(body).index(cap)
before = txt(cap)
after = CAPTION.read_text(encoding="utf8").strip()
assert after.startswith("图1：")

maths = cap.xpath("./m:oMath", namespaces=NS)
assert [''.join(x.xpath('.//m:t/text()', namespaces=NS)) for x in maths] == ["Ri", "ti", "Ri", "ti"]
assert not cap.xpath(".//w:commentRangeStart|.//w:bookmarkStart", namespaces=NS)
matches = list(re.finditer("Rᵢ|tᵢ", after))
assert [m.group() for m in matches] == ["Rᵢ", "tᵢ", "Rᵢ", "tᵢ"], [m.group() for m in matches]

for node in list(cap):
    if node.tag != "{" + NS["w"] + "}pPr":
        cap.remove(node)
start = 0
for match, math in zip(matches, maths):
    run = etree.SubElement(cap, "{" + NS["w"] + "}r")
    etree.SubElement(run, "{" + NS["w"] + "}t").text = after[start:match.start()]
    cap.append(deepcopy(math))
    start = match.end()
run = etree.SubElement(cap, "{" + NS["w"] + "}r")
etree.SubElement(run, "{" + NS["w"] + "}t").text = after[start:]

blips = fig.xpath(".//a:blip", namespaces=NS)
assert len(blips) == 1
rid = blips[0].get("{" + NS["r"] + "}embed")
rels = etree.fromstring(parts["word/_rels/document.xml.rels"])
media = "word/" + {n.get("Id"): n.get("Target") for n in rels}[rid]
assert media == "word/media/image1.png", media
inline = blips[0].xpath("ancestor::wp:inline", namespaces=NS)[0]
width = round(170.0 / 25.4 * 914400)
with Image.open(FIG) as im:
    height = round(width * im.height / im.width)
for e in inline.xpath("./wp:extent|.//a:xfrm/a:ext", namespaces=NS):
    e.set("cx", str(width))
    e.set("cy", str(height))
for pr in inline.xpath("./wp:docPr", namespaces=NS):
    pr.set("name", "Figure 1 FlexPPI model framework")
    pr.set("descr", "FlexPPI framework: static backbone input, residue/pair/frame features, four shared geometric attention blocks, shared readout MLP for RMSF and Hchi mean and log-variance, and variance-aware MD supervision")


def prop(p, name, value):
    pr = p.find("w:pPr", NS)
    if pr is None:
        pr = etree.Element("{" + NS["w"] + "}pPr")
        p.insert(0, pr)
    e = pr.find("w:" + name, NS)
    if e is None:
        e = etree.SubElement(pr, "{" + NS["w"] + "}" + name)
    e.set("{" + NS["w"] + "}val", str(value))


prop(fig, "keepNext", 1)
prop(cap, "keepLines", 1)

alltext = "\n".join(txt(p) for p in body)
assert len(root.xpath(".//w:tbl", namespaces=NS)) == 10
assert len(root.xpath(".//m:oMath", namespaces=NS)) == 36
assert len(root.xpath(".//a:blip", namespaces=NS)) == 6
assert len(re.findall(r"^\[\d+\]", alltext, re.M)) == 28

updates = {media: FIG.read_bytes(),
           "word/document.xml": etree.tostring(root, xml_declaration=True,
                                               encoding="UTF-8", standalone=True)}
with zipfile.ZipFile(OUT, "w") as z:
    for info in infos:
        z.writestr(info, updates.get(info.filename, parts[info.filename]))
with zipfile.ZipFile(OUT) as z:
    assert z.testzip() is None

changed = sorted(updates)
untouched = [i.filename for i in infos if i.filename not in updates]
audit = {"status": "built_pending_render_QA", "source": str(SOURCE),
         "source_sha256": sha(SOURCE), "document": str(OUT), "document_sha256": sha(OUT),
         "figure1": str(FIG), "figure1_sha256": sha(FIG),
         "caption_before": before, "caption_after": after,
         "caption_body_index": capidx, "figure_body_index": figidx,
         "changed_zip_parts": changed, "untouched_parts": len(untouched),
         "width_EMU": width, "height_EMU": height, "width_mm": round(width / 36000, 2),
         "height_mm": round(height / 36000, 2),
         "native_tables": 10, "native_equations": 36, "references": 28,
         "main_figures": 5, "supplementary_figures": 1,
         "prose_changed": False, "new_MD": 0, "new_training": 0, "new_inference": 0}
(QA / "integration_audit.json").write_text(json.dumps(audit, ensure_ascii=False, indent=2),
                                           encoding="utf8")
print(json.dumps({k: audit[k] for k in ("document", "document_sha256", "width_mm",
                                        "height_mm", "changed_zip_parts")}, ensure_ascii=False))
