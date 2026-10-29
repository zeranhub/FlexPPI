"""Read-only OOXML comparison for the SB27 to KRAS manuscript replacement.

Writes only audit JSON files. Never modifies either DOCX.
"""
import argparse
from difflib import SequenceMatcher
import hashlib
import json
import re
from pathlib import Path
from zipfile import ZipFile

from lxml import etree

NS = {
    "w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main",
    "m": "http://schemas.openxmlformats.org/officeDocument/2006/math",
    "a": "http://schemas.openxmlformats.org/drawingml/2006/main",
    "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
}
W = "{" + NS["w"] + "}"
OLD_CASE = re.compile(r"SB27|SB47|(?<![A-Za-z0-9])Q7(?![A-Za-z0-9])|(?<![A-Za-z0-9])R7(?![A-Za-z0-9])|2AK4|4JRY|LPEPLPQGQLTAY|447185", re.I)


def sha(data):
    return hashlib.sha256(data).hexdigest()


def nodehash(node):
    return sha(etree.tostring(node, method="c14n"))


def text(node):
    return "".join(node.xpath(".//w:t/text()", namespaces=NS))


def read_docx(path):
    with ZipFile(path) as z:
        parts = {n: z.read(n) for n in z.namelist()}
    root = etree.fromstring(parts["word/document.xml"])
    body = root.find("w:body", NS)
    relationships = etree.fromstring(parts["word/_rels/document.xml.rels"])
    relmap = {r.get("Id"): r.get("Target") for r in relationships}
    visible_images = []
    for pindex, paragraph in enumerate(root.xpath("//w:p", namespaces=NS)):
        for rid in paragraph.xpath(".//a:blip/@r:embed", namespaces=NS):
            target = "word/" + relmap[rid]
            visible_images.append({"p_index": pindex, "rid": rid, "part": target, "sha256": sha(parts[target])})
    all_xml_hits = []
    for name, data in parts.items():
        if not name.endswith(".xml"):
            continue
        try:
            r = etree.fromstring(data)
        except etree.XMLSyntaxError:
            continue
        for p in r.xpath("//w:p", namespaces=NS):
            t = text(p)
            if OLD_CASE.search(t):
                all_xml_hits.append({"part": name, "text": t})
        # Catch diagram, drawing metadata, custom properties and other string attributes.
        for node in r.iter():
            for key, value in node.attrib.items():
                if OLD_CASE.search(value):
                    all_xml_hits.append({"part": name, "attribute": key, "value": value})
    anchors = {}
    for tag in ["commentRangeStart", "commentRangeEnd", "commentReference", "bookmarkStart", "bookmarkEnd"]:
        anchors[tag] = [dict(n.attrib) for n in root.xpath("//w:" + tag, namespaces=NS)]
    table_hashes = [nodehash(n) for n in root.xpath("//w:tbl", namespaces=NS)]
    math_hashes = [nodehash(n) for n in root.xpath("//m:oMath", namespaces=NS)]
    return {
        "path": str(path), "sha256": sha(Path(path).read_bytes()), "parts": parts,
        "root": root, "body": body, "anchors": anchors, "table_hashes": table_hashes,
        "math_hashes": math_hashes, "visible_images": visible_images,
        "old_case_hits": all_xml_hits,
    }


def public_manifest(doc):
    return {
        "path": doc["path"], "sha256": doc["sha256"],
        "parts_sha256": {n: sha(b) for n, b in doc["parts"].items()},
        "body_element_count": len(doc["body"]),
        "table_count": len(doc["table_hashes"]), "table_hashes": doc["table_hashes"],
        "OMML_count": len(doc["math_hashes"]), "OMML_hashes": doc["math_hashes"],
        "visible_images": doc["visible_images"], "anchors": doc["anchors"],
        "old_case_hits": doc["old_case_hits"],
    }


def compare(old, new):
    op, np = old["parts"], new["parts"]
    unchanged = [n for n in op.keys() & np.keys() if op[n] == np[n]]
    changed = [n for n in op.keys() & np.keys() if op[n] != np[n]]
    removed = sorted(op.keys() - np.keys())
    added = sorted(np.keys() - op.keys())
    same_body_length = len(old["body"]) == len(new["body"])
    body_changes = []
    allowed_old_indices = {2, 9, 49, 51, 52, 94, 95, 96, 97, 98, 99, 100, 116, 118, 150}
    old_hashes = [nodehash(n) for n in old["body"]]
    new_hashes = [nodehash(n) for n in new["body"]]
    change_blocks = []
    for operation, a0, a1, b0, b1 in SequenceMatcher(a=old_hashes, b=new_hashes, autojunk=False).get_opcodes():
        if operation == "equal":
            continue
        change_blocks.append({
            "operation": operation, "old_indices": list(range(a0, a1)),
            "new_indices": list(range(b0, b1)),
            "old_texts": [text(n) for n in old["body"][a0:a1]],
            "new_texts": [text(n) for n in new["body"][b0:b1]],
        })
    body_changes_confined = all(
        block["old_indices"] and set(block["old_indices"]) <= allowed_old_indices
        for block in change_blocks
    )
    if same_body_length:
        for index, (a, b) in enumerate(zip(old["body"], new["body"])):
            if nodehash(a) != nodehash(b):
                body_changes.append({"index": index, "old": text(a), "new": text(b), "old_tag": etree.QName(a).localname, "new_tag": etree.QName(b).localname})
    else:
        body_changes.append({"error": "Body element counts differ; inspect text-aligned diff separately."})
    # These are the seven figures outside Figure 6 actually referenced by v40.
    retained_images = ["word/media/image1.png", "word/media/image2.png", "word/media/image3.png", "word/media/image21.png", "word/media/image5.png", "word/media/image40.png", "word/media/image7.png"]
    image_checks = {n: n in np and op[n] == np[n] for n in retained_images}
    comments = [n for n in op if "comment" in n.lower()]
    checks = {
        "all_36_native_equations_identical": old["math_hashes"] == new["math_hashes"] and len(new["math_hashes"]) == 36,
        "all_10_native_tables_identical": old["table_hashes"] == new["table_hashes"] and len(new["table_hashes"]) == 10,
        "all_original_comment_parts_identical": all(n in np and op[n] == np[n] for n in comments),
        "all_comment_and_bookmark_anchors_identical": old["anchors"] == new["anchors"],
        "all_7_retained_figure_files_identical": all(image_checks.values()),
        "8_visible_figures": len(new["visible_images"]) == 8,
        "no_old_case_text_or_attributes": not new["old_case_hits"],
        "body_changes_confined_to_authorized_case_locations": body_changes_confined,
    }
    expected_parts = {"word/document.xml", "word/media/image41.png"}
    unexpected_parts = sorted(set(changed) - expected_parts)
    # Core-properties edits are informational metadata; record rather than fail.
    substantive_unexpected = [n for n in unexpected_parts if n not in {"docProps/core.xml", "docProps/app.xml"}]
    checks["only_expected_substantive_package_parts_changed"] = not substantive_unexpected and not removed and not added
    return {
        "old": old["path"], "new": new["path"], "new_sha256": new["sha256"],
        "checks": checks, "all_structural_checks_pass": all(checks.values()),
        "changed_parts": sorted(changed), "added_parts": added, "removed_parts": removed,
        "unexpected_changed_parts": unexpected_parts,
        "retained_image_checks": image_checks, "old_case_hits": new["old_case_hits"],
        "body_changes": body_changes, "new_visible_images": new["visible_images"],
        "same_body_element_count": same_body_length, "body_change_blocks": change_blocks,
        "comment_parts": comments,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("old", type=Path)
    parser.add_argument("new", type=Path, nargs="?")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    old = read_docx(args.old)
    result = compare(old, read_docx(args.new)) if args.new else public_manifest(old)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({"report": str(args.output), "checks": result.get("checks"), "old_case_hit_count": len(result["old_case_hits"]), "changed_parts": result.get("changed_parts")}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
