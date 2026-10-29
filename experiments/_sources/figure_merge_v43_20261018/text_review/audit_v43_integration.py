"""Read-only package and cross-reference audit for the two merged figures.

Usage: bundled_python audit_v43_integration.py --new NEW.docx
       --figure2 FULL_FIG2.png --figure5 FULL_FIG5.png --output REPORT.json
Only the report is written; all manuscript and image inputs are read-only.
"""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import re
from zipfile import ZipFile
from lxml import etree

BASE = Path(r"E:\Project\FlexPPI\paper\revision_v14_20261010")
SOURCE = BASE / "paper_v42_20261017_中文修订_蛋白背景与实例动机_三线表.docx"
NS = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main", "m": "http://schemas.openxmlformats.org/officeDocument/2006/math", "a": "http://schemas.openxmlformats.org/drawingml/2006/main", "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships"}


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def nh(node):
    return sha(etree.tostring(node, method="c14n"))


def txt(node):
    return "".join(node.xpath(".//w:t/text()", namespaces=NS))


def read(path):
    with ZipFile(path) as z:
        parts = {n: z.read(n) for n in z.namelist()}
    root = etree.fromstring(parts["word/document.xml"])
    body = root.find("w:body", NS)
    rels = etree.fromstring(parts["word/_rels/document.xml.rels"])
    rm = {x.get("Id"): x.get("Target") for x in rels}
    active = []
    for index, node in enumerate(body):
        for rid in node.xpath(".//a:blip/@r:embed", namespaces=NS):
            media = "word/" + rm[rid]
            active.append({"body_index": index, "rid": rid, "media": media, "sha256": sha(parts[media])})
    return {"path": str(path), "sha256": sha(path.read_bytes()), "parts": parts, "root": root, "body": body, "active": active}


def hashes(doc, query):
    return [nh(n) for n in doc["root"].xpath(query, namespaces=NS)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--new", required=True, type=Path)
    ap.add_argument("--figure2", required=True, type=Path)
    ap.add_argument("--figure5", required=True, type=Path)
    ap.add_argument("--output", required=True, type=Path)
    args = ap.parse_args()
    old, new = read(SOURCE), read(args.new)
    op, np = old["parts"], new["parts"]
    comments = [n for n in op if "comment" in n.lower()]
    anchor_tags = ["commentRangeStart", "commentRangeEnd", "commentReference", "bookmarkStart", "bookmarkEnd"]
    anchors = lambda d: {tag: [dict(n.attrib) for n in d["root"].xpath("//w:" + tag, namespaces=NS)] for tag in anchor_tags}
    refs = lambda d: [nh(n) for n in d["body"] if re.match(r"^\[\d+\]", txt(n))]
    expected_media = [
        sha(op["word/media/image1.png"]), sha(args.figure2.read_bytes()),
        sha(op["word/media/image5.png"]), sha(op["word/media/image41.png"]),
        sha(args.figure5.read_bytes()), sha(op["word/media/image7.png"]),
    ]
    original_removed_hashes = [sha(op["word/media/" + n]) for n in ["image2.png", "image3.png", "image21.png", "image40.png"]]
    actual_media = [x["sha256"] for x in new["active"]]
    preserved_media = ["word/media/" + n for n in ["image1.png", "image5.png", "image41.png", "image7.png"]]

    # All substantive original content outside Results remains protected, except
    # the method and supplement paragraphs that require figure-number changes.
    results_start = next(i for i, n in enumerate(old["body"]) if txt(n).startswith("3.1 "))
    discussion_start = next(i for i, n in enumerate(old["body"]) if txt(n).startswith("4 "))
    protected = []
    for index, n in enumerate(old["body"]):
        if results_start <= index < discussion_start:
            continue
        t = txt(n)
        if t.startswith(("侧链实例使用泛素", "正文图4采用")):
            continue
        if not t and not n.xpath(".//w:drawing | .//m:oMath | .//w:commentReference | .//w:bookmarkStart", namespaces=NS):
            continue
        protected.append({"source_body_index": index, "text": t, "hash": nh(n)})
    current = Counter(nh(n) for n in new["body"])
    missing_protected = []
    for item in protected:
        if current[item["hash"]] == 0:
            missing_protected.append(item)
        else:
            current[item["hash"]] -= 1

    texts = [txt(n) for n in new["body"]]
    captions = [t for t in texts if re.match(r"^图[1-9]：", t)]
    expected_headings = ["3.1 两类动态预测与完整结构环境的贡献", "3.2 同类型残基的局部动态差异", "3.3 KIX结构域中主链活动区域的定位", "3.4 KRAS G12D单点抗原变化的主链动态响应", "3.5 联合读数的构象对应与群体验证", "3.5.1 侧链预测对低自由能取向范围与角度占据的验证", "3.5.2 泛素同类型侧链多样性的预测与轨迹对照"]
    observed_headings = [t for t in texts if re.match(r"^3\.[1-9](?:\.[1-9])?\s", t)]
    old_figure_hits = [{"index": i, "text": t} for i, t in enumerate(texts) if re.search(r"图\s*[67](?!\d)", t)]
    old_section_hits = [{"index": i, "text": t} for i, t in enumerate(texts) if re.search(r"^3\.6(?:\.|\s)", t)]
    fig2caption = next((t for t in texts if t.startswith("图2：")), "")
    fig5caption = next((t for t in texts if t.startswith("图5：")), "")
    changed_parts = sorted(n for n in op.keys() & np.keys() if op[n] != np[n])
    unexpected_parts = [n for n in changed_parts if n not in {"word/document.xml", "word/_rels/document.xml.rels", "[Content_Types].xml"} and not n.startswith("word/media/")]
    checks = {
        "source_v42_SHA_matches": old["sha256"] == "0d56dd11e3abab0b7e12bbd75ef33c352d949bc305212c01f9eb5509541df887",
        "10_native_tables_identical": len(hashes(new, "//w:tbl")) == 10 and hashes(old, "//w:tbl") == hashes(new, "//w:tbl"),
        "36_native_OMML_equations_identical": len(hashes(new, "//m:oMath")) == 36 and hashes(old, "//m:oMath") == hashes(new, "//m:oMath"),
        "28_references_identical": len(refs(new)) == 28 and refs(old) == refs(new),
        "all_original_comment_parts_identical": all(n in np and op[n] == np[n] for n in comments),
        "all_comment_and_bookmark_anchors_identical": anchors(old) == anchors(new),
        "four_retained_media_files_byte_identical": all(n in np and op[n] == np[n] for n in preserved_media),
        "exactly_six_active_images_in_correct_order": actual_media == expected_media,
        "old_2_3_4_7_images_not_active": all(h not in actual_media for h in original_removed_hashes),
        "outside_Results_protected_content_identical": not missing_protected,
        "five_main_captions_numbered_in_order": [re.match(r"^图([1-9])", t).group(1) for t in captions] == ["1", "2", "3", "4", "5"],
        "new_results_heading_sequence_correct": observed_headings == expected_headings,
        "no_legacy_Fig6_or_Fig7_text_reference": not old_figure_hits,
        "no_legacy_section36_heading": not old_section_hits,
        "Fig2_no_removed_ranking_distribution_caption": "逐记录Spearman相关分布" not in fig2caption,
        "Fig5_caption_has_83_records_and33142_positions": "83" in fig5caption and "33,142" in fig5caption,
        "only_expected_package_part_types_changed": not unexpected_parts,
    }
    report = {
        "status": "pass" if all(checks.values()) else "fail", "document_sha256": new["sha256"],
        "document": str(args.new), "source_document": str(SOURCE), "checks": checks,
        "active_images": new["active"], "expected_merged_image_paths": [str(args.figure2), str(args.figure5)],
        "changed_parts": changed_parts, "unexpected_changed_parts": unexpected_parts,
        "added_parts": sorted(np.keys() - op.keys()), "removed_parts": sorted(op.keys() - np.keys()),
        "missing_protected_content": missing_protected,
        "legacy_figure_reference_hits": old_figure_hits, "legacy_section_heading_hits": old_section_hits,
        "observed_results_headings": observed_headings,
        "manual_scientific_review_required": ["新图5b为配对相关差，不是条件独立贡献", "新图5e不是逐χ角预测", "新图5f及角宽来自MD，不是模型自由能预测", "四类召回完整保留于主文一句和S9", "83/33142与3063/824988分母分清"]
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({"status": report["status"], "document_sha256": new["sha256"], "checks": checks}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
