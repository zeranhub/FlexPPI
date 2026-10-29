"""Read-only OOXML audit for five protein-background additions in v42."""
import hashlib
import json
import re
from pathlib import Path
from zipfile import ZipFile

from lxml import etree

BASE = Path(r"E:\Project\FlexPPI\paper\revision_v14_20261010")
OLD = BASE / "paper_v41_20261017_中文修订_KRAS新抗原动态响应_三线表.docx"
NEW = BASE / "paper_v42_20261017_中文修订_蛋白背景与实例动机_三线表.docx"
OUT = BASE / "experiment_completion/protein_background_v42_20261017/manuscript_QA"
NS = {
    "w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main",
    "m": "http://schemas.openxmlformats.org/officeDocument/2006/math",
    "a": "http://schemas.openxmlformats.org/drawingml/2006/main",
}


def sha(value):
    return hashlib.sha256(value).hexdigest()


def nodehash(node):
    return sha(etree.tostring(node, method="c14n"))


def text(node):
    return "".join(node.xpath(".//w:t/text()", namespaces=NS))


def read(path):
    with ZipFile(path) as z:
        parts = {n: z.read(n) for n in z.namelist()}
    return parts, etree.fromstring(parts["word/document.xml"])


def main():
    old_parts, old = read(OLD)
    new_parts, new = read(NEW)
    old_body = old.find("w:body", NS)
    new_body = new.find("w:body", NS)
    old_hashes = [nodehash(n) for n in old_body]
    new_hashes = [nodehash(n) for n in new_body]
    changes = []
    for i, (a, b) in enumerate(zip(old_body, new_body)):
        if nodehash(a) != nodehash(b):
            changes.append({"OOXML_body_index": i, "old": text(a), "new": text(b)})
    hashes = lambda root, expr: [nodehash(n) for n in root.xpath(expr, namespaces=NS)]
    old_refs = [nodehash(n) for n in old.xpath("//w:body/w:p", namespaces=NS) if re.match(r"^\[\d+\]", text(n))]
    new_refs = [nodehash(n) for n in new.xpath("//w:body/w:p", namespaces=NS) if re.match(r"^\[\d+\]", text(n))]
    anchor_tags = ["commentRangeStart", "commentRangeEnd", "commentReference", "bookmarkStart", "bookmarkEnd"]
    anchors = lambda root: {tag: [dict(n.attrib) for n in root.xpath("//w:" + tag, namespaces=NS)] for tag in anchor_tags}
    changed_parts = sorted(n for n in old_parts.keys() & new_parts.keys() if old_parts[n] != new_parts[n])
    media = [n for n in old_parts if n.startswith("word/media/")]
    comments = [n for n in old_parts if "comment" in n.lower()]
    old_changed_indices = [49, 53, 90, 95, 103]
    new_background_prefixes = [
        "结构实例考察主链活动区域定位", "侧链实例使用泛素与EDD",
        "转录共激活因子CBP", "KRAS是一种", "泛素是由76个氨基酸",
    ]
    replacement_mapping = []
    normalized_new_hashes = new_hashes.copy()
    for old_index, prefix in zip(old_changed_indices, new_background_prefixes):
        matches = [i for i, n in enumerate(new_body) if text(n).startswith(prefix)]
        if len(matches) != 1:
            raise ValueError(f"Expected one revised paragraph for prefix {prefix!r}: {matches}")
        new_index = matches[0]
        normalized_new_hashes[new_index] = old_hashes[old_index]
        replacement_mapping.append({"old_index": old_index, "new_index": new_index, "new_text": text(new_body[new_index])})
    old_caption_index = next(i for i, n in enumerate(old_body) if text(n).startswith("图6："))
    old_figure_index = old_caption_index - 1
    fig_hash, cap_hash = old_hashes[old_figure_index], old_hashes[old_caption_index]
    figure_caption_preserved = new_hashes.count(fig_hash) == 1 and new_hashes.count(cap_hash) == 1
    if figure_caption_preserved:
        reduced_old = [h for h in old_hashes if h not in {fig_hash, cap_hash}]
        reduced_new = [h for h in normalized_new_hashes if h not in {fig_hash, cap_hash}]
        order_preserved_other_than_fig6 = reduced_old == reduced_new
        new_figure_index, new_caption_index = new_hashes.index(fig_hash), new_hashes.index(cap_hash)
    else:
        order_preserved_other_than_fig6 = False
        new_figure_index = new_caption_index = None
    checks = {
        "source_SHA_matches_final_v41": sha(OLD.read_bytes()) == "58a15ba884180814dc562d276a3fb1c0e321b95a73cc302e0b0e68a4dafc0a68",
        "only_document_XML_changed": changed_parts == ["word/document.xml"] and old_parts.keys() == new_parts.keys(),
        "body_element_count_unchanged": len(old_body) == len(new_body),
        "only_five_authorized_paragraphs_changed_and_Fig6_caption_moved": order_preserved_other_than_fig6,
        "Figure6_and_caption_original_XML_identical": figure_caption_preserved,
        "all_8_visible_figures_preserved": len(new.xpath("//w:drawing", namespaces=NS)) == 8 and hashes(old, "//w:drawing") == hashes(new, "//w:drawing"),
        "all_media_files_byte_identical": all(old_parts[n] == new_parts[n] for n in media),
        "all_10_tables_identical": len(new.xpath("//w:tbl", namespaces=NS)) == 10 and hashes(old, "//w:tbl") == hashes(new, "//w:tbl"),
        "all_36_native_equations_identical": len(new.xpath("//m:oMath", namespaces=NS)) == 36 and hashes(old, "//m:oMath") == hashes(new, "//m:oMath"),
        "all_comments_parts_byte_identical": all(old_parts[n] == new_parts[n] for n in comments),
        "all_comment_bookmark_anchors_identical": anchors(old) == anchors(new),
        "all_28_references_identical": len(new_refs) == 28 and old_refs == new_refs,
    }
    report = {
        "status": "pass" if all(checks.values()) else "fail",
        "document_sha256": sha(NEW.read_bytes()),
        "source": str(OLD), "candidate": str(NEW), "candidate_SHA256": sha(NEW.read_bytes()),
        "checks": checks, "all_package_checks_pass": all(checks.values()),
        "changed_parts": changed_parts, "changed_paragraphs": changes,
        "authorized_background_replacements": replacement_mapping,
        "Figure6_reordering": {"old_figure_index": old_figure_index, "old_caption_index": old_caption_index, "new_figure_index": new_figure_index, "new_caption_index": new_caption_index},
        "note_on_indices": "body_index为OOXML w:body直接子元素索引，与python-docx段落列表索引可能不同。",
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "v42_package_integrity_audit.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({"all_package_checks_pass": all(checks.values()), "checks": checks, "candidate_SHA256": report["candidate_SHA256"]}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
