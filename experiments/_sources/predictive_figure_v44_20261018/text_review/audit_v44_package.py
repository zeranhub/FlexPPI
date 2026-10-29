"""Read-only, source-bound OOXML preservation and reference audit for v44."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import re
from zipfile import ZipFile
from lxml import etree

ROOT = Path("E:/Project/FlexPPI/experiments/_sources/predictive_figure_v44_20261018")
SOURCE = ROOT.parent.parent / "paper_v43_20261018_中文修订_能力与构象合并图_三线表.docx"
SOURCE_SHA = "5d32c2af90a6aa8268ac0ef50ba46480370fc798bcc890d91fa71123dfb51c5a"
NS = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main", "a": "http://schemas.openxmlformats.org/drawingml/2006/main", "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships", "m": "http://schemas.openxmlformats.org/officeDocument/2006/math"}


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def nh(node):
    return sha(etree.tostring(node, method="c14n"))


def text(node):
    return "".join(node.xpath(".//w:t/text()", namespaces=NS))


def read(path):
    with ZipFile(path) as z:
        parts = {n: z.read(n) for n in z.namelist()}
    tree = etree.fromstring(parts["word/document.xml"])
    body = tree.find("w:body", NS)
    rels = etree.fromstring(parts["word/_rels/document.xml.rels"])
    relmap = {n.get("Id"): n.get("Target") for n in rels}
    images = []
    for i, node in enumerate(body):
        for rid in node.xpath(".//a:blip/@r:embed", namespaces=NS):
            part = "word/" + relmap[rid]
            images.append({"body_index": i, "rid": rid, "part": part, "sha256": sha(parts[part])})
    return {"path": str(path), "sha256": sha(path.read_bytes()), "parts": parts, "tree": tree, "body": body, "images": images}


def hashes(doc, query):
    return [nh(n) for n in doc["tree"].xpath(query, namespaces=NS)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--new", required=True, type=Path)
    ap.add_argument("--figure2", required=True, type=Path)
    ap.add_argument("--output", required=True, type=Path)
    args = ap.parse_args()
    old, new = read(SOURCE), read(args.new)
    proposal = json.loads((ROOT / "text_review/manuscript_proposal.json").read_text(encoding="utf-8"))
    prefix_edits = [p["prefix"] for p in proposal["replacements"]] + ["图2："]
    old_nodes = []
    for node in old["body"]:
        tx = text(node)
        if any(tx.startswith(p) for p in prefix_edits):
            continue
        # The root builder may change display dimensions for the replaced figure.
        if any(i["body_index"] == list(old["body"]).index(node) and i["part"] == "word/media/image2.png" for i in old["images"]):
            continue
        old_nodes.append({"text": tx, "hash": nh(node)})
    observed = Counter(nh(n) for n in new["body"])
    missing = []
    for item in old_nodes:
        if observed[item["hash"]] > 0:
            observed[item["hash"]] -= 1
        else:
            missing.append(item)
    anchors = lambda d: {k: [dict(n.attrib) for n in d["tree"].xpath("//w:" + k, namespaces=NS)] for k in ["commentRangeStart", "commentRangeEnd", "commentReference", "bookmarkStart", "bookmarkEnd"]}
    refs = lambda d: [nh(n) for n in d["body"] if re.match(r"^\[\d+\]", text(n))]
    op, np = old["parts"], new["parts"]
    changed = sorted(k for k in op.keys() & np.keys() if op[k] != np[k])
    new_texts = [text(n) for n in new["body"]]
    matches = []
    for rep in proposal["replacements"]:
        matches.append({"prefix": rep["prefix"], "count": new_texts.count(rep["new"]), "old_removed": rep["old"] not in new_texts if rep["old"] != rep["new"] else True})
    inserted = [{"prefix": ins["text"][:30], "count": new_texts.count(ins["text"])} for ins in proposal["insertions"]]
    fig2 = next((s for s in new_texts if s.startswith("图2：")), "")
    caption = (ROOT / "figure2/caption_zh.md").read_text(encoding="utf-8").strip()
    old_image_hashes = [x["sha256"] for x in old["images"]]
    expected_image_hashes = list(old_image_hashes)
    assert len(expected_image_hashes) == 6
    expected_image_hashes[1] = sha(args.figure2.read_bytes())
    checks = {
        "source_v43_SHA_matches": old["sha256"] == SOURCE_SHA,
        "10_native_three_line_tables_identical": len(hashes(new, "//w:tbl")) == 10 and hashes(old, "//w:tbl") == hashes(new, "//w:tbl"),
        "36_OMML_equations_identical": len(hashes(new, "//m:oMath")) == 36 and hashes(old, "//m:oMath") == hashes(new, "//m:oMath"),
        "28_references_identical": len(refs(new)) == 28 and refs(old) == refs(new),
        "comments_parts_identical": all(np.get(k) == raw for k, raw in op.items() if "comment" in k.lower()),
        "all_comment_bookmark_anchors_identical": anchors(old) == anchors(new),
        "only_document_and_Fig2_media_part_changed": set(changed) == {"word/document.xml", "word/media/image2.png"},
        "no_package_parts_added_or_removed": op.keys() == np.keys(),
        "all_other_body_nodes_exactly_preserved": not missing,
        "all_replacement_paragraphs_present_once": all(v["count"] == 1 and v["old_removed"] for v in matches),
        "all_insertion_paragraphs_present_once": all(v["count"] == 1 for v in inserted),
        "Fig2_caption_matches_reviewed_asset": fig2 == caption,
        "six_images_same_order_and_Fig5_unchanged": [x["sha256"] for x in new["images"]] == expected_image_hashes,
        "five_numbered_main_captions": [re.match(r"^图([1-9])：", s).group(1) for s in new_texts if re.match(r"^图[1-9]：", s)] == list("12345"),
        "removed_Fig2_panels_not_described_as_active": not any("图2e、f进一步按实际输入链数" in s or "MAE经验累积分布" in s for s in new_texts),
    }
    report = {"status": "passed" if all(checks.values()) else "failed", "document": str(args.new), "document_sha256": new["sha256"], "source_document_sha256": old["sha256"], "checks": checks, "changed_parts": changed, "missing_protected_body_nodes": missing, "replacement_matches": matches, "insertions": inserted, "active_images": new["images"], "figure2_sha256": sha(args.figure2.read_bytes())}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({"status": report["status"], "document_sha256": new["sha256"], "checks": checks}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
