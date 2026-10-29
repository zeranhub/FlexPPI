import json,hashlib
from pathlib import Path
root=Path(r"E:\Project\FlexPPI\experiments\_sources\mel5_compact_v55_20261022")
doc=Path(r"E:\Project\FlexPPI\paper\revision_v14_20261010\paper_v55_20261022_中文修订_MEL5精简联合图4_三线表.docx")
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
dhash=sha(doc)
assert dhash=="056c2500e95d3992cf25f31980751b9936587c3dbe9faee6cf9738c88388153c"
notes={
1:"Title, abstract and introduction readable; no clipping, overlapping text or missing characters.",
2:"Methods headings and text render normally; pagination whitespace intact.",
3:"Figure1 all panels and complete Chinese caption visible; subsequent text does not overlap figure.",
4:"Equations1/2, radicals, subscripts and formula numbers render normally.",
5:"Equations3/4 and training methods readable without right-margin clipping.",
6:"Methods2.5/2.6, named methods, units and inline mathematical notation readable.",
7:"Metrics and KIX MD methods render without missing characters or overlap.",
8:"MEL5 methods and equation5 visible; revised case geometry and comparison definitions remain readable.",
9:"Methods continuation and results headings render normally.",
10:"Table1 remains complete native three-line table with all columns, rows, caption and table notes visible.",
11:"Figure2 all panels, legends, colorbar, numeric table and caption visible without clipping.",
12:"Results subsections readable; pagination whitespace normal.",
13:"Figure3 complete multi-panel figure and caption visible; no cropping.",
14:"MEL5 text fits page. References correctly use Figure4d upper/lower instead of old e/f. No prose/formatting damage.",
15:"Figure4 has only a–d. Former ranked bars removed. New d contains full upper CDR delta-RMSF and lower pMHC delta-Hchi maps, all three contrast rows and separate colorbars. Complete updated caption is on same page and describes d upper/lower, units, gray cells and comparison baselines. Subsequent section3.5 starts below without collision."
}
pages=[]
for n in range(1,16):
    p=root/"manuscript_QA"/"pages"/f"page-{n}.png"
    assert p.exists()
    pages.append({"page":n,"path":str(p),"sha256":sha(p),"review_method":"individually opened with view_image","status":"passed","finding":notes[n]})
report={
"status":"passed",
"reviewer_task":"/root/v54_layout_review",
"scope":"Independent visual review of rendered v55 manuscript pages1–15, each opened individually.",
"document_path":str(doc),
"document_sha256":dhash,
"expected_document_sha256":dhash,
"pages_reviewed":list(range(1,16)),
"page_count_reviewed":15,
"checks":{"all_pages_opened_individually":True,"missing_characters_detected":False,"cropping_detected":False,"overlaps_detected":False,"three_line_table_page10_intact":True,"equations_intact":True,"figure4_exactly_four_panels":True,"ranked_bar_removed":True,"figure4_d_upper_lower_heatmaps":True,"figure4_caption_same_page":True,"old_e_f_references_absent_in_case_prose":True},
"required_changes":[],
"pages":pages
}
out=root/"manuscript_QA"/"independent_pages_1_15_review.json"
out.write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding="utf-8")
print(json.dumps({"status":"passed","doc_sha":dhash,"pages":len(pages),"report":str(out)},ensure_ascii=False))
