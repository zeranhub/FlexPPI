import json,hashlib
from pathlib import Path
root=Path(r"E:\Project\FlexPPI\experiments\_sources\mel5_composite_v54_20261022")
doc=Path(r"E:\Project\FlexPPI\paper\revision_v14_20261010\paper_v54_20261022_中文修订_MEL5六面板组合图4_三线表.docx")
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
dhash=sha(doc)
assert dhash=="90b66a4029e22b0f67e52d3176839e547a8265690e5ff21d385ab7b5fecd1027"
notes={
1:"Title, abstract and introduction readable; no clipping or missing characters.",
2:"Method section heading and text render normally; page-break whitespace intentional.",
3:"Figure1 all panels and complete Chinese caption visible; subsequent prose does not overlap.",
4:"Equations1/2, radicals, subscripts, numbering and surrounding Chinese prose render normally.",
5:"Equations3/4 and training-method prose render normally; no right-margin clipping.",
6:"Methods2.5/2.6 text, named methods, units and inline mathematical notation readable.",
7:"Metrics and KIX MD methods readable; no table/figure collision.",
8:"Revised MEL5 methods and equation5 readable; formula and figure-reference text not clipped.",
9:"Continuation of methods and result headings render normally.",
10:"Table1 uses intact native three-line layout; full columns, rows, caption and table notes visible.",
11:"Figure2 all panels, legends, colorbar, table and complete caption visible; no clipping.",
12:"Results subsections render normally; whitespace follows existing figure pagination.",
13:"Figure3 entire multi-panel composite and complete caption visible; no cropping or overlap.",
14:"Revised MEL5 background, model results and contact-mechanism paragraphs readable; no missing characters.",
15:"New Figure4 all six panels fully visible, including four-state peptide legend, all23 bars, all45 CDR positions and both heatmaps. Complete Chinese caption appears on same page and defines gray cells and contrast baselines. Following section3.5 starts clearly without collision."
}
pages=[]
for n in range(1,16):
    p=root/"manuscript_QA"/"pages"/f"page-{n}.png"
    assert p.exists()
    pages.append({"page":n,"path":str(p),"sha256":sha(p),"review_method":"individually opened with view_image","status":"passed","finding":notes[n]})
report={
"status":"passed",
"reviewer_task":"/root/v54_layout_review",
"scope":"Independent rendered layout review of manuscript pages1–15; each page opened individually.",
"document_path":str(doc),
"document_sha256":dhash,
"expected_document_sha256":dhash,
"pages_reviewed":list(range(1,16)),
"page_count_reviewed":15,
"checks":{"all_pages_opened_individually":True,"missing_characters_detected":False,"cropping_detected":False,"overlaps_detected":False,"three_line_table_page10_intact":True,"equations_intact":True,"figure4_six_panels_intact":True,"figure4_caption_same_page":True},
"required_changes":[],
"pages":pages
}
out=root/"manuscript_QA"/"independent_pages_1_15_review.json"
out.write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding="utf-8")
print(json.dumps({"status":"passed","doc_sha":dhash,"pages":len(pages),"report":str(out)},ensure_ascii=False))
