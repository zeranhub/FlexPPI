"""Record completed individual visual inspection, without changing rendered artifacts."""
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parent.parent
DOC=ROOT.parent.parent/'paper_v54_20261022_中文修订_MEL5六面板组合图4_三线表.docx'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
assert sha(DOC)=='90b66a4029e22b0f67e52d3176839e547a8265690e5ff21d385ab7b5fecd1027'
notes={
8:'Methods retain all four structures, three contrasts, 30-site union, 23/23/22 eligibility, common numbering, fixed atom pairs and MD geometric criteria. Formula5 and page bottom remain inside margins with no missing glyphs.',
15:'Full six-panel Figure4 is present; all a–f panels, local atomic labels, four-state peptide legend, CDR curves, ranked23 labels, both heatmap/colorbars and full Chinese caption fit the page with no overlap or clipping. Next section heading and opening lines fit below.',
16:'Body sections3.5.1/3.5.2 render with complete Greek/chi/degree symbols and text; no overflow.',
17:'Unchanged Figure5 all a–f panels and full caption fit on one page; charts, labels and equations visible with no truncation.',
18:'Discussion, conclusion and code/data sections are complete and inside page margins; no missing symbols or collisions.',
19:'References1–23 render in continuous readable entries within page boundaries.',
20:'References24–28 render correctly. Remaining whitespace is the existing separation before supplementary material.',
21:'Supplementary heading, data-source paragraph and three-line TableS1 are intact, aligned and unclipped.',
22:'Supplementary FigureS1 all a–f panels, legend, caption and accompanying paragraph fit without cropping.',
23:'Supplementary training text and three-line TablesS2/S3 fit; row labels and numeric columns remain aligned.',
24:'Supplementary baseline text and three-line TableS4 fit, including note beneath table.',
25:'Supplementary family-split text and three-line TableS5 render fully with aligned columns and intact lower rule.',
26:'Supplementary prediction-driven validation text and TablesS6/S7 fit; lower explanatory paragraph remains above footer.',
27:'Supplementary KIX text and three-line TableS8 render fully, including explanatory paragraph.',
28:'Supplementary paired-readout text and three-line TableS9 render fully with all rows, numbers and lower rule.',
}
rows=[]
for number,note in notes.items():
    p=ROOT/f'manuscript_QA/pages/page-{number}.png'
    assert p.exists()
    rows.append({'page':number,'path':str(p),'SHA256':sha(p),'individually_opened_with_view_image':True,'status':'passed','observations':note,'blocking_issues':[]})
report={'status':'passed','DOCX_path':str(DOC),'DOCX_SHA256':sha(DOC),'reviewed_pages':[8]+list(range(15,29)),'review_type':'Each page opened separately with view_image and visually inspected; no contact-sheet-only review.','read_only_review':True,'rendered_page_count_total':28,'all_reviewed_pages_inside_margins':True,'Figure4_and_caption_on_same_page':15,'formulas_figures_three_line_tables_preserved':True,'required_changes':[],'pages':rows}
(ROOT/'manuscript_QA/independent_pages_15_28_review.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({'status':'passed','pages':len(rows),'file':str(ROOT/'manuscript_QA/independent_pages_15_28_review.json')}))
