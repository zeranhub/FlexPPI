"""Record separately completed single-page v55 visual review; no source edits."""
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parent.parent
DOC=ROOT.parent.parent/'paper_v55_20261022_中文修订_MEL5精简联合图4_三线表.docx'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
assert sha(DOC)=='056c2500e95d3992cf25f31980751b9936587c3dbe9faee6cf9738c88388153c'
notes={
15:'New Figure4 contains only four labeled panels a–d. Original ranked23 Hchi bar is absent. Upper and lower heatmaps form single paneld with separate RMSF/Hchi colorbars and units. Local atom labels, peptide state legend, six-CDR curves and both matrix label rows fit. Full updated Chinese caption is on the same page, followed by section3.5 and opening paragraph; no overlap or clipping.',
16:'Sections3.5.1/3.5.2 render fully and legibly with chi, rho, superscripts and degree symbols; page ends with complete paragraph and no overflow.',
17:'Existing Figure5 all six subpanels and complete caption fit within margins with no missing labels, crop or overlap.',
18:'Discussion, conclusion and code/data sections render completely inside margins, with no glyph loss or collision.',
19:'References1–23 render as readable uninterrupted entries inside page boundaries.',
20:'References24–28 render correctly; existing whitespace separates the bibliography from supplementary material.',
21:'Supplementary heading and source text plus three-line TableS1 are intact, aligned and unclipped.',
22:'Supplementary FigureS1 all six subpanels, legends, complete caption and following text fit without clipping.',
23:'Supplementary training text and three-line TablesS2/S3 fit; numeric columns and lower explanatory text remain aligned.',
24:'Supplementary baseline text, three-line TableS4 and note render fully inside margins.',
25:'Supplementary family-split text and all rows of three-line TableS5 render completely.',
26:'Supplementary validation text, three-line TablesS6/S7 and lower explanatory paragraph fit without clipping or overlap.',
27:'Supplementary KIX text, three-line TableS8 and explanatory paragraph render fully.',
28:'Supplementary paired-readout text and three-line TableS9 render completely, including all numeric rows and lower rule.',
}
pages=[]
for n,note in notes.items():
    p=ROOT/f'manuscript_QA/pages/page-{n}.png'
    assert p.exists()
    pages.append({'page':n,'path':str(p),'PNG_SHA256':sha(p),'individually_opened_with_view_image':True,'status':'passed','observations':note,'blocking_issues':[]})
report={'status':'passed','DOCX_path':str(DOC),'DOCX_SHA256':sha(DOC),'reviewed_pages':list(range(15,29)),'review_type':'Each page opened separately with view_image and visually inspected; no contact-sheet-only review.','read_only_review':True,'rendered_page_count_total':28,'all_reviewed_pages_inside_margins':True,'Figure4_and_complete_caption_on_same_page':15,'required_changes':[],'pages':pages}
(ROOT/'manuscript_QA/independent_pages_15_28_review.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({'status':'passed','pages_reviewed':len(pages),'DOCX_SHA256':sha(DOC),'file':str(ROOT/'manuscript_QA/independent_pages_15_28_review.json')}))
