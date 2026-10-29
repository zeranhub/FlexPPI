from pathlib import Path
from datetime import datetime,timezone
from hashlib import sha256
from zipfile import ZipFile
import json
from pypdf import PdfReader
R=Path(r'E:\Project\FlexPPI\experiments\_sources\mel5_tcr_binding_20261028\manuscript\v128')
P=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010\paper_v128_20261028_中文修订_单点突变与残基动态叙事.docx')
source=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010\paper_v127_20261028_中文修订_TCR结合模型对照.docx')
change=json.loads((R/'revision_changes.json').read_text(encoding='utf-8'))
assert sha256(P.read_bytes()).hexdigest()==change['target_sha256']
assert sha256(source.read_bytes()).hexdigest()==change['source_sha256']
page_info=json.loads((R/'render/page_text.json').read_text(encoding='utf-8'))
changed=[x['page'] for x in page_info if not x['pixel_identical_to_previous_reviewed_page']]
assert changed==[13,15,34,35]
assert len(PdfReader(next((R/'render').glob('*.pdf'))).pages)==35
qa={'status':'passed','reviewed_utc':datetime.now(timezone.utc).isoformat(),'document':str(P),'document_sha256':change['target_sha256'],'pages':35,'new_or_changed_pages_visually_reviewed':[13,15,34,35],'remaining_31_pages':'Pixel-identical to individually reviewed v127 pages; prior QA lineage retained.','previous_QA':str(R.parent/'render_QA.json'),'visual_checks':{'paragraph_flow_and_page_breaks':'pass','Chinese_Greek_and_math_glyphs':'pass','supplement_G_table_and_figure_captions':'pass','clipping_or_overlapping_text':'none observed'},'content_review':'Independent agent verified proposed section3.4 paragraphs against existing experiments and model audits; passed. I7A supporting evidence and G4A disagreement retained.','structural_checks':{'only_original_package_part_changed':'word/document.xml','all_50_media_binary_unchanged':True,'all_existing_tables_fields_comments_sections_drawings_unchanged':True,'original_dynamic_and_contact_evidence_unchanged':True},'page_hashes':[{'page':x['page'],'sha256':sha256((R/'render'/f"page-{x['page']}.png").read_bytes()).hexdigest(),'review':'visually inspected' if x['page'] in changed else 'pixel-identical to previous inspected page'} for x in page_info]}
(R/'render_QA.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8')
change['render_validation']='passed: 35 pages; modified pages13,15,34,35 visually reviewed; other31 pixel-identical to prior reviewed pages'
change['independent_content_review']='passed; numerical evidence, predictor scope, training coverage caveat and hypothesis language verified'
(R/'revision_changes.json').write_text(json.dumps(change,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({'target':str(P),'pages':35,'review':'passed','original_media_preserved':50},ensure_ascii=False))
