from pathlib import Path
from zipfile import ZipFile,ZIP_DEFLATED
from pypdf import PdfReader
import json,hashlib
R=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010');V=R/'section33_expansion_20261025/model_first_story_search_20261025/revision_v87_general_story_compact'
fig=V/'figure3/fig3_joint_flexibility_compact';doc=R/'paper_v87_20261026_中文修订_互补柔性背景与图3凝练.docx'
qa=json.loads((V/'figure3/figure_QA.json').read_text());qa['visual_review']='passed: full-resolution composition, independent count/value and structure review, and final Word-scale page11 review'
qa['font_embedding']='ArialMT and Arial-BoldMT embedded TrueType subsets; pdffonts verified'
qa['manual_style_resolutions']={'AP-0':'Mandatory verbatim palette uses DIVERGING with three spaces, while the validator requires a single space; all baseline values match.','AP-4':'Legends occupy dedicated figure-coordinate white space outside data axes.','AP-6':'Contact bars represent exact per-residue counts; shielding bars represent one static descriptor per residue. Neither is a replicate group mean requiring a jitter layer.'}
(V/'figure3/figure_QA.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8')
pdf=PdfReader(str(fig)+'.pdf');w=float(pdf.pages[0].mediabox.width)*25.4/72;h=float(pdf.pages[0].mediabox.height)*25.4/72;assert abs(w-183)<.02 and abs(h-190)<.02
with ZipFile(doc) as z:assert z.read('word/media/image5.png')==Path(str(fig)+'.png').read_bytes()
finalpdf=list((V/'render/release_final').glob('*.pdf'))[0];assert len(PdfReader(finalpdf).pages)==29
plan=json.loads((V/'manuscript/v87_text_plan.json').read_text(encoding='utf-8'));assert plan['paragraph_replacements']['68']=='3.3 主链与侧链柔性的互补解析'
assert 'CaM' not in plan['paragraph_replacements']['68'] and plan['insert_before_original_paragraph_69']
audit={'status':'complete','source_docx':'paper_v86_20261026_中文修订_图3完整接触与侧链环境.docx','final_docx':str(doc),'docx_sha256':hashlib.sha256(doc.read_bytes()).hexdigest(),'pages':29,'section_title':plan['paragraph_replacements']['68'],'inserted_background_before_original_paragraph':69,'main_figure_mm':[w,h],'previous_height_mm':225,'height_reduction_percent':(225-h)/225*100,'same_fontsize_floor_as_v86':6.2,'contact_summary_uses_all3600_pairs':True,'all25_contact_rows_retained':True,'all21_sidechains_retained':True,'all92_profile_values_retained':True,'full_contact_matrix_preserved':'figure3/source_contact_map/all25x144_contacts.pdf + source_data/all3600_CaMKI_CaM_distances.csv','new_model_or_MD_runs':False,'new_structural_asset':'True same-view M316-L45/S99 local magnification from 1A06','changed_zip_parts':['word/document.xml','word/media/image5.png'],'all_other_images_and_six_tables_unchanged':True,'visual_QA':'First9 pages exact pixels of audited v86; independent pages10–19 and19–29 reviewed; final two caption precision edits affect page11 only, reviewed primary, remaining28 pages exact pixels of reviewed release01.','main_figure_sha256':qa['sha256'],'caption_precision_edits':['25 backbone positions and21 sidechain positions','1A06 residues317–318 unresolved, not biologically deleted']}
(V/'completion_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
zipout=R/'paper_v87_20261026_互补柔性与紧凑图3_源文件与PyMOL.zip'
folders=['figure3','pymol','pymol_context','pymol_context_local','data_audit','supplement_1MXE','manuscript']
files=[]
for folder in folders:files.extend(p for p in (V/folder).rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.suffix!='.pyc')
for name in ['README_交付与复绘.md','completion_audit.json','build_v87_manuscript.py','design_review.md','final_figure_review.md','narrative/3.3_定稿.md','narrative/paragraph_replacements.json','narrative/final_text_and_pages10_18_review.md','render/v86_v87_pixel_comparison.json','render/final_caption_precision_pixel_comparison.json','render/render_docx_word_wrapper.py','render/export_word_readonly.ps1']:files.append(V/name)
with ZipFile(zipout,'w',ZIP_DEFLATED) as z:
 for p in sorted(set(files)):z.write(p,Path(V.name)/p.relative_to(V))
with ZipFile(zipout) as z:assert z.testzip() is None
print(json.dumps({'docx_sha256':audit['docx_sha256'],'source_zip':str(zipout),'source_files':len(set(files)),'zip_bytes':zipout.stat().st_size,'zip_sha256':hashlib.sha256(zipout.read_bytes()).hexdigest()},ensure_ascii=False))