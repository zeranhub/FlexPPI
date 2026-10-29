from pathlib import Path
from zipfile import ZipFile,ZIP_DEFLATED
from pypdf import PdfReader
import json,hashlib
R=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010');V=R/'section33_expansion_20261025/model_first_story_search_20261025/revision_v88_two_panel_story'
fig=V/'figure3/fig3_two_panel_story';doc=R/'paper_v88_20261026_中文修订_图3双面板主线.docx'
qa=json.loads((V/'figure3/figure_QA.json').read_text());qa['visual_review']='passed: primary full figure and Word page11; independent main figure,92 values and pages10–19; independent supplementary pages19–29'
qa['font_embedding']='ArialMT and Arial-BoldMT TrueType subsets embedded, pdffonts verified'
qa['fontsize_note']='7.2pt minimum at183mm master; approximately6.45pt after164mm placement inWord, visually checked'
qa['manual_style_resolutions']={'AP-0':'Mandatory palette literal has three spaces after DIVERGING; validator falsely requires one. All required baseline values unchanged.','AP-4':'Figure-coordinate numerical and structural legends are placed in dedicated white space outside data axes.'}
(V/'figure3/figure_QA.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8')
pdf=PdfReader(str(fig)+'.pdf');w=float(pdf.pages[0].mediabox.width)*25.4/72;h=float(pdf.pages[0].mediabox.height)*25.4/72;assert abs(w-183)<.02 and abs(h-100)<.02
with ZipFile(doc) as z:
 assert z.testzip() is None
 assert z.read('word/media/image5.png')==Path(str(fig)+'.png').read_bytes()
renderpdf=list((V/'render/release_01').glob('*.pdf'))[0];assert len(PdfReader(renderpdf).pages)==29
text=(V/'manuscript/v88_final_text.txt').read_text(encoding='utf-8');assert all(x not in text for x in ['图3c','图3d','图3e','图3f','图3a、c'])
assert '3.3 主链与侧链柔性的互补解析' in text
audit={'status':'complete','final_docx':str(doc),'docx_sha256':hashlib.sha256(doc.read_bytes()).hexdigest(),'pages':29,'main_figure_mm':[w,h],'previous_main_height_mm':190,'height_and_area_reduction_percent':(190-h)/190*100,'main_panels':['a structural location','b complete predicted RMSF and Hchi profiles'],'retired_independent_panels':['d lobe-contact counts','e Hchi paired sidechain display with shielding','f autoinhibitory context','intermediate M316 environment comparison'],'retired_data_preserved':True,'all25_backbone_and21_sidechain_positions':True,'all92_profile_values_match_original_family_fold1':True,'all_axes_start_at_zero':True,'full_intra_same_bound_input_not_biological_states':True,'changed_zip_parts':['word/document.xml','word/media/image5.png'],'all_six_tables_and_other_images_identical_to_v87':True,'new_model_or_MD_runs':False,'all_pages_QA':'1–9 exact pixel identity to visually audited v87; actual primary and independent10–19; independent19–29','main_png_sha256':qa['sha256']['png'],'final_references':'Figure3 a/b only throughout'}
(V/'completion_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
zipout=R/'paper_v88_20261026_图3双面板主线_源文件与PyMOL.zip'
files=[]
for folder in ['figure3','pymol','data_audit','supplement_1MXE','manuscript','source_analysis']:
 files.extend(p for p in (V/folder).rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.suffix!='.pyc')
for name in ['README_交付与复绘.md','completion_audit.json','build_v88_manuscript.py','narrative/3.3_定稿.md','narrative/paragraph_replacements.json','narrative/final_main_text_figure_and_pages_review.md','render/render_docx_word_wrapper.py','render/export_word_readonly.ps1']:files.append(V/name)
with ZipFile(zipout,'w',ZIP_DEFLATED) as z:
 for p in sorted(set(files)):z.write(p,Path(V.name)/p.relative_to(V))
with ZipFile(zipout) as z:assert z.testzip() is None
print(json.dumps({'word':str(doc),'figure':str(fig)+'.pdf','zip':str(zipout),'source_files':len(set(files)),'zip_bytes':zipout.stat().st_size,'docx_sha256':audit['docx_sha256'],'zip_sha256':hashlib.sha256(zipout.read_bytes()).hexdigest()},ensure_ascii=False))