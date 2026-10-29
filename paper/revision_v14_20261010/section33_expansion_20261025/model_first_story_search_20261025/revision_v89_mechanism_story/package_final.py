from pathlib import Path
from zipfile import ZipFile,ZIP_DEFLATED
from pypdf import PdfReader
import json,hashlib
R=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010');V=R/'section33_expansion_20261025/model_first_story_search_20261025/revision_v89_mechanism_story';doc=R/'paper_v89_20261026_中文修订_图3调控机制假说.docx';fig=V/'figure3/fig3_mechanism_story'
qa=json.loads((V/'figure3/figure_QA.json').read_text(encoding='utf-8'))
qa['visual_review']='Passed: final PNG independently viewed; Word pages10–13 primary,14–18 independent,19–29 independent; final corrected page24 primary; pages1–9 exact PNG identity to prior audited release.'
qa['font_embedding']='ArialMT and Arial-BoldMT embedded TrueType subsets, verified using pdffonts.'
qa['font_size_in_Word164mm_pt']=7.9*164/183
qa['manual_style_resolutions']={'AP-0':'Mandatory palette baseline preserved verbatim; validator spacing false positive.','AP-4':'Figure-coordinate legend outside both data axes.','CL-2':'PDF mediabox182.999999999995x178.000000000015mm.'}
(V/'figure3/figure_QA.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8')
with ZipFile(doc) as z:
 assert z.testzip() is None
 assert z.read('word/media/image5.png')==Path(str(fig)+'.png').read_bytes()
sourceqa=json.loads((V/'figure3/independent_data_QA.json').read_text(encoding='utf-8'));assert sourceqa['all92_match_frozen_exactly'] and sourceqa['PNG_sha256_reviewed']==qa['sha256']['png']
renderpdf=next((V/'render/release_02').glob('*.pdf'));assert len(PdfReader(renderpdf).pages)==29
assert all((V/'render/release_02'/f'page-{n}.png').exists() for n in range(1,30))
text=(V/'manuscript/v89_final_text.txt').read_text(encoding='utf-8');assert 'Y215' in text and '3.3 主链与侧链柔性的互补解析' in text
build=json.loads((V/'manuscript/v89_build_audit.json').read_text(encoding='utf-8'))
audit={'status':'complete','final_docx':str(doc),'docx_sha256':hashlib.sha256(doc.read_bytes()).hexdigest(),'pages':29,'main_figure_mm':[183,178],'main_figure_png_sha256':qa['sha256']['png'],'panels':['a observed regulatory contacts and CaM anchor context','b full absolute RMSF and Hchi predictions','c explicit mechanistic hypothesis and occupancy-matched test'],'all92_values_exactly_match_frozen_family_fold1':True,'source_data_sha256':sourceqa['local_source_sha256'],'same_bound_input_not_free_vs_bound':True,'mechanism_not_directly_measured':True,'new_model_or_MD_runs':False,'protected_other_images_and6tables':True,'changed_zip_parts':build['changed_zip_parts'],'changed_paragraphs':build['modified_paragraphs'],'all_pages_QA':'1–9 exact prior pixel identity;10–13 primary;14–18 independent;19–29 independent; corrected24 primary; final02 28/29 pages exact reviewedrelease01 identity','closed_review_issue':'F298 4.5A neighbour Y215 added to paragraph163; complete six-neighbour distance CSV/JSON retained; inset highlights five nearer neighbours.'}
(V/'completion_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
# Include short provenance and primary-source links, not copies of third-party full articles/SI.
claims='''# 机制证据对应关系

真实实验与结构背景：Clapperton等2002的摘要明确K296/F298的自抑制/CaM结合序列重叠、W303接C叶瓣、M316接N叶瓣；C叶瓣结合紧、N叶瓣对激活必需。官方SI图S3显示孤立调控肽可抑制CaMKI293。这支持调控段具有两类接触能力，未测量回接动力学。来源：https://doi.org/10.1021/bi026660t ；官方SI https://acs.figshare.com/articles/journal_contribution/3586317 。

自抑制催化核心背景来自Goldberg等1996及1A06：https://doi.org/10.1016/S0092-8674(00)81066-1 ；https://www.rcsb.org/structure/1A06 。结合调控肽结构1MXE：https://www.rcsb.org/structure/1MXE 。

模型观测：Full在W303–M316定义的连续14位点平均RMSF0.513Å，21个有效侧链Hχ呈位置差异。Full/Intra是同一结合态输入上的两种训练模型。完整读数及MD复核见源数据。

新机制假说：双锚点约束可能减少调控段抑制性回接，尚无直接回接事件或激活测量；图3c明确标Hypothesis。占据率匹配后扰动N叶–M316并联测局部运动和活性，是未来验证设计，不是本次已做实验。残余螺旋、先C后N、尾部牵绳功能和构象选择路径没有写成证据。
'''
(V/'narrative/primary_sources_and_claims.md').write_text(claims,encoding='utf-8')
zipout=R/'paper_v89_20261026_图3调控机制假说_源文件与PyMOL.zip'
files=[]
for folder in ['figure3','pymol_mechanism','data_audit','data_audit_only','supplement_1MXE','source_analysis','manuscript']:
 files.extend(p for p in (V/folder).rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.suffix!='.pyc')
for name in ['README_交付与复绘.md','completion_audit.json','figure_contract.md','build_v89_manuscript.py','prepare_manuscript.py','prepare_release_audit.py','audit_last_render.py','narrative/3.3_定稿.md','narrative/paragraph_replacements.json','narrative/primary_sources_and_claims.md','narrative/independent_mechanism_story_review.md','narrative/final_figure_and_text_review.md','narrative/final_supplement_pages_QA.md','render/render_docx_word_wrapper.py','render/export_word_readonly.ps1']:
 files.append(V/name)
for p in files:assert p.is_file(),p
with ZipFile(zipout,'w',ZIP_DEFLATED) as z:
 for p in sorted(set(files)):z.write(p,Path(V.name)/p.relative_to(V))
with ZipFile(zipout) as z:assert z.testzip() is None
print(json.dumps({'word':str(doc),'png':str(fig)+'.png','pdf':str(fig)+'.pdf','zip':str(zipout),'source_files':len(set(files)),'zip_bytes':zipout.stat().st_size,'word_sha256':audit['docx_sha256'],'zip_sha256':hashlib.sha256(zipout.read_bytes()).hexdigest()},ensure_ascii=False))