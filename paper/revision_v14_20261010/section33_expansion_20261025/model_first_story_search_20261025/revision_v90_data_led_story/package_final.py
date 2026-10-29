from pathlib import Path
from zipfile import ZipFile,ZIP_DEFLATED
from pypdf import PdfReader
from dataclasses import asdict
import json,hashlib,shutil,importlib.util,sys
R=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010');S=R/'section33_expansion_20261025/model_first_story_search_20261025';O=S/'revision_v89_mechanism_story';V=S/'revision_v90_data_led_story';doc=R/'paper_v90_20261026_中文修订_图3数据呈现与机制推断.docx';fig=V/'figure3/fig3_data_led_story'
for folder in ['data_audit','data_audit_only','supplement_1MXE']:
 shutil.copytree(O/folder,V/folder,dirs_exist_ok=True)
# Final text is limited to observed structural context and predictions; there is no drawn mechanism or activity outcome.
text=(V/'manuscript/v90_final_text.txt').read_text(encoding='utf-8');assert '图3c' not in text
q=json.loads((V/'figure3/figure_QA.json').read_text(encoding='utf-8'));q['visual_review']='Passed: actual final PNG and Word pages10–12 reviewed; remaining26 pages exact PNG identity to previously verified V89.';q['font_embedding']='ArialMT and Arial-BoldMT embedded TrueType subsets, pdffonts verified.';q['minimum_font_in_Word164mm_pt']=7.9*164/183
(V/'figure3/figure_QA.json').write_text(json.dumps(q,ensure_ascii=False,indent=2),encoding='utf-8')
qa_path=Path(r'C:\Users\DELL\.codex\skills\academic-figure-skill\scripts\qa_validator.py');sp=importlib.util.spec_from_file_location('figqa',qa_path);mod=importlib.util.module_from_spec(sp);sys.modules['figqa']=mod;sp.loader.exec_module(mod);findings=[];source=(V/'figure3/compose_fig3_1MXE.py').read_text(encoding='utf-8')
for name in dir(mod):
 if name.startswith('check_'):
  f=getattr(mod,name)(source);findings.extend(f if isinstance(f,list) else [f])
(V/'figure3/style_QA.json').write_text(json.dumps({'findings':[asdict(f) for f in findings],'manual_resolutions':{'AP-0':'Mandatory baseline preserved verbatim; validator whitespace literal false positive.','AP-4':'Global legends are figure-coordinate artists outside axes.','CL-2':'Actual PDF mediabox183x154mm verified.'}},ensure_ascii=False,indent=2),encoding='utf-8')
p=PdfReader(str(fig)+'.pdf').pages[0];dim=[float(p.mediabox.width)*25.4/72,float(p.mediabox.height)*25.4/72];assert abs(dim[0]-183)<.01 and abs(dim[1]-154)<.01
assert len(PdfReader(next((V/'render/release_01').glob('*.pdf'))).pages)==29
with ZipFile(doc) as z:
 assert z.testzip() is None
 assert z.read('word/media/image5.png')==Path(str(fig)+'.png').read_bytes()
audit={'status':'complete','removed':'V89 c mechanism schematic and activity-prediction arrows','panels':['a actual structural context and contact enlargement','b all absolute model profiles plus repeat Full RMSF expanded y scale'],'canvas_mm':dim,'original_unique_profile_values':92,'Full_RMSF_repeated_values':25,'repeat_not_independent_data':True,'new_model_or_MD_runs':False,'no_MD_reference_in_main_figure':True,'all25_RMSF_and21_Hchi_positions':True,'mechanistic_interpretation_in_body_only':True,'pages':29,'all_pages_QA':'Actual10–12; remaining26 exact identity to verifiedV89','all6tables_and_other_images_preserved':True,'source_docx':'paper_v89_20261026_中文修订_图3调控机制假说.docx','final_docx':str(doc),'docx_sha256':hashlib.sha256(doc.read_bytes()).hexdigest(),'figure_png_sha256':q['sha256']['png']}
(V/'completion_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
readme='''# V90 图3数据呈现与机制推断

本版彻底撤去V89的c机制示意、回接虚线和预测活性箭头。主图仅保留真实结构背景和数据；机制解释留在正文，并明确属于依据结构、柔性与既有功能实验提出的假说。

figure3/fig3_data_led_story.pdf为183×154mm矢量主文件，SVG文字可编辑，PNG450dpi。相比V89高度178mm减少24mm。主图a保留1A06/1MXE真实结构及F298放大；b完整保留25位RMSF和21个有效Hχ的Full/Intra曲线。左下Full: expanded y scale只是将同一Full的全部25个RMSF值在0–2.3Å纵轴范围再画一次，以看清中央低值变化，不是新数据或独立证据。原双曲线全范围0–8.2Å也保留，无断轴、无裁值。

source_data/all_plotted_dynamics.csv含92个唯一读数，repeated_Full_RMSF_expanded.csv含25个重复放大读数。独立核查均与V89冻结数据精确相等。无新训练、推理或MD；Intra仍为同一结合态坐标上的链内模型预测，不代表游离肽。

figure3/compose_fig3_1MXE.py从脚本目录读取相对路径，依赖numpy/pandas/matplotlib/Pillow，可重导PDF/SVG/PNG。结构灰色、调控段紫色、功能位点青色；Full蓝圆实线、Intra暖赭方点虚线。最小7.9pt，在Word164mm宽下约7.08pt。

pymol/保存真实结构输入、渲染PNG、PSE、脚本、全部距离与投影。图a继承V89已验证资产，未改变坐标。1A06仅观测至316，1MXE不含激酶催化域；F298局部图突出较近5个邻居，完整6对象（含Y215）均在距离CSV/JSON。

narrative/3.3_定稿.md和figure3/caption_zh.md为本版最终正文和图注。build_v90_manuscript.py以V89为源，局部更新p71/p73/p75/p76与p72的图像，其他段落、6表和其他图片保持相同。Word仍29页；10–12页实际看图核查，其余26页与已完成视觉检查的V89逐像素相同。render/的包装脚本使用本机Word只读导出与内置render_docx.py生成内部QA页面。

data_audit、data_audit_only、supplement_1MXE为继承的模型和MD证据，用于复核主图及图S3。它们不属于新模型结果。核心W303–M316按结构锚点预先定义，其14位Full RMSF全部0.447–0.689Å；原文只据连续受约束和侧链位置差异结合文献[29,30]作机制推断，没有直接测得激活或抑制性回接。数据审计中的5个时间窗口来自同一条轨迹，不是独立重复。
'''
(V/'README_交付与复绘.md').write_text(readme,encoding='utf-8')
zipout=R/'paper_v90_20261026_图3数据呈现_源文件与PyMOL.zip';files=[]
for folder in ['figure3','pymol','data_audit','data_audit_only','supplement_1MXE','manuscript']:
 files.extend(p for p in (V/folder).rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.suffix!='.pyc')
for name in ['README_交付与复绘.md','completion_audit.json','figure_contract.md','build_v90_manuscript.py','prepare_manuscript.py','prepare_figure.py','audit_render.py','narrative/3.3_定稿.md','narrative/paragraph_replacements.json','render/render_docx_word_wrapper.py','render/export_word_readonly.ps1']:
 files.append(V/name)
with ZipFile(zipout,'w',ZIP_DEFLATED) as z:
 for p in sorted(set(files)):z.write(p,Path(V.name)/p.relative_to(V))
with ZipFile(zipout) as z:assert z.testzip() is None
print(json.dumps({'word':str(doc),'pdf':str(fig)+'.pdf','png':str(fig)+'.png','zip':str(zipout),'files':len(set(files)),'docx_sha256':audit['docx_sha256'],'zip_bytes':zipout.stat().st_size},ensure_ascii=False))