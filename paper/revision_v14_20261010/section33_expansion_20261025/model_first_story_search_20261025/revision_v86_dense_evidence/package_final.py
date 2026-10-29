from pathlib import Path
from zipfile import ZipFile,ZIP_DEFLATED
from pypdf import PdfReader
from PIL import Image
import json,hashlib
R=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010');V=R/'section33_expansion_20261025/model_first_story_search_20261025/revision_v86_dense_evidence'
fig=V/'figure3/fig3_1MXE_dense_evidence';doc=R/'paper_v86_20261026_中文修订_图3完整接触与侧链环境.docx'
qa=json.loads((V/'figure3/figure_QA.json').read_text());qa['visual_review']='passed: primary full figure and Word-scale review; independent d/e/c data audit and f scientific/visual audit'
qa['font_embedding']='ArialMT and Arial-BoldMT embedded TrueType subsets (pdffonts verified)'
qa['manual_style_resolutions']={'AP-0':'Palette baseline is copied verbatim; validator falsely requires single space in DIVERGING = whereas mandatory source uses three spaces. Exact colour values unchanged.','AP-4':'Both legends are figure-coordinate artists in dedicated white space outside data axes.','AP-6':'Each shielding bar is one observed static residue descriptor, not a group mean; no missing replicate points.'}
(V/'figure3/figure_QA.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8')
pdf=PdfReader(str(fig)+'.pdf');w=float(pdf.pages[0].mediabox.width)*25.4/72;h=float(pdf.pages[0].mediabox.height)*25.4/72
assert abs(w-183)<.02 and abs(h-225)<.02
with ZipFile(doc) as z:assert z.read('word/media/image5.png')==Path(str(fig)+'.png').read_bytes()
render=V/'render/release_final';finalpdf=list(render.glob('*.pdf'))[0];assert len(PdfReader(finalpdf).pages)==29
pixel=json.loads((V/'render/v85_v86_pixel_comparison.json').read_text());assert all(pixel[str(i)]['pixel_identical_to_audited_v85'] for i in range(1,10))
refs='''# 图3期刊图参考与具体改图依据

本版实际查阅了下列正式论文中的具体图。仅借鉴证据组织、区域标记和版式；所有新图数值及结构来自本项目，不移用参考图的实验、统计或机制结论。

1. Weber et al., 2024, Briefings in Bioinformatics 25(1), bbad504，Figure 6，Unsupervised and supervised AI on molecular dynamics simulations reveals complex characteristics of HLA-A2-peptide immunogenicity。正式论文：https://academic.oup.com/bib/article/25/1/bbad504/7560312 。图中把完整结构/锚点、接触距离、主链/侧链残基曲线串联；受体区域用浅透明带表示。用于本版a/b与c/e的衔接及c的浅区域带，不引入该文的免疫原性或时间轨迹结论。

2. Zhu et al., 2022, Nature Communications，Neural relational inference to learn long-range allosteric interactions in proteins from molecular dynamics simulations，Figure 2。正式图页：https://www.nature.com/articles/s41467-022-29331-3/figures/2 。结构—完整残基矩阵—区域分区互相对应，分区条在矩阵边界。用于d的完整25×144静态接触图和轴外叶瓣条，不推断本项目未计算的变构网络。

3. Genome Biology，2024，FoldPAthreader论文，Figure 4。正式论文：https://link.springer.com/article/10.1186/s13059-024-03291-x 。图中结构功能区配色与逐残基曲线上的细色带对应。用于c的完整序列与W303–M316细区域条，不用厚黑框覆盖数据。

4. Nguyen et al., 2025, Nature Communications，A domain-swapped CaMKII conformation facilitates linker-mediated allosteric regulation，Figure 7。正式图页：https://www.nature.com/articles/s41467-025-63249-w/figures/7 。局部接触、完整序列轨道及相邻功能定量彼此定位。此文对象为CaMKII，和本项目CaMKI不同，仅作图形组织参考。

本版落实为：d保留全部3600个残基对；e保留21个有χ角侧链，按类型排列并对齐静态遮蔽；f展示真实1A06与1MXE结构和M316的不同接触对象。静态埋藏与Hχ并非普遍单调，正文保留K304/K302等对照。c继续使用绝对RMSF与Hχ，不画主图MD参考，MD核查仍在S3。
'''
(V/'图3期刊参考与改图依据.md').write_text(refs,encoding='utf-8')
readme='''# 1MXE图3完整接触与侧链环境

最终主图为figure3/fig3_1MXE_dense_evidence.pdf（183×225 mm，Arial嵌入），同名SVG保留可编辑文本，PNG为450 dpi。图3已替换进paper_v86_20261026_中文修订_图3完整接触与侧链环境.docx。图注在figure3/caption_zh.md，正文在narrative/3.3_定稿.md。完整排版29页，经渲染核查；内部页图不作为交付图稿。

本版d为全部25×144的起始结构接触图：≤4 Å深青，4–4.5 Å浅青，>4.5 Å白；53+26+3521格，无位置抽样。e为全部21个有定义Hχ的侧链，按类型分组，以蓝圆/暖赭方点对比Full/Intra，右侧青条表示静态CaM遮蔽比例。f为1A06真实自抑制CaMKI与1MXE结合短肽结构，对比共有294–316片段，1A06缺失317–318不补造；标出M316邻近L45/S99。虚线问号为机制假说。c是完整25/21位置的绝对预测曲线，使用结构锚点定义的区域带与同坐标序列，未按预测高低筛区。

模型与MD结果沿用已完成、经审计的全组成链家族留出第1折预测；本轮没有重新训练或声称对1A06新增推理、MD或激活实验。Full与Intra是同一结合态输入的独立训练模型，不能解释为游离/结合实验状态。1MXE不含激酶催化域，f不表示模拟观察到激活或重建抑制接触。S3和其他原图保留，图3在主文中不绘MD参考。

复绘主图：用包含numpy/pandas/matplotlib/Pillow的Python运行figure3/compose_fig3_1MXE.py，脚本使用包内相对路径。当前可用解释器为E:/Project/FlexPPI/experiments/_sources/local_model_benchmark_20261023/runtime_py312/Scripts/python.exe。结构重绘运行pymol/render_1MXE_v85.py与pymol_context/render_structural_contexts.py；后者一次生成最终透明度、接触侧链棒、PSE、裁边图和landmarks，运行环境D:/Anaconda/envs/pvb_cpu/python.exe。两结构全部原子只作刚体变换，不扭曲片段；PSE验证见pymol_context/structure_QA.json。

data_audit保留完整模型/MD与结构来源，data_design保留完整距离与静态环境字段；figure3/source_data是主图的冻结直接输入。narrative/paragraph_replacements.json是Word构建的最终文本来源；manuscript/v86_text_plan.json与v86_build_audit.json记录最终内容及最小替换审计。源论文保持不变。

验证通过：主图全部3600格、21侧链、92曲线值逐项精确；所有柔性坐标从零开始；字体嵌入；PNG无越界文字；Word图文同页，无新增孤页。29页中25页与已审定v85逐像素相同，4页仅为授权图3、正文及方法变化。自动样式检查的三个模式识别误报及人工核查理由见figure3/figure_QA.json。

参考图与具体借鉴见图3期刊参考与改图依据.md。
'''
(V/'README_交付与复绘.md').write_text(readme,encoding='utf-8')
audit={'status':'complete','docx':str(doc),'docx_pages':29,'docx_sha256':hashlib.sha256(doc.read_bytes()).hexdigest(),'figure_pdf_mm':[w,h],'figure_png_sha256':hashlib.sha256(Path(str(fig)+'.png').read_bytes()).hexdigest(),'all_25x144_pairs':3600,'all_valid_sidechains':21,'all_profile_values':92,'modified_zip_parts':['word/document.xml','word/media/image5.png'],'unchanged_pages_pixel_exact':25,'changed_render_pages':[10,11,12,24],'all_pages_visual_QA':'1–9 exact pixel identity with previously visually audited release; 10 primary image review; 11–18 independent primary scientific reviewer; 19–29 independent data reviewer','new_model_runs':False,'new_MD_runs':False,'1A06_role':'static experimental structural context, no missing atoms rebuilt','style_and_data_QA':'passed with documented literal-parser false positives resolved manually'}
(V/'completion_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
zipout=R/'paper_v86_20261026_图3完整证据_源文件与PyMOL.zip'
folders=['figure3','pymol','pymol_context','data_audit','data_design','supplement_1MXE','manuscript']
files=[]
for folder in folders:
 files.extend(p for p in (V/folder).rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.suffix not in ['.pyc'])
for name in ['README_交付与复绘.md','图3期刊参考与改图依据.md','completion_audit.json','build_v86_manuscript.py','narrative/3.3_定稿.md','narrative/paragraph_replacements.json','narrative/final_scientific_visual_review.md','render/v85_v86_pixel_comparison.json','render/render_docx_word_wrapper.py','render/export_word_readonly.ps1']:
 files.append(V/name)
with ZipFile(zipout,'w',ZIP_DEFLATED) as z:
 for p in sorted(set(files)):z.write(p,Path(V.name)/p.relative_to(V))
with ZipFile(zipout) as z:assert z.testzip() is None
print(json.dumps({'docx_sha256':audit['docx_sha256'],'source_zip':str(zipout),'files':len(set(files)),'zip_bytes':zipout.stat().st_size,'zip_sha256':hashlib.sha256(zipout.read_bytes()).hexdigest()},ensure_ascii=False))