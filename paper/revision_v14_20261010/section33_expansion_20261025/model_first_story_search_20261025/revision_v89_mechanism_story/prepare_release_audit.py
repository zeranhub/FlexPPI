from pathlib import Path
import importlib.util,sys,json,hashlib,shutil
from dataclasses import asdict
from pypdf import PdfReader
R=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010');S=R/'section33_expansion_20261025/model_first_story_search_20261025';V=S/'revision_v89_mechanism_story';O=S/'revision_v88_two_panel_story'
# Read-only style checks against the complete composer, then document manual resolutions.
f=Path(r'C:\Users\DELL\.codex\skills\academic-figure-skill\scripts\qa_validator.py');sp=importlib.util.spec_from_file_location('figure_qa',f);mod=importlib.util.module_from_spec(sp);sys.modules['figure_qa']=mod;sp.loader.exec_module(mod)
source=(V/'figure3/compose_fig3_1MXE.py').read_text(encoding='utf-8');findings=[]
for name in dir(mod):
 if name.startswith('check_'):
  result=getattr(mod,name)(source)
  findings.extend(result if isinstance(result,list) else [result])
qa={'automatic_findings':[asdict(f) for f in findings],'manual_resolutions':{'AP-0':'The verbatim mandatory palette has DIVERGING followed by three spaces. Literal validator demands one. All baseline values preserved.','AP-4':'Full/Intra legend and structural component key are figure-coordinate artists in dedicated white space, not inside data axes.','CL-2':'Actual PDF mediabox independently checked at183x178mm.'}}
(V/'figure3/style_validator_with_manual_resolutions.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8')
# Frozen source files remain available for a self-contained figure release.
for folder in ['data_audit','supplement_1MXE','source_analysis']:
 shutil.copytree(O/folder,V/folder,dirs_exist_ok=True)
# First nine pages have no edits; compare exact pixels rather than assuming layout identity.
identical=[]
for i in range(1,10):
 a=V/'render/release_01'/f'page-{i}.png';b=O/'render/release_01'/f'page-{i}.png'
 assert a.read_bytes()==b.read_bytes(),i
 identical.append(i)
figpdf=PdfReader(V/'figure3/fig3_mechanism_story.pdf');p=figpdf.pages[0];dim=[float(p.mediabox.width)*25.4/72,float(p.mediabox.height)*25.4/72]
assert abs(dim[0]-183)<.02 and abs(dim[1]-178)<.02
pdf=next((V/'render/release_01').glob('*.pdf'));reader=PdfReader(pdf);assert len(reader.pages)==29
text='\n'.join(p.extract_text() for p in reader.pages)
assert '3.3 主链与侧链柔性的互补解析' in text
assert all(w not in text for w in ['图3d','图3e','图3f'])
(V/'manuscript/render_audit.json').write_text(json.dumps({'pages':len(reader.pages),'unchanged_pages_exact_png_identity':identical,'figure_pdf_dimensions_mm':dim,'modified_text_indices':[70,71,74,75,76,73,163],'main_figure_page':11,'minimum_font_at_Word164mm_pt':7.9*164/183,'all_six_tables_and_other_images_unchanged':True,'new_model_or_MD_runs':False},ensure_ascii=False,indent=2),encoding='utf-8')
readme='''# 图3调控机制假说版 V89

新版将结构、预测与功能假说连接为：调控片段可形成抑制接触 → CaM双锚点捕获 → 约束可能减少抑制接触重建。图3c是假说示意，未生成完整激酶活化结构或新的模型/MD轨迹。

## 图与正文

figure3/fig3_mechanism_story.pdf：183 × 178 mm矢量主文件，字体已嵌入；SVG保留可编辑文字；PNG为450 dpi。源数据含全部25个RMSF和21个有效Hχ位置，Full/Intra共92点，与冻结家族留出第1折结果精确一致。主图不展示MD参考，原核查在supplement_1MXE中保留。

figure3/compose_fig3_1MXE.py：以脚本所在目录为基准读取相对路径，依赖numpy/pandas/matplotlib/Pillow。修改排版后重运行可生成PDF/SVG/PNG。字色语义为结构灰色、调控段紫色、功能位点青色；Full蓝圆实线、Intra暖赭方点虚线。图内没有表格或灰色解释性小字。Word以164 mm宽插入，最小文字约7.08 pt。

narrative/3.3_定稿.md与figure3/caption_zh.md为最终正文和图注；paragraph_replacements.json列出对V88的局部文本修改。build_v89_manuscript.py以V88为输入，替换图3及授权段落，保留全部6表和其他图片。

## 结构与机制

pymol_mechanism/保存真实1A06与1MXE原子输入、透明结构图、可编辑PSE、渲染脚本、投影landmark和接触距离。全景使用相同相机旋转，按版面分别缩放；F298接触区在自抑制全景右侧局部放大。加载坐标未发生内部形变。1A06只观测至316，不补建317–318；1MXE仅含CaM与294–318调控肽，不含激酶催化域。

F298在1A06中邻近G214/I210/F104/Y184/P216/Y215，K296不能画成同样的远端抑制锁扣。1MXE的K296/F298也接触CaM，图3c保留这条关系。图3b中的Intra是链内编码模型对结合态坐标的预测，不是游离肽，不支持残余螺旋或先C后N的结合时间顺序。图3c形状表示功能关系，不能用于读取距离或角度。

文献支持：调控肽可抑制CaMKI293；W303/C叶瓣和M316/N叶瓣为双锚点；C叶瓣结合紧，N叶瓣对激活必需。模型提出的具体联系是：局部约束可能减少抑制接触重建。在CaM结合占据率相近时，N叶瓣–M316扰动若同时降低局部约束与激酶活化，可作为未来检验；该实验本次未做。

## 原始文献

Clapperton等，Biochemistry 2002，doi:10.1021/bi026660t；官方SI doi:10.1021/bi026660t.s001，https://acs.figshare.com/articles/journal_contribution/3586317 。已核读摘要与全部6页SI，图S3为调控肽抑制CaMKI293实验。Goldberg等，Cell 1996，doi:10.1016/S0092-8674(00)81066-1 。官方结构：https://www.rcsb.org/structure/1A06 和 https://www.rcsb.org/structure/1MXE 。主文全文未取得的实验细节未补写。

## 数据与检查

data_audit为继承的原始模型/MD证据；data_audit_only包含新机制审计、所有位点和时间窗口复算。5个100 ns窗口来自同一条轨迹，并非独立重复。source_analysis保留先前静态环境与全接触矩阵以供复核；这些不是恢复到主图的独立面板。完整Word已渲染29页检查，内部渲染PDF和页面PNG不列为投稿交付。
'''
(V/'README_交付与复绘.md').write_text(readme,encoding='utf-8')
print(json.dumps({'pages':29,'first_nine_pixel_identity':identical,'PDF_mm':dim,'automatic_categories':{c:sum(f.category==c for f in findings) for c in ['PASS','WARN','FAIL']}},ensure_ascii=False))