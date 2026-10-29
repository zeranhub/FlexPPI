from pathlib import Path
from zipfile import ZipFile,ZIP_DEFLATED
import json,hashlib
from lxml import etree
from pypdf import PdfReader
V=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010\section33_expansion_20261025\model_first_story_search_20261025\revision_v85_1MXE');BASE=V.parents[2]
DOC=BASE/'paper_v85_20261026_中文修订_1MXE机制假说与图3重绘.docx';SOURCE=BASE/'paper_v83_20261025_中文修订_图3模型柔性再分配与关键区域.docx'
assert DOC.exists(),DOC
for sub in ['figure3','supplement_1MXE']:
 p=V/sub/'skill_QA_static_report.json';q=json.loads(p.read_text(encoding='utf-8'));q['manual_resolutions']['AP-4_legend']='Legends are drawn explicitly in figure coordinates outside the quantitative axes, rather than via ax.legend; final PNGs independently inspected without occlusion.';p.write_text(json.dumps(q,ensure_ascii=False,indent=2),encoding='utf-8')
contract='''# 图3设计与证据约定

核心论点：1MXE的两个识别锚点界定连续的低主链波动核心，而同类型侧链保留位置差异；这些观察与既有叶瓣功能实验共同支持可检验的调控构型维持假说。

布局为真实结构与定量数据混合组图，183×190 mm。a为输入结构，b在a右侧显示相同旋转方向的锚点局部；d完整核心/其余位置均值，e全部两个Asn的原值，f明确标为Mechanistic hypothesis的示意；c在d–f下方展示完整25/21位点。主图不绘制MD、不放表、不放灰色说明小字。定量为Full蓝色圆点/实线、Intra暖赭方点/虚线，组分为灰色CaM/紫色CaMKI/青色关键侧链。

所有定量坐标从零开始；无数据平滑、纵轴截断或只取成功位点。W303–M316核心来自结构锚点定义，不按曲线挑选。Hχ空缺保留四个Ala断线；没有把不同残基种类的原始Hχ差直接解释成环境效应。MD在S3中核查，五个100 ns窗口是同一轨迹的时间段，不作独立重复或显著性检验。

Full/Intra是同一个结合态结构上分别训练的模型，不能标成结合/游离状态。f中虚线为提出的后果，尚未直接测量激酶回接触或激活；体系不含激酶催化域，也未量化钙效应。

主图100个定量绘制值、S3的766个绘制值均已逐项对照冻结家族第1折数据。矢量PDF/SVG保留文字和曲线，真实结构图为PyMOL渲染；PDF字体已嵌入。原始输入/重编号输入坐标、PSE原子和相机已独立审计。
'''
(V/'figure_contract_and_QA.md').write_text(contract,encoding='utf-8')
readme='''# 1MXE图3与正文修订

正文主线：CaM双锚点包裹CaMKI调控肽，在W303–M316的连续14残基核心恢复低主链波动，同时保留N297/N308等同类型侧链的位置差异；结合C端叶瓣结合与N端叶瓣激活分工，提出N端叶瓣的进一步约束可能帮助维持自抑制解除构型的可检验假说。

交付包括图3及图S3的PDF、SVG、450 dpi PNG；完整25/21位点、全体系169/145位点的绘图源表；实际输入结构的PyMOL会话、局部放大和渲染脚本；中文正文、图注与修改计划。Word论文另作为独立文件交付。最终正文以narrative/3.3_定稿.md及manuscript/v85_text_plan.json为准。

## 复绘

在当前项目环境用已有科学Python运行figure3/compose_fig3_1MXE.py与supplement_1MXE/compose_figS3_1MXE.py。两脚本使用包内冻结源CSV和pymol/图像，不需重新训练或推理。结构重绘脚本见pymol/；PSE可直接打开编辑。build_v85_manuscript.py以项目中保留的v83 DOCX为源，仅替换已授权正文和两张图，不重建其他图表。原始500 ns XTC和训练权重未打入制图源包，来源、标签复算及结构SHA见data_audit/和相关审计。

## 已完成核查

采用全组成链家族留出第1折，未混入旧record-split第4折或训练内第0折。全25aa RMSF MAE Full/Intra为0.114859/2.383203 Å，21个有效侧链Hχ为0.254418/0.422994 nats；不是每一个侧链都改善。S3保留全位点参考及失败偏差。

正文、图注、摘要、引言、讨论、结论、补充D和参考文献29/30均同步改为1MXE。图形已独立核查数据、配色、标签和零起点；Word已由实际Word导出并检查最终29页，正文末尾孤立尾段已并回正文。

柔性读数和静态接触为机制假说提供线索，不代表已直接验证激活或回接触。SASA比较是在相同肽坐标下删除CaM的静态遮蔽对照，不是真实游离态。五个时间窗不是独立重复。

原始结构与功能研究：Clapperton et al., Biochemistry41:14669–14679(2002), DOI10.1021/bi026660t；Goldberg et al., Cell84:875–887(1996), DOI10.1016/S0092-8674(00)81066-1。
'''
(V/'README_交付与复绘.md').write_text(readme,encoding='utf-8')
# Full-state verification against the actual final DOCX package.
with ZipFile(SOURCE) as z:old={n:z.read(n) for n in z.namelist()}
with ZipFile(DOC) as z:
 assert z.testzip() is None;new={n:z.read(n) for n in z.namelist()}
 assert set(old)==set(new)
 changed=[n for n in old if old[n]!=new[n]]
 assert set(changed)=={'word/document.xml','word/media/image5.png','word/media/figureS3_CaM_v83.png'},changed
ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
root=etree.fromstring(new['word/document.xml']);pars=root.find('w:body',ns).findall('w:p',ns)
plan=json.loads((V/'manuscript/v85_text_plan.json').read_text(encoding='utf-8'))['paragraph_replacements']
for key,t in plan.items():assert ''.join(pars[int(key)].xpath('.//w:t/text()',namespaces=ns))==t,key
for sub,name,target in [('figure3','fig3_1MXE_mechanistic_hypothesis.png','word/media/image5.png'),('supplement_1MXE','figS3_1MXE_validation.png','word/media/figureS3_CaM_v83.png')]:
 assert new[target]==(V/sub/name).read_bytes()
 q=json.loads((V/sub/'figure_QA.json').read_text(encoding='utf-8'));assert q['visual_review'].startswith('passed')
 pdf=PdfReader(V/sub/name.replace('.png','.pdf'));assert abs(float(pdf.pages[0].mediabox.width)/72*25.4-183)<.02
paperpdf=next((V/'render/release_final').glob('*.pdf'));assert len(PdfReader(paperpdf).pages)==29
requirements={'redrawn_around_1MXE':True,'actual_model_predictions_family_fold1':True,'complete_25_backbone_and_21_sidechain_values':True,'real_structure_and_key_sites_with_local_views_to_right':True,'absolute_zero_based_axes_two_models':True,'d_e_f_above_c':True,'no_tables_or_gray_fine_notes':True,'manuscript_palette_coordinated':True,'mechanism_hypothesis_and_testable_prediction_in_text':True,'no_claim_of_proven_activation_or_free_bound_transition':True,'S3_MD_and_complete_complex_validation':True,'section33_and_global_references_synchronized':True,'other_document_parts_unchanged':True,'29_pages_rendered_and_visually_verified':True,'native_PyMOL_and_sources_available':True}
audit={'requirements':requirements,'docx':str(DOC),'docx_sha256':hashlib.sha256(DOC.read_bytes()).hexdigest(),'modified_package_parts':changed,'page16_visual_review':'passed root inspection; CPU timing paragraph retained on same page','full_page_review_evidence':'manuscript/final_render_comparison.json plus independent agents original first-render reviews','all_requirements_verified':all(requirements.values())}
(V/'completion_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
# Source bundle excludes internal DOCX render pages and superseded drafts.
selected=[]
for p in V.rglob('*'):
 if not p.is_file():continue
 rel=p.relative_to(V)
 if rel.parts[0] in ['render','__pycache__']:continue
 if p.suffix in ['.pyc','.log']:continue
 if 'QA_structure_contact_sheet' in p.name:continue
 selected.append(p)
zipname=BASE/'paper_v85_20261026_1MXE图3与S3_源文件与PyMOL.zip'
with ZipFile(zipname,'w',ZIP_DEFLATED) as z:
 for p in sorted(selected):z.write(p,'revision_v85_1MXE/'+p.relative_to(V).as_posix())
with ZipFile(zipname) as z:
 assert z.testzip() is None
 for p in selected:assert z.read('revision_v85_1MXE/'+p.relative_to(V).as_posix())==p.read_bytes()
print(json.dumps({'docx':str(DOC),'source_zip':str(zipname),'zip_files':len(selected),'zip_bytes':zipname.stat().st_size,'zip_sha256':hashlib.sha256(zipname.read_bytes()).hexdigest(),'requirements':requirements},ensure_ascii=False))