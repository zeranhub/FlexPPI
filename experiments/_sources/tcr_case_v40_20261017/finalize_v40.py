from pathlib import Path
from datetime import datetime,timezone,timedelta
import hashlib,json,zipfile,shutil
ROOT=Path(__file__).resolve().parent;BASE=ROOT.parent.parent;QA=ROOT/'manuscript_QA'
DOC=BASE/'paper_v40_20261017_中文修订_TCR识别环境变化实例_三线表.docx'
TMP=Path(r'C:\Users\DELL\AppData\Local\Temp\flexppi_completion_20261010')
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
readj=lambda p:json.loads(Path(p).read_text(encoding='utf8'))
write=lambda p,x:Path(p).write_text(json.dumps(x,ensure_ascii=False,indent=2),encoding='utf8')
audit=readj(QA/'numerical_structural_audit.json');integration=readj(QA/'integration_audit.json')
assert audit['status']=='passed' and sha(DOC)==audit['document_sha256']==integration['document_sha256']
pages=[{'page':n,'status':'passed','sha256':sha(TMP/'render_v40'/f'page-{n}.png')} for n in range(1,31)]
write(QA/'final_page_visual_review.json',{'status':'passed','document_sha256':sha(DOC),'method':'Root directly opened every final150dpi original page PNG, pages1–30; figure6 also independently inspected at450dpi.','pages':pages,'required_fixes':[],'notes':'All captions remain with their figures, native tables and equations readable; inherited page breaks reserve space for complete figures.'})
write(ROOT/'figure6/figure6_visual_QA.json',{'status':'passed','image_sha256':sha(ROOT/'figure6/fig6_tcr_case.png'),'passes':{'data_integrity':'all13 peptide positions and11valid Hchi positions; original predictions and recomputed3replicate labels checked','typography':'English labels; embedded TrueType PDF fonts;7–10pt base design; no overlap after legend revision','molecular_views':'MHC aligned, shared camera/scale/crop; true coordinates; pMHC-only intervention labels explicit','visual_verification':'450dpi standalone PNG and Word150dpi page17 directly viewed'},'limitations_in_result':'Q9 Hchi increase is visible and described as not reproduced; directional agreement refers to MD means.'})
audit['visual_review']='passed_all30pages_directly_viewed';write(QA/'numerical_structural_audit.json',audit)
integration.update({'status':'complete','render_status':'passed','pages':30,'visual_review':'passed_all30pages_directly_viewed','remaining_required_work':[]});write(QA/'integration_audit.json',integration)
note=ROOT/'修改说明_zh.md'
note.write_text('''v40新增TCR–pMHC实例，2026-10-17

在244条ITrace记录的完整家族留出预测中，筛查124对相同肽/MHC序列、不同TCR的结构组合，并另核对已报道的KRAS、HBV和A6对照。选用识别方式差异有原始结构研究支持、且主链与Q7侧链变化均有预测响应的SB27(2AK4)/SB47(4JRY)。完整配对与指标见local_pair_screen及change_screen；实例为筛查后选择，不作为未筛选全队列的平均成绩。

该对照中，13肽LPEPLPQGQLTAY与MHC输入序列完全相同，全部ITrace记录不在所用第0折模型的训练集合内。重新分析每种复合物的3条200ns轨迹，全部保存帧为1001–1003帧。原始文件未改动；逐肽位置的标签均值及有效位置方差在3e-6容差内重现。模型是论文当前联合RMSF/Hchi的2输出模型，检查点SHA256为ff2974406d5d1ee07e6dc7118e43e57c0acc9849cc8feba7f5290476fca15ec5。

全13个位点的RMSF变化Spearman相关0.945，变化MAE0.093Å，13/13方向与MD均值相同。Q7联合Hchi的MD变化1.074nats、模型变化1.099nats；SB47三条轨迹Q7熵均高于SB27全部三条。Hchi的全肽变化没有主链同等的吻合程度，Q9参考增幅未重现，图6d及正文保留这一事实。图显示全部13个RMSF位置和11个Hchi有效位置，无χ角的G8/A12仅作为无值标记。

在hit-gpu上以同一模型补做4次CPU推理：两个完整输入复算、两个移除TCR的输入。复算与缓存最大差异低于2e-5；移除TCR时肽、MHC、β2微球蛋白的原始原子坐标行完全保留。Q7预测差由1.099变为−0.018nats。此计算量化预测对TCR输入环境的依赖，并未运行游离pMHC的MD。

原文献为Liu等JBC2013，doi:10.1074/jbc.M112.447185。4JRY的PDB标题与原论文在HLA亚型命名上存在差异，但原始RCSB实体序列及实际输入的MHC序列逐位一致，正文使用HLA-B35。文献中Q7/G8在SB47晶体的暴露描述不等同于轨迹中持续无接触；复算的接触数据完整保留，但不作简单接触消失的机制结论。

新增正文3.5、英文图6、方法2.7中的实验定义及参考文献28；摘要、引言、讨论和结论同步融入实际发现。原泛素/自由能章节变为3.6，图6变为图7；原有7张活动图像的字节及10张原生三线表、36项原生数学均保持。批注内容及锚点保留。30页逐页渲染检查通过。本轮没有重新训练或新跑MD。
''',encoding='utf8')
(ROOT/'README_zh.md').write_text('''v40 TCR识别环境变化实例复现材料

最终稿位于本目录上两级，文件名paper_v40_20261017_中文修订_TCR识别环境变化实例_三线表.docx。正文新增3.5/图6；原泛素与自由能分析为3.6/图7。原稿v39的图1–5与S1及其源数据包仍适用。

figure6/：最终英文图PNG/PDF/SVG、所有13个肽位置的源表、中文图注、PyMOL会话及数据/视觉核验。plot_case.py为可直接运行的绘图脚本，使用既有plot_runtime；render_structures.py使用现有headless PyMOL。
sb27_sb47/：六条原始轨迹的逐位标签与Q7角度/占据/接触摘要、原始文件SHA256、模型输入和完整/移除TCR推理。analyze_sb27_sb47.py从E:/data/ITrace原始文件重新计算；原始XTC依ITrace许可获取。
context_inference.py：hit-gpu上的同一联合模型输入对照；运行目录/root/flexppi_completion_20261010/tcr_case_v40_20261017，环境/root/flexppi_completion_20261010/env/bin/python，原源代码/root/flexppi_revision_20261002/corrected_training_v1/src。模型检查点和原始预测已由既有完成实验备份及本包SHA记录追溯。
family_predictions/：244条ITrace当前家族留出Full/Intra预测及清单。正文实例仅用Full，未使用历史3输出模型的预测。
local_pair_screen/、change_screen/、change_screen_focus/：完整序列匹配、配对筛查及文献候选对比。biology_sources/保留所用原始文献/数据库信息；unrelated_PMC3668746_excluded.md是错误地址返回的无关论文，不用于分析，也不纳入交付压缩包。
manuscript_QA/：原稿/新稿SHA、逐数值/逐结构核验、30页页面检查记录与全文提取。可用build_v40.py重建文稿，export_v40_word.ps1与render_v40.py以隐藏Word导出后交由官方render_docx.py渲染，audit_v40.py复核数值与文档结构。
''',encoding='utf8')
mapping=readj(ROOT.parent/'reader_revision_v39_20261016/CURRENT_FIGURE_MAP.json')
mapping['7']=mapping.pop('6');mapping['7']['source_revision']='v39'
mapping['6']={'source_folder':'figure6','png':'fig6_tcr_case.png','pdf':'fig6_tcr_case.pdf','svg':'fig6_tcr_case.svg','new_in_revision':'v40'}
for key,item in mapping.items():
    item['absolute_source_folder']=str((ROOT if key=='6' else ROOT.parent/'reader_revision_v39_20261016')/item['source_folder'])
write(ROOT/'CURRENT_FIGURE_MAP.json',mapping)
pointer=ROOT.parent/'completed_results_20261012/delivery_status.json'
snapshot=QA/'delivery_status_v39_snapshot.json'
if not snapshot.exists():shutil.copy2(pointer,snapshot)
prior=readj(snapshot);assert prior['document_sha256']==integration['source_sha256']
files=sorted(p for p in ROOT.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.name not in ('v40_package_manifest.json','delivery_status.json','unrelated_PMC3668746_excluded.md'))
manifest={'document':str(DOC),'document_sha256':sha(DOC),'previous_document':integration['source'],'files':{p.relative_to(ROOT).as_posix():sha(p) for p in files}}
write(ROOT/'v40_package_manifest.json',manifest)
package=BASE/'paper_v40_20261017_TCR实例图与实验源数据.zip'
with zipfile.ZipFile(package,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    for p in files:z.write(p,p.relative_to(ROOT).as_posix())
    z.write(ROOT/'v40_package_manifest.json','v40_package_manifest.json')
with zipfile.ZipFile(package) as z:
    assert z.testzip() is None
    for name,digest in manifest['files'].items():assert hashlib.sha256(z.read(name)).hexdigest()==digest
prior.update({'status':'chinese_manuscript_v40_TCR_case_delivered','completed_beijing':datetime.now(timezone(timedelta(hours=8))).isoformat(),'previous_document':integration['source'],'previous_document_sha256':integration['source_sha256'],'document':str(DOC),'document_sha256':sha(DOC),'pages':30,'main_figures':7,'supplementary_figures':1,'figures':8,'references':28,'render_and_visual_review':'passed_all30pages_directly_reviewed','content_review':'passed_source_numerical_and_scope_audit','new_figure_package':str(package),'new_figure_package_sha256':sha(package),'new_experiment_note':str(note),'changes_this_round':'Added same-pMHC SB27/SB47 example in main3.5 and EnglishFigure6;3x200ns reference recomputation each;complex input perturbation with same model;ubiquitin/freeenergy now3.6/Figure7','ubiquitin_figure_number':7,'TCR_case_figure_number':6,'TCR_case_pdbs':['2ak4','4jry'],'TCR_case_delta_RMSF_rho':.945054945054945,'TCR_case_delta_RMSF_MAE_A':.09310372059161846,'TCR_case_Q7_delta_Hchi_MD':1.074399700209652,'TCR_case_Q7_delta_Hchi_model_full':audit['Q7_predicted_delta_full'],'TCR_case_Q7_delta_Hchi_model_TCR_removed':audit['Q7_predicted_delta_TCR_removed'],'new_MD_production_this_round':0,'new_model_training_this_round':0,'new_model_inference_this_round':4,'new_MD_analysis_this_round':6,'remaining_required_work':[]})
write(ROOT/'delivery_status.json',prior);write(pointer,prior)
(ROOT.parent/'completion_status_zh.md').write_text(f'当前交付 {DOC}，SHA256 {sha(DOC)}，30页。正文3.5与英文图6加入同一肽/MHC序列下的SB27/SB47动态变化实例；六条200ns轨迹标签复算，全部13个主链及11个侧链有效位置保留；4次同模型输入推理。原泛素与自由能章节为3.6/图7。10张原生三线表、36项数学、原批注和原有7幅活动图像均保留。30页逐页检查通过。源包 {package}。本轮无新增训练或MD生产。\n',encoding='utf8')
print(json.dumps({'status':'complete','document':str(DOC),'document_sha256':sha(DOC),'pages':30,'package':str(package),'package_MB':round(package.stat().st_size/1e6,2)},ensure_ascii=False))
