from pathlib import Path
from datetime import datetime,timezone,timedelta
import hashlib,json,zipfile,shutil
ROOT=Path(__file__).resolve().parent;BASE=ROOT.parent.parent;QA=ROOT/'manuscript_QA'
DOC=BASE/'paper_v39_20261017_中文修订_读者表述与英文图表_三线表.docx'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
readj=lambda p:json.loads(Path(p).read_text(encoding='utf8'))
write=lambda p,x:Path(p).write_text(json.dumps(x,ensure_ascii=False,indent=2),encoding='utf8')
numeric=readj(QA/'numerical_structural_audit.json');integration=readj(QA/'integration_audit.json')
assert numeric['status']=='passed' and sha(DOC)==numeric['document_sha256']==integration['output_sha256']
reviewed=[]
for f in ['visual_pages_01_08.json','visual_pages_09_17.json','visual_pages_18_27.json']:
    v=readj(QA/f);assert v['status'] in ('pass','passed','PASS')
    assert sha(DOC) in json.dumps(v)
    for page in v['pages']:
        n=page['page'];assert page.get('status',page.get('verdict')) in ('pass','passed','PASS');reviewed.append(n)
        digest=page.get('sha256',page.get('image_sha256',page.get('png_sha256')))
        if digest:assert sha(Path(r'C:\Users\DELL\AppData\Local\Temp\flexppi_completion_20261010\render_v39')/f'page-{n}.png')==digest
assert sorted(reviewed)==list(range(1,28))
assert readj(QA/'independent_content_review.json')['status']=='pass'
integration.update({'render_status':'passed','pages':27,'all27pages_directly_reviewed':True,'independent_content_review':'pass'})
write(QA/'integration_audit.json',integration)
write(QA/'final_page_visual_review.json',{'status':'passed','document_sha256':sha(DOC),'pages':27,'all_pages_directly_viewed':reviewed,'required_fixes':[]})
note=ROOT/'修改说明_zh.md'
note.write_text('''本轮修订（v39，2026-10-17）

表1改为两级指标表头，方法名称左对齐，基础对照、已有RMSF方法与FlexPPI链内/完整模型按比较关系排列，统一短名称。A为原五折3063条记录，B合并原表2全部5方法/4指标的ThermoPCD结果。所有52项原始结果保留；不将不同训练语料的方法按一个总名次混排，不把本文非最优数值加粗。

取消原3.4外部验证独立节、图5和表2，3.1保留三句话并引用表1B。外部Hχ相对链内改善9.4%，相对静态特征改善1.5%，RMSF幅值与两对照接近；摘要不再单独突出外部9.4%。训练排除633条、保留2430条、32条属于一个关联家族及单模型的事实集中在方法中说明。完整逐记录结果仍在原实验源数据中，没有把被删图表搬入补充材料。

正文去除“冻结”“发布权重”“发布模型”等重复工程术语，已有方法采用作者提供模型的方案集中说明一次。保留χ角分箱、训练预算、固定阈值和单位换算等理解实验所需的条件。指定a子图第二行说明已删除；原始ρ=.979旁的正文引用与图中内容一致，类型内ρ=.868留在正文解释。

全部活动图像的可见文字均为英文，图注仍为中文。图1、原图6及S1原本已为英文，保留其图像；图2/3/4及原图7改为英文。原图6改为图5，原图7改为图6；原3.5改为3.4，原3.6及子节改为3.5及相应子节。源目录保留原图号以便追溯，当前图号与文件映射见CURRENT_FIGURE_MAP.json。

全文27页逐页检查通过；10张表均为Word原生可编辑三线表；27条参考文献、36项原生数学与原批注/锚点保留。原表1/2和9张补表数值、14份图2/3/4的数值CSV及原图7全部715个绘图数据均核对未变。本轮没有新增训练、推理或MD。
''',encoding='utf8')
mapping={
 '1':{'source_folder':'english_figures/fig1','png':'fig1_english.png','pdf':'fig1_english.pdf','svg':'fig1_english.svg'},
 '2':{'source_folder':'english_figures/fig2','png':'fig2_accuracy_and_chain_context.png','pdf':'fig2_accuracy_and_chain_context.pdf','svg':'fig2_accuracy_and_chain_context.svg'},
 '3':{'source_folder':'english_figures/fig3','png':'fig3_amplitude_and_residue_types.png','pdf':'fig3_amplitude_and_residue_types.pdf','svg':'fig3_amplitude_and_residue_types.svg'},
 '4':{'source_folder':'english_figures/fig4','png':'fig4_joint_complementarity.png','pdf':'fig4_joint_complementarity.pdf','svg':'fig4_joint_complementarity.svg'},
 '5':{'source_folder':'english_figures/fig6','png':'fig6_english.png','pdf':'fig6_english.pdf','svg':'fig6_english.svg','previous_figure_number':6},
 '6':{'source_folder':'english_figures/fig7','png':'fig7_english.png','pdf':'fig7_english.pdf','svg':'fig7_english.svg','previous_figure_number':7},
 'S1':{'source_folder':'english_figures/figS1','png':'figS1_observed_label_variation.png','pdf':'figS1_observed_label_variation.pdf','svg':'figS1_observed_label_variation.svg'}}
for item in mapping.values():
    for ext in ('png','pdf','svg'):assert (ROOT/item['source_folder']/item[ext]).exists()
write(ROOT/'CURRENT_FIGURE_MAP.json',mapping)
(ROOT/'README_zh.md').write_text('本包为v39英文图、绘图源数据及修订核验材料。CURRENT_FIGURE_MAP.json按新稿图号1–6、S1列出文件路径。english_figures子目录沿用原图号：fig6为当前图5（KIX/泛素结构实例），fig7为当前图6（全残基/逐角散点/自由能曲线）。其他图号不变。修改说明及表1参考论文链接在根目录；文稿与版面核验材料在manuscript_QA。复现仍使用原训练/MD源材料，本包没有新增实验。\n',encoding='utf8')
pointer=ROOT.parent/'completed_results_20261012/delivery_status.json'
snapshot=QA/'delivery_status_v38_snapshot.json'
if not snapshot.exists():shutil.copy2(pointer,snapshot)
prior=readj(snapshot)
assert prior['document_sha256']==integration['source_sha256']
files=[p for p in ROOT.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.name not in ('v39_package_manifest.json','delivery_status.json')]
files=sorted(files)
manifest={'document':str(DOC),'document_sha256':sha(DOC),'previous_document':integration['source'],'files':{p.relative_to(ROOT).as_posix():sha(p) for p in files}}
write(ROOT/'v39_package_manifest.json',manifest)
package=BASE/'paper_v39_20261017_英文图与修订源文件.zip'
with zipfile.ZipFile(package,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    for p in files:z.write(p,p.relative_to(ROOT).as_posix())
    z.write(ROOT/'v39_package_manifest.json','v39_package_manifest.json')
with zipfile.ZipFile(package) as z:
    assert z.testzip() is None
    for name,digest in manifest['files'].items():assert hashlib.sha256(z.read(name)).hexdigest()==digest
for key in list(prior):
    if key.startswith('figure7_'):prior.pop(key)
prior.update({'status':'chinese_manuscript_v39_reader_revision_and_English_figures_delivered','completed_beijing':datetime.now(timezone(timedelta(hours=8))).isoformat(),'previous_document':integration['source'],'previous_document_sha256':integration['source_sha256'],'document':str(DOC),'document_sha256':sha(DOC),'pages':27,'main_figures':6,'supplementary_figures':1,'figures':7,'main_tables':1,'supplementary_tables':9,'three_line_tables':10,'render_and_visual_review':'passed_all27pages_directly_reviewed','content_review':'passed_independent_all_data_retained_and_cross_references_updated','scope':'Chinese manuscript with English visible figure labels; no full English text translation or external submission.','new_figure_package':str(package),'new_figure_package_sha256':sha(package),'new_experiment_note':str(note),'changes_this_round':'RebuiltTable1A_B;integratedThermoPCDinto3.1;removedold3.4Figure5Table2;Englishallfigures;removedredundantfrozen_release_jargon;renumberedfiguresandsections;allnumericdatapreserved','external_standalone_section_removed':True,'external_results_table':'1B','ubiquitin_figure_number':6,'ubiquitin_figure_panel_letters':['a','b','c'],'ubiquitin_figure_source_points':715,'ubiquitin_figure_scatter_points':151,'ubiquitin_a_type_mean_label_displayed':False,'ubiquitin_a_legend_upper_right':True,'new_MD_production_this_round':0,'new_model_training_this_round':0,'new_model_inference_this_round':0,'new_MD_analysis_this_round':0,'remaining_required_work':[]})
write(ROOT/'delivery_status.json',prior);write(pointer,prior)
(ROOT.parent/'completion_status_zh.md').write_text(f'当前交付 {DOC}，SHA256 {sha(DOC)}，27页。表1重新组织为A五折与B外部结果，全部原表数值保留；取消旧3.4、图5和表2，外部结果融入3.1三句话及表1B。全文清理冻结/发布权重/发布模型等重复用语，图内可见文字统一英文，指定a子图第二行已删。旧图6/7改为图5/6，章节与交叉引用同步。10张原生三线表、36项数学、27条文献及原批注保留；27页逐页检查及独立内容审查通过。源包 {package}。本轮无新增训练、推理或MD。\n',encoding='utf8')
print(json.dumps({'status':'complete','document':str(DOC),'document_sha256':sha(DOC),'pages':27,'package':str(package),'package_MB':round(package.stat().st_size/1e6,2)},ensure_ascii=False))
