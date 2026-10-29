from pathlib import Path
from datetime import datetime,timezone,timedelta
import hashlib,json,zipfile,shutil
ROOT=Path(__file__).resolve().parent;BASE=ROOT.parent.parent
QA=ROOT/'manuscript_QA';FIG=ROOT/'figures';EXP=ROOT/'experiment';SRC=ROOT/'source'
STEM='fig7_clean_chi_scatter_and_MD_curves'
DOC=BASE/'paper_v37_20261015_中文修订_图7删去汇总表与布局精简_三线表.docx'
SOURCE=BASE/'paper_v36_20261015_中文修订_图7逐角散点与三子图布局_三线表.docx'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
readj=lambda p:json.loads(Path(p).read_text(encoding='utf8'))
write=lambda p,x:Path(p).write_text(json.dumps(x,ensure_ascii=False,indent=2),encoding='utf8')
pages=readj(QA/'page_pixel_comparison.json')
changed=[p['page'] for p in pages if not p['pixel_identical_to_reviewed_v36']]
assert len(pages)==28
for p in pages:
    p['review']='final_page_PNG_directly_visually_reviewed' if p['page'] in changed else 'pixel_identical_to_already_reviewed_final_v36'
visual=readj(QA/'final_page_visual_review.json')
assert visual['status']=='passed' and visual['changed_pages_directly_reviewed']==changed
integration=readj(QA/'integration_audit.json');numeric=readj(QA/'numerical_content_audit.json');scientific=readj(QA/'independent_scientific_review.json');experiment=readj(EXP/'audit.json');figure=readj(FIG/(STEM+'.audit.json'))
assert scientific['status']==numeric['status']==experiment['status']=='passed'
assert sha(DOC)==integration['document_sha256']==numeric['document_sha256']
assert sha(DOC) in json.dumps(scientific)
for ext in ['png','pdf','svg']:assert sha(FIG/(STEM+'.'+ext))==figure['output_sha256'][ext]
assert experiment['checks_passed']==1149
for name,digest in experiment['output_sha256'].items():assert sha(EXP/name)==digest
for name,item in experiment['original_source_provenance'].items():assert sha(item['read_path'])==item['original_sha256']
integration.update({'review':'passed_all28pages;changedpages_directlyreviewed;unchangedpages_pixelidentical;all564originalvaluesretained;151scatterpoints_exact;8associationsrecomputed;3paragraphsupdated','pages':28})
write(QA/'integration_audit.json',integration)
note=ROOT/'图7删去汇总表与布局精简说明.md'
note.write_text('v37按要求删除图7b下方整个2×4相关系数矩阵，保留χ1–χ4全部151个逐角散点及各n/ρ，统一坐标和原点重合位置保留。画布收紧为183×147mm，a/b/c布局保持，a仅标Q2/Q41，c原432条曲线数据及低自由能角区宽度保留。图7b现在仅展示联合Hχ预测与MD单角细箱熵的关联，方法、结果、图注同步删去该矩阵的宽度相关说明。四χ细箱熵ρ0.497505/0.447186/0.524704/0.127273均保留；每χ宽度相关及全部八个原系数仍保存在experiment来源表内。原83记录群体低自由能宽度验证和c的物理宽度值保留。新715条图中数据为132个全序列值+151散点+432曲线点，与v36对应数值严格相同；未新增MD、训练、推断或统计分析。原生三线表11张、数学元素36个、参考文献27条及批注保留。\n',encoding='utf8')
pointer=ROOT.parent/'completed_results_20261012/delivery_status.json';snapshot=QA/'delivery_status_v36_snapshot.json'
assert not snapshot.exists()
shutil.copy2(pointer,snapshot)
prior=readj(snapshot)
assert prior['document']==str(SOURCE) and sha(SOURCE)==prior['document_sha256']
oldarchive=Path(prior['new_figure_package']);assert sha(oldarchive)==prior['new_figure_package_sha256']
(ROOT/'README_zh.md').write_text('本包包含v37图7全矢量PDF/SVG、300dpi PNG、715条完整绘图数据、全部151角明细与原始缓存、作图和文稿集成/核验脚本。figures/plot_clean_chi_scatter_and_MD_curves.py优先读取包内source及experiment。b下汇总表已删除，原8个系数仍保存于来源统计文件。旧源包 '+str(oldarchive)+'，SHA256 '+sha(oldarchive)+'。本轮为图形精简与对应文字修订，没有新增实验或分析。\n',encoding='utf8')
files=[]
for directory,prefix in [(FIG,'figures'),(QA,'QA'),(EXP,'experiment'),(SRC,'source')]:
    files.extend((p,prefix+'/'+p.name) for p in sorted(directory.iterdir()) if p.is_file())
files.extend((p,p.name) for p in sorted(ROOT.iterdir()) if p.is_file() and p.name not in ['v37_package_manifest.json','delivery_status.json'])
assert len({name for p,name in files})==len(files)
manifest={'document_sha256':sha(DOC),'source_document_sha256':sha(SOURCE),'previous_source_archive':str(oldarchive),'previous_source_archive_sha256':sha(oldarchive),'files':{name:sha(p) for p,name in files}}
write(ROOT/'v37_package_manifest.json',manifest);files.append((ROOT/'v37_package_manifest.json','v37_package_manifest.json'))
package=BASE/'paper_v37_20261015_图7简化散点与源数据.zip'
with zipfile.ZipFile(package,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    for path,name in files:z.write(path,name)
with zipfile.ZipFile(package) as z:
    assert z.testzip() is None
    assert all(hashlib.sha256(z.read(name)).hexdigest()==digest for name,digest in manifest['files'].items())
prior.update({'status':'chinese_manuscript_v37_figure7_clean_scatter_delivered','completed_beijing':datetime.now(timezone(timedelta(hours=8))).isoformat(),'previous_document':str(SOURCE),'previous_document_sha256':sha(SOURCE),'document':str(DOC),'document_sha256':sha(DOC),'pages':28,'render_and_visual_review':f'passed_all28pages_{len(changed)}changedpages_directview_{28-len(changed)}exactpixels','content_review':'passed_abc_panel_layout;all151scatterpoints;summarymatrixremoved;original432curves','new_figure_package':str(package),'new_figure_package_sha256':sha(package),'new_experiment_note':str(note),'changes_this_round':'Removedwhole2x4summarymatrix;retainedall151scatterpoints;compactedlayout;all715retainedvaluessame;3textparagraphsupdated.','new_numerical_review':'passed_root715points_exact;all151scattervalues_verified;all8sourcecoefficients_preserved','new_MD_production_this_round':0,'new_model_training_this_round':0,'new_model_inference_this_round':0,'new_MD_analysis_this_round':0,'figure7_panel_letters':['a','b','c'],'figure7_scatter_points':151,'figure7_total_source_points':715,'figure7_summary_matrix_displayed':False,'remaining_required_work':[]})
write(ROOT/'delivery_status.json',prior);write(pointer,prior)
(ROOT.parent/'completion_status_zh.md').write_text(f'当前交付 {DOC}，SHA256 {sha(DOC)}，28页。图7b下整个汇总表删除，保留χ1–χ4全部151个散点及n/ρ，紧凑a/b/c布局，a只标Q2/Q41。715条保留数值严格等同v36。方法、结果和图注同步，原每χ宽度系数保存在源数据，c曲线及83记录主验证仍在。三线表、原数学、批注及其他图均保留，28页排版核验通过。新版源包 {package}。本轮无新增分析、MD、训练或推断。\n',encoding='utf8')
print(json.dumps({'status':'complete','document':str(DOC),'document_sha256':sha(DOC),'pages':28,'package':str(package),'package_MB':round(package.stat().st_size/1e6,2),'checks_experiment':1149,'checks_figure':figure.get('independent_check_count',figure.get('checks_passed','see_figure_audit'))},ensure_ascii=False))
