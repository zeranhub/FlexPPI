from pathlib import Path
from datetime import datetime,timezone,timedelta
import hashlib,json,zipfile,shutil
ROOT=Path(__file__).resolve().parent;BASE=ROOT.parent.parent
QA=ROOT/'manuscript_QA';FIG=ROOT/'figures';EXP=ROOT/'experiment';SRC=ROOT/'source'
STEM='fig7_rich_chi_scatter_and_MD_curves'
DOC=BASE/'paper_v36_20261015_中文修订_图7逐角散点与三子图布局_三线表.docx'
SOURCE=BASE/'paper_v35_20261015_中文修订_图7全残基逐χ占据关联_三线表.docx'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
readj=lambda p:json.loads(Path(p).read_text(encoding='utf8'))
write=lambda p,x:Path(p).write_text(json.dumps(x,ensure_ascii=False,indent=2),encoding='utf8')
pages=readj(QA/'page_pixel_comparison.json')
changed=[p['page'] for p in pages if not p['pixel_identical_to_reviewed_v35']]
assert len(pages)==28
for p in pages:
    p['review']='final_page_PNG_directly_visually_reviewed' if p['page'] in changed else 'pixel_identical_to_already_reviewed_final_v35'
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
note=ROOT/'图7逐角散点与三子图布局说明.md'
note.write_text('v36按要求将图7上排拆为a、b，下排三幅自由能轴统一为c。a保留全66有效残基的联合Hχ预测/参考及原统计，仅标Q2/Q41。b由四幅χ1–χ4同类型中心化散点及2×4相关系数概览构成，完整显示151个有效角，n=66/52/23/10；散点横轴为联合Hχ预测相对同类型均值的偏差，纵轴为MD单角细箱占据熵相对同类型均值的偏差，使用共同线性范围，重叠原点不抖动不删点。八个原相关系数全部保留，χ4弱/负关联不变。c的全部432原MD曲线点、角区总宽度与坐标范围保留；Q41/Q2仍为泛素六GLN内冻结预测低/高端点。全文仅更新方法、结果与图注三段对应关系，没有新增实验、MD生产、推断或训练。原83记录主验证数值、11张三线表、36个数学元素、27条参考文献及原批注均保留。新723条绘图数据为132全序列值、151散点、8摘要及432曲线值；所有数据来源与v35相同。\n',encoding='utf8')
pointer=ROOT.parent/'completed_results_20261012/delivery_status.json';snapshot=QA/'delivery_status_v35_snapshot.json'
assert not snapshot.exists()
shutil.copy2(pointer,snapshot)
prior=readj(snapshot)
assert prior['document']==str(SOURCE) and sha(SOURCE)==prior['document_sha256']
oldarchive=Path(prior['new_figure_package']);assert sha(oldarchive)==prior['new_figure_package_sha256']
(ROOT/'README_zh.md').write_text('本包包含v36图7PNG、全矢量PDF/SVG、723条完整绘图数据、151个有效角与264个角槽明细、四角双指标关联、原角度缓存、作图/分析/文稿集成脚本及核验。figures/plot_rich_chi_scatter_and_MD_curves.py使用包内source与experiment可重绘。a/b/c的统计定义与v35相同，新增展示已有151条逐角散点，没有新训练或MD。旧版源包 '+str(oldarchive)+'，SHA256 '+sha(oldarchive)+'；可追溯已有角度复算与83记录群体验证。\n',encoding='utf8')
files=[]
for directory,prefix in [(FIG,'figures'),(QA,'QA'),(EXP,'experiment'),(SRC,'source')]:
    files.extend((p,prefix+'/'+p.name) for p in sorted(directory.iterdir()) if p.is_file())
files.extend((p,p.name) for p in sorted(ROOT.iterdir()) if p.is_file() and p.name not in ['v36_package_manifest.json','delivery_status.json'])
assert len({name for p,name in files})==len(files)
manifest={'document_sha256':sha(DOC),'source_document_sha256':sha(SOURCE),'previous_source_archive':str(oldarchive),'previous_source_archive_sha256':sha(oldarchive),'files':{name:sha(p) for p,name in files}}
write(ROOT/'v36_package_manifest.json',manifest);files.append((ROOT/'v36_package_manifest.json','v36_package_manifest.json'))
package=BASE/'paper_v36_20261015_图7逐角散点与源数据.zip'
with zipfile.ZipFile(package,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    for path,name in files:z.write(path,name)
with zipfile.ZipFile(package) as z:
    assert z.testzip() is None
    assert all(hashlib.sha256(z.read(name)).hexdigest()==digest for name,digest in manifest['files'].items())
prior.update({'status':'chinese_manuscript_v36_figure7_rich_chi_scatter_delivered','completed_beijing':datetime.now(timezone(timedelta(hours=8))).isoformat(),'previous_document':str(SOURCE),'previous_document_sha256':sha(SOURCE),'document':str(DOC),'document_sha256':sha(DOC),'pages':28,'render_and_visual_review':f'passed_all28pages_{len(changed)}changedpages_directview_{28-len(changed)}exactpixels','content_review':'passed_abc_panel_layout;all151scatterpoints;all8coefficients;original432curves','new_figure_package':str(package),'new_figure_package_sha256':sha(package),'new_experiment_note':str(note),'changes_this_round':'Splitupperrowa_andb;4chi_scatterplots151points_plus8coefficientmatrix;3MDcurveaxes_sharepanelc;3textparagraphscrossrefsupdated.','new_numerical_review':'passed_root723points_exact;all151scattervalues_verified;all8coefficients_recomputed','new_MD_production_this_round':0,'new_model_training_this_round':0,'new_model_inference_this_round':0,'new_MD_analysis_this_round':0,'figure7_panel_letters':['a','b','c'],'figure7_scatter_points':151,'figure7_total_source_points':723,'remaining_required_work':[]})
write(ROOT/'delivery_status.json',prior);write(pointer,prior)
(ROOT.parent/'completion_status_zh.md').write_text(f'当前交付 {DOC}，SHA256 {sha(DOC)}，28页。图7上排a=全66位点联合Hχ剖面、b=全部151角的四χ散点加八系数矩阵；下排三幅MD自由能轴共用c。a只标Q2/Q41，b保留全四角/两指标/重合点/负值，c432原数据保持。方法、结果与图注三段同步。全部723绘图源值核验、28页排版检查通过，原三线表、数学、批注与83记录统计保留。新版源包 {package}。本轮仅改图形和文稿，未新增分析、MD、训练或推断。\n',encoding='utf8')
print(json.dumps({'status':'complete','document':str(DOC),'document_sha256':sha(DOC),'pages':28,'package':str(package),'package_MB':round(package.stat().st_size/1e6,2),'checks_experiment':1149,'checks_figure':figure.get('independent_check_count',figure.get('checks_passed','see_figure_audit'))},ensure_ascii=False))
