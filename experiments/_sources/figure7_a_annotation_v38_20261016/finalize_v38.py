from pathlib import Path
from datetime import datetime,timezone,timedelta
import hashlib,json,zipfile,shutil
ROOT=Path(__file__).resolve().parent;BASE=ROOT.parent.parent
QA=ROOT/'manuscript_QA';FIG=ROOT/'figures';EXP=ROOT/'experiment';SRC=ROOT/'source'
STEM='fig7_clear_a_annotations_and_MD_curves'
DOC=BASE/'paper_v38_20261016_中文修订_图7标注说明与图例调整_三线表.docx'
SOURCE=BASE/'paper_v37_20261015_中文修订_图7删去汇总表与布局精简_三线表.docx'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
readj=lambda p:json.loads(Path(p).read_text(encoding='utf8'))
write=lambda p,x:Path(p).write_text(json.dumps(x,ensure_ascii=False,indent=2),encoding='utf8')
pages=readj(QA/'page_pixel_comparison.json')
changed=[p['page'] for p in pages if not p['pixel_identical_to_reviewed_v37']]
assert len(pages)==28
for p in pages:
    p['review']='final_page_PNG_directly_visually_reviewed' if p['page'] in changed else 'pixel_identical_to_already_reviewed_final_v37'
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
integration.update({'review':'passed_all28pages;changedpages_directlyreviewed;unchangedpages_pixelidentical;all564originalvaluesretained;151scatterpoints_exact;8associationsrecomputed;caption_updated_only','pages':28})
write(QA/'integration_audit.json',integration)
note=ROOT/'图7标注图例与指标定义说明.md'
note.write_text('v38按要求删除a图横轴标题“泛素残基编号”，保留数字刻度与Q2/Q41；图例移至a折线图右上方外侧空白，术语改为“扣除氨基酸类型均值后 ρ 0.868”。全部715条绘图数据与v37相同，b/c完整区域逐像素一致。a中的1–72号位置均属泛素A链，G10/A28/G35/A46/G47/G53无可定义χ角，Hχ不适用，因此没有绘点并在此断线；66个有效位置包括三处有效零值Pro均保留。相关校正从预测和MD参考中减去相同训练集氨基酸类型参考均值μ_a，再对两组偏差计算Spearmanρ=0.868；均值来自第0折以外四折训练有效残基的普通平均。该均值不同于b每χ层内的案例类型均值，a折线本身仍为原始Hχ。旧逐残基profile与原始冻结输出可在v29源包追溯。文稿只更新图注一段，其余正文、全部数学、原批注、11张原生三线表及其他7幅图均保留。本轮无新增MD、训练或推断。\n',encoding='utf8')
pointer=ROOT.parent/'completed_results_20261012/delivery_status.json';snapshot=QA/'delivery_status_v37_snapshot.json'
assert not snapshot.exists()
shutil.copy2(pointer,snapshot)
prior=readj(snapshot)
assert prior['document']==str(SOURCE) and sha(SOURCE)==prior['document_sha256']
oldarchive=Path(prior['new_figure_package']);assert sha(oldarchive)==prior['new_figure_package_sha256']
(ROOT/'README_zh.md').write_text('本包包含v38图7全矢量PDF/SVG、300dpi PNG、715条完整绘图数据、全部151角源表、原始角度缓存、作图与文稿核验脚本。运行figures/plot_clear_a_annotations_and_MD_curves.py可从包内source/experiment重绘。本轮只改a文字与图例，图注说明类型均值与无χ断线，数值和b/c保持。旧源包 '+str(oldarchive)+'，SHA256 '+sha(oldarchive)+'。\n',encoding='utf8')
files=[]
for directory,prefix in [(FIG,'figures'),(QA,'QA'),(EXP,'experiment'),(SRC,'source')]:
    files.extend((p,prefix+'/'+p.name) for p in sorted(directory.iterdir()) if p.is_file())
files.extend((p,p.name) for p in sorted(ROOT.iterdir()) if p.is_file() and p.name not in ['v38_package_manifest.json','delivery_status.json'])
assert len({name for p,name in files})==len(files)
manifest={'document_sha256':sha(DOC),'source_document_sha256':sha(SOURCE),'previous_source_archive':str(oldarchive),'previous_source_archive_sha256':sha(oldarchive),'files':{name:sha(p) for p,name in files}}
write(ROOT/'v38_package_manifest.json',manifest);files.append((ROOT/'v38_package_manifest.json','v38_package_manifest.json'))
package=BASE/'paper_v38_20261016_图7标注图例与源数据.zip'
with zipfile.ZipFile(package,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    for path,name in files:z.write(path,name)
with zipfile.ZipFile(package) as z:
    assert z.testzip() is None
    assert all(hashlib.sha256(z.read(name)).hexdigest()==digest for name,digest in manifest['files'].items())
prior.update({'status':'chinese_manuscript_v38_figure7_clear_a_annotations_delivered','completed_beijing':datetime.now(timezone(timedelta(hours=8))).isoformat(),'previous_document':str(SOURCE),'previous_document_sha256':sha(SOURCE),'document':str(DOC),'document_sha256':sha(DOC),'pages':28,'render_and_visual_review':f'passed_all28pages_{len(changed)}changedpages_directview_{28-len(changed)}exactpixels','content_review':'passed_a_legend_upperright;x_title_removed;type_mean_definition_clear;gaps_6AlaGly_same_chain;715valuesexact','new_figure_package':str(package),'new_figure_package_sha256':sha(package),'new_experiment_note':str(note),'changes_this_round':'Removeda_xlabel;movedlegend_upperright;renamedtypecorrelation;captionexplainsame_training_mean_and_nochi_gaps;715valuesexact;bc_pixelidentical.','new_numerical_review':'passed_root715points_exact;all151scattervalues_verified;all8sourcecoefficients_preserved','new_MD_production_this_round':0,'new_model_training_this_round':0,'new_model_inference_this_round':0,'new_MD_analysis_this_round':0,'figure7_panel_letters':['a','b','c'],'figure7_scatter_points':151,'figure7_total_source_points':715,'figure7_summary_matrix_displayed':False,'figure7_a_legend_upper_right':True,'figure7_a_x_title_displayed':False,'figure7_a_gaps':[10,28,35,46,47,53],'remaining_required_work':[]})
write(ROOT/'delivery_status.json',prior);write(pointer,prior)
(ROOT.parent/'completion_status_zh.md').write_text(f'当前交付 {DOC}，SHA256 {sha(DOC)}，28页。图7a横轴标题删除，图例移右上方；类型均值校正改为“扣除氨基酸类型均值后”，图注解释同一训练参考均值与A链六个无χ的Gly/Ala位置造成断线。全部715值与v37一致，b/c逐像素相同。只改图注一个正文段落，其余正文、三线表、数学、批注与83记录主验证保留，全部28页排版核验。新版源包 {package}。本轮无新实验。\n',encoding='utf8')
print(json.dumps({'status':'complete','document':str(DOC),'document_sha256':sha(DOC),'pages':28,'package':str(package),'package_MB':round(package.stat().st_size/1e6,2),'checks_experiment':1149,'checks_figure':figure.get('independent_check_count',figure.get('checks_passed','see_figure_audit'))},ensure_ascii=False))
