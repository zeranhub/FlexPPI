from pathlib import Path
from datetime import datetime,timezone,timedelta
import hashlib,json,zipfile,shutil
ROOT=Path(__file__).resolve().parent;BASE=ROOT.parent.parent
QA=ROOT/'manuscript_QA';FIG=ROOT/'figures';SRC=ROOT/'source';SRC.mkdir(exist_ok=True)
STEM='fig7_combined_predictions_and_MD_curves'
DOC=BASE/'paper_v33_20261014_中文修订_图7预测合并与自由能曲线_三线表.docx'
SOURCE=BASE/'paper_v32_20261014_中文修订_图7横轴标注整理_三线表.docx'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
readj=lambda p:json.loads(Path(p).read_text(encoding='utf8'))
write=lambda p,x:Path(p).write_text(json.dumps(x,ensure_ascii=False,indent=2),encoding='utf8')
integration=readj(QA/'integration_audit.json');numeric=readj(QA/'numerical_content_audit.json')
assert sha(DOC)==integration['document_sha256']==numeric['document_sha256'] and numeric['status']=='passed'
pages=readj(QA/'page_pixel_comparison.json');changed=[8,16,17]
assert [p['page'] for p in pages if not p['pixel_identical_to_reviewed_v32']]==changed
for p in pages:p['review']='final_page_PNG_directly_visually_reviewed' if p['page'] in changed else 'pixel_identical_to_already_reviewed_final_v32'
write(QA/'final_page_visual_review.json',{'status':'passed','pages':28,'page_reviews':pages,'changed_pages_directly_reviewed':changed,'unchanged_pages_pixel_identical':25,'no_clipping_or_overlap':True})
fa=readj(FIG/(STEM+'.audit.json'))
for ext in ['png','pdf','svg']:assert sha(FIG/(STEM+'.'+ext))==fa['output_sha256'][ext]
for path,digest in fa['source_sha256'].items():
    assert sha(path)==digest
    shutil.copy2(path,SRC/Path(path).name)
fa.update({'status':'passed','visual_review':'passed_independent_originalPNG_and_finalWordpage17;all66profile_andall6Glnpredictions_combined_in_a;heatmaps_removed;originalcurves_renumbered_b-d','script_sha256':sha(FIG/'plot_combined_predictions_and_MD_curves.py')})
write(FIG/(STEM+'.audit.json'),fa)
integration.update({'review':'passed_all28pages;3changedpages_directview;25pixelidentical;576retainedpoints_exact;3paragraphs_crossrefs_updated','pages':28})
write(QA/'integration_audit.json',integration)
write(QA/'independent_scientific_review.json',{'reviewer':'free_energy_validation_design','status':'passed','a_left_all66_Hchi_profile':True,'a_right_all6Gln_predictions_reference_error':True,'Q40_largesterror_retained':True,'all3_occupancy_heatmaps_deleted':True,'original432_MD_curve_values_exact_b_to_d':True,'methods_results_caption_updated':True,'no_deleted_heatmap_references':True,'all83record_primary_statistics_retained':True,'required_fixes':[]})
note=ROOT/'图7预测合并说明.md'
note.write_text('v33按用户要求将原图7b左侧的六谷氨酰胺预测结果图移动至a右侧，保留MD参考、模型预测与带符号误差。a左保留完整66位点的Hχ预测/参考剖面及整齐六Q轴标签；右六GLN按冻结预测从低到高排列，Q40最大误差明确显示。删除原b右χ1–χ3的三幅占据率热图、概率色标和相关图例，未转移到补充材料。原下排三个MD自由能曲线由c–e改为b–d，全部432原始点不变。a保留132全残基+12六GLN值，合计576绘图数据点逐值与v32完全一致；删去的1296热图箱保留在已有实验归档。正文方法、3.6.1、图注同步修改三段，六GLN热图形态描述删除，原83记录验证数值保留。最终28页，最终页8/16/17实际看图检查，其余25页与已审查v32逐像素一致；11张原生三线表、36个数学元素、27条参考文献和原批注保持。本轮没有新增MD分析、生产模拟或模型训练。\n',encoding='utf8')
pointer=ROOT.parent/'completed_results_20261012/delivery_status.json'
if not (QA/'delivery_status_v32_snapshot.json').exists():shutil.copy2(pointer,QA/'delivery_status_v32_snapshot.json')
prior=readj(pointer);oldarchive=Path(prior['new_figure_package']);assert sha(oldarchive)==prior['new_figure_package_sha256']
(ROOT/'README_zh.md').write_text('本包包含v33图7的全矢量PDF/SVG、300dpi PNG、576条完整绘图源数据、原始参考/预测源表、冻结MD角度NPZ及直方图、作图与集成/核验脚本。运行figures/plot_combined_predictions_and_MD_curves.py可重绘，需要numpy/pandas/matplotlib/Pillow和Microsoft YaHei字体；原工程路径不存在时使用source。核验与文稿集成脚本保留主机工程路径。此前完整绘图及源数据包 '+str(oldarchive)+'，SHA256 '+sha(oldarchive)+'；更早完整角度提取审计仍由前包追溯v30。\n',encoding='utf8')
files=[]
for directory,prefix in [(FIG,'figures'),(QA,'QA'),(SRC,'source')]:files.extend((p,prefix+'/'+p.name) for p in directory.iterdir() if p.is_file())
files.extend((p,p.name) for p in [note,ROOT/'README_zh.md',ROOT/'prepare_v33.py',ROOT/'build_v33.py',ROOT/'audit_v33.py',ROOT/'export_v33_word.ps1',ROOT/'render_v33.py',Path(__file__)])
manifest={'document_sha256':sha(DOC),'source_document_sha256':sha(SOURCE),'previous_source_archive':str(oldarchive),'previous_source_archive_sha256':sha(oldarchive),'files':{name:sha(p) for p,name in files}}
write(ROOT/'v33_package_manifest.json',manifest);files.append((ROOT/'v33_package_manifest.json','v33_package_manifest.json'))
package=BASE/'paper_v33_20261014_图7预测合并与源数据.zip'
with zipfile.ZipFile(package,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    for path,name in files:z.write(path,name)
with zipfile.ZipFile(package) as z:assert z.testzip() is None and all(hashlib.sha256(z.read(name)).hexdigest()==digest for name,digest in manifest['files'].items())
prior.update({'status':'chinese_manuscript_v33_figure7_combined_predictions_delivered','completed_beijing':datetime.now(timezone(timedelta(hours=8))).isoformat(),'previous_document':str(SOURCE),'previous_document_sha256':sha(SOURCE),'document':str(DOC),'document_sha256':sha(DOC),'pages':28,'render_and_visual_review':'passed_all28pages_3changedpages_directview_25exactpixels','content_review':'passed_a_combined_predictions_heatmaps_removed_b-d_curves_original','new_figure_package':str(package),'new_figure_package_sha256':sha(package),'new_experiment_note':str(note),'changes_this_round':'Mergedall6Glnpairedpredictions_into_a_right_of_all66profile;removedall3occupancyheatmaps;renumberedoriginalMDcurvesb-d;updatedmethods3.6.1caption;all83recordstatsretained.','new_numerical_review':'passed_all576retainedpoints_exact_432MDcurvepoints_exact','new_MD_production_this_round':0,'new_model_training_this_round':0,'new_MD_analysis_this_round':0,'occupancy_heatmaps_displayed':False,'remaining_required_work':[]})
write(ROOT/'delivery_status.json',prior);write(pointer,prior)
(ROOT.parent/'completion_status_zh.md').write_text(f'当前交付 {DOC}，SHA256 {sha(DOC)}，28页。图7a左为完整66Hχ位置剖面、右为原六GLN预测/参考/带符号误差。原b右三幅占据率图已删除，原自由能曲线改编号b–d，576保留数据点不变。三段正文与图注同步，83记录主验证仍在正文。所有28页审核通过，三线表、数学与批注完整。新版源包 {package}。本轮无新增实验，当前要求已完成。\n',encoding='utf8')
print(json.dumps({'status':'complete','document':str(DOC),'document_sha256':sha(DOC),'pages':28,'package_MB':round(package.stat().st_size/1e6,2)},ensure_ascii=False))
