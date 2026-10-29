from pathlib import Path
from datetime import datetime,timezone,timedelta
import hashlib,json,zipfile,shutil
ROOT=Path(__file__).resolve().parent;BASE=ROOT.parent.parent
QA=ROOT/'manuscript_QA';FIG=ROOT/'figures';DATA=ROOT/'experiment';SRC=ROOT/'source';SRC.mkdir(exist_ok=True)
TMP=Path(r'C:\Users\DELL\AppData\Local\Temp\flexppi_completion_20261010')
OLD=ROOT.parent/'all_residues_figure_v29_20261014'
STEM='fig7_full_residue_and_chi_occupancy'
DOC=BASE/'paper_v30_20261014_中文修订_全残基逐角占据与自由能验证_三线表.docx'
SOURCE=BASE/'paper_v29_20261014_中文修订_全残基与同类型预测及自由能曲线_三线表.docx'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
readj=lambda p:json.loads(Path(p).read_text(encoding='utf8'))
write=lambda p,x:Path(p).write_text(json.dumps(x,ensure_ascii=False,indent=2),encoding='utf8')
integration=readj(QA/'integration_audit.json');numeric=readj(QA/'numerical_content_audit.json')
assert sha(DOC)==integration['document_sha256']==numeric['document_sha256']
assert numeric['status']=='passed' and numeric['pages']==28
pages=readj(QA/'page_pixel_comparison.json')
changed=[8,9,10,11,16,17]
assert [p['page'] for p in pages if not p['pixel_identical_to_reviewed_v29']]==changed
for p in pages:
    assert sha(TMP/'render_v30'/f"page-{p['page']}.png")==p['PNG_sha256']
    p['review']='final_page_PNG_directly_visually_reviewed' if p['page'] in changed else 'pixel_identical_to_already_reviewed_final_v29'
write(QA/'final_page_visual_review.json',{'status':'passed','pages':28,'changed_pages_directly_reviewed':changed,'unchanged_pages_pixel_identical':22,'page_reviews':pages,'figure_and_caption_on_same_page':17,'no_clipping_or_overlap':True})
fa=readj(FIG/(STEM+'.audit.json'))
for ext in ['png','pdf','svg']:assert sha(FIG/(STEM+'.'+ext))==fa['output_sha256'][ext]
for path,digest in fa['source_sha256'].items():assert sha(path)==digest
fa.update({'status':'passed','visual_review':'passed_original_PNG_independent_review_and_final_Word_page17;66completepositions151angles6GLN18angles;no_clipping_or_legend_overlap','production_asset_review':'LineTrend/plot_trend.py and heatmap/plot_composition.py fullscriptsandcompanionPNGsreviewed;cross_type_parameter_inheritance;no_application_specific_native_run','script_sha256':sha(FIG/'plot_full_residue_and_chi_occupancy.py')})
write(FIG/(STEM+'.audit.json'),fa)
integration.update({'review':'passed_all28pages;6changedpagesdirectview;22pixelidentical;all2202figurevaluesand432originalcurvepointschecked','pages':28})
write(QA/'integration_audit.json',integration)
write(QA/'independent_scientific_review.json',{'reviewer':'free_energy_validation_design','status':'passed','checks':['aHchipredictions_errors_and_all151angle_ranges_aligned_by_position','b_all6Gln_Hchi_errors_3anglefulloccupancy_aligned_by_row','Q40_largestoverestimate_explicit','Q49_broadoccupancy_included','MD_curve_vs_scalar_model_scope_clear','common_probability_color_scale','no_actual_legend_or_text_overlap'],'caption_precision_fixes_applied':['a_Fstar_le1_explicit','b_origin0_explicit'],'required_remaining_fixes':[]})
for source,target in [(OLD/'source/A003R.npz','A003R_frozen_model.npz'),(OLD/'figures/fig7_all_residues_Gln_and_MD_free_energy.statistics.json','fig7_all_residues_Gln_and_MD_free_energy.statistics.json'),(ROOT.parent/'free_energy_curves_v27_20261014/case_curves/histograms_501frames.csv','histograms_501frames.csv')]:shutil.copy2(source,SRC/target)
note=ROOT/'图7逐角信息扩展说明.md'
note.write_text('v30重构图7a/b。a将泛素链66个有效联合熵预测与原MD五窗口均值沿作者编号展开，并在相同位置对齐带符号误差和χ1–χ4的MD低F*≤1角区总宽度。全部151个可定义角的位置保留：χ1、χ2、χ3、χ4分别有66、52、23、10个；66个位置中的113个不存在角槽使用缺失，不填0。三处有效零联合熵Pro保留。b全部6个GLN按冻结预测从低到高为Q41、Q31、Q40、Q62、Q49、Q2；同一行显示联合Hχ参考/预测、带符号误差及全部3角的72个5°箱帧占比，并标出F*≤1的实际角区。Q40高估0.368 nats、Q49宽占据均明确呈现。c–e的432个原MD曲线点和低F标记逐值不变。模型输出联合熵标量，各角占据及曲线是实际MD对照，未构造逐角预测。全部角度数据来自已有轨迹0–100 ns/200 ps/501帧，没有新增MD生产或训练。独立核验全部132264条含缺失角记录、10872个有效箱、604逐角逐原点指标；六GLN9018旧角值完全同值，432旧箱匹配。新增151角源NPZ、原角长CSV、直方图、宽度/细箱熵、提取及独立核验脚本随本包提供。统计是本案例的描述性比较，没有虚构重复或置信区间。a/b联合Hχ仍采用原五个1,000帧窗口均值，其平均绝对误差分别为0.080及0.105 nats；a排序相关0.979、训练类型均值校正后0.868。正文方法、3.6.1及图注同步改写；83条记录主验证不变。11张原生三线表、36个数学元素、27条参考文献、原批注及其他7幅活动图内容均保持。最终28页，修改页8、9、10、11、16、17实际看图核查，其余22页与已审查v29逐像素一致。\n',encoding='utf8')
oldarchive=BASE/'paper_v27_20261014_自由能曲线与模型验证源数据.zip'
assert sha(oldarchive)=='6729cb49b9fb30cc4e347f38761e6df717877e957a5bae740d45c9a186dec9d5'
(ROOT/'README_zh.md').write_text('本包包含新版图7的300dpi PNG与全矢量PDF/SVG、2202条作图源数据，66位点全部151角的501帧MD角度、完整5°直方图与指标、冻结预测和原标签参考、提取与独立核验脚本。experiment/README_数据定义与复算.md详细定义每个数组、位置索引和缺失角；图7逐角信息扩展说明.md解释图形层次及模型论证。可在具有numpy、pandas、matplotlib、Pillow和Microsoft YaHei中文字体的Python环境运行figures/plot_full_residue_and_chi_occupancy.py，原始工程路径不存在时自动使用source与experiment目录。build/render/audit脚本保留原主机路径用于本轮文稿审计。原完整MD轨迹没有重复打包，提取审计保留其哈希及位置映射；完整既有83记录与原案例验证源包仍为 '+str(oldarchive)+'，SHA256 '+sha(oldarchive)+'。\n',encoding='utf8')
pointer=ROOT.parent/'completed_results_20261012/delivery_status.json'
if not (QA/'delivery_status_v29_snapshot.json').exists():shutil.copy2(pointer,QA/'delivery_status_v29_snapshot.json')
files=[]
for directory,prefix in [(FIG,'figures'),(QA,'QA'),(SRC,'source'),(DATA,'experiment')]:files.extend((p,prefix+'/'+p.name) for p in directory.iterdir() if p.is_file())
files.extend((p,p.name) for p in [note,ROOT/'README_zh.md',ROOT/'build_v30.py',ROOT/'audit_v30.py',ROOT/'export_v30_word.ps1',ROOT/'render_v30.py',ROOT/'prepare_plot_v30.py',Path(__file__)])
manifest={'document':str(DOC),'document_sha256':sha(DOC),'source_document_sha256':sha(SOURCE),'retained_model_validation_source_archive':str(oldarchive),'retained_model_validation_source_archive_sha256':sha(oldarchive),'files':{name:sha(path) for path,name in files}}
write(ROOT/'v30_package_manifest.json',manifest);files.append((ROOT/'v30_package_manifest.json','v30_package_manifest.json'))
package=BASE/'paper_v30_20261014_图7逐角占据与源数据.zip'
with zipfile.ZipFile(package,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    for path,name in files:z.write(path,name)
with zipfile.ZipFile(package) as z:assert z.testzip() is None and all(hashlib.sha256(z.read(name)).hexdigest()==digest for name,digest in manifest['files'].items())
status=readj(pointer)
status.update({'status':'chinese_manuscript_v30_figure7_rich_angle_occupancy_delivered','completed_beijing':datetime.now(timezone(timedelta(hours=8))).isoformat(),'previous_document':str(SOURCE),'previous_document_sha256':sha(SOURCE),'document':str(DOC),'document_sha256':sha(DOC),'pages':28,'render_and_visual_review':'passed_all28pages_6changedpagesdirectview_22exactpixels','content_review':'passed_aligned_model_prediction_error_all151angles_and6Gln18occupancies_original432curvedata_exact','new_figure_package':str(package),'new_figure_package_sha256':sha(package),'new_experiment_note':str(note),'changes_this_round':'RebuiltFigure7a_as_all66positionprediction_and_error_profiles_plus151chiangle_MDwidths;b_all6GlnjointHchi/error_aligned_with_fullchi1-3probability_maps;c-e_originalcurves_retained;updatedmethods3.6.1caption.','new_MD_production_this_round':0,'new_model_training_this_round':0,'new_numerical_review':'passed_132264angle_records_10872bins_604metrics_2202figurevalues_432oldMDcurvepoints_exact','new_all_residue_validation_positions':66,'new_all_Gln_validation_positions':6,'new_all_defined_angle_positions':151,'new_same_type_angle_positions':18,'remaining_required_work':[]})
write(ROOT/'delivery_status.json',status);write(pointer,status)
(ROOT.parent/'completion_status_zh.md').write_text(f'当前交付 {DOC}，SHA256 {sha(DOC)}，28页。图7a按位置对齐66个联合熵预测/参考、误差及全部151个χ角MD低F角区宽度；b按预测排序对齐全部六GLN的联合熵/误差与18个角的完整占据；c–e保留原432曲线数据。正文方法、3.6.1与图注已更新，83记录主验证仍在正文。11张原生三线表、36个数学元素、27条参考文献及原批注保持，全部28页审核通过。新版完整图7源包 {package}，已有83记录原实验源包沿用 {oldarchive}。本轮分析既有MD，没有新增生产模拟或训练，当前要求已完成。\n',encoding='utf8')
print(json.dumps({'status':'complete','document':str(DOC),'document_sha256':sha(DOC),'pages':28,'package':str(package),'package_MB':round(package.stat().st_size/1e6,2)},ensure_ascii=False))
