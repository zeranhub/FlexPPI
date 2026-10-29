from pathlib import Path
from datetime import datetime,timezone,timedelta
import hashlib,json,zipfile,shutil
ROOT=Path(__file__).resolve().parent;BASE=ROOT.parent.parent
QA=ROOT/'manuscript_QA';FIG=ROOT/'figures';SRC=ROOT/'source';SRC.mkdir(exist_ok=True)
STEM='fig7_clean_axis_six_Gln_and_occupancy'
DOC=BASE/'paper_v32_20261014_中文修订_图7横轴标注整理_三线表.docx'
SOURCE=BASE/'paper_v31_20261014_中文修订_图7简化与六谷氨酰胺轴标注_三线表.docx'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
readj=lambda p:json.loads(Path(p).read_text(encoding='utf8'))
write=lambda p,x:Path(p).write_text(json.dumps(x,ensure_ascii=False,indent=2),encoding='utf8')
integration=readj(QA/'integration_audit.json');numeric=readj(QA/'numerical_content_audit.json')
assert sha(DOC)==integration['document_sha256']==numeric['document_sha256'] and numeric['status']=='passed'
pages=readj(QA/'page_pixel_comparison.json')
assert [p['page'] for p in pages if not p['pixel_identical_to_reviewed_v31']]==[17]
for p in pages:p['review']='final_page_PNG_directly_visually_reviewed' if p['page']==17 else 'pixel_identical_to_already_reviewed_final_v31'
write(QA/'final_page_visual_review.json',{'status':'passed','pages':28,'page_reviews':pages,'changed_page_directly_reviewed':17,'unchanged_pages_pixel_identical':27,'no_clipping_or_overlap':True})
fa=readj(FIG/(STEM+'.audit.json'))
for ext in ['png','pdf','svg']:assert sha(FIG/(STEM+'.'+ext))==fa['output_sha256'][ext]
for path,digest in fa['source_sha256'].items():
    assert sha(path)==digest
    shutil.copy2(path,SRC/Path(path).name)
fa.update({'status':'passed','visual_review':'passed_independent_originalPNG_and_finalWordpage17;numeric_ticks_Glnlabels_title_separate_rows;Q40Q41_leaderanchors_correct_no_overlap','script_sha256':sha(FIG/'plot_clean_axis_six_Gln_and_occupancy.py')})
write(FIG/(STEM+'.audit.json'),fa)
integration.update({'review':'passed_all28pages;1changedpage_directview;27pixelidentical;all1872pointsunchanged;axisanchors_correct','pages':28})
write(QA/'integration_audit.json',integration)
write(QA/'independent_scientific_review.json',{'reviewer':'free_energy_validation_design','status':'passed','numeric_ticks_one_row':True,'all6_Gln_labels_same_second_row':True,'centered_axis_title_third_row':True,'Q40_Q41_correct_leader_anchors':True,'no_label_line_or_btitle_overlap':True,'all1872_points_exact_v31':True,'required_fixes':[]})
note=ROOT/'图7横轴整理说明.md'
note.write_text('v32仅整理图7a横轴：恢复1/10/20/30/40/50/60/72常规刻度在第一行，全部六个GLN标签在第二行，轴标题居中在第三行。Q40与Q41文字分别向左/右分开，短引线实际锚点仍为残基40/41，不改变数据位置。引线安排在数字行下方，没有跨数字或标签。全部1872个图数据点、b–e、正文及图注与v31保持一致。最终28页，仅页17有像素变化并实际检查，其余27页与已审查v31逐像素一致；三线表、数学元素、批注完整保持。本轮无新增实验。\n',encoding='utf8')
pointer=ROOT.parent/'completed_results_20261012/delivery_status.json'
if not (QA/'delivery_status_v31_snapshot.json').exists():shutil.copy2(pointer,QA/'delivery_status_v31_snapshot.json')
prior=readj(pointer)
oldarchive=Path(prior['new_figure_package']);assert sha(oldarchive)==prior['new_figure_package_sha256']
(ROOT/'README_zh.md').write_text('本包包含v32图7的矢量PDF/SVG、300dpi PNG、1872点源表与作图/集成/核验脚本；source包含所需MD概率NPZ、残基源表、统计和曲线直方图。可运行figures/plot_clean_axis_six_Gln_and_occupancy.py重绘，需numpy/pandas/matplotlib/Pillow与Microsoft YaHei字体，原工程路径不存在时使用source。完整前版图7及源数据归档 '+str(oldarchive)+'，SHA256 '+sha(oldarchive)+'。\n',encoding='utf8')
files=[]
for directory,prefix in [(FIG,'figures'),(QA,'QA'),(SRC,'source')]:files.extend((p,prefix+'/'+p.name) for p in directory.iterdir() if p.is_file())
files.extend((p,p.name) for p in [note,ROOT/'README_zh.md',ROOT/'prepare_v32.py',ROOT/'build_v32.py',ROOT/'audit_v32.py',ROOT/'export_v32_word.ps1',ROOT/'render_v32.py',Path(__file__)])
manifest={'document_sha256':sha(DOC),'source_document_sha256':sha(SOURCE),'previous_source_archive':str(oldarchive),'previous_source_archive_sha256':sha(oldarchive),'files':{name:sha(p) for p,name in files}}
write(ROOT/'v32_package_manifest.json',manifest);files.append((ROOT/'v32_package_manifest.json','v32_package_manifest.json'))
package=BASE/'paper_v32_20261014_图7横轴整理与源数据.zip'
with zipfile.ZipFile(package,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    for path,name in files:z.write(path,name)
with zipfile.ZipFile(package) as z:assert z.testzip() is None and all(hashlib.sha256(z.read(name)).hexdigest()==digest for name,digest in manifest['files'].items())
prior.update({'status':'chinese_manuscript_v32_figure7_clean_axis_delivered','completed_beijing':datetime.now(timezone(timedelta(hours=8))).isoformat(),'previous_document':str(SOURCE),'previous_document_sha256':sha(SOURCE),'document':str(DOC),'document_sha256':sha(DOC),'pages':28,'render_and_visual_review':'passed_all28pages_1changedpage_directview_27exactpixels','content_review':'passed_clean_axis_all6Glnlabels_correct_1872points_text_unchanged','new_figure_package':str(package),'new_figure_package_sha256':sha(package),'new_experiment_note':str(note),'changes_this_round':'Figure7a_axis_regularnumeric_ticks_one_row_all6Gln_labels_secondrow_centeredtitle_thirdrow_Q40Q41_leaderanchors_exact;allfiguredataandtextunchanged.','new_numerical_review':'passed_1872points_exact_and_correctaxisanchors','new_MD_production_this_round':0,'new_model_training_this_round':0,'new_MD_analysis_this_round':0,'remaining_required_work':[]})
write(ROOT/'delivery_status.json',prior);write(pointer,prior)
(ROOT.parent/'completion_status_zh.md').write_text(f'当前交付 {DOC}，SHA256 {sha(DOC)}，28页。图7a横轴数字刻度、六GLN标签和居中轴标题分别独立成行，Q40/Q41文字分开且短引线保持实际锚点。全部1872个绘图数据点与正文不变。仅页17调整，其余27页与已审查v31逐像素一致，所有页面检查通过。源包 {package}。本轮无新增实验，当前要求已完成。\n',encoding='utf8')
print(json.dumps({'status':'complete','document':str(DOC),'document_sha256':sha(DOC),'pages':28,'package_MB':round(package.stat().st_size/1e6,2)},ensure_ascii=False))
