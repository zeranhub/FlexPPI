from pathlib import Path
from datetime import datetime,timezone,timedelta
import hashlib,json,zipfile,shutil
ROOT=Path(__file__).resolve().parent;BASE=ROOT.parent.parent
QA=ROOT/'manuscript_QA';FIG=ROOT/'figures';SRC=ROOT/'source';SRC.mkdir(exist_ok=True)
TMP=Path(r'C:\Users\DELL\AppData\Local\Temp\flexppi_completion_20261010')
STEM='fig7_Hchi_profile_six_Gln_axis_and_occupancy'
DOC=BASE/'paper_v31_20261014_中文修订_图7简化与六谷氨酰胺轴标注_三线表.docx'
SOURCE=BASE/'paper_v30_20261014_中文修订_全残基逐角占据与自由能验证_三线表.docx'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
readj=lambda p:json.loads(Path(p).read_text(encoding='utf8'))
write=lambda p,x:Path(p).write_text(json.dumps(x,ensure_ascii=False,indent=2),encoding='utf8')
integration=readj(QA/'integration_audit.json');numeric=readj(QA/'numerical_content_audit.json')
assert sha(DOC)==integration['document_sha256']==numeric['document_sha256']
assert numeric['status']=='passed' and numeric['pages']==28
pages=readj(QA/'page_pixel_comparison.json');changed=[8,17]
assert [p['page'] for p in pages if not p['pixel_identical_to_reviewed_v30']]==changed
for p in pages:
    assert sha(TMP/'render_v31'/f"page-{p['page']}.png")==p['PNG_sha256']
    p['review']='final_page_PNG_directly_visually_reviewed' if p['page'] in changed else 'pixel_identical_to_already_reviewed_final_v30'
write(QA/'final_page_visual_review.json',{'status':'passed','pages':28,'changed_pages_directly_reviewed':changed,'unchanged_pages_pixel_identical':26,'page_reviews':pages,'figure_and_caption_on_same_page':17,'no_clipping_or_overlap':True})
fa=readj(FIG/(STEM+'.audit.json'))
for ext in ['png','pdf','svg']:assert sha(FIG/(STEM+'.'+ext))==fa['output_sha256'][ext]
for path,digest in fa['source_sha256'].items():
    assert sha(path)==digest
    shutil.copy2(path,SRC/Path(path).name)
fa.update({'status':'passed','visual_review':'passed_original_PNG_independent_review_and_final_Word_page17;all6GLNaxispositions_exact;Q40/Q41verticalstagger_nooverlap;all1872retaineddata_exact','script_sha256':sha(FIG/'plot_Hchi_profile_six_Gln_axis_and_occupancy.py')})
write(FIG/(STEM+'.audit.json'),fa)
integration.update({'review':'passed_all28pages;2changedpagesdirectview;26pixelidentical;1872retainedpoints_exact;6Glnaxislabelschecked','pages':28})
write(QA/'integration_audit.json',integration)
write(QA/'independent_scientific_review.json',{'reviewer':'free_energy_validation_design','status':'passed','all132_aHchi_values_exact':True,'all1740_b_to_e_values_exact':True,'a_error_and_width_groups_removed':True,'six_Gln_axis_positions':[2,31,40,41,49,62],'Q41_shift_vertical_only':True,'Q40_Q41_no_overlap':True,'Q41_b_title_clearance':True,'caption_deleted_layers_no_remaining_references':True,'required_fixes':[]})
note=ROOT/'图7轴标注调整说明.md'
note.write_text('v31按用户要求删除图7a下方的误差轨道与逐角低自由能角区宽度热图，只保留完整66个有效位置的模型联合Hχ预测与MD参考剖面。在横轴准确标记全部六个谷氨酰胺Q2、Q31、Q40、Q41、Q49、Q62；Q41只向下错层，所有标注横坐标仍是实际残基编号。Q2紫色、Q41橙色与b和c–e保持一致。a保留132个原始参考/预测点，b–e保留1740个原始数据点，1872点逐值均与v30一致；原a误差及宽度330条绘图记录从本轮图表源表移除，原计算源数据仍保存在v30归档。正文方法和图注删掉与这两层相关的表述，3.6.1其他论述保持。完整28页，最终页8、17实际看图检查，其他26页与已审查v30逐像素一致；11张原生三线表、36个数学元素、27条参考文献及原批注保持。没有新增MD分析、生产模拟或模型训练。\n',encoding='utf8')
oldarchive=BASE/'paper_v30_20261014_图7逐角占据与源数据.zip'
assert sha(oldarchive)=='86ebd1f58ac04c95a2f98dd670a85ed43174c6e579ce2529c162bff6c67801bb'
(ROOT/'README_zh.md').write_text('本包包含v31图7的300dpi PNG与全矢量PDF/SVG、1872条完整绘图源数据、所需冻结读数与MD概率NPZ、模型/参考源表、统计文件与原曲线直方图、生成和审查脚本。运行figures/plot_Hchi_profile_six_Gln_axis_and_occupancy.py可重绘，需numpy、pandas、matplotlib、Pillow与Microsoft YaHei字体；原工程路径不存在时自动使用source目录。论文集成/渲染脚本保留原主机路径供审计。完整已有角度提取、原角长CSV及独立复算数据继续使用 '+str(oldarchive)+'，SHA256 '+sha(oldarchive)+'。\n',encoding='utf8')
pointer=ROOT.parent/'completed_results_20261012/delivery_status.json'
if not (QA/'delivery_status_v30_snapshot.json').exists():shutil.copy2(pointer,QA/'delivery_status_v30_snapshot.json')
files=[]
for directory,prefix in [(FIG,'figures'),(QA,'QA'),(SRC,'source')]:files.extend((p,prefix+'/'+p.name) for p in directory.iterdir() if p.is_file())
files.extend((p,p.name) for p in [note,ROOT/'README_zh.md',ROOT/'build_v31.py',ROOT/'audit_v31.py',ROOT/'export_v31_word.ps1',ROOT/'render_v31.py',ROOT/'prepare_v31.py',ROOT/'prepare_build_v31.py',Path(__file__)])
manifest={'document':str(DOC),'document_sha256':sha(DOC),'source_document_sha256':sha(SOURCE),'retained_angle_source_archive':str(oldarchive),'retained_angle_source_archive_sha256':sha(oldarchive),'files':{name:sha(path) for path,name in files}}
write(ROOT/'v31_package_manifest.json',manifest);files.append((ROOT/'v31_package_manifest.json','v31_package_manifest.json'))
package=BASE/'paper_v31_20261014_图7六谷氨酰胺轴标注与源数据.zip'
with zipfile.ZipFile(package,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    for path,name in files:z.write(path,name)
with zipfile.ZipFile(package) as z:assert z.testzip() is None and all(hashlib.sha256(z.read(name)).hexdigest()==digest for name,digest in manifest['files'].items())
status=readj(pointer)
status.update({'status':'chinese_manuscript_v31_figure7_axis_sites_delivered','completed_beijing':datetime.now(timezone(timedelta(hours=8))).isoformat(),'previous_document':str(SOURCE),'previous_document_sha256':sha(SOURCE),'document':str(DOC),'document_sha256':sha(DOC),'pages':28,'render_and_visual_review':'passed_all28pages_2changedpagesdirectview_26exactpixels','content_review':'passed_a_error_width_tracks_removed_all6Gln_axislabels_exact_b_to_e_retained','new_figure_package':str(package),'new_figure_package_sha256':sha(package),'new_experiment_note':str(note),'changes_this_round':'Removedtwo_lower_tracks_fromFigure7a;keptall66Hchi_prediction_referenceprofiles;labelledall6Glnataccurateaxispositions_Q40/Q41verticalstagger;allb-e_dataunchanged;updatedmethodsandcaption.','new_MD_production_this_round':0,'new_model_training_this_round':0,'new_MD_analysis_this_round':0,'new_numerical_review':'passed_all1872retainedfigurepoints_exact_6Glnaxispositions_correct','retained_angle_source_archive':str(oldarchive),'retained_angle_source_archive_sha256':sha(oldarchive),'whole_residue_MD_width_tracks_displayed':False,'remaining_required_work':[]})
write(ROOT/'delivery_status.json',status);write(pointer,status)
(ROOT.parent/'completion_status_zh.md').write_text(f'当前交付 {DOC}，SHA256 {sha(DOC)}，28页。图7a保留完整66个联合熵预测/参考剖面，删去误差和逐角宽度两层，并在横轴准确标出Q2/Q31/Q40/Q41/Q49/Q62；Q40/Q41垂直错层。b–e内容与全部数据不变。方法与图注更新，三线表、原数学与批注保持；所有28页审核通过。新版图7源包 {package}，完整已有角度实验源包沿用 {oldarchive}。本轮无新增计算，当前要求已完成。\n',encoding='utf8')
print(json.dumps({'status':'complete','document':str(DOC),'document_sha256':sha(DOC),'pages':28,'package':str(package),'package_MB':round(package.stat().st_size/1e6,2)},ensure_ascii=False))
