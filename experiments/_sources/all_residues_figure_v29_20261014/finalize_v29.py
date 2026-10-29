from pathlib import Path
from datetime import datetime,timezone,timedelta
import hashlib,json,zipfile,shutil
ROOT=Path(__file__).resolve().parent;BASE=ROOT.parent.parent
QA=ROOT/'manuscript_QA';FIG=ROOT/'figures';SRC=ROOT/'source';SRC.mkdir(exist_ok=True)
TMP=Path(r'C:\Users\DELL\AppData\Local\Temp\flexppi_completion_20261010')
STEM='fig7_all_residues_Gln_and_MD_free_energy'
DOC=BASE/'paper_v29_20261014_中文修订_全残基与同类型预测及自由能曲线_三线表.docx'
SOURCE=BASE/'paper_v28_20261014_中文修订_图7三条自由能曲线_三线表.docx'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
readj=lambda p:json.loads(Path(p).read_text(encoding='utf8'))
write=lambda p,x:Path(p).write_text(json.dumps(x,ensure_ascii=False,indent=2),encoding='utf8')
integration=readj(QA/'integration_audit.json');numeric=readj(QA/'numerical_content_audit.json')
assert sha(DOC)==integration['document_sha256']==numeric['document_sha256']
assert numeric['status']=='passed' and numeric['pages']==28
pages=readj(QA/'page_pixel_comparison.json')
assert [p['page'] for p in pages if not p['pixel_identical_to_reviewed_v28']]==[8,16,17]
for p in pages:
    assert sha(TMP/'render_v29'/f"page-{p['page']}.png")==p['PNG_sha256']
    p['review']='final_page_PNG_directly_visually_reviewed' if p['page'] in [8,16,17] else 'pixel_identical_to_already_reviewed_final_v28'
write(QA/'final_page_visual_review.json',{'status':'passed','pages':28,'changed_pages_directly_reviewed':[8,16,17],'unchanged_pages_pixel_identical':25,'page_reviews':pages,'figure_and_caption_on_same_page':17,'no_clipping_or_overlap':True})
fa=readj(FIG/(STEM+'.audit.json'))
for ext in ['png','pdf','svg']:assert sha(FIG/(STEM+'.'+ext))==fa['output_sha256'][ext]
for path,digest in fa['source_sha256'].items():assert sha(path)==digest
fa['visual_review']='passed_original_PNG_independent_review_and_final_Word_page17; complete66/6candidatepoints; original432curvepoints_exact; no_clipping_or_legend_overlap'
write(FIG/(STEM+'.audit.json'),fa)
integration.update({'review':'passed_all28pages;3changedpagesdirectview;25pixelidentical;66/6predictionvaluesand432curvedatachecked','pages':28})
write(QA/'integration_audit.json',integration)
write(QA/'independent_scientific_review.json',{'reviewer':'free_energy_validation_design','status':'passed','checks':['66uniquevalidresidues','all3validzeroreferencesretained','6GLNexactlyallQwithin66','largesterrorQ40explicit','Q41/Q2predictionextrema','all66andall6metricsindependentlyrecalculated','shapesdenotesource_colorsdenotesite','captiondistinguishesfivewindowlabelsand501framecurves'],'required_fixes':[]})
for path in fa['source_sha256']:shutil.copy2(path,SRC/Path(path).name)
shutil.copy2(ROOT.parent/'curve_only_v28_20261014/figures/fig7_model_selected_three_MD_free_energy_curves.display_widths.csv',SRC/'original_curve_widths_four_origins.csv')
note=ROOT/'图7全位点对照说明.md'
note.write_text('v29在原三幅自由能曲线上方增加两幅预测评价图：a展示泛素链72个输入残基中全部66个具有Hχ参考值的位置；包含三个有效零值，六个无Hχ参考值的Gly/Ala不作为熵评价位置。b展示全部六个谷氨酰胺Q2、Q31、Q40、Q41、Q49、Q62的MD参考与冻结模型预测，逐位标注绝对误差，包括所有66个位置中误差最大的Q40（0.368 nats）。a的MAE为0.0795305065 nats、Spearman相关0.9793749878，减去训练折类型均值后的相关0.8681919720；b的MAE为0.1050252914 nats。上排数据使用原五个连续1,000帧标签窗口，来自A003R第0折冻结完整模型，未重新训练或筛选误差较小位置。Q41/Q2仍按六个位点预测最小/最大选择，下排c–e来自原0–100 ns、200 ps间隔501帧的MD曲线，原432个曲线数据与底部低自由能角区逐值未变；模型预测标量Hχ，曲线本身是MD参考。上排明确量化完整预测表现，下排检验预测端点对应的实际取向范围。正文3.6.1、方法与图注同步修改。83条记录主结果保留；原已删除的群体范围扫描图没有恢复。11张原生三线表、36个原生数学元素、27条参考文献及原批注保留，其他7幅活动图片内容不变。最终28页，修改页8、16、17实际看图检查，其余25页与已审查v28逐像素一致。本轮使用既有预测和MD源数据，无新增MD生产或训练。\n',encoding='utf8')
oldarchive=BASE/'paper_v27_20261014_自由能曲线与模型验证源数据.zip'
assert sha(oldarchive)=='6729cb49b9fb30cc4e347f38761e6df717877e957a5bae740d45c9a186dec9d5'
(ROOT/'README_zh.md').write_text('本包包含图7的PNG、矢量PDF/SVG、全部66个侧链熵有效位置和全部六个谷氨酰胺的原始预测/参考、逐位误差、全部绘图点、原MD直方图及冻结预测NPZ，以及作图、集成与核验脚本。图7全位点对照说明.md记录本轮变化与准确数值。具有numpy、pandas、matplotlib、Pillow的Python环境可运行figures/plot_all_predictions_and_curves.py重绘；中文字体使用Microsoft YaHei。脚本优先读取原始工程路径，原路径不存在时读取本包source目录。Word集成和排版核验脚本保留原主机路径，用于审计本轮文稿，不是跨平台排版工具。完整既有MD角度提取、83记录检验和相关源数据继续见 '+str(oldarchive)+'，SHA256 '+sha(oldarchive)+'。\n',encoding='utf8')
pointer=ROOT.parent/'completed_results_20261012/delivery_status.json'
if not (QA/'delivery_status_v28_snapshot.json').exists():shutil.copy2(pointer,QA/'delivery_status_v28_snapshot.json')
files=[]
for directory,prefix in [(FIG,'figures'),(QA,'QA'),(SRC,'source')]:files.extend((p,prefix+'/'+p.name) for p in directory.iterdir() if p.is_file())
files.extend((p,p.name) for p in [note,ROOT/'README_zh.md',ROOT/'build_v29.py',ROOT/'audit_v29.py',ROOT/'export_v29_word.ps1',ROOT/'render_v29.py',Path(__file__)])
manifest={'document':str(DOC),'document_sha256':sha(DOC),'source_document_sha256':sha(SOURCE),'retained_model_validation_source_archive':str(oldarchive),'retained_model_validation_source_archive_sha256':sha(oldarchive),'files':{name:sha(p) for p,name in files}}
write(ROOT/'v29_package_manifest.json',manifest);files.append((ROOT/'v29_package_manifest.json','v29_package_manifest.json'))
package=BASE/'paper_v29_20261014_图7全位点预测与曲线源数据.zip'
with zipfile.ZipFile(package,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    for path,name in files:z.write(path,name)
with zipfile.ZipFile(package) as z:assert z.testzip() is None and all(hashlib.sha256(z.read(name)).hexdigest()==digest for name,digest in manifest['files'].items())
status=readj(pointer)
status.update({'status':'chinese_manuscript_v29_figure7_all_residues_and_Gln_delivered','completed_beijing':datetime.now(timezone(timedelta(hours=8))).isoformat(),'previous_document':str(SOURCE),'previous_document_sha256':sha(SOURCE),'document':str(DOC),'document_sha256':sha(DOC),'pages':28,'render_and_visual_review':'passed_all28pages_3changedpagesdirectview_25exactpixels','content_review':'passed_all66_and_all6_predictions_displayed_original432MDcurvedata_exact','new_figure_package':str(package),'new_figure_package_sha256':sha(package),'new_experiment_note':str(note),'changes_this_round':'Addedall66Hchi_validresidue_predictionvalidationandall6Glnpairedpredictionsaboveoriginal3MDcurves; addedcompleteMAE/rho/errorsandselectionbasis; updatedmethods3.6.1caption/crossrefs; all83recordresultsretained.','new_MD_production_this_round':0,'new_model_training_this_round':0,'new_numerical_review':'passed_66valid3zero6Gln_510plottedpoints_432originalMDcurvepointsexact','new_all_residue_validation_positions':66,'new_all_Gln_validation_positions':6,'remaining_required_work':[]})
write(ROOT/'delivery_status.json',status);write(pointer,status)
(ROOT.parent/'completion_status_zh.md').write_text(f'当前交付 {DOC}，SHA256 {sha(DOC)}，28页。图7上排完整展示66个Hχ有效残基与全部六个谷氨酰胺，逐位误差包含Q40最大误差；下排c–e保留原432个MD自由能曲线数据点。正文、方法及图注更新，83记录主结果保留。11张原生三线表、36个数学元素、27条参考文献及原批注保持，所有页面检查通过。新版源数据包 {package}，完整既有MD与群体验证源数据沿用 {oldarchive}。本轮无新增MD生产或训练，当前要求已完成。\n',encoding='utf8')
print(json.dumps({'status':'complete','document':str(DOC),'document_sha256':sha(DOC),'pages':28,'package':str(package),'package_MB':round(package.stat().st_size/1e6,2)},ensure_ascii=False))
