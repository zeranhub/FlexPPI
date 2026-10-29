from pathlib import Path
from datetime import datetime,timezone,timedelta
import hashlib,json,zipfile,shutil
ROOT=Path(r'E:\Project\FlexPPI\experiments\_sources\figure7_all_chi_v35_20261015');BASE=ROOT.parent.parent
QA=ROOT/'manuscript_QA';FIG=ROOT/'figures';EXP=ROOT/'experiment';SRC=ROOT/'source'
STEM='fig7_all_chi_associations_and_MD_curves'
DOC=BASE/'paper_v35_20261015_中文修订_图7全残基逐χ占据关联_三线表.docx'
SOURCE=BASE/'paper_v34_20261015_中文修订_图7六谷氨酰胺顺序对照_三线表.docx'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
readj=lambda p:json.loads(Path(p).read_text(encoding='utf8'))
write=lambda p,x:Path(p).write_text(json.dumps(x,ensure_ascii=False,indent=2),encoding='utf8')
pages=readj(QA/'page_pixel_comparison.json');changed=[7,8,9,10,11,16,17,18]
assert [p['page'] for p in pages if not p['pixel_identical_to_reviewed_v34']]==changed
for p in pages:
    p['review']='final_page_PNG_directly_visually_reviewed' if p['page'] in changed else 'pixel_identical_to_already_reviewed_final_v34'
write(QA/'final_page_visual_review.json',{'status':'passed','pages':28,'page_reviews':pages,'changed_pages_directly_reviewed':changed,'unchanged_pages_pixel_identical':20,'no_new_clipping_overlap_or_pagination_defects':True,'figure7_PNG_directly_reviewed':True,'three_line_tables_preserved':True,'reviewer':'root'})
integration=readj(QA/'integration_audit.json');numeric=readj(QA/'numerical_content_audit.json');scientific=readj(QA/'independent_scientific_review.json');experiment=readj(EXP/'audit.json');figure=readj(FIG/(STEM+'.audit.json'))
assert scientific['status']==numeric['status']==experiment['status']=='passed'
assert sha(DOC)==integration['document_sha256']==numeric['document_sha256']
assert sha(DOC) in json.dumps(scientific)
for ext in ['png','pdf','svg']:assert sha(FIG/(STEM+'.'+ext))==figure['output_sha256'][ext]
assert experiment['checks_passed']==1149
for name,digest in experiment['output_sha256'].items():assert sha(EXP/name)==digest
for name,item in experiment['original_source_provenance'].items():assert sha(item['read_path'])==item['original_sha256']
integration.update({'review':'passed_all28pages;8changedpages_directlyreviewed;20pixelidentical;all564originalvaluesretained;8newassociationsrecomputed;4paragraphsupdated','pages':28})
write(QA/'integration_audit.json',integration)
note=ROOT/'图7全残基逐χ占据关联说明.md'
note.write_text('v35将图7a右侧六谷氨酰胺配对图替换为泛素全部可定义χ角的分层占据关联。χ1–χ4分别含66、52、23、10个位点，共151个有效角；113个未定义槽保持缺失，未补零，三个Pro保留。每角层内，将冻结联合Hχ预测及单角MD指标分别减去该层同氨基酸类型均值，再计算全体残差的Spearman相关；单例类型及常数组保留为零残差并计入并列秩。细箱占据熵相关为0.497505、0.447186、0.524704、0.127273；低F*角区总宽度相关为0.195150、0.395550、0.639980、−0.182930，全部四角与两指标显示。该图评价联合读数与逐角MD实际占据的对应，未声称模型直接输出χ1–χ4预测。\n\na左保留全部66个有效残基Hχ曲线及原统计，仅标注Q2与Q41。b–d的原432条MD曲线点及坐标保留，端点选择仍为泛素六个谷氨酰胺内联合Hχ预测最低Q41与最高Q2。MD逐角指标由已有0–100 ns、200 ps间隔501帧提取，5°箱、四原点平均协议与正文一致；曲线仍显示原点0°。原主评价联合熵参考仍采用5,000帧分成五窗后平均。\n\n新图572条绘图记录由原132个全序列预测/参考值、432曲线值及8项新分层相关组成；原六GLN的12个重复配对点从图中移除，其源值仍保存在source与实验表内，全部66位点与Q40完整保留。方法、结果、图注与讨论四段同步修改，83记录主验证与其全部数值保留。原生三线表11张、数学元素36个、参考文献27条及原批注完整保留；28页全部核验，8变化页直接查看、20页与审查通过v34逐像素一致。未新增MD生产模拟、训练或推断；本轮只对已有角度缓存复算。\n',encoding='utf8')
pointer=ROOT.parent/'completed_results_20261012/delivery_status.json';snapshot=QA/'delivery_status_v34_snapshot.json'
assert not snapshot.exists()
shutil.copy2(pointer,snapshot)
prior=readj(snapshot)
assert prior['document']==str(SOURCE) and sha(SOURCE)==prior['document_sha256']
oldarchive=Path(prior['new_figure_package']);assert sha(oldarchive)==prior['new_figure_package_sha256']
(ROOT/'README_zh.md').write_text('本包包含v35图7PNG、全矢量PDF/SVG、572条完整绘图数据、151个有效角与264个角槽明细、四角双指标关联、原角度缓存、作图/分析/文稿集成脚本及核验。experiment/analyze_all_chi_associations.py优先使用随包source，可从原501帧角度重算全部604个角与箱原点组合的细箱熵及低F*宽度。figures/plot_all_chi_associations_and_MD_curves.py读取包内source与experiment统计。所需numpy/pandas/matplotlib/Pillow与Microsoft YaHei字体。原始Hχ参考与预测未变，全部66残基保留；本轮没有新训练或生产MD。旧版源包 '+str(oldarchive)+'，SHA256 '+sha(oldarchive)+'；可追溯早期角度提取与83记录群体验证。\n',encoding='utf8')
files=[]
for directory,prefix in [(FIG,'figures'),(QA,'QA'),(EXP,'experiment'),(SRC,'source')]:
    files.extend((p,prefix+'/'+p.name) for p in sorted(directory.iterdir()) if p.is_file())
files.extend((p,p.name) for p in sorted(ROOT.iterdir()) if p.is_file() and p.name not in ['v35_package_manifest.json','delivery_status.json'])
assert len({name for p,name in files})==len(files)
manifest={'document_sha256':sha(DOC),'source_document_sha256':sha(SOURCE),'previous_source_archive':str(oldarchive),'previous_source_archive_sha256':sha(oldarchive),'files':{name:sha(p) for p,name in files}}
write(ROOT/'v35_package_manifest.json',manifest);files.append((ROOT/'v35_package_manifest.json','v35_package_manifest.json'))
package=BASE/'paper_v35_20261015_图7逐χ占据关联与源数据.zip'
with zipfile.ZipFile(package,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    for path,name in files:z.write(path,name)
with zipfile.ZipFile(package) as z:
    assert z.testzip() is None
    assert all(hashlib.sha256(z.read(name)).hexdigest()==digest for name,digest in manifest['files'].items())
prior.update({'status':'chinese_manuscript_v35_figure7_allchi_associations_delivered','completed_beijing':datetime.now(timezone(timedelta(hours=8))).isoformat(),'previous_document':str(SOURCE),'previous_document_sha256':sha(SOURCE),'document':str(DOC),'document_sha256':sha(DOC),'pages':28,'render_and_visual_review':'passed_all28pages_8changedpages_directview_20exactpixels','content_review':'passed_all66_and151definedangles;fourchi_twoMDmetrics_associations;onlyQ2Q41labels;original_b-d_curves','new_figure_package':str(package),'new_figure_package_sha256':sha(package),'new_experiment_note':str(note),'changes_this_round':'Replacedall6GLNpairedplotwithchi1-chi4MDoccupancyassociationfor151angles;onlyQ2Q41axislabels;all564retainedoriginalvaluesunchanged;8newassociations;methodsresultscaptiondiscussionupdated.','new_numerical_review':'passed_1149experiment_checks_55figure_checks_root8coefficientsrecomputed_564originalvaluesexact','new_MD_production_this_round':0,'new_model_training_this_round':0,'new_model_inference_this_round':0,'new_MD_analysis_this_round':'cached_501frame_angles_recomputed_for151angles;no_rawtrajectory_extraction','six_Gln_display_order':None,'six_Gln_paired_predictions_displayed':False,'signed_error_numbers_displayed':False,'Q40_retained_in_all_stats':True,'new_all_chi_layers':[1,2,3,4],'new_all_chi_layer_counts':[66,52,23,10],'new_all_chi_defined_angle_count':151,'new_all_chi_unavailable_slot_count':113,'remaining_required_work':[]})
write(ROOT/'delivery_status.json',prior);write(pointer,prior)
(ROOT.parent/'completion_status_zh.md').write_text(f'当前交付 {DOC}，SHA256 {sha(DOC)}，28页。图7a左完整66位点剖面保留，仅标Q2/Q41；a右改为全部151有效χ角的四层双MD指标关联，n=66/52/23/10。全部8项关联系数均展示，χ4弱/负值保留；模型仍只预测联合Hχ。564条原图保留数值完全一致，b–d原自由能曲线保持。方法、3.6.1结果、图注与讨论同步；83记录主验证与三线表保留。所有28页排版核验通过。新版源包 {package}。本轮仅复算原角度缓存，未新增MD生产或训练。\n',encoding='utf8')
print(json.dumps({'status':'complete','document':str(DOC),'document_sha256':sha(DOC),'pages':28,'package':str(package),'package_MB':round(package.stat().st_size/1e6,2),'checks_experiment':1149,'checks_figure':55},ensure_ascii=False))
