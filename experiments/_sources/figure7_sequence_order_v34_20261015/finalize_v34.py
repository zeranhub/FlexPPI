from pathlib import Path
from datetime import datetime,timezone,timedelta
import hashlib,json,zipfile,shutil
ROOT=Path(__file__).resolve().parent;BASE=ROOT.parent.parent
QA=ROOT/'manuscript_QA';FIG=ROOT/'figures';SRC=ROOT/'source'
STEM='fig7_sequence_order_predictions_and_MD_curves'
DOC=BASE/'paper_v34_20261015_中文修订_图7六谷氨酰胺顺序对照_三线表.docx'
SOURCE=BASE/'paper_v33_20261014_中文修订_图7预测合并与自由能曲线_三线表.docx'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
readj=lambda p:json.loads(Path(p).read_text(encoding='utf8'))
write=lambda p,x:Path(p).write_text(json.dumps(x,ensure_ascii=False,indent=2),encoding='utf8')
integration=readj(QA/'integration_audit.json');numeric=readj(QA/'numerical_content_audit.json')
visual=readj(QA/'final_page_visual_review.json');scientific=readj(QA/'independent_scientific_review.json')
assert visual['status']==scientific['status']==numeric['status']=='passed'
assert sha(DOC)==integration['document_sha256']==numeric['document_sha256']
fa=readj(FIG/(STEM+'.audit.json'))
for ext in ['png','pdf','svg']:assert sha(FIG/(STEM+'.'+ext))==fa['output_sha256'][ext]
integration.update({'review':'passed_all28pages_with_changedpages_directly_reviewed;576values_retained;all6sequence_order;3paragraphs_updated','pages':28})
write(QA/'integration_audit.json',integration)
note=ROOT/'图7六谷氨酰胺顺序对照说明.md'
note.write_text('v34将图7a右侧六个谷氨酰胺按残基编号Q2、Q31、Q40、Q41、Q49、Q62排列，删除全部位点的带符号误差数字列，保留MD参考与模型预测的配对点。Q2和Q41沿用紫色与橙色，对应下方冻结模型预测最高/最低端点的原逐角MD自由能曲线。a左66位点及轴标签、a右Hχ数值范围、b–d坐标与数据保持。全部576条源数据保留；12个六GLN点仅行位置改变，Hχ及原signed_error仍在CSV，其他564点数值完全一致；0.105 nats六GLN MAE包含Q40。方法、结果及图注三段同步，删除Q40误差独句，未加入机制解释。11张原生三线表、36个数学元素、27条参考文献及原批注保留。未新增轨迹分析、MD生产或训练。\n',encoding='utf8')
pointer=ROOT.parent/'completed_results_20261012/delivery_status.json'
snapshot=QA/'delivery_status_v33_snapshot.json'
assert not snapshot.exists()
shutil.copy2(pointer,snapshot)
prior=readj(snapshot)
assert prior['document']==str(SOURCE) and sha(SOURCE)==prior['document_sha256']
oldarchive=Path(prior['new_figure_package']);assert sha(oldarchive)==prior['new_figure_package_sha256']
(ROOT/'README_zh.md').write_text('本包包含v34图7PNG、全矢量PDF/SVG、576条完整绘图数据、源数据、作图与Word集成/核验脚本。作图脚本优先使用随包source，所需numpy/pandas/matplotlib/Pillow及Microsoft YaHei字体。六GLN所有参考值与预测值未变，signed_error仍保留在源CSV但不重复绘制为文字列。六GLN按序列位置展示，b–d仍选原冻结模型预测最低Q41与最高Q2。旧版源包 '+str(oldarchive)+'，SHA256 '+sha(oldarchive)+'。旧包可追溯全部角度提取与83记录群体验证来源。\n',encoding='utf8')
files=[]
for directory,prefix in [(FIG,'figures'),(QA,'QA'),(SRC,'source')]:
    files.extend((p,prefix+'/'+p.name) for p in sorted(directory.iterdir()) if p.is_file())
files.extend((p,p.name) for p in sorted(ROOT.iterdir()) if p.is_file() and p.name not in ['v34_package_manifest.json','delivery_status.json'])
assert len({name for p,name in files})==len(files)
manifest={'document_sha256':sha(DOC),'source_document_sha256':sha(SOURCE),'previous_source_archive':str(oldarchive),'previous_source_archive_sha256':sha(oldarchive),'files':{name:sha(p) for p,name in files}}
write(ROOT/'v34_package_manifest.json',manifest);files.append((ROOT/'v34_package_manifest.json','v34_package_manifest.json'))
package=BASE/'paper_v34_20261015_图7顺序对照与源数据.zip'
with zipfile.ZipFile(package,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    for path,name in files:z.write(path,name)
with zipfile.ZipFile(package) as z:
    assert z.testzip() is None
    assert all(hashlib.sha256(z.read(name)).hexdigest()==digest for name,digest in manifest['files'].items())
changed=visual['changed_pages_directly_reviewed']
prior.update({'status':'chinese_manuscript_v34_figure7_sequence_order_delivered','completed_beijing':datetime.now(timezone(timedelta(hours=8))).isoformat(),'previous_document':str(SOURCE),'previous_document_sha256':sha(SOURCE),'document':str(DOC),'document_sha256':sha(DOC),'pages':28,'render_and_visual_review':f'passed_all28pages_{len(changed)}changedpages_directview_{28-len(changed)}exactpixels','content_review':'passed_all66_andall6_retained;sequence_order;redundant_error_numbers_removed;original_b-d_curves','new_figure_package':str(package),'new_figure_package_sha256':sha(package),'new_experiment_note':str(note),'changes_this_round':'Retainedall6GLN_in_sequence_order;removedredundant_signed_error_numbercolumn;pairedreferencepredictionretained;all576valuesandall432MDcurvesunchanged;methodsresultscaptionupdated.','new_numerical_review':'passed_all576valuesretained_12Gln_row_positions_only_changed_432MDcurvepoints_exact','new_MD_production_this_round':0,'new_model_training_this_round':0,'new_MD_analysis_this_round':0,'six_Gln_display_order':[2,31,40,41,49,62],'signed_error_numbers_displayed':False,'Q40_retained_in_all_stats':True,'remaining_required_work':[]})
write(ROOT/'delivery_status.json',prior);write(pointer,prior)
(ROOT.parent/'completion_status_zh.md').write_text(f'当前交付 {DOC}，SHA256 {sha(DOC)}，28页。图7a左完整66位点剖面保留，a右全六GLN按残基编号排列，删除全部带符号误差数字列，保留参考与预测配对点。所有576原绘图数值未变，六点仅改变行序；Q40及全体统计保留。b–d自由能曲线不变。方法、3.6.1结果与图注同步修改，83记录群体验证与三线表保留。所有28页排版核验通过。新版源包 {package}。本轮无新增实验，当前要求已完成。\n',encoding='utf8')
print(json.dumps({'status':'complete','document':str(DOC),'document_sha256':sha(DOC),'pages':28,'package_MB':round(package.stat().st_size/1e6,2)},ensure_ascii=False))
