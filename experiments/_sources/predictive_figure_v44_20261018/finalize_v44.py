"""Publish the audited version and exact editable figure/data package."""
from pathlib import Path
from datetime import datetime,timezone,timedelta
import hashlib,json,shutil,zipfile

ROOT=Path(__file__).resolve().parent;BASE=ROOT.parent.parent;QA=ROOT/'manuscript_QA'
DOC=BASE/'paper_v44_20261018_中文修订_预测一致性与高动态识别_三线表.docx'
TMP=Path(r'C:\Users\DELL\AppData\Local\Temp\flexppi_completion_20261010\render_v44')
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
read=lambda p:json.loads(Path(p).read_text(encoding='utf8'))
def write(p,value):Path(p).write_text(json.dumps(value,ensure_ascii=False,indent=2),encoding='utf8')
digest=sha(DOC);integration=read(QA/'integration_audit.json')
assert digest==integration['document_sha256']
for name in ['independent_package_audit.json','independent_scientific_audit.json']:
    audit=read(QA/name)
    assert audit['document_sha256']==digest and audit['status'] in ('pass','passed')
for name in ['data_audit.json','visual_QA.json']:
    assert read(ROOT/'figure2'/name)['status'] in ('pass','passed')
page_map=read(QA/'v44_page_map.json');pages=len(page_map)
reviews=[]
for p in sorted(QA.glob('visual_pages_*.json')):
    report=read(p)
    assert report['status'] in ('pass','passed') and report['document_sha256']==digest
    for row in report['pages']:
        n=row['page'];assert row['status'] in ('pass','passed')
        reviews.append({'page':n,'status':'passed','reviewer':report.get('reviewer',p.stem),
         'sha256':sha(TMP/f'page-{n}.png')})
reviews.sort(key=lambda row:row['page'])
assert [row['page'] for row in reviews]==list(range(1,pages+1))
write(QA/'final_page_visual_review.json',{'status':'passed','document_sha256':digest,'pages':reviews,
 'method':'Every final150dpi page PNG individually opened at original resolution by assigned reviewers; root inspected the complete new Figure2 and its manuscript page.',
 'required_fixes':[]})
integration.update(status='complete',pages=pages,render_status='passed',
 visual_review=f'passed_all{pages}pages_individually_reviewed',remaining_required_work=[])
write(QA/'integration_audit.json',integration)
old_map=read(ROOT.parent/'figure_merge_v43_20261018/CURRENT_FIGURE_MAP.json')
old_map['2']={'source_folder':'figure2','absolute_source_folder':str(ROOT/'figure2'),
 'png':'fig2_predictive_capability.png','pdf':'fig2_predictive_capability.pdf','svg':'fig2_predictive_capability.svg',
 'new_in_revision':'v44','panels':['a prediction agreement RMSF','b prediction agreement Hchi',
 'c high-dynamic PR RMSF','d high-dynamic PR Hchi','e within-type RMSF gain','f within-type Hchi gain']}
for entry in old_map.values():
    for ext in ['png','pdf','svg']:assert (Path(entry['absolute_source_folder'])/entry[ext]).is_file()
write(ROOT/'CURRENT_FIGURE_MAP.json',old_map)
note=ROOT/'修改说明_zh.md'
note.write_text(f'''v44 预测一致性与高动态位点识别，2026-10-18

图2按数值预测、位点筛选、同类型增益三层重绘。a/b使用五折折外预测中全部RMSF/Hχ有效位置，按每条记录总权相同绘制二维密度，加入理想一致线和参考值分箱下的预测趋势；MAE及Spearman沿用逐记录计算再等权平均。c/d把每条记录MD第90百分位及以上位置作为高值集合，保留边界全部并列值，以连续预测分数形成PR曲线。各记录先算阶梯AP和AUROC，再等权汇总；宏PR曲线与平均AP分开定义。并列预测按相同阈值整体处理，不用残基排列强行拆开。无两类参考标签的记录明确标记不定义，并保存覆盖与原因。同类型Hχ筛选另在共同位置上核对。

e/f精确继承v43的同类型MAE降幅，包含全部20种RMSF与18种Hχ类型、A/G不适用位置及P/V负收益。旧MAE累积分布和链数分组图不再占用主图，完整统计和源数据保留于v43源包；表1及正文保留整体误差和完整/链内对照。本文不把AUROC替代连续回归误差，不新增ACC。

方法2.6、结果3.1/3.2及图注已同步，并保留图1/3/4/5/S1、KIX/KRAS/泛素背景、全部10张原生三线表、36项数学、28条文献和批注锚点。统计复算使用现有逐残基预测，没有新模型训练、推理或MD。数值与独立审计见analysis及text_review；包级与{pages}页逐张视觉检查见manuscript_QA。
''',encoding='utf8')
(ROOT/'README_zh.md').write_text('''v44 可复算预测评价与图2

当前中文稿位于上两级，文件名paper_v44_20261018_中文修订_预测一致性与高动态识别_三线表.docx。CURRENT_FIGURE_MAP.json索引全部有效图源。

analysis保存预先确定的评价协议、全部输入身份核对、折外逐残基预测/参考、各记录PR阈值、AP/AUROC、实际阳性比例、记录等权宏曲线及分箱统计。图例中的AP是记录AP的平均，不是宏PR曲线下面积。Hχ同类型评价采用记录内类型等权再记录等权。图像中不裁掉尾部数据，也不下采样。

figure2提供PNG、PDF、SVG、作图脚本、中文图注和数据/视觉核对。text_review提供独立科学复核与实际采用的逐段文本。build_v44.py仅替换原图2媒体部件和列明的正文；manuscript_QA保存文档结构、全部页面视觉检查、正文和页码。export_v44_word.ps1及render_v44.py记录实际渲染过程。其余图及旧图2的完整诊断来源见CURRENT_FIGURE_MAP.json和v43源数据包，未新增训练或MD。
''',encoding='utf8')
pointer=ROOT.parent/'completed_results_20261012/delivery_status.json'
snapshot=QA/'delivery_status_v43_snapshot.json'
if not snapshot.exists():shutil.copy2(pointer,snapshot)
prior=read(snapshot);assert prior['document_sha256']==integration['source_sha256']
files=sorted(p for p in ROOT.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.name not in ('v44_edit_manifest.json','delivery_status.json'))
manifest={'document':str(DOC),'document_sha256':digest,'previous_document':integration['source'],
 'files':{p.relative_to(ROOT).as_posix():sha(p) for p in files}}
write(ROOT/'v44_edit_manifest.json',manifest)
package=BASE/'paper_v44_20261018_预测评估图与源数据.zip'
with zipfile.ZipFile(package,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=5) as z:
    for p in files:z.write(p,p.relative_to(ROOT).as_posix())
    z.write(ROOT/'v44_edit_manifest.json','v44_edit_manifest.json')
with zipfile.ZipFile(package) as z:
    assert z.testzip() is None
    for name,d in manifest['files'].items():assert hashlib.sha256(z.read(name)).hexdigest()==d
prior.update(status='chinese_manuscript_v44_predictive_evaluation_delivered',
 completed_beijing=datetime.now(timezone(timedelta(hours=8))).isoformat(),
 document=str(DOC),document_sha256=digest,pages=pages,
 previous_document=integration['source'],previous_document_sha256=integration['source_sha256'],
 changes_this_round='Replaced Figure2 with prediction/reference density, high-dynamic-site macro precision-recall and within-type MAE gain; added AP/AUROC evaluation and synchronized Methods/Results/caption.',
 content_review='passed_independent_numerical_scientific_and_native_package_preservation_audits',
 render_and_visual_review=f'passed_all{pages}pages_individually_reviewed',
 document_edit_package=str(package),document_edit_package_sha256=sha(package),
 new_figure_package=str(package),new_figure_package_sha256=sha(package),
 new_experiment_note=str(note),predictive_evaluation_revision_note=str(note),
 current_figure_map=str(ROOT/'CURRENT_FIGURE_MAP.json'),
 new_predictive_evaluation_records=3063,new_high_dynamic_threshold='within-record MD90th percentile, inclusive ties',
 new_predictive_evaluation_AP='record-equal step average precision',new_predictive_evaluation_AUROC='record-equal AUROC',
 new_MD_production_this_round=0,new_model_training_this_round=0,new_model_inference_this_round=0,new_MD_analysis_this_round=0,
 remaining_required_work=[])
write(ROOT/'delivery_status.json',prior);write(pointer,prior)
(ROOT.parent/'completion_status_zh.md').write_text(f'当前交付 {DOC}，SHA256 {digest}，{pages}页。图2包含预测一致性、PR曲线及同类型MAE收益；方法和结果同步，AP/AUROC基于全部现有折外预测复算。全部10张原生三线表、36项数学、28参考及其余图保持；全部页面逐张检查通过。\n',encoding='utf8')
print(json.dumps({'status':'complete','document':str(DOC),'document_sha256':digest,'pages':pages,'figure_sources':str(package)},ensure_ascii=False))
