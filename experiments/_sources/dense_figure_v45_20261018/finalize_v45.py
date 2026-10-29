"""Publish the audited version and exact editable figure/data package."""
from pathlib import Path
from datetime import datetime,timezone,timedelta
import hashlib,json,shutil,zipfile

ROOT=Path(__file__).resolve().parent;BASE=ROOT.parent.parent;QA=ROOT/'manuscript_QA'
DOC=BASE/'paper_v45_20261018_中文修订_紧凑多层模型评价图_三线表.docx'
TMP=Path(r'C:\Users\DELL\AppData\Local\Temp\flexppi_completion_20261010\render_v45')
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
assert read(ROOT/'analysis/addition_audit.json')['status']=='passed'
assert read(ROOT/'text_review/independent_metrics_audit.json')['status'] in ('pass','passed')
page_map=read(QA/'v45_page_map.json');pages=len(page_map)
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
old_map=read(ROOT.parent/'predictive_figure_v44_20261018/CURRENT_FIGURE_MAP.json')
old_map['2']={'source_folder':'figure2','absolute_source_folder':str(ROOT/'figure2'),
 'png':'fig2_dense_capability.png','pdf':'fig2_dense_capability.pdf','svg':'fig2_dense_capability.svg',
 'new_in_revision':'v45','panels':['a prediction agreement RMSF','b prediction agreement Hchi',
 'c high-dynamic PR RMSF','d high-dynamic PR Hchi','e within-type Hchi PR','f source-stratified AP']}
for entry in old_map.values():
    for ext in ['png','pdf','svg']:assert (Path(entry['absolute_source_folder'])/entry[ext]).is_file()
write(ROOT/'CURRENT_FIGURE_MAP.json',old_map)
note=ROOT/'修改说明_zh.md'
note.write_text(f'''v45 紧凑多层模型评价图，完成于2026-10-19

原图2e/f的两幅类型MAE降幅散点图退出主图，相关数值结论和完整源数据保留。新图2按连续数值比较、三个筛选任务及三个来源内的表现重排。a/b保留完整模型全部961,321个RMSF与824,988个Hχ有效位置的记录等权密度，并叠加Full/Intra/Static/AA四种方法按MD参考分箱的原始预测均值；Full的20箱趋势和1,000次来源分层重抽记录区间与v44逐值一致。RMSF完整高值尾保留，预测未重新校准。

c/d/e并列RMSF、Hχ及同类型Hχ的PR曲线。c/d沿用每记录MD第90百分位及以上的定义；e在同记录同类型内重定义阈值、先对42,108个有效记录–类型组形成曲线，再在每记录内对类型等权、最后对3,063条记录等权。图例AP仍为相同层级下的原始阶梯AP均值，分别为完整模型0.730、0.926和0.623，全部完整入选端点保留，未用宏曲线面积替代。f为四方法×三任务×三来源的36格AP比较，显示来源样本量及原始AP，只回答来源内的表现而非来源留出验证。

方法2.6、结果3.1/3.2及图注已同步，并保留图1/3/4/5/S1、KIX/KRAS/泛素背景、全部10张原生三线表、36项数学、28条文献和批注锚点。统计复算使用现有逐残基预测，没有新模型训练、推理或MD。数值与独立审计见analysis及text_review；包级与{pages}页逐张视觉检查见manuscript_QA。
''',encoding='utf8')
(ROOT/'README_zh.md').write_text('''v45 可复算预测评价与图2

当前中文稿位于上两级，文件名paper_v45_20261018_中文修订_紧凑多层模型评价图_三线表.docx。CURRENT_FIGURE_MAP.json索引全部有效图源。

analysis保存新增评价协议、四方法分箱预测趋势、同类型Hχ完整阈值曲线与两级等权宏曲线、36格来源AP，以及v44主PR及指标副本。完整Full密度和既有全部有效点/阈值数据继承v44源包；新增结果的每一项输入身份、掩码和位置顺序均重新核对。图例AP是同级逐记录AP的平均，不是宏PR曲线下面积。图像不裁尾、不下采样，也不重新校准预测。

figure2提供PNG、PDF、SVG、作图脚本、中文图注和数据/视觉核对。text_review提供独立科学复核与实际采用的逐段文本。build_v45.py仅替换图2媒体部件和列明的正文；manuscript_QA保存文档结构、全部页面视觉检查、正文和页码。export_v45_word.ps1及render_v45.py记录实际渲染过程。其余图及旧图2的完整诊断来源见CURRENT_FIGURE_MAP.json和v44/v43源数据包，未新增训练、推理或MD。
''',encoding='utf8')
pointer=ROOT.parent/'completed_results_20261012/delivery_status.json'
snapshot=QA/'delivery_status_v44_snapshot.json'
if not snapshot.exists():shutil.copy2(pointer,snapshot)
prior=read(snapshot);assert prior['document_sha256']==integration['source_sha256']
files=sorted(p for p in ROOT.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.name not in ('v45_edit_manifest.json','delivery_status.json'))
manifest={'document':str(DOC),'document_sha256':digest,'previous_document':integration['source'],
 'files':{p.relative_to(ROOT).as_posix():sha(p) for p in files}}
write(ROOT/'v45_edit_manifest.json',manifest)
package=BASE/'paper_v45_20261018_紧凑模型评价图与源数据.zip'
with zipfile.ZipFile(package,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=5) as z:
    for p in files:z.write(p,p.relative_to(ROOT).as_posix())
    z.write(ROOT/'v45_edit_manifest.json','v45_edit_manifest.json')
with zipfile.ZipFile(package) as z:
    assert z.testzip() is None
    for name,d in manifest['files'].items():assert hashlib.sha256(z.read(name)).hexdigest()==d
prior.update(status='chinese_manuscript_v45_compact_model_evaluation_delivered',
 completed_beijing=datetime.now(timezone(timedelta(hours=8))).isoformat(),
 document=str(DOC),document_sha256=digest,pages=pages,
 previous_document=integration['source'],previous_document_sha256=integration['source_sha256'],
 changes_this_round='Removed sparse type-MAE gain panels from Figure2; enlarged Full prediction/reference densities with four-method binned trends; aligned three PR tasks including within-type Hchi; added36 source-stratified AP cells and synchronized Methods/Results/caption.',
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
(ROOT.parent/'completion_status_zh.md').write_text(f'当前交付 {DOC}，SHA256 {digest}，{pages}页。图2删除旧类型MAE点图，含放大的两项预测密度与四方法趋势、三任务PR及36格来源AP；方法、结果和图注同步。全部10张原生三线表、36项数学、28参考及其余图保持；全部页面逐张检查通过。\n',encoding='utf8')
print(json.dumps({'status':'complete','document':str(DOC),'document_sha256':digest,'pages':pages,'figure_sources':str(package)},ensure_ascii=False))
