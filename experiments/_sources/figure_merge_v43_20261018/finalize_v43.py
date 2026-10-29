"""Finalize delivery only after numerical, package and all-page visual gates pass."""
from pathlib import Path
from datetime import datetime,timezone,timedelta
import hashlib,json,shutil,zipfile

ROOT=Path(__file__).resolve().parent;BASE=ROOT.parent.parent;QA=ROOT/'manuscript_QA'
DOC=BASE/'paper_v43_20261018_中文修订_能力与构象合并图_三线表.docx'
TMP=Path(r'C:\Users\DELL\AppData\Local\Temp\flexppi_completion_20261010\render_v43')
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
read=lambda p:json.loads(Path(p).read_text(encoding='utf8'))
def write(p,value):Path(p).write_text(json.dumps(value,ensure_ascii=False,indent=2),encoding='utf8')
digest=sha(DOC)
assert digest=='5d32c2af90a6aa8268ac0ef50ba46480370fc798bcc890d91fa71123dfb51c5a'
integration=read(QA/'integration_audit.json');assert integration['document_sha256']==digest
package_audit=read(QA/'independent_package_audit.json')
assert package_audit['document_sha256']==digest and all(package_audit['checks'].values())
science=read(QA/'independent_scientific_audit.json')
assert science['document_sha256']==digest and science['status'] in ('pass','passed')
for folder in ['figure2','figure5']:
    assert read(ROOT/folder/'data_audit.json')['status'] in ('pass','passed')
assert read(ROOT/'text_review/cohort83_source_audit.json')['status'] in ('pass','passed')
page_map=read(QA/'v43_page_map.json');pages=len(page_map);assert pages==28
records=[]
for filename,first,last,reviewer in [
 ('visual_pages_1_10.json',1,10,'kras_figure6_build'),
 ('visual_pages_11_17.json',11,17,'root'),
 ('visual_pages_18_28.json',18,28,'kras_manuscript_text_review')]:
    report=read(QA/filename)
    assert report['status'] in ('pass','passed') and report['document_sha256']==digest
    for n in range(first,last+1):
        records.append({'page':n,'status':'passed','reviewer':reviewer,'sha256':sha(TMP/f'page-{n}.png')})
assert [r['page'] for r in records]==list(range(1,pages+1))
write(QA/'final_page_visual_review.json',{'status':'passed','document_sha256':digest,'pages':records,
 'method':'All28 final150dpi original page PNGs individually opened by assigned reviewers; root directly inspected every Results figure and reordered Results page11–17.',
 'required_fixes':[],'figure_caption_layout':'Both new merged figures and their complete captions fit together on page11 and page17.'})
integration.update(status='complete',pages=pages,render_status='passed',visual_review='passed_all28pages_individually_reviewed',remaining_required_work=[])
write(QA/'integration_audit.json',integration)

old_map=read(ROOT.parent/'protein_background_v42_20261017/CURRENT_FIGURE_MAP.json')
new_map={}
for new,old in [('1','1'),('3','5'),('4','6'),('S1','S1')]:
    new_map[new]=dict(old_map[old]);new_map[new]['previous_figure_number']=old
for n,stem,old in [('2','fig2_model_capability',['2','3']),('5','fig5_joint_readouts_conformations',['4','7'])]:
    new_map[n]={'source_folder':f'figure{n}','absolute_source_folder':str(ROOT/f'figure{n}'),
      'png':stem+'.png','pdf':stem+'.pdf','svg':stem+'.svg','merged_previous_figures':old,'new_in_revision':'v43'}
new_map={k:new_map[k] for k in ['1','2','3','4','5','S1']}
for entry in new_map.values():
    for kind in ['png','pdf','svg']:assert (Path(entry['absolute_source_folder'])/entry[kind]).is_file()
write(ROOT/'CURRENT_FIGURE_MAP.json',new_map)
note=ROOT/'修改说明_zh.md'
note.write_text('''v43 图2–3及图4与7合并，2026-10-18

原图2与3重组为新图2，六个子图按左列RMSF、右列Hχ排列：第一行呈现全部3,063条记录的MAE累积分布；第二行呈现同氨基酸类型内完整模型相对链内模型及静态特征对照的MAE降幅；第三行呈现按实际输入链数分组的配对误差差值。各记录数据与负收益保留，原图3的幅度诊断源数据完整归档。整体排名成绩由表1承担，删除重复的排名分布图。

原图4与7重组为新图5，先以两个读数的全队列关系引入，再展示83条记录、33,142个同类型有效位置中侧链与主链输出的配对关联差，以及仅按预测分组所得的真实MD低自由能角区宽度差；接着展示泛素全部66个有效位置的联合熵参考与预测、全部151个单角观测的关联，以及模型选出Q41/Q2后的三条逐角MD自由能曲线。曲线的原始箱值、零计数断线、阈值与色条不变。两个输出的四类组合召回保留在正文一句及完整表S9中。

正文按图的证据顺序调整为3.5.1全队列验证、3.5.2泛素具体案例。原图5（KIX及泛素结构）改为图3，原图6（KRAS G12D）改为图4，图像内容不变；全部正文、方法与补充材料交叉引用同步。保留v42新增的蛋白背景。表1及九张补充表共10张原生三线表、36项原生数学、28条文献及原批注/锚点完整保持。未新增MD、模型训练或推理。

两张新图均提供PNG、纯矢量SVG、字体嵌入PDF及作图脚本/逐点数据。源字号最低9 pt，在Word中以167.64 mm宽插入后最低约8.24 pt。当前稿共28页，已由独立隐藏Word实例只读导出，再用文档技能render_docx.py生成150 dpi逐页PNG，全部逐张检查通过。图2及完整图注位于第11页，图5及完整图注位于第17页。
''',encoding='utf8')
(ROOT/'README_zh.md').write_text('''v43 合并图与文稿整合

当前中文稿位于上两级，文件名paper_v43_20261018_中文修订_能力与构象合并图_三线表.docx。CURRENT_FIGURE_MAP.json为五张正文图及S1的有效源图索引。

figure2和figure5目录分别提供合并图PNG/PDF/SVG、中文图注、作图脚本、源数据与数值/视觉核对。新图2来源于原图2/3，新图5来源于原图4/7及已完成的83记录构象占据验证；本轮没有新MD、训练和推理。图5b展示配对相关差，不代表条件独立贡献；图5e关联同一联合Hχ输出与逐角MD指标；图5f曲线来自MD。完整背景、训练划分与此前科学结果继承v42。

build_v43.py执行实际OOXML整合；text_review给出采用的逐段正文与独立科学审计；manuscript_QA保存包级完整性、逐页视觉记录、正文和页码。export_v43_word.ps1及render_v43.py记录实际渲染。所有10张表、36项数学、28条文献及原批注锚点均保持。
''',encoding='utf8')

pointer=ROOT.parent/'completed_results_20261012/delivery_status.json'
snapshot=QA/'delivery_status_v42_snapshot.json'
if not snapshot.exists():shutil.copy2(pointer,snapshot)
prior=read(snapshot);assert prior['document_sha256']==integration['source_sha256']
files=sorted(p for p in ROOT.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.name not in ('v43_edit_manifest.json','delivery_status.json'))
manifest={'document':str(DOC),'document_sha256':digest,'previous_document':integration['source'],
 'files':{p.relative_to(ROOT).as_posix():sha(p) for p in files}}
write(ROOT/'v43_edit_manifest.json',manifest)
package=BASE/'paper_v43_20261018_合并图与文稿源数据.zip'
with zipfile.ZipFile(package,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    for p in files:z.write(p,p.relative_to(ROOT).as_posix())
    z.write(ROOT/'v43_edit_manifest.json','v43_edit_manifest.json')
with zipfile.ZipFile(package) as z:
    assert z.testzip() is None
    for name,d in manifest['files'].items():assert hashlib.sha256(z.read(name)).hexdigest()==d
prior.update(status='chinese_manuscript_v43_merged_figures_delivered',
 completed_beijing=datetime.now(timezone(timedelta(hours=8))).isoformat(),
 document=str(DOC),document_sha256=digest,pages=pages,figures=6,main_figures=5,supplementary_figures=1,
 previous_document=integration['source'],previous_document_sha256=integration['source_sha256'],
 changes_this_round='Merged previous Figures2/3 into newFigure2 and Figures4/7 into newFigure5; retained KIX/ubiquitin structuralFigure3 and KRASFigure4; reordered population-to-ubiquitin Results and synchronized all references.',
 content_review='passed_independent_source_numerical_scientific_and_package_preservation_audits',
 render_and_visual_review='passed_all28pages_individually_reviewed',
 document_edit_package=str(package),document_edit_package_sha256=sha(package),
 new_figure_package=str(package),new_figure_package_sha256=sha(package),new_experiment_note=str(note),
 merged_figures_revision_note=str(note),current_figure_map=str(ROOT/'CURRENT_FIGURE_MAP.json'),
 TCR_case_figure_number=4,ubiquitin_figure_number=5,ubiquitin_figure_panel_letters=['d','e','f'],
 joint_readouts_figure_number=5,joint_readouts_figure_panel_letters=['a','b','c'],
 model_capability_figure_number=2,model_capability_figure_panel_letters=['a','b','c','d','e','f'],
 new_MD_production_this_round=0,new_model_training_this_round=0,new_model_inference_this_round=0,new_MD_analysis_this_round=0,
 remaining_required_work=[])
write(ROOT/'delivery_status.json',prior);write(pointer,prior)
(ROOT.parent/'completion_status_zh.md').write_text(f'当前交付 {DOC}，SHA256 {digest}，28页。原图2–3合为新图2，原图4与7合为新图5；KIX及泛素结构图改为图3，KRAS实例改为图4。保留全部10张原生三线表、36项数学、28参考、蛋白背景及原批注和锚点。全部28页逐张检查通过。\n',encoding='utf8')
print(json.dumps({'status':'complete','document':str(DOC),'document_sha256':digest,'pages':pages,'figure_sources':str(package)},ensure_ascii=False))
