"""Publish only the final reviewed manuscript and the reproducible figure masters."""
from pathlib import Path
from datetime import datetime,timezone,timedelta
import hashlib,json,zipfile,shutil
ROOT=Path(__file__).resolve().parent;BASE=ROOT.parent.parent;QA=ROOT/'manuscript_QA'
assert not (ROOT/'USER_REVIEW_20261019.md').exists(), 'Visual direction rejected by user; do not promote this draft.'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
read=lambda p:json.loads(Path(p).read_text(encoding='utf8'))
def write(p,v):Path(p).write_text(json.dumps(v,ensure_ascii=False,indent=2),encoding='utf8')
integration=read(QA/'integration_audit.json');doc=Path(integration['document']);digest=sha(doc)
assert digest==integration['document_sha256']
package_audit=read(QA/'package_preservation_audit.json')
assert package_audit['status']=='passed' and package_audit['document_sha256']==digest
science=read(ROOT/'audit/final_scientific_review.json')
assert science['status']=='passed' and science['document_sha256']==digest
figure=Path(integration['figure1']);assert science['figure1_sha256']==sha(figure)
review=read(ROOT/'audit/final_visual_review.json')
assert review['status']=='passed' and review['figure1_sha256']==sha(figure)
render=read(ROOT/'figure1/render_audit.json')
for ext in ['png','pdf','svg']:
    assert render['output_sha256'][ext]==sha(ROOT/f'figure1/fig1_framework.{ext}')
assert render['min_word_font_pt']>=8 and not render['out_of_canvas_text']
render.update(status='passed',independent_visual_review='../audit/final_visual_review.json',
 independent_scientific_review='../audit/final_scientific_review.json')
write(ROOT/'figure1/render_audit.json',render)
write(ROOT/'figure1/scientific_manifest.json',{
 'status':'passed','document_sha256':digest,'figure1_sha256':sha(figure),
 'source_hashes':render['source_hashes'],'output_sha256':render['output_sha256'],
 'drawing_script_sha256':sha(ROOT/'figure1/draw_framework.py'),
 'caption_sha256':sha(ROOT/'figure1/caption_zh.md'),
 'code_and_config_audit':'../audit/method_content.json',
 'protected_predictions_audit':'../audit/protected_840_data.json',
 'selected_atomic_frames_audit':'../audit/selected_frame_audit.json',
 'new_training':0,'new_inference':0,'new_MD':0})
page_map=read(QA/'v46_page_map.json');pages=len(page_map)
tmp=Path(r'C:\Users\DELL\AppData\Local\Temp\flexppi_completion_20261010\render_v46')
reviews=[]
for p in sorted(QA.glob('visual_pages_*.json')):
    item=read(p);assert item['status']=='passed' and item['document_sha256']==digest
    for row in item['pages']:
        assert row['status']=='passed'
        n=row['page'];reviews.append({'page':n,'status':'passed','reviewer':item['reviewer'],'sha256':sha(tmp/f'page-{n}.png')})
reviews.sort(key=lambda r:r['page'])
assert [r['page'] for r in reviews]==list(range(1,pages+1))
write(QA/'final_page_visual_review.json',{'status':'passed','document_sha256':digest,'pages':reviews,'required_fixes':[]})
integration.update(status='complete',pages=pages,render_status='passed',remaining_required_work=[])
write(QA/'integration_audit.json',integration)
fm=read(ROOT.parent/'dense_figure_v45_20261018/CURRENT_FIGURE_MAP.json')
fm['1']={'source_folder':'figure1','absolute_source_folder':str(ROOT/'figure1'),
 'png':'fig1_framework.png','svg':'fig1_framework.svg','pdf':'fig1_framework.pdf','new_in_revision':'v46',
 'panels':['a model overview','b residue local geometry','c fixed-frame invariant point attention','d MD supervision']}
for entry in fm.values():
    for ext in ['png','svg','pdf']:assert (Path(entry['absolute_source_folder'])/entry[ext]).is_file()
write(ROOT/'CURRENT_FIGURE_MAP.json',fm)
(ROOT/'修改说明_zh.md').write_text(f'''v46 模型框架图重绘，完成于2026-10-19

图1沿用静态结构输入、几何编码、联合动态读出及MD监督的整体叙事，重画各组件。b子图用主链原子与局部坐标的放大关系说明几何，不再将多组坐标轴和残基名拥挤地叠在蛋白带状结构上。其余子图统一字体、语义配色、线条、箭头和间距，以明确输入、几何注意力、模型输出与监督标签的关系。图内全部可见文字为英文，中文图注同步更新。

源码核对确认：模型使用一个共享StateHead中间网络，最终为RMSF和Hχ分别输出均值及对数方差通道。引言和方法中三处读出层措辞已据此澄清，未改变模型或训练。MD标签观测方差与模型残差方差在损失中相加；局部坐标在编码中保持固定。原4g9f示例观测及预测来源保留，细节见figure1源文件和audit。

其余正文、结论、图2–5/S1、10张原生三线表、36个公式、28条文献及批注均保留。最终稿{pages}页，经文档包完整性、科学内容、图形及逐页视觉检查。没有新增训练、推理或MD。
''',encoding='utf8')
(ROOT/'README_zh.md').write_text('''v46 图1源文件与核对记录

figure1/fig1_framework.svg和.pdf为矢量母版，.png为Word中实际使用的高分辨率版本；同目录含可复现脚本、素材来源和中文图注。CURRENT_FIGURE_MAP.json索引全文现用图源。audit包含实现核对、期刊示例设计参考及最终独立科学/视觉审查；manuscript_QA记录文档包和最终全部页面检查。build_v46.py为局部OOXML替换脚本，不会重建无关文档部件。

原始v45文稿和旧图1保持不变。最新DOCX在上两级目录，文件名paper_v46_20261019_中文修订_模型框架图重绘_三线表.docx。没有新增模型训练、推理或MD。
''',encoding='utf8')
pointer=ROOT.parent/'completed_results_20261012/delivery_status.json'
snapshot=QA/'delivery_status_v45_snapshot.json'
if not snapshot.exists():shutil.copy2(pointer,snapshot)
prior=read(snapshot);assert prior['document_sha256']==integration['source_sha256']
files=sorted(p for p in ROOT.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.name not in ('v46_edit_manifest.json','delivery_status.json'))
manifest={'document':str(doc),'document_sha256':digest,'previous_document':integration['source'],
 'files':{p.relative_to(ROOT).as_posix():sha(p) for p in files}}
write(ROOT/'v46_edit_manifest.json',manifest)
archive=BASE/'paper_v46_20261019_模型框架图与源文件.zip'
with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=5) as z:
    for p in files:z.write(p,p.relative_to(ROOT).as_posix())
    z.write(ROOT/'v46_edit_manifest.json','v46_edit_manifest.json')
with zipfile.ZipFile(archive) as z:
    assert z.testzip() is None
    for name,d in manifest['files'].items():assert hashlib.sha256(z.read(name)).hexdigest()==d
prior.update(status='chinese_manuscript_v46_redrawn_framework_delivered',
 completed_beijing=datetime.now(timezone(timedelta(hours=8))).isoformat(),document=str(doc),document_sha256=digest,pages=pages,
 previous_document=integration['source'],previous_document_sha256=integration['source_sha256'],
 changes_this_round='Redrew Figure1 components, especially local-frame geometry; unified typography, semantic colors and attention/training flows; clarified shared readout implementation in three prose locations and synchronized caption.',
 content_review='passed_source_code_scientific_review_and_native_package_preservation',render_and_visual_review=f'passed_all{pages}pages_individually_reviewed',
 document_edit_package=str(archive),document_edit_package_sha256=sha(archive),new_figure_package=str(archive),new_figure_package_sha256=sha(archive),
 new_experiment_note=str(ROOT/'修改说明_zh.md'),current_figure_map=str(ROOT/'CURRENT_FIGURE_MAP.json'),
 new_MD_production_this_round=0,new_model_training_this_round=0,new_model_inference_this_round=0,new_MD_analysis_this_round=0,remaining_required_work=[])
write(ROOT/'delivery_status.json',prior);write(pointer,prior)
(ROOT.parent/'completion_status_zh.md').write_text(f'当前交付 {doc}，SHA256 {digest}，{pages}页。图1各组成重绘，b局部几何重点重构，澄清共享读出层并同步图注。图2及其余科学内容保留；10张原生三线表、36公式、28参考与批注完整；全部页面逐张检查通过。\n',encoding='utf8')
print(json.dumps({'status':'complete','document':str(doc),'figure1':str(figure),'pages':pages,'source_archive':str(archive)},ensure_ascii=False))
