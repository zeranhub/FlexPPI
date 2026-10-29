"""Promote only the exact manuscript, figure and pages that passed review."""
from pathlib import Path
from datetime import datetime,timezone,timedelta
import hashlib,json,zipfile,shutil
ROOT=Path(__file__).resolve().parent;BASE=ROOT.parent.parent;QA=ROOT/'manuscript_QA'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
read=lambda p:json.loads(Path(p).read_text(encoding='utf8'))
def write(p,v):Path(p).write_text(json.dumps(v,ensure_ascii=False,indent=2),encoding='utf8')
integration=read(QA/'integration_audit.json');doc=Path(integration['document']);digest=sha(doc)
assert digest==integration['document_sha256']
for name in ['package_preservation_audit.json']:
    report=read(QA/name);assert report['status']=='passed' and report['document_sha256']==digest
science=read(ROOT/'audit/final_scientific_review.json')
assert science['status']=='passed' and science['document_sha256']==digest
assert science['figure4_sha256']==sha(integration['figure4'])
narrative=read(ROOT/'literature/final_narrative_review.json')
assert narrative['status']=='passed' and narrative['document_sha256']==digest
figure=read(ROOT/'figure4/visual_QA.json')
assert figure['status']=='passed' and figure['PNG_SHA256']==sha(integration['figure4'])
data=read(ROOT/'figure4/data_audit.json');assert data['status']=='passed'
for ext in ['png','pdf','svg']:
    assert data['output_SHA256']['fig4_kras_affinity.'+ext]==sha(ROOT/('figure4/fig4_kras_affinity.'+ext))
page_map=read(QA/'v48_page_map.json');pages=len(page_map)
tmp=Path(r'C:\Users\DELL\AppData\Local\Temp\flexppi_completion_20261010\render_v48')
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
fm=read(ROOT.parent/'framework_v47_20261019/CURRENT_FIGURE_MAP.json')
fm['4']={'source_folder':'figure4','absolute_source_folder':str(ROOT/'figure4'),
 'png':'fig4_kras_affinity.png','svg':'fig4_kras_affinity.svg','pdf':'fig4_kras_affinity.pdf','new_in_revision':'v48',
 'panels':['a WT and G12D input structures and published SPR affinities','b model absolute RMSF at all ten peptide positions']}
for entry in fm.values():
    for ext in ['png','svg','pdf']:assert (Path(entry['absolute_source_folder'])/entry[ext]).is_file()
write(ROOT/'CURRENT_FIGURE_MAP.json',fm)
(ROOT/'修改说明_zh.md').write_text(f'''v48 KRAS亲和力机制与绝对RMSF，完成于2026-10-19

图4b改为WT与G12D两条模型绝对RMSF曲线，展示全部10个肽位置的原始预测，去除MD曲线、误差棒及相关/MAE等比较量。a保留输入结构、序列、叠合及文献SPR亲和力；字体调整至Word幅宽下至少8 pt。英文图中文字和中文图注同步。

摘要、引言、方法、3.4、讨论和结论同步改写。案例从模型预测的区域响应出发：N端近似重合，G6/D6及邻近未替换G7、A5、V8波动下降。与原研究SPR、晶体结构及25×4 ns MD/MM-PBSA发现的D6–HLA相互作用和静电/焓机制对应，提供逐位置动态支持。ITrace三条200 ns重复轨迹仅作补充佐证。删除本案例“低估变化幅度”、零变化基线和以MD相似性为主的论述。

图2–3、图5/S1、10张原生三线表、36个公式、28条文献、批注和其余文字保留。最终稿{pages}页，经过数据、独立科学审查、文档包保存与全部页面视觉检查。没有新增训练、推理或MD。Hχ是否扩展此案例另见只读审查记录，不自动向正文添加未经筛选的解释。
''',encoding='utf8')
(ROOT/'README_zh.md').write_text('''v48源文件和核查记录。figure4目录含Word实际使用PNG、PDF/SVG导出、绘图脚本、全部20个绘图值及旧图源档案；literature核实原研究实验与机制；audit核实原预测和结构近接触及最终科学结论；manuscript_QA核查原生文档包与最终28页。build_v48.py仅局部修改原生OOXML，不重建其他文档部件。CURRENT_FIGURE_MAP.json索引现用全部图源。原始v47稿和旧图不改动。
''',encoding='utf8')
pointer=ROOT.parent/'completed_results_20261012/delivery_status.json'
snapshot=QA/'delivery_status_v47_snapshot.json'
if not snapshot.exists():shutil.copy2(pointer,snapshot)
prior=read(snapshot);assert prior['document_sha256']==integration['source_sha256']
files=sorted(p for p in ROOT.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.name not in ('v48_edit_manifest.json','delivery_status.json'))
manifest={'document':str(doc),'document_sha256':digest,'previous_document':integration['source'],
 'files':{p.relative_to(ROOT).as_posix():sha(p) for p in files}}
write(ROOT/'v48_edit_manifest.json',manifest)
archive=BASE/'paper_v48_20261019_KRAS图文与源文件.zip'
with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=5) as z:
    for p in files:z.write(p,p.relative_to(ROOT).as_posix())
    z.write(ROOT/'v48_edit_manifest.json','v48_edit_manifest.json')
with zipfile.ZipFile(archive) as z:
    assert z.testzip() is None
    for name,d in manifest['files'].items():assert hashlib.sha256(z.read(name)).hexdigest()==d
prior.update(status='chinese_manuscript_v48_KRAS_absolute_RMSF_and_affinity_mechanism_delivered',
 completed_beijing=datetime.now(timezone(timedelta(hours=8))).isoformat(),document=str(doc),document_sha256=digest,pages=pages,
 previous_document=integration['source'],previous_document_sha256=integration['source_sha256'],
 changes_this_round='Figure 4b uses two absolute model RMSF curves over all ten peptide positions; KRAS case centers on local model response and primary SPR/crystallography/MD evidence for HLA-mediated antigen recognition, with MD as corroboration.',
 content_review='passed_primary_literature_model_data_scientific_review_and_native_package_preservation',render_and_visual_review=f'passed_all{pages}pages_individually_reviewed',
 document_edit_package=str(archive),document_edit_package_sha256=sha(archive),new_figure_package=str(archive),new_figure_package_sha256=sha(archive),
 new_experiment_note=str(ROOT/'修改说明_zh.md'),current_figure_map=str(ROOT/'CURRENT_FIGURE_MAP.json'),
 new_MD_production_this_round=0,new_model_training_this_round=0,new_model_inference_this_round=0,new_MD_analysis_this_round=0,remaining_required_work=[])
prior.update(TCR_case_MD_curves_displayed=False,TCR_case_MD_similarity_metrics_displayed=False,
 TCR_case_plot_readout='absolute model RMSF for WT and G12D at all 10 peptide positions',
 TCR_case_MD_evidence_role='supplementary corroboration',TCR_case_MD_corroboration_trajectories=6,
 TCR_case_affinity_evidence='published SPR, same JDIa41b1 construct',
 TCR_case_mechanism_evidence='published X-ray HLA network, temperature-dependent SPR and MD/MM-PBSA; model local backbone dynamics')
write(ROOT/'delivery_status.json',prior);write(pointer,prior)
(ROOT.parent/'completion_status_zh.md').write_text(f'当前交付 {doc}，SHA256 {digest}，{pages}页。图4b改为WT/G12D两条绝对预测RMSF曲线，KRAS实例围绕模型局部动态响应与SPR/晶体结构/MD机制改写。其余图及10张原生三线表、36公式、28参考与批注完整；全部页面逐张检查通过。\n',encoding='utf8')
print(json.dumps({'status':'complete','document':str(doc),'figure4':integration['figure4'],'pages':pages,'source_archive':str(archive)},ensure_ascii=False))
