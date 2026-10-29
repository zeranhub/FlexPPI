"""Promote only the reviewed manuscript, figure, native package and rendered pages."""
from pathlib import Path
from datetime import datetime,timezone,timedelta
import hashlib,json,zipfile,shutil
ROOT=Path(__file__).resolve().parent;BASE=ROOT.parent.parent;QA=ROOT/'manuscript_QA'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
read=lambda p:json.loads(Path(p).read_text(encoding='utf8'))
def write(p,v):Path(p).write_text(json.dumps(v,ensure_ascii=False,indent=2),encoding='utf8')
integration=read(QA/'integration_audit.json');doc=Path(integration['document']);digest=sha(doc)
assert digest==integration['document_sha256']
assert sha(integration['source'])==integration['source_sha256']
package=read(QA/'package_preservation_audit.json')
assert package['status']=='passed' and package['document_sha256']==digest
science=read(ROOT/'scientific_review/final_scientific_review.json')
assert science['status']=='passed' and science['document_sha256']==digest
assert science['figure4_sha256']==sha(integration['figure4'])
narrative=read(ROOT/'literature_review/final_narrative_review.json')
assert narrative['status']=='passed' and narrative['document_sha256']==digest
figure=read(ROOT/'figure4/visual_QA.json')
assert figure['status']=='passed' and figure['PNG_SHA256']==sha(integration['figure4'])
data=read(ROOT/'figure4/data_audit.json');assert data['status']=='passed'
assert data['PNG_SHA256']==sha(integration['figure4'])
inventory=read(ROOT/'case_source_inventory.json');assert inventory['status']=='passed'
for row in inventory['files']:
    assert sha(ROOT/row['archive_path'])==row['sha256']==sha(row['original_path'])
page_map=read(QA/'v49_page_map.json');pages=len(page_map)
tmp=Path(r'C:\Users\DELL\AppData\Local\Temp\flexppi_completion_20261010\render_v49')
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
fm=read(ROOT.parent/'kras_affinity_v48_20261019/CURRENT_FIGURE_MAP.json')
fm['4']={'source_folder':'figure4','absolute_source_folder':str(ROOT/'figure4'),
 'png':'fig4_mel5_joint.png','svg':'fig4_mel5_joint.svg','pdf':'fig4_mel5_joint.pdf','new_in_revision':'v49',
 'panels':['a alpha24beta17 structure, Q72-beta53/beta54 contact network and published SPR affinities',
           'b absolute model RMSF across all 45 positions in six CDRs',
           'c absolute model Hchi at all 23 valid pMHC interface positions',
           'd contact geometry occupancy in three existing MD replicas per state']}
for entry in fm.values():
    for ext in ['png','svg','pdf']:assert (Path(entry['absolute_source_folder'])/entry[ext]).is_file()
write(ROOT/'CURRENT_FIGURE_MAP.json',fm)
(ROOT/'修改说明_zh.md').write_text(f'''v49 MEL5–α24β17双动态机制实例，完成于2026-10-19

用MEL5与其亲和成熟衍生物α24β17替换原KRAS G12D案例，正文3.4和图4同步重写。补充MART-1、HLA-A*02:01呈递、天然MEL5与亲和成熟背景，以及原研究SPR和晶体结构证据。第28条参考文献改为Madura等2013年原研究，其他27条保留。

新图4包含真实4JFF结构与Q72–βPhe53/βGly54放大、六个完整CDR的45位模型RMSF、全部23个有效pMHC界面位点的模型Hχ、已有ITrace三重复轨迹的接触占据佐证。全部可见标签为英文，模型曲线使用绝对读数，未画MD曲线或将Hχ定向为RMSF的原因。共同编号与工程输入+1映射已核对。

机制论证围绕未改变类型的两个残基：HLA Q72模型Hχ由1.792降至1.375 nats，βGly54模型RMSF由1.118降至0.747 Å。具体原子接触与三重复占据支持亲和成熟重排接触网络、局部侧链与骨架受到协同约束的解释，并与原研究解离减慢相呼应。亲和力来自原文SPR；已有ITrace轨迹的补充接触分析与原研究实验分别注明。

摘要、引言、结构实例方法、讨论、结论同步修改。图1–3、5/S1、10张原生三线表、36个公式、批注和不相关正文保持。最终稿{pages}页，数据、文献、文档包及全部页面视觉核查通过。本轮无新增训练、推理或MD生产；已有六条轨迹的再分析来源和完整位点数据已归档。
''',encoding='utf8')
(ROOT/'README_zh.md').write_text('''v49源文件与核查记录。figure4为Word使用的PNG及可编辑SVG、PDF、真实结构渲染与绘图脚本；case_sources为原始预测、完整区域及已有轨迹再分析的轻量备份，原始大轨迹的路径及SHA保留于analysis_audit.json；literature_review核查原文事实，scientific_review核查数字与机制；manuscript_QA核查原生文档包及所有渲染页。build_v49.py对原v48仅局部编辑原生OOXML。CURRENT_FIGURE_MAP.json索引现用全部图源。原v48与旧案例归档保持不变。
''',encoding='utf8')
pointer=ROOT.parent/'completed_results_20261012/delivery_status.json'
snapshot=QA/'delivery_status_v48_snapshot.json'
if not snapshot.exists():shutil.copy2(pointer,snapshot)
prior=read(snapshot);assert prior['document_sha256']==integration['source_sha256']
files=sorted(p for p in ROOT.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.name not in ('v49_edit_manifest.json','delivery_status.json'))
manifest={'document':str(doc),'document_sha256':digest,'previous_document':integration['source'],
 'files':{p.relative_to(ROOT).as_posix():sha(p) for p in files}}
write(ROOT/'v49_edit_manifest.json',manifest)
archive=BASE/'paper_v49_20261019_MEL5图文与源文件.zip'
with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=5) as z:
    for p in files:z.write(p,p.relative_to(ROOT).as_posix())
    z.write(ROOT/'v49_edit_manifest.json','v49_edit_manifest.json')
with zipfile.ZipFile(archive) as z:
    assert z.testzip() is None
    for name,d in manifest['files'].items():assert hashlib.sha256(z.read(name)).hexdigest()==d
prior={k:v for k,v in prior.items() if not k.startswith('TCR_case_')}
prior.update(status='chinese_manuscript_v49_MEL5_joint_residue_mechanism_delivered',
 completed_beijing=datetime.now(timezone(timedelta(hours=8))).isoformat(),document=str(doc),document_sha256=digest,pages=pages,
 previous_document=integration['source'],previous_document_sha256=integration['source_sha256'],
 changes_this_round='Replace KRAS with MEL5 to alpha24beta17; model absolute RMSF and Hchi locate Q72-beta53/beta54 interface reorganization, with primary SPR/crystallography and existing ITrace contact-occupancy corroboration.',
 content_review='passed_primary_literature_model_data_scientific_review_and_native_package_preservation',
 render_and_visual_review=f'passed_all{pages}pages_individually_reviewed',
 document_edit_package=str(archive),document_edit_package_sha256=sha(archive),new_figure_package=str(archive),new_figure_package_sha256=sha(archive),
 new_experiment_note=str(ROOT/'修改说明_zh.md'),current_figure_map=str(ROOT/'CURRENT_FIGURE_MAP.json'),
 new_MD_production_this_round=0,new_model_training_this_round=0,new_model_inference_this_round=0,
 new_MD_analysis_this_round=6,new_MD_analysis_frames=6012,remaining_required_work=[],
 TCR_case_name='MEL5 affinity maturation to alpha24beta17',TCR_case_pdbs=['3hg1','4jff'],TCR_case_figure_number=4,
 TCR_case_training_ITrace_records=0,TCR_case_all_CDR_RMSF_positions=45,TCR_case_valid_pMHC_interface_Hchi_positions=23,
 TCR_case_MD_curves_displayed=False,TCR_case_MD_similarity_metrics_displayed=False,TCR_case_contact_occupancies_displayed=True,
 TCR_case_plot_readouts='absolute model RMSF and Hchi in complete predefined regions; local contact occupancy from existing MD',
 TCR_case_MD_evidence_role='local contact-network corroboration',TCR_case_MD_corroboration_trajectories=6,
 TCR_case_affinity_evidence='published SPR, MEL5 KD18 micromolar and alpha24beta17 KD0.6 nanomolar, same HLA-A*02:01/ELA pMHC',
 TCR_case_mechanism_evidence='published crystal contact reorganization and SPR dissociation; paired model Q72 Hchi/Gly54 RMSF and existing ITrace contact geometry',
 references=28,three_line_tables=10,figures=6)
write(ROOT/'delivery_status.json',prior);write(pointer,prior)
(ROOT.parent/'completion_status_zh.md').write_text(f'当前交付 {doc}，SHA256 {digest}，{pages}页。MEL5–α24β17替换KRAS，图4展示完整区域模型双读数及残基接触机制；摘要、方法、讨论、结论和文献同步。10原生三线表、36公式、28参考与批注完整；全部页面逐张检查通过。\n',encoding='utf8')
print(json.dumps({'status':'complete','document':str(doc),'figure4':integration['figure4'],'pages':pages,'source_archive':str(archive)},ensure_ascii=False))
