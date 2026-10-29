"""Finalize the audited KRAS manuscript and its local reproducibility record."""
from pathlib import Path
from datetime import datetime,timezone,timedelta
import hashlib,json,shutil,zipfile

ROOT=Path(__file__).resolve().parent
BASE=ROOT.parent.parent
QA=ROOT/'manuscript_QA'
DOC=BASE/'paper_v41_20261017_中文修订_KRAS新抗原动态响应_三线表.docx'
SEARCH=ROOT.parent/'itrace_case_search_only_20261017'
TMP=Path(r'C:\Users\DELL\AppData\Local\Temp\flexppi_completion_20261010')
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
read=lambda p:json.loads(Path(p).read_text(encoding='utf8'))
def write(p,x):Path(p).write_text(json.dumps(x,ensure_ascii=False,indent=2),encoding='utf8')
digest=sha(DOC)
assert digest=='58a15ba884180814dc562d276a3fb1c0e321b95a73cc302e0b0e68a4dafc0a68'
comparison=read(ROOT/'integration_audit/v41_package_comparison.json')
assert comparison['new_sha256']==digest and comparison['all_structural_checks_pass']
groups=[('visual_pages_1_10.json',range(1,11),'kras_figure6_build'),
        ('visual_pages_11_20.json',range(11,21),'kras_manuscript_text_review'),
        ('visual_pages_21_29.json',range(21,30),'kras_integration_audit')]
pages=[]
for filename,numbers,reviewer in groups:
    report=read(QA/filename)
    assert (report.get('status')=='passed' or report.get('overall_status')=='passed' or report.get('all_pass') is True),filename
    assert digest in json.dumps(report),filename
    for n in numbers:
        pages.append({'page':n,'status':'passed','reviewer':reviewer,
                      'sha256':sha(TMP/'render_v41'/f'page-{n}.png')})
assert len(read(QA/'v41_page_map.json'))==29
write(QA/'final_page_visual_review.json',{
 'status':'passed','document_sha256':digest,'pages':pages,
 'method':'All final150dpi original page PNGs individually opened by assigned reviewers; root independently inspected final pages1,7,8,16,17,18,19,21 and Figure6 standalone.',
 'root_additional_pages':[1,7,8,16,17,18,19,21],'required_fixes':[],
 'layout_fix':'Complete Figure6 caption kept together; final case text, figure and caption share page16.'})
integration=read(QA/'integration_audit.json')
assert integration['document_sha256']==digest
integration.update(status='complete',render_status='passed',pages=29,
                   visual_review='passed_all29pages_individually_reviewed',remaining_required_work=[])
write(QA/'integration_audit.json',integration)

data=ROOT/'case_sources';data.mkdir(exist_ok=True)
source_files=[
 SEARCH/'promising_cases_replica_validation/7ow5_7ow6_all_positions_three_replicas.csv',
 SEARCH/'KRAS_story_static_and_dynamic_checks.json',
 SEARCH/'case_story_static_check.py',
 SEARCH/'model_audit/verified_summary.json',
 SEARCH/'biology_first/sources/kras_g12d_2022_exa.txt',
 SEARCH/'replica_recheck/7ow5/7ow5_model_input.pdb',
 SEARCH/'replica_recheck/7ow6/7ow6_model_input.pdb']
source_manifest=[]
for source in source_files:
    assert source.is_file(),source
    target=data/source.name
    shutil.copy2(source,target)
    source_manifest.append({'original':str(source),'local_copy':target.relative_to(ROOT).as_posix(),'sha256':sha(target)})
write(data/'source_manifest.json',source_manifest)

note=ROOT/'修改说明_zh.md'
note.write_text('''v41 KRAS G12D实例替换，2026-10-17

将v40中的SB27/SB47及Q7侧链实例完整替换为KRAS野生型/G12D对照，正文3.5和图6编号保留。摘要、引言、方法2.7、讨论、结论及参考[28]同步替换，不再保留旧实例筛选、侧链变化或移除TCR推理的论述。

实例为同一JDIa41b1 TCR及HLA-A*11:01识别VVVGAGGVGK与VVVGADGVGK，仅肽第6位G→D，对应蛋白G12D。生物学背景引用Poole等Nature Communications2022，13:5333，doi:10.1038/s41467-022-32811-1。已报道KD为野生型3.0μM、G12D743±18pM；此亲和差异作为研究背景，论文模型分析直接验证主链动态响应。

主模型为论文当前完整家族第0折联合模型，训练ATLAS1,866条与DynaRepo564条，ITrace0条；检查点SHA256 ff2974406d5d1ee07e6dc7118e43e57c0acc9849cc8feba7f5290476fca15ec5。正文没有用较高分的旧strict模型替换主模型。MD使用ITrace每状态3条200ns既有轨迹，保留全部1001–1003保存帧，标签已在先前实例筛查中复算并核对原缓存。本轮不新增训练、MD生产或模型推理。

全10肽位点ΔRMSF Spearman0.806、MAE0.229Å；其余9个未替换位点Spearman0.733、MAE0.182Å，较零变化MAE0.303Å降低40.2%，8/9方向相符，7/9位置的三条G12D轨迹RMSF均低于三条野生型轨迹。全部10位及三重复均值差SEM完整绘图。正文描述中央A5/G7/V8响应及模型幅度低估。

初稿中考虑过链内模型对照，但该缓存对短肽的绝对值严重失准，未能清楚解释环境贡献，最终删除相关结果段及本例方法句；原始检查与数值保留于内部text_review审计，不在正文加入防御性说明。KRAS实例围绕完整主模型对静态输入下肽段动态变化的预测证据展开，不声称预测结合亲和力或重新发现原研究机制。

原文档的其余7幅可见图、10张原生三线表、36项原生数学、所有原批注内容及锚点均逐项相同；包内仅document.xml与原图6的image41.png改变。29页以隐藏的独立Word实例只读导出，再由文档技能render_docx.py生成150dpi逐页PNG，全部逐张检查通过。3.5三段正文、英文图6及完整图注同在第16页。
''',encoding='utf8')
(ROOT/'README_zh.md').write_text('''v41 KRAS G12D实例复现材料

主文档位于本目录上两级，文件名paper_v41_20261017_中文修订_KRAS新抗原动态响应_三线表.docx。该版本替换v40中的旧TCR实例，正文3.5、图6；泛素与自由能仍为3.6和图7。交付文稿为中文，图中可见文字为英文，所有10表仍为原生三线表。

figure6包含最终PNG/PDF/SVG、逐点预测与三重复MD源数据、统计定义、原始结构渲染、PyMOL会话和脚本。主图只用当前Full模型。case_sources保留所用结构输入、模型划分核验、全残基参考表、静态结构检查及原始文献摘录，source_manifest逐文件记录来源和SHA256。原始XTC轨迹仍在ITrace数据目录，未重复打包。

build_v41.py仅修改原稿包内document.xml和原图6的PNG，其余部件保持字节相同。text_review包含科学叙述与私有控制有效性核验；integration_audit含原包/新包逐部件比较；manuscript_QA记录29页逐张视觉检查和文稿SHA256。export_v41_word.ps1和render_v41.py为实际渲染路径。

CURRENT_FIGURE_MAP.json指向本稿全部图源。图1–5、图7及S1沿用v39图源，本目录只重建图6。所有分析、图源和文稿均可追溯到模型及既有MD参考；本次替换不新增MD、训练或推理。
''',encoding='utf8')
mapping=read(ROOT.parent/'tcr_case_v40_20261017/CURRENT_FIGURE_MAP.json')
mapping['6']={'source_folder':'figure6','png':'fig6_kras_case.png','pdf':'fig6_kras_case.pdf',
              'svg':'fig6_kras_case.svg','new_in_revision':'v41','absolute_source_folder':str(ROOT/'figure6')}
write(ROOT/'CURRENT_FIGURE_MAP.json',mapping)
pointer=ROOT.parent/'completed_results_20261012/delivery_status.json'
snapshot=QA/'delivery_status_v40_snapshot.json'
if not snapshot.exists():shutil.copy2(pointer,snapshot)
prior=read(snapshot)
assert prior['document_sha256']==integration['source_sha256']
# Superseded case-specific fields must not identify the previous case as current.
for k in list(prior):
    if k.startswith('TCR_case_') or k in {'new_figure_package','new_figure_package_sha256'}:prior.pop(k)
files=sorted(p for p in ROOT.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.name not in {'delivery_status.json','v41_package_manifest.json'})
manifest={'document':str(DOC),'document_sha256':digest,'previous_document':integration['source'],
          'files':{p.relative_to(ROOT).as_posix():sha(p) for p in files}}
write(ROOT/'v41_package_manifest.json',manifest)
package=BASE/'paper_v41_20261017_KRAS实例图与实验源数据.zip'
with zipfile.ZipFile(package,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    for p in files:z.write(p,p.relative_to(ROOT).as_posix())
    z.write(ROOT/'v41_package_manifest.json','v41_package_manifest.json')
with zipfile.ZipFile(package) as z:
    assert z.testzip() is None
    for filename,digest_file in manifest['files'].items():assert hashlib.sha256(z.read(filename)).hexdigest()==digest_file
prior.update(status='chinese_manuscript_v41_KRAS_G12D_case_delivered',
 completed_beijing=datetime.now(timezone(timedelta(hours=8))).isoformat(),
 previous_document=integration['source'],previous_document_sha256=integration['source_sha256'],
 document=str(DOC),document_sha256=digest,pages=29,
 content_review='passed_KRAS_source_numerical_scope_and_package_audit',
 render_and_visual_review='passed_all29pages_individually_reviewed',
 changes_this_round='Replaced SB27/SB47 and Q7 with KRAS G12D WT/mutant dynamic-response case; main3.5 and EnglishFigure6; abstract, methods, discussion, conclusion and reference28 synchronized.',
 new_figure_package=str(package),new_figure_package_sha256=sha(package),new_experiment_note=str(note),
 TCR_case_figure_number=6,TCR_case_pdbs=['7ow5','7ow6'],TCR_case_name='KRAS G12D',
 TCR_case_delta_RMSF_all10_rho=.806060606060606,TCR_case_delta_RMSF_all10_MAE_A=.22888186163487184,
 TCR_case_delta_RMSF_unchanged9_rho=.7333333333333333,TCR_case_delta_RMSF_unchanged9_MAE_A=.18150934934249716,
 TCR_case_unchanged9_improvement_over_zero=.4016165813029179,
 TCR_case_MD_reference_trajectories=6,TCR_case_training_ITrace_records=0,
 new_MD_production_this_round=0,new_model_training_this_round=0,new_model_inference_this_round=0,
 new_MD_analysis_this_round=0,remaining_required_work=[])
write(ROOT/'delivery_status.json',prior)
write(pointer,prior)
(ROOT.parent/'completion_status_zh.md').write_text(f'当前交付 {DOC}，SHA256 {digest}，29页。正文3.5与英文图6已替换为KRAS G12D；同一TCR/HLA的十肽G6D对照，全10位与未替换9位MD/预测比较。无新增训练、MD生产或推理。其余7幅可见图、10张原生三线表、36项数学和原批注/锚点全部保持。29页逐张视觉核查通过。源材料包 {package}。\n',encoding='utf8')
print(json.dumps({'status':'complete','document':str(DOC),'document_sha256':digest,'pages':29,'source_package':str(package)},ensure_ascii=False))
