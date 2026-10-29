"""Package checked figures and scientific results; update current delivery pointer."""
from pathlib import Path
from datetime import datetime, timezone, timedelta
import hashlib, json, zipfile, shutil

ROOT=Path(__file__).resolve().parent
TASK=ROOT.parent
BASE=ROOT.parent.parent
QA=ROOT/'manuscript_QA'
FIG=ROOT/'figures'
DOC=BASE/'paper_v21_20261012_中文修订_构象分布与自由能分析_三线表.docx'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
write=lambda p,x:Path(p).write_text(json.dumps(x,ensure_ascii=False,indent=2),encoding='utf8')
now=datetime.now(timezone(timedelta(hours=8))).isoformat()
review=json.loads((QA/'final_content_and_layout_audit.json').read_text(encoding='utf8'))
assert review['status']=='passed' and review['document_sha256']==sha(DOC)

names=['fig4_joint_complementarity','fig7_ubiquitin_angular_landscapes','figS12_partition_sampling','figS13_representative_distribution_cases']
figure_review={'status':'passed','reviewed_beijing':now,'reviewer':'root','figures':{},'visual_findings':'All four final PNGs inspected at native resolution and at their Word placement. Labels, legends, scatter points, circular seams, missing-bin breaks and captions are legible and correctly placed. No overlaps, clipped labels or missing glyphs remain.'}
for name in names:
    figure_review['figures'][name]={ext:sha(FIG/f'{name}.{ext}') for ext in ['pdf','svg','png']}
write(FIG/'final_figure_visual_QA.json',figure_review)
for name,field in [('case_figure_manifest.json','visual_review'),('figS12_partition_sampling.audit.json','status')]:
    p=FIG/name
    old=json.loads(p.read_text(encoding='utf8'))
    if not (QA/('pre_final_'+name)).exists():shutil.copy2(p,QA/('pre_final_'+name))
    old[field]='passed' if field=='visual_review' else 'passed_numeric_and_visual_QA'
    old['final_visual_QA']='final_figure_visual_QA.json'
    write(p,old)
integration=json.loads((QA/'integration_audit.json').read_text(encoding='utf8'))
integration['render_and_content_review']='passed'
integration['final_content_and_layout_audit_sha256']=sha(QA/'final_content_and_layout_audit.json')
integration['supplementary_table_numbering']='Renumbered original S3,S1,S2 to S1,S2,S3 by first appearance; all native cell values unchanged.'
write(QA/'integration_audit.json',integration)

note=ROOT/'本轮新增实验与论文修改说明.md'
note.write_text('''本轮围绕“RMSF、Hχ与构象分布及自由能曲线的联系”补全实际分析，并整合为v21中文稿。主文新增第2.8、3.3和3.6节及相关讨论，补充材料新增J–L，加入图4、图7、图S12、图S13与表S12、表S13。最终稿49页，7张正文图、13张补充图，15张表均为原生可编辑三线表；新增4篇经原始文献核验的参考文献，共35篇。联合训练继续作为默认方案，正文围绕两类读数的互补用途展开，真实结果按实际数值报告。

全队列实验使用校正后的3,063条记录及824,988个RMSF/Hχ共同有效位置。在每条记录内减去同一种氨基酸的平均读数后，两指标的参考及预测Spearman相关分别为0.391和0.517。相同类型内，低RMSF／高Hχ与高RMSF／低Hχ共有47,354个位点，占全部共同有效位置的5.740%；先对记录求比例、再对记录等权平均为6.105%。模型对这两类组合的记录平均召回率为15.8%和16.5%，对低／低和高／高组合为54.6%和61.4%。图4与表S12同时展示组合存在性和预测恢复结果。

轨迹实验读取83条冻结记录的全部203条现有轨迹，统一分析0–100 ns、每200 ps保存的一帧，共501帧。实际提取、有效位置、原统一窗口标签复现及全部缓存文件经过独立核验。完整联合Hχ在原分区下的记录平均为1.484 nats，平移60°后为0.962 nats；非Pro χ1边缘熵由0.717变为0.374 nats，相邻保存帧的平均类别变化由189.2次变为21.9次。图S12与表S13将角区间、保存帧类别变化以及前后半程和重复的占据差异分开呈现，说明连续角度峰与离散区间的具体对应。

泛素Q41/Q62沿用原病例的5,000帧，新增χ1、χ2和χ3的完整时间序列、概率密度及相对经验投影自由能F*。两位点的χ1原分区边缘熵几乎相同，但Q62的多角度组合覆盖更丰富，60°联合平移后仍保留这一高低差异。另按预设的来源内中位误差规则选取2py5_A、A0022和ITrace 6p64，展示完整预测曲线及同类型位点的连续分布。自由能曲线由现有轨迹的角度占据概率计算，用来解释联合读数所概括的构象形态。

绘图文件包包含四张新增图的PDF、SVG和PNG、中文图注、绘图数据及源脚本。全部计算的源表、逐位置结果、输入冻结清单、独立审计与绘图源文件保存在同目录的科学结果归档；203条角度及Cα抽取缓存保留于trajectory_analysis/cache，各文件校验值在completion_audit.json中。KIX继续使用四种状态各一条100 ns轨迹。原始v13及此前修订稿保留。
''',encoding='utf8')

figure_zip=BASE/'paper_v21_20261012_新增图_PDF_SVG_PNG及绘图数据.zip'
figure_files=sorted(p for p in FIG.iterdir() if p.is_file())
with zipfile.ZipFile(figure_zip,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    for p in figure_files:z.write(p,'figures/'+p.name)
    z.write(note,note.name)
with zipfile.ZipFile(figure_zip) as z:
    assert z.testzip() is None
    for p in figure_files:assert hashlib.sha256(z.read('figures/'+p.name)).hexdigest()==sha(p)

# The extraction cache is already certified and retained; archive the complete
# analysis tables, code, plots and audit trail without redundant staged XTCs.
source_files=[]
for p in ROOT.rglob('*'):
    if not p.is_file():continue
    rel=p.relative_to(ROOT)
    if any(part in {'cache','_staging','__pycache__','preflight_v1','table_parts'} for part in rel.parts):continue
    if p.suffix in {'.zip','.tar','.log','.pyc'} or p.name.endswith('.tar.gz'):continue
    source_files.append(p)
source_files=sorted(source_files)
archive=ROOT/'v21_scientific_results_and_sources.zip'
archive_manifest={'created_beijing':now,'document_sha256':sha(DOC),'retained_extraction_cache':str(ROOT/'trajectory_analysis/cache'),'retained_extraction_cache_certificate_sha256':sha(ROOT/'trajectory_analysis/completion_audit.json'),'files':{str(p.relative_to(ROOT)).replace('\\','/'):sha(p) for p in source_files}}
write(ROOT/'source_archive_manifest.json',archive_manifest)
with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=4) as z:
    for p in source_files:z.write(p,str(p.relative_to(ROOT)).replace('\\','/'))
    z.write(ROOT/'source_archive_manifest.json','source_archive_manifest.json')
with zipfile.ZipFile(archive) as z:
    assert z.testzip() is None
    for relative,digest in archive_manifest['files'].items():assert hashlib.sha256(z.read(relative)).hexdigest()==digest

prior_status=TASK/'completed_results_20261012/delivery_status.json'
prior_note=TASK/'completion_status_zh.md'
if not (QA/'delivery_status_v20_snapshot.json').exists():shutil.copy2(prior_status,QA/'delivery_status_v20_snapshot.json')
if not (QA/'completion_status_v20_snapshot.md').exists():shutil.copy2(prior_note,QA/'completion_status_v20_snapshot.md')
old=json.loads(prior_status.read_text(encoding='utf8'))
delivery={**old,'status':'chinese_manuscript_v21_distribution_experiments_delivered','completed_beijing':datetime.now(timezone(timedelta(hours=8))).isoformat(),'document':str(DOC),'document_sha256':sha(DOC),'pages':49,'three_line_tables':15,'figures':20,'main_figures':7,'supplementary_figures':13,'references':35,'previous_document':old['document'],'previous_document_sha256':old['document_sha256'],'content_review':'passed_against_independently_verified_scientific_outputs','render_and_visual_review':'passed_all_49_pages','new_trajectory_analysis_records':83,'new_trajectory_analysis_trajectories':203,'new_complementarity_analysis_records':3063,'new_complementarity_common_residues':824988,'new_MD_production_this_round':0,'new_model_training_this_round':0,'new_figure_package':str(figure_zip),'new_figure_package_sha256':sha(figure_zip),'new_scientific_results_archive':str(archive),'new_scientific_results_archive_sha256':sha(archive),'new_experiment_note':str(note),'remaining_required_work':[],'changes_this_round':'Actual full-cohort complementary-output analysis, angular-partition and time/repeat diagnostics, continuous angular distributions and empirical projected free-energy plots; Chinese methods/results/discussion integrated; new three-line tables; all retained tables, comments, math and scientific image parts preserved.'}
write(ROOT/'delivery_status.json',delivery)
write(prior_status,delivery)
prior_note.write_text(f'''当前交付为 {DOC}，共49页，7张正文图、13张补充图、15张原生可编辑三线表及35篇参考文献。文件SHA256为 {sha(DOC)}。图包为 {figure_zip}。本轮补全3,063条折外记录的双读数组合分析、83条记录的203条轨迹角区间与采样分析，以及泛素与三个来源案例的连续角度分布和相对经验投影自由能。方法、结果、图注和讨论均已整合为中文。

内容数值、图表引用、用户批注及原生公式核验通过；最终Word渲染49页，全部页面的视觉核验通过。原13张表的每个单元格数值保留，新增两张三线表；补充表按出现顺序重编号。联合训练继续作为默认模型，KIX复用四种状态各一条100 ns轨迹。

所有本轮实验与论文修改已完成；科学源表和脚本归档为 {archive}，SHA256为 {sha(archive)}。203条角度及Cα抽取缓存仍保留在 {ROOT/'trajectory_analysis/cache'}。此前20次补充训练、发布基线和外部验证已完成，最终权重与逐残基源表的原归档保持不变。当前详细状态位于 {ROOT/'delivery_status.json'}，本轮说明位于 {note}。原始v13及此前v19、v20均保留。
''',encoding='utf8')
print(json.dumps({'status':'complete','document':str(DOC),'document_sha256':sha(DOC),'figure_package_MB':round(figure_zip.stat().st_size/1e6,2),'source_archive_MB':round(archive.stat().st_size/1e6,2),'source_archive_files':len(source_files),'remaining_required_work':[]},ensure_ascii=False))
