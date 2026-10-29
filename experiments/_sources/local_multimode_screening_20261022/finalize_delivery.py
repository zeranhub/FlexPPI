"""Requirement-level artifact checks and recoverable full result archive."""
from pathlib import Path
from datetime import datetime,timezone
import json,csv,hashlib,zipfile,os

ROOT=Path(__file__).resolve().parent
BASE=ROOT.parent
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p,x):p.write_text(json.dumps(x,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf8')

def main():
    inference=json.loads((ROOT/'predictions/full_pool/progress.json').read_text())
    assert inference['status']=='passed' and inference['completed_this_invocation']==1211
    assert not (ROOT/'predictions/full_pool/.worker.lock').exists()
    qa=json.loads((ROOT/'result_figures/figure_QA.json').read_text())
    with (ROOT/'result_figures/figure_source_values.csv').open(encoding='utf-8-sig',newline='') as h:
        values=list(csv.DictReader(h))
    assert len(values)==588
    for fig in ['three_modes_yield','feedback_all_contexts','profile_specificity_controls']:
        for ext in ['pdf','png','svg']:assert (ROOT/'result_figures'/f'{fig}.{ext}').stat().st_size>1000
    qa.update(status='passed_personal_visual_and_numeric_QA',personal_visual_review_all3_final_PNGs=True,
        fonts_embedded_all3_PDFs=True,PDF_fonts='ArialMT/Arial-BoldMT Type0 with embedded TrueType',
        vector_PDF_pages_each=1,PDF_widths_mm=[179.8,180.5,177.0],all_labels_English=True,
        no_clipped_or_overlapping_labels_observed=True,source_value_rows=588,
        figures_show_strong_sequence_controls=True,post_initial_result_baseline_extension_disclosed_in_report=True)
    dump(ROOT/'result_figures/figure_QA.json',qa)
    audit=json.loads((ROOT/'evaluation/independent_result_audit.json').read_text());assert audit['total_metric_rows_recomputed']==104413
    results=[]
    for setting in ['primary','top_decile','proteinmpnn002','three_dose_AUC','controls','stronger_baselines']:
        p=ROOT/'evaluation'/('results_'+setting)/'execution_audit.json'
        a=json.loads(p.read_text());assert a['status']=='completed' and a['n_panels']==7
        results.append(dict(setting=setting,audit_sha256=sha(p)))
    mainpaper=BASE.parent/'paper_v49_20261019_中文修订_MEL5双动态机制实例_三线表.docx'
    paper_sha=sha(mainpaper);assert paper_sha=='5f79ddf7d657dc9b782c9fb432e373a1e7c170c1c88d019e528fbf73c366c222'
    cache_checks=0
    for p in (ROOT/'predictions/full_pool').glob('*/candidate_audit.json'):
        a=json.loads(p.read_text());assert a['status']=='passed'
        c=a['geometry_cache_validation']
        if c is not None:assert c['all_bit_identical'];cache_checks+=1
        for name,digest in a['output_SHA256'].items():assert sha(p.parent/name)==digest
    assert cache_checks==14
    deliverables=['本地多模式筛选实验_完整结果_20261022.md','三种模式_全部结果摘要.csv',
        'evaluation/candidate_rankings_seed20261022.csv','evaluation/independent_result_audit.json',
        'evaluation/decision_summary.json','result_figures/figure_QA.json']
    status=dict(status='completed',timestamp_utc=datetime.now(timezone.utc).isoformat(),
        requested_modes_tested=['parent-only site nomination','common-pool one-shot prioritization','batch feedback prioritization'],
        complete_antigen_variants=1197,receptor_contexts=7,antigen_groups=3,
        all_forward_inputs=1211,cpu_only=True,neural_training=False,MD=False,
        all14_template_cache_checks_bit_identical=True,all1211_output_commits_SHA_verified=True,
        all6_result_sets_completed=results,metric_rows_independently_recomputed=104413,
        manuscript_unchanged_sha256=paper_sha,historical_2431_model_not_manuscript2430=True,
        result='No consistent extra selection yield over strong sequence baselines; limited H-only feedback signal in TCR6, not validated directed engineering.',
        deliverables={p:sha(ROOT/p) for p in deliverables},
        DMS_receptor_branch='Not used: full matrix/construct mapping unavailable; antigen-site mode1 covered by full9x19 peptide saturation data.')
    dump(ROOT/'delivery_status.json',status)
    # Sources and all numeric outputs are retained, including unfavorable results.
    paths=[p for p in ROOT.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.suffix not in {'.zip','.tmp'}]
    archive=ROOT/'本地多模式实验_完整结果与源表_20261022.zip'
    temporary=archive.with_suffix('.zip.tmp')
    extras=[(BASE/'local_antigen_screening_20261021/model/epoch_040.pt','local_antigen_screening_20261021/model/epoch_040.pt'),
        (BASE/'local_antigen_screening_20261021/local_runtime/local_cpu_inference.py','local_antigen_screening_20261021/local_runtime/local_cpu_inference.py')]
    src=Path(r'E:\Project\FlexPPI\src\flexppi')
    extras.extend((p,'project_src/flexppi/'+p.relative_to(src).as_posix()) for p in src.rglob('*.py'))
    with zipfile.ZipFile(temporary,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=5) as z:
        for i,p in enumerate(sorted(paths),1):
            z.write(p,ROOT.name+'/'+p.relative_to(ROOT).as_posix())
            if i%1000==0:print(json.dumps(dict(packaged_files=i,total=len(paths))),flush=True)
        for p,name in extras:z.write(p,name)
        z.writestr('README_结果包.txt','本包包含三种模式的完整本地测试、全部1197候选及逐残基预测、全部阴性结果、6组主分析/敏感性/对照结果、图和源表。\n入口为local_multimode_screening_20261022/本地多模式筛选实验_完整结果_20261022.md。\n额外保留了实际使用的历史模型、复用加载器和项目源代码快照。运行脚本使用原工作目录中的绝对路径；迁移机器时需调整路径，保留已归档结果不覆盖。\n这不是已验证的抗原改造方案，也未修改论文稿。\n')
    os.replace(temporary,archive)
    with zipfile.ZipFile(archive) as z:
        assert z.testzip() is None
        for relative,digest in status['deliverables'].items():assert hashlib.sha256(z.read(ROOT.name+'/'+relative)).hexdigest()==digest
        assert hashlib.sha256(z.read('local_antigen_screening_20261021/model/epoch_040.pt')).hexdigest()=='b3eed060919df15f075329bc892f959d9355232be5feb561417edd49ea82b39f'
        n=len(z.namelist())
    package=dict(status='passed_full_archive_CRC_and_key_SHA_verification',path=str(archive),bytes=archive.stat().st_size,
        files=n,sha256=sha(archive),all_reported_artifacts_verified=True)
    dump(ROOT/'package_audit.json',package)
    print(json.dumps(package,ensure_ascii=False),flush=True)

if __name__=='__main__':main()
