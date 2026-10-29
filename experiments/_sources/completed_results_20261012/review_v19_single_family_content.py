"""Independent read-only scientific/numeric DOCX review, no document edits."""
from pathlib import Path
import hashlib,json,zipfile,sys
from datetime import datetime,timezone
from docx import Document
import pandas as pd

ROOT=Path('E:/Project/FlexPPI/experiments/_sources')
OUT=ROOT/'completed_results_20261012'
DOC=ROOT.parent/'paper_v19_20261012_中文修订_补充实验完成_三线表.docx'
EXPECTED='354c969a523f1f481caf0bf18bf7123f8405ff9f0c52d20e69a2214f48403032'
OLD=Path('C:/Users/DELL/.codex/worktrees/4666/FlexPPI/paper/revision_20261002/figures/corrected_external_results_preparation/completed/20261003T163520Z_03a3b308/all_32_per_record_metrics.csv')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def require(condition,message):
    if not condition:
        raise ValueError(message)


def main():
    sys.stdout.reconfigure(encoding='utf-8')
    doc=Document(DOC)
    text=[paragraph.text for paragraph in doc.paragraphs]
    tables=[[[cell.text for cell in row.cells] for row in table.rows] for table in doc.tables]
    checks,errors=[],[]
    def check(name,actual,expected):
        okay=bool(actual==expected)
        item=dict(name=name,passed=okay,actual=actual,expected=expected)
        checks.append(item)
        if not okay:errors.append(item)
    source_hashes={}
    def csv(path):
        source_hashes[str(path)]=sha(path)
        return pd.read_csv(path)
    single=csv(OUT/'summary/single_target_overall_record_macro.csv').set_index(['variant','metric'])
    family=csv(OUT/'summary/family_overall_record_macro.csv').set_index(['variant','metric'])
    source=csv(OUT/'summary/family_source_record_macro.csv').set_index(['variant','dataset','metric'])
    external=csv(OUT/'summary/external_2430_record_macro.csv').set_index(['variant','metric'])
    original_baseline=csv(OUT/'baseline_original_overall_values_verified.csv').set_index('method')
    family_baseline=csv(OUT/'baseline_family_overall_values_verified.csv').set_index('method')
    bf=json.loads((OUT/'audit/external_figure_and_backflip_reuse.json').read_text())['backflip_reuse']
    def row(value):
        return [f'{float(value):.3f}' if pd.notna(value) else '不适用' for value in value]
    expected=[]
    for name,variant in [('完整模型','full'),('仅链内模型','intrachain'),('静态特征模型','static_sequence_complex'),('氨基酸类型均值','aa_only')]:
        a,b=external.loc[(variant,'rmsf')],external.loc[(variant,'chi_entropy')]
        expected.append([name]+row([a.mae,a.spearman,b.mae,b.spearman]))
    expected.append(['BackFlip2发布模型']+row([bf['record_macro_mae'],bf['record_macro_spearman']])+['无','无'])
    check('Table2 complete new2430 cells',tables[1][1:],expected)
    expected=[]
    for metric,label,variant,training in [('rmsf','RMSF','rmsf_only','单目标'),('rmsf','RMSF','full','联合目标'),
        ('chi_entropy','Hχ','chi_entropy_only','单目标'),('chi_entropy','Hχ','full','联合目标')]:
        value=single.loc[(variant,metric)]
        expected.append([label,training]+row([value.mae,value.rmse,value.spearman]))
    check('TableS9 recordmacro all cells',tables[10][1:],expected)
    expected=[]
    for dataset,label in [('all','整体'),('atlas','ATLAS'),('dynarepo','DynaRepo'),('itrace','ITrace')]:
        for variant,name in [('full','完整'),('intrachain','链内')]:
            a,b=(family.loc[(variant,'rmsf')],family.loc[(variant,'chi_entropy')]) if dataset=='all' else (source.loc[(variant,dataset,'rmsf')],source.loc[(variant,dataset,'chi_entropy')])
            expected.append([label,name]+row([a.mae,a.spearman,b.mae,b.spearman]))
    check('TableS10 recordmacro source/all cells',tables[11][1:],expected)
    expected=[]
    for method,name in [('flexppi','完整FlexPPI'),('intrachain','链内FlexPPI'),('planet_full','PLANET full'),
        ('planet_seq','PLANET seq'),('rmsf_net_pdb','RMSF-net PDB'),('anm','ANM'),('gnm','GNM')]:
        a,b=original_baseline.loc[method],family_baseline.loc[method]
        expected.append([name]+row([a.mae,a.spearman,b.mae,b.spearman,a.z_mae]))
    check('TableS11 original/family/relative all cells',tables[12][1:],expected)
    # New Table1 released rows also match the independent baseline review.
    expected=[]
    for method,name in [('planet_full','PLANET full 发布权重'),('planet_seq','PLANET seq 发布权重'),('rmsf_net_pdb','RMSF-net PDB-only')]:
        a=original_baseline.loc[method]
        expected.append([name]+row([a.mae,a.spearman])+['无','无'])
    check('Table1 new released rows',tables[0][-3:],expected)
    # Historical TableS12 is checked directly against the old32 source, not
    # against any new strict2430 predictions.
    old=csv(OLD)
    prior=old[old.protocol.eq('strict_single_seed7')].groupby(['variant','metric'])[['mae','spearman']].mean()
    old_names={'完整模型':'full','仅链内模型':'intrachain','静态特征模型':'matched_static','氨基酸类型均值':'matched_aa'}
    expected=[]
    for name,variant in old_names.items():
        a,b=prior.loc[(variant,'rmsf')],prior.loc[(variant,'chi_entropy')]
        expected.append([name]+row([a.mae,a.spearman,b.mae,b.spearman]))
    expected.append(['BackFlip2发布模型']+row([bf['record_macro_mae'],bf['record_macro_spearman']])+['无','无'])
    check('TableS12 old2431 historical cells',tables[13][1:],expected)
    # Previously audited KIX sources should stay intact in the new document.
    kix=csv(ROOT/'kix_existing_reuse_20261010/summary/state_summary.csv').set_index('state')
    expected=[]
    for state,name in [('apo','游离KIX'),('myb','cMyb'),('mll','MLL'),('both','双伙伴')]:
        a=kix.loc[state]
        expected.append([name]+row([a.RMSF_MAE_angstrom,a.RMSF_spearman])+[str(int(a.top20_overlap))+'/18']+row([a.L12_MD_mean_RMSF_angstrom,a.L12_fold1_mean_RMSF_angstrom]))
    check('TableS7 retained KIX cells',tables[8][1:],expected)
    kix_delta=csv(ROOT/'kix_existing_reuse_20261010/summary/l12_delta_summary.csv').set_index('contrast')
    expected=[]
    for contrast,label in [('myb-apo','cMyb−游离KIX'),('mll-apo','MLL−游离KIX'),('both-apo','双伙伴−游离KIX'),
        ('mll-myb','MLL−cMyb'),('both-myb','双伙伴−cMyb'),('both-mll','双伙伴−MLL')]:
        a=kix_delta.loc[contrast]
        expected.append([label,f'{a.MD_delta_RMSF_angstrom:+.3f}',f'{a.fold1_delta_RMSF_angstrom:+.3f}',f'{a.absolute_delta_error_angstrom:.3f}'])
    check('TableS8 retained KIX all six deltas',tables[9][1:],expected)
    # Whole original trajectory cert and figure itself are not rewritten.
    figure=OUT/'figures/fig4_strict_external_paired.png'
    figure_sha=sha(figure)
    with zipfile.ZipFile(DOC) as archive:
        embedded=[name for name in archive.namelist() if name.startswith('word/media/') and hashlib.sha256(archive.read(name)).hexdigest()==figure_sha]
    check('Fig4 exact new2430 source image embedded',len(embedded),1)
    stats=csv(OUT/'figures/fig4_strict_external_paired_statistics.csv')
    for item in stats.itertuples(index=False):
        value=external.loc[(item.variant,item.metric)]
        check('Fig4 stats '+item.variant+' '+item.metric,round(item.mean_mae,12),round(value.mae,12))
    main='\n'.join(text[:106])
    check('Main new protocol training2430 stated', '2,430' in main,True)
    check('Main no stale all-target external claim',any('两项目标的平均误差' in t and '严格外部' in t and '降低' in t for t in text[:106]),False)
    check('Joint contribution no longer untested',any('联合学习相对于独立单目标模型的贡献仍需单独检验' in t for t in text),False)
    check('ITrace only fold0 explained','仅在第0折测试' in text[40] and '该折训练不含ITrace' in text[40],True)
    check('Long MMseqs and short exact protocol explicit','至少20' in text[40] and '短' in text[40] and '精确' in text[40],True)
    check('New external reference no longer points to old S1','覆盖情况见图S1' in text[78],False)
    check('Family no fabricated source SD',any('±' in cell for row in tables[11] for cell in row),False)
    historical_mentions=[]
    for index,t in enumerate(text):
        if '632' in t or '2,431' in t or '2431' in t:
            historical_mentions.append(dict(paragraph=index,text=t,historical_explicit=('旧' in t or '历史' in t)))
    check('Every old632/2431 paragraph historical explicit',all(row['historical_explicit'] for row in historical_mentions),True)
    issues=[]
    if '对每条输入链执行' in text[40] and 'MMseqs2' in text[40]:
        issues.append(dict(severity='minor',paragraph=40,topic='short-chain search wording',
            text=text[40],recommendation='将“对每条输入链执行…MMseqs2序列搜索”改为“对长度至少20个残基的输入链执行…搜索，短肽采用精确序列连接”，与实际3770长唯一序列搜索和275条训练短链exact处理一致。'))
    if '覆盖情况见图S1' in text[78]:
        issues.append(dict(severity='major_cross_reference',paragraph=78,topic='new protocol points to historical figure',
            text=text[78],recommendation='新2430方案的序列与移除覆盖引用方法和补充H；图S1仅说明旧632/2431历史关联，勿用它作为新633/2430主方案的关系证明。'))
    if '0.004621' in text[208] or '0.002989' in text[208]:
        issues.append(dict(severity='optional_precision',paragraph=208,topic='display precision',
            recommendation='按用户统一3位小数偏好可显示MAE下降约0.005 Å/0.003 nats，或仅显示0.8%/1.8%；原精度留在源CSV。数值本身正确，不构成科学阻塞。'))
    report=dict(status='numerics_passed_with_text_clarifications' if not errors else 'failed_numeric_checks',
        document=str(DOC),sha256=sha(DOC),expected_input_sha256=EXPECTED,input_hash_matches=sha(DOC)==EXPECTED,
        checked_utc=datetime.now(timezone.utc).isoformat(),scope='Independent whole-text scientific review; table2/S7/S9/S10/S11/S12 source comparison; exact newFig4 media/statistics. Visual QA owned by root.',
        read_only=True,document_modified=False,pipeline_modified=False,inference_or_certificate_rerun=False,
        paragraph_count=len(text),table_count=len(tables),numeric_checks=checks,numeric_failure_count=len(errors),
        text_issues=issues,historical_old_protocol_occurrences=historical_mentions,
        scientific_conclusions=dict(joint_small_mean_gain_not_uniform=True,family_performance_declines=True,
            new_external_full_RMSF_worse_than_intra_static=True,new_external_full_Hchi_better_than_static_about1p5percent=True,
            ITrace_one_component_fold0_no_fivefold_SD=True,old_protocol_not_mixed_into_new_table=True),
        source_sha256=source_hashes,external_Fig4_media_sha256=figure_sha,external_Fig4_embedded_name=embedded)
    report['initial_document_sha256']='d5104be36fa74196b336d0518729daeb386a7da9f5a8e3d9a4910af96e2e3dab'
    report['required_text_clarifications_resolved_by_root']=['Short sequences handled by exact sequence connections; MMseqs applies to length>=20',
        'New633/2430 primary external analysis references methods and supplementH, while S1 is explicitly historical']
    report['required_unresolved_scientific_issues']=[] if not errors and all(issue['severity']=='optional_precision' for issue in issues) else errors+issues
    (OUT/'v19_single_family_content_review.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status=report['status'],sha256=report['sha256'],numeric_checks=len(checks),numeric_failures=len(errors),text_issues=issues),ensure_ascii=False,indent=2))
    if errors:raise SystemExit(1)


if __name__=='__main__':main()
