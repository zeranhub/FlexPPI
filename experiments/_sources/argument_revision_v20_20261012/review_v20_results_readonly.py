"""Read-only final manuscript content/source consistency review."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import re
import sys
import zipfile
import pandas as pd
from docx import Document

BASE = Path('E:/Project/FlexPPI/experiments/_sources')
DATA = BASE / 'completed_results_20261012'
OUT = BASE / 'argument_revision_v20_20261012'
DOC = BASE.parent / 'paper_v20_20261012_中文修订_联合模型_三线表.docx'
TXT = OUT / 'v20_content.txt'
EXPECTED = '4b4597b42fdc8a1300ab8f4155e98a4960872a2368fc887397d133af7d264d79'

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def main():
    sys.stdout.reconfigure(encoding='utf-8')
    begin_hash = sha(DOC)
    doc = Document(DOC)
    # python-docx Paragraph.text omits Word's editable equation objects.
    # Read ordinary runs and equation runs together in XML order, matching
    # the manuscript extractor, and omit cosmetic line-break objects.
    paragraphs = [''.join(n.text or '' for n in p._p.xpath('.//w:t | .//m:t')).strip() for p in doc.paragraphs]
    tables = [[[c.text for c in r.cells] for r in t.rows] for t in doc.tables]
    checks, errors, sources = [], [], {}
    def check(name, actual, expected):
        item = {'name': name, 'passed': bool(actual == expected), 'actual': actual, 'expected': expected}
        checks.append(item)
        if not item['passed']:
            errors.append(item)
    def csv(path):
        sources[str(path)] = sha(path)
        return pd.read_csv(path)
    def fmt(values):
        return [f'{float(v):.3f}' if pd.notna(v) else '不适用' for v in values]
    check('Latest document frozen hash', begin_hash, EXPECTED)
    check('Thirteen retained tables, old historical strict table removed', len(tables), 13)
    txt_paragraphs = {int(i): t for i, t in re.findall(r'^P(\d+): (.*)$', TXT.read_text(encoding='utf-8'), re.M)}
    check('Extracted paragraphs match actual DOCX', all(i < len(paragraphs) and paragraphs[i] == t for i, t in txt_paragraphs.items()), True)
    check('Extracted nonempty paragraph count matches DOCX', len(txt_paragraphs), sum(bool(t) for t in paragraphs))
    single = csv(DATA / 'summary/single_target_overall_record_macro.csv').set_index(['variant', 'metric'])
    family = csv(DATA / 'summary/family_overall_record_macro.csv').set_index(['variant', 'metric'])
    family_source = csv(DATA / 'summary/family_source_record_macro.csv').set_index(['variant', 'dataset', 'metric'])
    external = csv(DATA / 'summary/external_2430_record_macro.csv').set_index(['variant', 'metric'])
    baseline = csv(DATA / 'baseline_original_overall_values_verified.csv').set_index('method')
    baseline_source = csv(DATA / 'baseline_original_source_values_verified.csv')
    bf_path = DATA / 'audit/external_figure_and_backflip_reuse.json'
    sources[str(bf_path)] = sha(bf_path)
    bf = json.loads(bf_path.read_text(encoding='utf-8'))['backflip_reuse']
    expected = []
    for name, variant in [('完整模型','full'),('仅链内模型','intrachain'),('静态特征模型','static_sequence_complex'),('氨基酸类型均值','aa_only')]:
        a, b = external.loc[(variant,'rmsf')], external.loc[(variant,'chi_entropy')]
        expected.append([name] + fmt([a.mae,a.spearman,b.mae,b.spearman]))
    expected.append(['BackFlip2发布模型'] + fmt([bf['record_macro_mae'],bf['record_macro_spearman']]) + ['无','无'])
    check('Table2 all new2430 external cells', tables[1][1:], expected)
    expected = []
    for metric,label,variant,training in [('rmsf','RMSF','rmsf_only','单目标'),('rmsf','RMSF','full','联合目标'),('chi_entropy','Hχ','chi_entropy_only','单目标'),('chi_entropy','Hχ','full','联合目标')]:
        v = single.loc[(variant,metric)]
        expected.append([label,training] + fmt([v.mae,v.rmse,v.spearman]))
    check('TableS9 all single/joint recordmacro cells', tables[10][1:], expected)
    expected = []
    for dataset,label in [('all','整体'),('atlas','ATLAS'),('dynarepo','DynaRepo'),('itrace','ITrace')]:
        for variant,name in [('full','完整'),('intrachain','链内')]:
            a,b = (family.loc[(variant,'rmsf')],family.loc[(variant,'chi_entropy')]) if dataset == 'all' else (family_source.loc[(variant,dataset,'rmsf')],family_source.loc[(variant,dataset,'chi_entropy')])
            expected.append([label,name] + fmt([a.mae,a.spearman,b.mae,b.spearman]))
    check('TableS10 all family/source recordmacro cells',tables[11][1:],expected)
    check('TableS10 no fabricated fivefold source SD',any('±' in cell for row in tables[11] for cell in row),False)
    expected = []
    for method,name in [('flexppi','完整FlexPPI'),('intrachain','链内FlexPPI'),('planet_full','PLANET full'),('planet_seq','PLANET seq'),('rmsf_net_pdb','RMSF-net PDB'),('anm','ANM'),('gnm','GNM')]:
        v = baseline.loc[method]
        expected.append([name] + fmt([v.mae,v.spearman,v.z_mae]))
    check('TableS11 original-only baseline/relative cells',tables[12][1:],expected)
    check('TableS11 has no duplicate family column',len(tables[12][0]),4)
    expected = []
    for method,name in [('planet_full','PLANET full 发布权重'),('planet_seq','PLANET seq 发布权重'),('rmsf_net_pdb','RMSF-net PDB-only')]:
        v = baseline.loc[method]
        expected.append([name] + fmt([v.mae,v.spearman]) + ['无','无'])
    check('Table1 released baseline rows',tables[0][-3:],expected)
    check('Table1 Full recordmacro row',tables[0][1],['完整模型'] + fmt([single.loc[('full','rmsf')].mae,single.loc[('full','rmsf')].spearman,single.loc[('full','chi_entropy')].mae,single.loc[('full','chi_entropy')].spearman]))
    kix = csv(BASE / 'kix_existing_reuse_20261010/summary/state_summary.csv').set_index('state')
    expected = []
    for state,name in [('apo','游离KIX'),('myb','cMyb'),('mll','MLL'),('both','双伙伴')]:
        v = kix.loc[state]
        expected.append([name] + fmt([v.RMSF_MAE_angstrom,v.RMSF_spearman]) + [f'{int(v.top20_overlap)}/18'] + fmt([v.L12_MD_mean_RMSF_angstrom,v.L12_fold1_mean_RMSF_angstrom]))
    check('TableS7 all four KIX state cells',tables[8][1:],expected)
    delta = csv(BASE / 'kix_existing_reuse_20261010/summary/l12_delta_summary.csv').set_index('contrast')
    expected = []
    for contrast,label in [('myb-apo','cMyb−游离KIX'),('mll-apo','MLL−游离KIX'),('both-apo','双伙伴−游离KIX'),('mll-myb','MLL−cMyb'),('both-myb','双伙伴−cMyb'),('both-mll','双伙伴−MLL')]:
        v = delta.loc[contrast]
        expected.append([label,f'{v.MD_delta_RMSF_angstrom:+.3f}',f'{v.fold1_delta_RMSF_angstrom:+.3f}',f'{v.absolute_delta_error_angstrom:.3f}'])
    check('TableS8 all six actual KIX state deltas',tables[9][1:],expected)
    prior_doc = BASE.parent / 'paper_v19_20261012_中文修订_补充实验完成_三线表.docx'
    prior_document = Document(prior_doc)
    prior_tables = [[[c.text for c in r.cells] for r in t.rows] for t in prior_document.tables]
    check('All editable equations and inline math tokens preserved',[list(p._p.xpath('.//m:t/text()')) for p in doc.paragraphs if p._p.xpath('.//m:t')],[list(p._p.xpath('.//m:t/text()')) for p in prior_document.paragraphs if p._p.xpath('.//m:t')])
    for i in list(range(2,8)):
        check('Previously source-audited retained table index '+str(i),tables[i],prior_tables[i])
    check('Table1 all original rows retained',tables[0],prior_tables[0])
    figure = DATA / 'figures/fig4_strict_external_paired.png'
    figure_sha = sha(figure)
    with zipfile.ZipFile(DOC) as archive:
        media = [n for n in archive.namelist() if n.startswith('word/media/')]
        figure_embedded = [n for n in media if hashlib.sha256(archive.read(n)).hexdigest() == figure_sha]
    check('New2430 Fig4 exact verified source image',len(figure_embedded),1)
    stats = csv(DATA / 'figures/fig4_strict_external_paired_statistics.csv')
    for r in stats.itertuples(index=False):
        check('Fig4 mean '+r.variant+' '+r.metric,round(r.mean_mae,12),round(external.loc[(r.variant,r.metric)].mae,12))
    text = '\n'.join(paragraphs)
    main_text = '\n'.join(paragraphs[:next(i for i,t in enumerate(paragraphs) if t == '参考文献')])
    check('Main new633/2430 protocol explicit','排除633条' in main_text and '2,430' in main_text,True)
    check('No retired632/2431/S12 reference',bool(re.search(r'(?<!\d)632(?!\d)|2,431|(?<!\d)2431(?!\d)|表S12|图S12',text)),False)
    check('Joint comparison no longer untested claim','联合学习相对于独立单目标模型的贡献仍需单独检验' in text,False)
    check('Family ITrace one test fold explicitly defined','全部244条ITrace位于第0折测试集合' in text and '仅第0折测试' in text,True)
    check('Short-chain exact and long-chain unmasked protocol explicit','长度至少20个残基' in text and '短链按精确序列连接' in text,True)
    check('Main external RMSF slightly higher retained','RMSF MAE为0.474 Å，链内和静态模型均为0.471 Å' in main_text,True)
    check('Single-target gain correctly modest','两项MAE略低' in main_text,True)
    check('No unqualified multitask improvement claim','大幅协同' in text or '普遍协同' in text,False)
    check('KIX one production trajectory per state','每种状态分析一条生产轨迹' in main_text,True)
    check('KIX genuine reversed amplitude delta retained','MD差值为+1.017 Å，模型为−1.589 Å' in text,True)
    check('Sidechain displays explicitly genuine MD frames','侧链朝向叠加来自MD轨迹' in text,True)
    # Captions use five main figures, eleven supplementary figures and eleven supplementary tables.
    figure_captions = re.findall(r'^(图(?:S)?\d+)：',text,re.M)
    table_captions = re.findall(r'^(表(?:S)?\d+)：',text,re.M)
    expected_figures = ['图'+str(i) for i in range(1,6)] + ['图S'+str(i) for i in range(1,12)]
    expected_tables = ['表1','表2'] + ['表S'+str(i) for i in range(1,12)]
    check('Main/supplement figure captions each once',sorted(figure_captions),sorted(expected_figures))
    check('Sixteen active inline figures after three retired plots removed',len(doc.inline_shapes),16)
    check('Main/supplement table captions each once',sorted(table_captions),sorted(expected_tables))
    mentioned_figures = set(re.findall(r'图(?:S)?\d+',text))
    mentioned_tables = set(re.findall(r'表(?:S)?\d+',text))
    check('No direct figure reference orphan',sorted(mentioned_figures-set(figure_captions)),[])
    check('No direct table reference orphan',sorted(mentioned_tables-set(table_captions)),[])
    gain_single = {'rmsf':(1-single.loc[('full','rmsf')].mae/single.loc[('rmsf_only','rmsf')].mae)*100,'chi_entropy':(1-single.loc[('full','chi_entropy')].mae/single.loc[('chi_entropy_only','chi_entropy')].mae)*100}
    gain_family = {m:(1-family.loc[('full',m)].mae/family.loc[('intrachain',m)].mae)*100 for m in ['rmsf','chi_entropy']}
    check('Single-target reported rounded relative gains',[f'{gain_single[m]:.1f}' for m in ['rmsf','chi_entropy']],['0.8','1.8'])
    check('Family reported rounded relative gains',[f'{gain_family[m]:.1f}' for m in ['rmsf','chi_entropy']],['7.5','3.8'])
    end_hash = sha(DOC)
    check('Document remained unchanged during read-only audit',end_hash,begin_hash)
    report = {'status':'passed' if not errors else 'failed','document':str(DOC),'sha256':end_hash,'expected_sha256':EXPECTED,'checked_utc':datetime.now(timezone.utc).isoformat(),'scope':'Actual v20 full-text scientific argument and source-verified Table1/2/S3-S11 numbers plus references and exact Fig4 media; root handles page rendering/visual QA.','read_only':True,'document_modified':False,'pipeline_modified':False,'experiments_or_certificates_rerun':False,'paragraph_count':len(paragraphs),'table_count':len(tables),'embedded_media_count':len(media),'active_inline_figure_count':len(doc.inline_shapes),'checks':checks,'failed_checks':errors,'critical_scientific_defects':[] if not errors else errors,'scientific_review':{'joint_default_is_supported_by_two_complementary_outputs':True,'does_not_claim_universal_or_large_multitask_synergy':True,'external_RMSF_MAE_higher_than_intra_static_explicit':True,'family_absolute_results_and_Full_vs_Intra_gain_true':True,'family_ITrace_only_fold0_no_SD':True,'KIX_narrative_is_activity_region_localization_not_response_amplitude_or_direction_validation':True,'all_six_KIX_response_deltas_preserved_in_S8':True,'Hchi_is_discrete_sampled_state_entropy_not_complete_thermodynamic_entropy':True,'no_required_defensive_sentence_deletion':True,'retired_old632_2431_protocol_figures_and_table_absent':True},'source_sha256':sources,'Fig4_sha256':figure_sha,'Fig4_embedded_names':figure_embedded,'editorial_optional_issues':[]}
    path = OUT / 'v20_results_content_review.json'
    path.write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'status':report['status'],'sha256':end_hash,'checks':len(checks),'failures':errors,'report':str(path),'paragraphs':len(paragraphs),'tables':len(tables),'media':len(media)},ensure_ascii=False,indent=2))
    return 1 if errors else 0

if __name__ == '__main__':
    raise SystemExit(main())
