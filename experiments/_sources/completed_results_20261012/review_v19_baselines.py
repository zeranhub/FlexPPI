"""Read-only v19 baseline, reference connectivity and numerical-content review."""
from pathlib import Path
from collections import Counter
import hashlib
import json
import re
import zipfile
import datetime
import pandas as pd
from docx import Document
from lxml import etree

TASK=Path(r'E:\Project\FlexPPI\experiments\_sources')
OUT=TASK/'completed_results_20261012'
DOC=TASK.parent/'paper_v19_20261012_中文修订_补充实验完成_三线表.docx'
PRIOR=TASK.parent/'paper_v18_20261010_中文修订_复用KIX单轨迹_三线表.docx'
FROZEN=Path(r'E:\Project\FlexPPI\experiments\core_corrected_pipeline')
EXPECTED_SHA='d5104be36fa74196b336d0518729daeb386a7da9f5a8e3d9a4910af96e2e3dab'
NS={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}


def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def fmt(value): return f'{float(value):.3f}'
def contents(table): return [[c.text for c in row.cells] for row in table.rows]


def main():
    reviewed_sha=sha(DOC)
    document=Document(DOC)
    paragraphs={i:p.text for i,p in enumerate(document.paragraphs)}
    text='\n'.join(paragraphs.values())
    assert len(document.tables)==14
    original=pd.read_csv(TASK/'baselines/additional_original_oof/record_macro_metrics.csv').set_index(['method','dataset'])
    family=pd.read_csv(TASK/'baselines/additional_family_oof/record_macro_metrics.csv').set_index(['method','dataset'])
    numerical=json.loads((OUT/'baseline_numeric_audit.json').read_text())
    reference=json.loads((OUT/'planet_reference_metadata.json').read_text())
    neural=pd.read_csv(FROZEN/'neural_exports/evaluation/per_record_metrics.csv')
    static=pd.read_csv(FROZEN/'static_corrected/oof_per_record_metrics.csv')
    neural_macro=neural.groupby(['variant','metric'])[['mae','spearman']].mean()
    static_macro=static.groupby(['variant','metric'])[['mae','spearman']].mean()
    table1=contents(document.tables[0])
    checked_cells=0
    for index,label,variant,source in [(1,'完整模型','full',neural_macro),
                                      (2,'仅链内模型','intrachain',neural_macro),
                                      (3,'静态特征模型','static_sequence_complex',static_macro),
                                      (4,'氨基酸类型均值','aa_only',static_macro)]:
        expected=[label,fmt(source.loc[(variant,'rmsf'),'mae']),fmt(source.loc[(variant,'rmsf'),'spearman']),
                  fmt(source.loc[(variant,'chi_entropy'),'mae']),fmt(source.loc[(variant,'chi_entropy'),'spearman'])]
        assert table1[index]==expected,(index,table1[index],expected)
        checked_cells+=5
    bf=numerical['backflip2_prior_scoring']['overall']
    assert table1[5]==['BackFlip2发布模型',fmt(bf['mae']),fmt(bf['spearman']),'无','无']
    checked_cells+=5
    for index,method,label in [(6,'planet_full','PLANET full 发布权重'),
                              (7,'planet_seq','PLANET seq 发布权重'),
                              (8,'rmsf_net_pdb','RMSF-net PDB-only')]:
        row=original.loc[(method,'all')]
        expected=[label,'不适用' if method=='rmsf_net_pdb' else fmt(row.mae),fmt(row.spearman),'无','无']
        assert table1[index]==expected,(index,table1[index],expected)
        checked_cells+=5
    table11=contents(document.tables[12])
    names=['完整FlexPPI','链内FlexPPI','PLANET full','PLANET seq','RMSF-net PDB','ANM','GNM']
    methods=['flexppi','intrachain','planet_full','planet_seq','rmsf_net_pdb','anm','gnm']
    for index,(method,name) in enumerate(zip(methods,names),1):
        first=original.loc[(method,'all')]; second=family.loc[(method,'all')]
        relative=method in {'rmsf_net_pdb','anm','gnm'}
        expected=[name,'不适用' if relative else fmt(first.mae),fmt(first.spearman),
                  '不适用' if relative else fmt(second.mae),fmt(second.spearman),fmt(first.z_mae)]
        assert table11[index]==expected,(index,table11[index],expected)
        checked_cells+=6
    assert f'家族完整/链内标准化MAE为{fmt(family.loc[("flexppi","all"),"z_mae"])}/{fmt(family.loc[("intrachain","all"),"z_mae"])}' in text
    full=original.loc[('flexppi','all')]; planet=original.loc[('planet_full','all')]
    assert planet.mae<full.mae and planet.spearman>full.spearman and planet.rmse>full.rmse
    assert f'整体RMSF MAE为{fmt(planet.mae)} Å、ρ为{fmt(planet.spearman)}' in text
    assert f'RMSE为{fmt(planet.rmse)} Å，高于FlexPPI的{fmt(full.rmse)} Å' in text
    for dataset in ['atlas','dynarepo','itrace']:
        row=original.loc[('planet_full',dataset)]
        assert f'{fmt(row.mae)} Å/{fmt(row.spearman)}' in text
    assert '+1.381 Å' in text and '+2.150 Å' in text
    # Check the other changed primary numerical tables against their saved source summaries.
    single=pd.read_csv(OUT/'summary/single_target_overall_record_macro.csv').set_index(['variant','metric'])
    for index,variant,metric,label in [(1,'rmsf_only','rmsf','单目标'),(2,'full','rmsf','联合目标'),
                                       (3,'chi_entropy_only','chi_entropy','单目标'),(4,'full','chi_entropy','联合目标')]:
        row=single.loc[(variant,metric)]
        assert contents(document.tables[10])[index]==['RMSF' if metric=='rmsf' else 'Hχ',label,fmt(row.mae),fmt(row.rmse),fmt(row.spearman)]
        checked_cells+=5
    rmsf_gain=single.loc[('rmsf_only','rmsf'),'mae']-single.loc[('full','rmsf'),'mae']
    hchi_gain=single.loc[('chi_entropy_only','chi_entropy'),'mae']-single.loc[('full','chi_entropy'),'mae']
    assert f'{rmsf_gain:.6f} Å' in text and f'{hchi_gain:.6f} nats' in text
    assert single.loc[('full','rmsf'),'rmse']>single.loc[('rmsf_only','rmsf'),'rmse']
    assert single.loc[('full','rmsf'),'spearman']<single.loc[('rmsf_only','rmsf'),'spearman']
    family_overall=pd.read_csv(OUT/'summary/family_overall_record_macro.csv').assign(dataset='all')
    family_source=pd.read_csv(OUT/'summary/family_source_record_macro.csv')
    families=pd.concat([family_overall,family_source],ignore_index=True).set_index(['variant','dataset','metric'])
    table10=contents(document.tables[11])
    for index,dataset,variant in [(1,'all','full'),(2,'all','intrachain'),(3,'atlas','full'),(4,'atlas','intrachain'),
                                  (5,'dynarepo','full'),(6,'dynarepo','intrachain'),(7,'itrace','full'),(8,'itrace','intrachain')]:
        rmsf=families.loc[(variant,dataset,'rmsf')]; hchi=families.loc[(variant,dataset,'chi_entropy')]
        source_name={'all':'整体','atlas':'ATLAS','dynarepo':'DynaRepo','itrace':'ITrace'}[dataset]
        expected=[source_name,'完整' if variant=='full' else '链内',fmt(rmsf.mae),fmt(rmsf.spearman),fmt(hchi.mae),fmt(hchi.spearman)]
        assert table10[index]==expected
        checked_cells+=6
    external=pd.read_csv(OUT/'summary/external_2430_record_macro.csv').set_index(['variant','metric'])
    table2=contents(document.tables[1])
    for index,variant,label in [(1,'full','完整模型'),(2,'intrachain','仅链内模型'),
                                (3,'static_sequence_complex','静态特征模型'),(4,'aa_only','氨基酸类型均值')]:
        r=external.loc[(variant,'rmsf')]; h=external.loc[(variant,'chi_entropy')]
        assert table2[index]==[label,fmt(r.mae),fmt(r.spearman),fmt(h.mae),fmt(h.spearman)]
        checked_cells+=5
    for comparator in ['intrachain','static_sequence_complex']:
        assert external.loc[('full','rmsf'),'mae']>external.loc[(comparator,'rmsf'),'mae']
        assert external.loc[('full','chi_entropy'),'mae']<external.loc[(comparator,'chi_entropy'),'mae']
        for target in ['rmsf','chi_entropy']:
            assert external.loc[('full',target),'spearman']>external.loc[(comparator,target),'spearman']
    assert '但RMSF幅值误差略高' in text
    rmsf_static=static_macro.loc[('static_sequence_complex','rmsf'),'mae']
    hchi_static=static_macro.loc[('static_sequence_complex','chi_entropy'),'mae']
    relative_static=[100*(rmsf_static-neural_macro.loc[('full','rmsf'),'mae'])/rmsf_static,
                     100*(hchi_static-neural_macro.loc[('full','chi_entropy'),'mae'])/hchi_static]
    assert '原体系留出的五折评价' in paragraphs[2]
    assert f'{relative_static[0]:.1f}%' in paragraphs[2] and f'{relative_static[1]:.1f}%' in paragraphs[2]
    # Citation expansion handles ranged citations, leaving mathematical zero-based intervals aside.
    refs={}; mentions={}
    for index,paragraph in paragraphs.items():
        match=re.match(r'^\[(\d+)\]\s',paragraph)
        if match:
            number=int(match.group(1)); assert number not in refs
            refs[number]=dict(paragraph=index,text=paragraph)
            continue
        for group in re.findall(r'\[([1-9][0-9]*(?:[–,-][1-9][0-9]*)*)\]',paragraph):
            for item in group.split(','):
                ends=re.split(r'[–-]',item)
                values=range(int(ends[0]),int(ends[-1])+1) if len(ends)>1 else [int(item)]
                for number in values: mentions.setdefault(number,[]).append(index)
    assert set(refs)==set(range(1,32)) and set(mentions)==set(refs)
    assert all('PLANET' in paragraphs[index] for index in mentions[31])
    assert 'ProDy 2.0' in refs[3]['text'] and 'ProDy: protein dynamics' in refs[30]['text']
    assert reference['title'] in refs[31]['text'] and reference['doi'] in refs[31]['text']
    assert 'Samuel Sledzieski and Sonya M. Hanson' in refs[31]['text']
    assert '2026' in refs[31]['text'] and 'version 2' in refs[31]['text']
    assert '[30,31]' not in text
    prody_context=next((dict(paragraph=index,text=p) for index,p in paragraphs.items() if 'ProDy 2.6.1' in p),None)
    assert prody_context and '[3,30]' in prody_context['text']
    # Restore-byte safeguard independently compared with the prior annotated v18.
    with zipfile.ZipFile(PRIOR) as a,zipfile.ZipFile(DOC) as b:
        old=a.read('word/comments.xml'); new=b.read('word/comments.xml')
        texts=lambda raw:[(c.get('{'+NS['w']+'}id'),''.join(c.xpath('.//w:t/text()',namespaces=NS)))
                          for c in etree.fromstring(raw).xpath('/w:comments/w:comment',namespaces=NS)]
        anchors=lambda z:Counter((e.tag,e.get('{'+NS['w']+'}id')) for e in etree.fromstring(z.read('word/document.xml')).iter()
                                 if e.tag.endswith(('commentRangeStart','commentRangeEnd','commentReference')))
        assert old==new and texts(old)==texts(new) and anchors(a)==anchors(b)
        comments=dict(comparison_source=str(PRIOR),byte_identical=True,texts_and_ids_identical=True,
                      anchors_identical=True,comments=len(texts(old)),
                      comments_xml_sha256=hashlib.sha256(new).hexdigest())
    limit_phrases=['encoded结构编码','默认300 K及float32','nm乘以10转换为Å',
                   '只评价排序及记录内标准化剖面误差，不赋予Å单位',
                   '未与原Chimera操作逐体素验证','未使用真实冷冻电镜密度',
                   '训练交叠未逐条核实','不意味着发布模型已按本文家族划分重训',
                   '不能作为相同家族划分下从头训练的结果','不应被扩大为所有柔性指标上的领先']
    assert all(phrase in text for phrase in limit_phrases)
    assert '末列对应原体系留出' in text
    report=dict(status='passed_no_required_baseline_or_reference_corrections',
                created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                document=str(DOC),document_sha256=reviewed_sha,initial_request_sha256=EXPECTED_SHA,
                document_updated_by_root_during_review=(reviewed_sha!=EXPECTED_SHA),DOCX_edited=False,
                scope='Baseline/reference and changed numerical-content logic, not visual layout or complete scientific reference metadata re-verification.',
                tables=14,table_cells_checked=checked_cells,
                checks=dict(Table1_all_eight_method_rows_match=True,TableS11_all_seven_methods_both_protocols_match=True,
                            S11_zMAE_original_column_and_family_caption_correct=True,
                            PLANET_MAE_rho_better_RMSE_worse_correct=True,PLANET_source_differences_correct=True,
                            TableS9_single_target_values_and_gain_direction_correct=True,
                            TableS10_family_values_correct=True,Table2_new_2430_values_and_metric_directions_correct=True,
                            abstract_original_static_MAE_reductions=relative_static,
                            source_or_family_scope_not_mixed=True,published_models_not_claimed_same_training_corpus=True,
                            no_RMSF_net_Angstrom_MAE=True,no_SHP_as_Hchi=True,
                            encoded_300K_float32_nm_times10_and_simulated_density_limits_explicit=True,
                            no_universal_RMSF_superiority_claim_remaining=True,
                            citation_31_PLANET_only=True,ProDy2_reference3_and_ProDy1_reference30_correct=True,
                            old_30_31_citation_removed=True),
                references=dict(numbers=sorted(refs),missing_references=[],uncited_references=[],
                                PLANET_31=refs[31],PLANET_31_citation_paragraphs=mentions[31],
                                ProDy2_3=refs[3],ProDy1_30=refs[30],ProDy_supplement_context=prody_context),
                comments=comments,issues=[],must_fix_issues=[],
                source_sha256={str(p):sha(p) for p in [TASK/'baselines/additional_original_oof/record_macro_metrics.csv',
                    TASK/'baselines/additional_family_oof/record_macro_metrics.csv',OUT/'summary/single_target_overall_record_macro.csv',
                    OUT/'summary/family_overall_record_macro.csv',OUT/'summary/family_source_record_macro.csv',
                    OUT/'summary/external_2430_record_macro.csv',OUT/'planet_reference_metadata.json',OUT/'baseline_numeric_audit.json',
                    FROZEN/'neural_exports/evaluation/per_record_metrics.csv',FROZEN/'static_corrected/oof_per_record_metrics.csv']},
                conclusion='All checked table cells and body interpretations match the frozen saved results. References 1–31 are connected, reference 31 correctly names PLANET and does not serve ProDy contexts. Six inherited annotations and their anchor inventory remain intact, including byte-identical comments XML. No baseline/reference content correction is required before root visual QA.')
    assert sha(DOC)==reviewed_sha, 'The document changed while this audit was running; rerun on the new saved document.'
    (OUT/'v19_baseline_content_review.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'status':report['status'],'document_sha256':report['document_sha256'],
                      'table_cells_checked':checked_cells,'references':len(refs),'comments':comments,'issues':[]},indent=2))


if __name__=='__main__': main()
