"""Independent checks of delivered tables against raw experimental CSVs."""
from pathlib import Path
from collections import Counter,defaultdict
import csv,json,math,statistics,hashlib,datetime
P=Path(__file__).parent
def read(path):return list(csv.DictReader(Path(path).open(encoding='utf-8-sig',newline='')))
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
primary=read(P/'SPR_labels_observations.csv')
candidates=read(P/'SPR_labels_candidates_primary.csv')
assert len(primary)==36 and len(candidates)==30
assert len(set(x['observation_id'] for x in primary))==36
assert len(set(x['candidate_id'] for x in candidates))==30
expected={}
for pid,source in [('1G4_NYESO','elife-67092-fig1-data2-v3.csv'),('A6_TAX','elife-67092-fig1-data3-v3.csv')]:
    raw=list(csv.reader((P/'raw'/source).open(encoding='utf-8-sig',newline='')))[2:]
    for r in raw:
        if r:
            expected[(pid,r[1])]=(float(r[6]) if float(r[6])>20 else float(r[2]))
source_checks=0
for r in primary:
    if r['study_id']=='pettmann_et_al':
        assert float(r['KD_uM'])==expected[(r['panel_id'],r['peptide_sequence'])]
        source_checks+=1
    assert abs(float(r['pKD_M'])-(6-math.log10(float(r['KD_uM']))))<1e-12
    assert r['temperature_C']=='37' and r['assay']=='SPR' and r['KD_unit']=='uM'
    assert r['KD_relation']=='=' and r['censoring']=='none_reported'
    assert r['parent_same_measurement_batch']=='not_reported'
index={r['observation_id']:r for r in primary}
g=defaultdict(list)
for r in primary:
    parent=index[r['parent_sample']]
    assert parent['study_id']==r['study_id'] and parent['panel_id']==r['panel_id']
    assert abs(float(r['delta_pKD_vs_parent'])-(float(r['pKD_M'])-float(parent['pKD_M'])))<1e-12
    g[r['candidate_id']].append(float(r['delta_pKD_vs_parent']))
for r in candidates:
    assert abs(float(r['delta_pKD_vs_parent'])-statistics.median(g[r['candidate_id']]))<1e-12
    assert int(r['n_study_observations'])==len(g[r['candidate_id']])
file_report=json.loads((P/'label_validation_report.json').read_text(encoding='utf8'))
for name,record in file_report['files'].items():assert sha(record['path'])==record['sha256'],name
source_report=json.loads((P/'source_manifest.json').read_text(encoding='utf8'))
for f in source_report['files']:assert sha(f['path'])==f['sha256'],f['path']
report={'status':'passed','checked_at_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'primary_rule':'reported_hybrid','raw_Pettmann_primary_rule_checks':source_checks,'original_source_value_checks_in_build_report':36,'parent_delta_checks':36,'candidate_aggregation_checks':30,'independent_TCR_series':2,'study_groups':{':'.join(k):n for k,n in Counter((r['panel_id'],r['study_id']) for r in primary).items()},'n_censored':0,'all_label_and_raw_source_hashes_match':True,'no_model_outputs_read':True,'raw_source_files_checked':len(source_report['files'])}
(P/'independent_labels_audit.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf8')
deliver=['README.md','SPR_candidates.csv','SPR_labels_candidates_primary.csv','SPR_labels_observations.csv','SPR_labels.json','templates.json','label_rule_freeze.json','label_validation_report.json','independent_labels_audit.json','source_manifest.json']
deliver += [f'SPR_labels_{kind}_{variant}.csv' for kind in ['candidates','observations'] for variant in ['reported_hybrid','all_constrained','all_fitted']]
(P/'delivery_manifest.json').write_text(json.dumps({'status':'passed','main_candidate_table':'SPR_candidates.csv','main_label_table':'SPR_labels_candidates_primary.csv','main_observation_table':'SPR_labels_observations.csv','files':{n:{'path':str(P/n),'sha256':sha(P/n)} for n in deliver}},ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps(report,ensure_ascii=False,indent=2))
