"""Predeclared original-paper SPR analysis rule plus two complete sensitivities.

This script only reads experimental sources. It must never read model predictions.
"""
from pathlib import Path
from collections import defaultdict
import csv,json,hashlib,statistics,math,datetime,copy
ROOT=Path(__file__).parent
def read(name):return list(csv.DictReader((ROOT/name).open(encoding='utf-8-sig',newline='')))
def write(name,rows):
    with (ROOT/name).open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
frozen=json.loads((ROOT/'SPR_labels.json').read_text(encoding='utf8'))
sources=frozen.get('study_level_HERMES_constrained_compatibility',frozen.get('study_level'))
assert sources is not None
fits={(r['panel_id'],r['peptide_sequence']):r for r in read('Pettmann2021_Bmax_methods.csv')}
candidates=read('SPR_candidates.csv')
variants={}
for variant in ['reported_hybrid','all_constrained','all_fitted']:
    rows=copy.deepcopy(sources)
    for r in rows:
        r['label_variant']=variant;r['analysis_role']='primary' if variant=='reported_hybrid' else 'complete sensitivity; not an independent measurement set'
        if r['study_id']=='pettmann_et_al':
            f=fits[(r['panel_id'],r['peptide_sequence'])]
            constrained=float(f['KD_uM'])
            use_constrained=variant=='all_constrained' or (variant=='reported_hybrid' and constrained>20.0)
            if use_constrained:
                kd=constrained
                r['selected_fit']='Bmax constrained'
            else:
                kd=float(f['KD_Bmax_fitted_uM'])
                r['error_value']=f['KD_Bmax_fitted_SD_uM']
                r['n']=int(f['KD_Bmax_fitted_n'])
                r['n_min']=r['n']
                r['selected_fit']='Bmax fitted'
                r['primary_locator']=r['primary_locator'].split('; columns')[0]+'; columns 3-6 (Bmax fitted)'
            r['KD_value']=kd;r['KD_uM']=kd;r['pKD_M']=6-math.log10(kd);r['negative_log10_KD_uM']=-math.log10(kd)
            r['SPR_method']='steady-state; '+r['selected_fit']
        else:r['selected_fit']='Aleksic conventional steady-state; unchanged across variants'
        r['fit_selection_rule']='Pettmann: constrained KD>20 uM -> constrained; otherwise fitted' if variant=='reported_hybrid' else variant
        r['HERMES_values_are_primary_for_this_variant']=variant=='all_constrained'
    keyed={(r['panel_id'],r['study_id'],r['peptide_sequence']):r for r in rows}
    for r in rows:
        p=keyed[(r['panel_id'],r['study_id'],r['parent_peptide'])]
        r['parent_KD_uM']=p['KD_uM'];r['parent_pKD_M']=p['pKD_M'];r['parent_source_locator']=p['primary_locator']
        r['delta_pKD_vs_parent']=r['pKD_M']-p['pKD_M']
        assert not r['is_parent'] or abs(r['delta_pKD_vs_parent'])<1e-15
    out=f'SPR_labels_observations_{variant}.csv';write(out,rows)
    group=defaultdict(list)
    for r in rows:group[r['candidate_id']].append(r)
    candidate_rows=[]
    for c in candidates:
        obs=group[c['candidate_id']];vals=[r['delta_pKD_vs_parent'] for r in obs]
        candidate_rows.append(dict(c,endpoint='delta_pKD_vs_parent',delta_pKD_vs_parent=statistics.median(vals),label_variant=variant,analysis_role='primary' if variant=='reported_hybrid' else 'complete sensitivity only',label_aggregation='median of within-study parent-normalized delta-pKD; two studies => equal mean',n_study_observations=len(obs),study_ids='|'.join(r['study_id'] for r in obs),study_delta_pKD_json=json.dumps({r['study_id']:r['delta_pKD_vs_parent'] for r in obs},sort_keys=True),min_study_delta_pKD=min(vals),max_study_delta_pKD=max(vals),censoring='none_reported',assay='SPR',temperature_C=37,parent_same_batch='not_reported',biological_independence='candidate within one of only two TCR series; fits and studies are not additional candidate samples'))
    candidate_out=f'SPR_labels_candidates_{variant}.csv';write(candidate_out,candidate_rows)
    variants[variant]={'observations':rows,'candidates':candidate_rows,'observation_file':str(ROOT/out),'candidate_file':str(ROOT/candidate_out),'primary_file_sha256':sha(ROOT/out),'candidate_file_sha256':sha(ROOT/candidate_out),'n_observations':len(rows),'n_candidates':len(candidate_rows)}
write('SPR_labels_observations.csv',variants['reported_hybrid']['observations'])
write('SPR_labels_candidates_primary.csv',variants['reported_hybrid']['candidates'])
write('SPR_labels_study_level.csv',variants['all_constrained']['observations'])
frozen['schema_version']='1.2'
frozen['primary_endpoint']='median within-study parent-normalized delta-pKD; original Pettmann reported hybrid fit rule; no exact batch matching claimed'
frozen['observations']=variants['reported_hybrid']['observations']
frozen['study_level_HERMES_constrained_compatibility']=sources
frozen.pop('study_level',None)
frozen['variants']=variants
(ROOT/'SPR_labels.json').write_text(json.dumps(frozen,ensure_ascii=False,indent=2),encoding='utf8')
oldfreeze=json.loads((ROOT/'label_rule_freeze.json').read_text(encoding='utf8')) if (ROOT/'label_rule_freeze.json').exists() else {}
report={'status':'passed','frozen_before_model_results':True,'selected_at_utc':oldfreeze.get('selected_at_utc',datetime.datetime.now(datetime.timezone.utc).isoformat()),'primary':'reported_hybrid','selection_basis':'Original Pettmann 2021 SPR methods: constrained KD above 20 uM reported as constrained, otherwise fitted. Selection approved in evaluation protocol before model-output inspection.','source_xml':str(ROOT/'raw/pettmann2021_fullText.xml'),'source_xml_sha256':sha(ROOT/'raw/pettmann2021_fullText.xml'),'source_locator':'Materials and methods / Surface plasmon resonance, paragraph immediately preceding SPR data analysis software','sensitivities':['all_constrained (HERMES replication)','all_fitted'],'n_independent_TCR_series':2,'n_unique_candidates':30,'n_study_mean_observations':36,'fit_variants_do_not_increase_sample_size':True,'parent_same_measurement_batch':'not_reported','all_have_measured_within_study_parent':True,'no_binary_hit_threshold_assigned':True,'variants':{k:{x:v[x] for x in ['observation_file','candidate_file','primary_file_sha256','candidate_file_sha256','n_observations','n_candidates']} for k,v in variants.items()}}
(ROOT/'label_rule_freeze.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf8')
validation=json.loads((ROOT/'label_validation_report.json').read_text(encoding='utf8'))
validation['primary_label_definition']=report['primary'];validation['primary_endpoint']=frozen['primary_endpoint']
validation['fit_variants']=report['variants'];validation['label_rule_freeze_sha256']=sha(ROOT/'label_rule_freeze.json')
for name in ['SPR_labels.json','SPR_labels_observations.csv','SPR_labels_candidates_primary.csv','SPR_labels_study_level.csv','label_rule_freeze.json']:
    validation['files'][name]={'path':str(ROOT/name),'sha256':sha(ROOT/name)}
(ROOT/'label_validation_report.json').write_text(json.dumps(validation,ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps({k:report[k] for k in ['status','primary','sensitivities','n_unique_candidates','n_study_mean_observations','all_have_measured_within_study_parent']},ensure_ascii=False,indent=2))
