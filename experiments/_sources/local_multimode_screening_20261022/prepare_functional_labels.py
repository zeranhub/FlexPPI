"""Source-only curation and endpoint binding, before new FlexPPI inference.

Run using the bundled document Python (openpyxl available). Source CSV fit
values are retained for provenance but never used as the primary phenotype.
"""
from pathlib import Path
from datetime import datetime, timezone
import csv, json, math, hashlib, re
import openpyxl

ROOT=Path(__file__).resolve().parent
RAW=ROOT/'antigen_labels/raw'
OUT=ROOT/'antigen_labels/curated'
OUT.mkdir(exist_ok=True)
AA=set('ACDEFGHIKLMNPQRSTVWY')
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p,x):p.write_text(json.dumps(x,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf8')
def csvout(p,rr):
    with p.open('w',encoding='utf-8-sig',newline='') as h:
        w=csv.DictWriter(h,fieldnames=list(rr[0]));w.writeheader();w.writerows(rr)
def num(x):
    try:v=float(x)
    except (ValueError,TypeError):return None
    return v if math.isfinite(v) else None

assert not list((ROOT/'predictions').glob('**/prediction.npz')), 'Source amendment must precede inference.'
source=RAW/'data_preprocessing/avidity_curves.xlsx'
rows=list(openpyxl.load_workbook(source,read_only=True,data_only=True).active.values)
assert len(rows)==183 and len(rows[0])==34
allrows=[];panels=[];parents=[];qc=[]
for tcr in range(1,8):
    panel=f'MSKCC_TCR{tcr}';block=(tcr-1)*5
    assert rows[2][block+1:block+4]==(0.01,1,100)
    with (RAW/f'mskcc_tcr{tcr}_ec50_sat_mut_af3_no_template.csv').open() as h:
        src=list(csv.DictReader(h))
    parent=src[0]['wt_seq'];assert len(parent)==9 and len(src)==180
    lookup={r['mutant']:r for r in src};assert len(lookup)==180
    wt=[];rr=[]
    for excel_row,r in enumerate(rows[3:],4):
        mutation=str(r[block]);m=re.fullmatch(r'([A-Z])(\d+)([A-Z])',mutation);assert m,mutation
        a,p,b=m[1],int(m[2]),m[3];assert parent[p-1]==a
        sequence=parent[:p-1]+b+parent[p:]
        s=lookup[mutation];assert s['sequence']==sequence
        v=[num(x) for x in r[block+1:block+4]]
        assert v[2] is not None
        if a==b:wt.append(v);continue
        rr.append(dict(panel_id=panel,candidate_id=f'{panel}__{mutation}',family_id='CMV_NLV' if tcr<=3 else 'gp100_IMD' if tcr<=6 else 'PDAC_GRL',
            sequence=sequence,parent_candidate_id=f'{panel}__parent',mutation=mutation,position_key=f'B:{p}',
            position=p,parent_aa=a,mutant_aa=b,endpoint='normalized_CD8_4_1BB_response_at_highest_tested_dose',
            value_oriented_higher_better=v[2],label_qc_pass=True,censoring='none_reported',source_id='Luksza2022_via_Visani2025_raw_xlsx',
            is_hit=int(v[2]>=.3),response_0_01=v[0],response_1=v[1],response_100=v[2],
            complete_three_dose=all(x is not None for x in v),log_dose_response_area=(v[0]+2*v[1]+v[2])/4 if all(x is not None for x in v) else None,
            source_excel_row=excel_row,source_excel_block=block,source_reported_fit_value=s['-log10(EC50)'],
            source_fit_is_reliable=s['is_reliable'],source_fit_used_as_primary=False))
    assert len(rr)==171 and len({x['sequence'] for x in rr})==171
    assert len(wt)==9 and all(v==wt[0] for v in wt)
    assert all(sum(r['position']==p for r in rr)==19 for p in range(1,10))
    rr.sort(key=lambda r:r['candidate_id']);allrows+=rr
    parents.append(dict(panel_id=panel,parent_candidate_id=f'{panel}__parent',sequence=parent,
        response_0_01=wt[0][0],response_1=wt[0][1],response_100=wt[0][2],
        repeated_WT_rows=9,independent_measurements_inferred=1))
    p=dict(panel_id=panel,family_id=rr[0]['family_id'],task_object='antigen_peptide_mutation',
        endpoint=rr[0]['endpoint'],units='source-normalized activation response; dose unit pending original SI verification',
        complete_candidate_pool=True,positive_only_selection=False,higher_is_better=True,parent_available=True,
        declared_hit_threshold=.3,hit_kind='fixed_threshold',hit_definition='Operational recognition retained: normalized highest-dose response>=0.3; source fitting heuristic, not a validated biological threshold',
        threshold_fixed_without_model_outputs=True,hit_is_above_parent_improvement=False,
        endpoint_is_proxy=True,mode1_eligible=True,mode2_eligible=True,mode3_eligible=True,
        biological_context=rr[0]['family_id'],candidate_count=171,parent_sequence=parent,
        modes2_3_initial_labels_n=18,mode3_initial_labels_n=20,mode2_extra_budget=30,mode3_extra_budget=40)
    panels.append(p)
    qc.append(dict(panel_id=panel,mutants=171,positions=9,substitutions_per_position=19,
        complete_highest_dose=171,complete_three_dose=sum(r['complete_three_dose'] for r in rr),
        hits=sum(r['is_hit'] for r in rr),WT_highest_dose=wt[0][2],duplicate_WT_rows=9,
        responses_above_one=sum(r['response_100']>1 for r in rr),no_response_candidates_retained=True))
csvout(OUT/'curated_candidates.csv',allrows);csvout(OUT/'parents.csv',parents)
dump(OUT/'panels.json',dict(panels=panels))
dump(OUT/'label_qc.json',dict(status='passed_source_only',source_sha256=sha(source),panels=qc,
    curves_normalized_not_probabilities=True,all_1197_mutants_retained=True,
    missing_low_dose_not_imputed=True,source_EC50_not_primary=True,
    source_EC50_issue='calc_obs_C converts log10(K) to ln(K), CSV later named -log10(EC50); sparse/outside-range fits and nonresponder sentinels unsuitable as exact affinity labels.',
    above_parent_improvement_not_claimed='Parent responses are near saturation and no experimental error supplied; small positive differences cannot establish improved potency.',
    independent_units='7 receptor contexts within3 antigen families; not7 independent antigen systems',
    predicted_outputs_seen=False))

ep=ROOT/'evaluation';protocol=ep/'multimode_protocol.json'
archive=ep/'multimode_protocol_pre_functional_endpoint_amendment.json'
assert not archive.exists();archive.write_bytes(protocol.read_bytes())
a=json.loads(protocol.read_text(encoding='utf8'))
a['status']='frozen_source_bound_before_new_inference'
a['functional_endpoint_amendment']={
 'timestamp_utc':datetime.now(timezone.utc).isoformat(),'new_predictions_seen':False,'raw_functional_labels_seen_for_source_QC':True,
 'reason':'Original three-dose curves available for all7 saturated peptide panels. Fitted CSV EC50 has log-base/coverage ambiguity; use complete observed highest-dose response for every candidate.',
 'primary_endpoint':'Source-normalized4-1BB-positive CD8 response at highest tested dose100 in original XLSX. Not affinity; normalization can exceed1.',
 'primary_hit':'response>=0.3, an operational source-fitting low-response heuristic; not a validated biological/clinical threshold; never above-parent improvement.',
 'secondary_top_decile':'Repeat all three modes against pool90th-percentile threshold as offline relative retrieval only, keeping identical selected sequences. Does not imply improved potency.',
 'secondary_three_dose_AUC':'Full three-dose log-concentration normalized trapezoid area, only complete source curves. Missing lower doses never imputed. Label-only cohort frozen before predictions; all three modes repeated, descriptive sensitivity only.',
 'mode1_scope':'Antigen mutation-site prioritization:9 positions with19 nonparent substitutions each, two selected sites=38 assays. Recognition-retaining variant yield, not proven affinity enhancement.',
 'replication':'7 receptor contexts and3 antigen groups; receptor1/2 are related and share very similar alpha chains.',
 'static_primary':'HERMES py000 fixed and ProteinMPNN v48_020 full_pep_masked, AF3_no_template published scores, both oriented higher pnlogp better. Parent mean across9 WT technical rows.',
 'static_sensitivity':'ProteinMPNN v48_002 full_pep_masked AF3_no_template is a predeclared additional source-only comparator setting; repeat all modes2/3, report both settings without selecting a winning setting.',
 'curated_labels_sha256':sha(OUT/'curated_candidates.csv'),'label_QC_sha256':sha(OUT/'label_qc.json')}
a['mode1_site_nomination']['scope']=a['functional_endpoint_amendment']['mode1_scope']
a['mode1_site_nomination']['success']='Primary joint-high two-site recognition-retaining yield must exceed both static strategies and mean random in every context. A companion win cannot replace it. This tests antigen experimental prioritization, not above-parent affinity improvement.'
a['mode1_site_nomination']['metrics']=['Number of recognition-retaining variants in equal38-assay budget','best and mean measured normalized response','fraction of all pool recognition-retaining variants recovered','secondary relative-top-decile recovery']
a['data_contract']['EC50']='Source fitted EC50 values archived only. Primary observed functional response source-bound by functional_endpoint_amendment.'
dump(protocol,a)
b=json.loads((ep/'protocol_freeze_binding.json').read_text())
b.update(current_protocol_sha256=sha(protocol),status=a['status'],data_eligibility_binding_pending=False,
    no_new_model_predictions_seen=True,new_functional_label_values_seen_for_QC=True,
    functional_endpoint_amendment_utc=a['functional_endpoint_amendment']['timestamp_utc'],
    labels_sha256=sha(OUT/'curated_candidates.csv'),label_QC_sha256=sha(OUT/'label_qc.json'))
dump(ep/'protocol_freeze_binding.json',b)
print(json.dumps(dict(status='frozen',panels=qc,protocol_sha256=sha(protocol)),ensure_ascii=False))
