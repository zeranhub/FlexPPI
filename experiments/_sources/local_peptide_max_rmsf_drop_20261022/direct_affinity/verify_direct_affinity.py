"""Independent standard-library numerical verification of exported scientific results."""
from pathlib import Path
from collections import defaultdict
import csv
import json
import math
import hashlib
from statistics import mean, median

ROOT=Path(__file__).resolve().parent
OLD=ROOT.parents[1]/'local_antigen_screening_20261021'
def rows(p):
    with p.open(encoding='utf-8-sig',newline='') as f:
        return list(csv.DictReader(f))
def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()
def rank(x):
    order=sorted(range(len(x)),key=x.__getitem__)
    result=[0.]*len(x)
    start=0
    while start<len(order):
        end=start+1
        while end<len(order) and x[order[end]]==x[order[start]]:
            end+=1
        value=(start+1+end)/2
        for j in order[start:end]:
            result[j]=value
        start=end
    return result
def rho(x,y):
    a,b=rank(x),rank(y)
    am,bm=mean(a),mean(b)
    return sum((u-am)*(v-bm) for u,v in zip(a,b))/math.sqrt(sum((u-am)**2 for u in a)*sum((v-bm)**2 for v in b))
sources=json.loads((ROOT/'source_only_protocol.json').read_text(encoding='utf-8'))['source_files_sha256']
assert all(digest(Path(p))==h for p,h in sources.items())
export=rows(ROOT/'candidate_scores_all_templates.csv')
manifest=json.loads((OLD/'structures/fixed_backbone_manifest.json').read_text())
meta={r['candidate_id']:r for r in manifest['candidates']}
verified_scores={}
max_error=0.
for r in export:
    cid=r['inference_candidate_id'];m=meta[cid]
    p=OLD/'predictions/fixed_backbone'/cid/'prediction.csv'
    parent=OLD/'predictions/fixed_backbone'/m['parent_inference_candidate_id']/'prediction.csv'
    pr=[a for a in rows(p) if a['chain']==m['peptide_chain']]
    wr=[a for a in rows(parent) if a['chain']==m['peptide_chain']]
    assert len(pr)==len(wr)==9
    assert ''.join(a['amino_acid'] for a in pr)==m['peptide_sequence']
    delta=[float(w['rmsf'])-float(a['rmsf']) for a,w in zip(pr,wr)]
    score=max(delta)
    assert abs(score-float(r['score_max_rmsf_drop_A']))<1e-12
    for i,(a,w,d) in enumerate(zip(pr,wr,delta),1):
        for expected,key in [(float(a['rmsf']),f'mutant_RMSF_P{i}_A'),(float(w['rmsf']),f'parent_RMSF_P{i}_A'),(d,f'drop_RMSF_P{i}_A')]:
            error=abs(expected-float(r[key]));max_error=max(max_error,error)
            assert error<1e-12
    verified_scores[cid]=score

candidate_labels={};study_labels={}
for version in ['reported_hybrid','all_constrained','all_fitted']:
    cand={r['candidate_id']:float(r['delta_pKD_vs_parent']) for r in rows(OLD/f'experimental_labels/SPR_labels_candidates_{version}.csv')}
    obs=rows(OLD/f'experimental_labels/SPR_labels_observations_{version}.csv')
    bycid=defaultdict(list)
    for r in obs:
        pk=-math.log10(float(r['KD_uM'])*1e-6)
        parent=-math.log10(float(r['parent_KD_uM'])*1e-6)
        assert abs(pk-float(r['pKD_M']))<1e-12
        difference=pk-parent
        assert abs(difference-float(r['delta_pKD_vs_parent']))<1e-12
        study_labels[version,r['study_id'],r['candidate_id']]=difference
        bycid[r['candidate_id']].append(difference)
    for cid,value in cand.items():
        assert abs(value-median(bycid[cid]))<1e-12
        candidate_labels[version,cid]=value

checked=0;max_rho_error=0.
for r in rows(ROOT/'spearman_and_pair_order_all_sensitivities.csv'):
    group=[]
    for score in export:
        if int(score['n_mutations'])!=1:
            continue
        if r['analysis_level']=='pooled_primary_templates_descriptive':
            include=score['primary_template']=='True'
        else:
            include=(score['panel_id']==r['panel_id'] and score['template']==r['template'])
        if not include:
            continue
        if r['analysis_level']=='within_original_study':
            key=(r['SPR_variant'],r['study'],score['candidate_id'])
            if key not in study_labels:
                continue
            label=study_labels[key]
        else:
            label=candidate_labels[r['SPR_variant'],score['candidate_id']]
        group.append((verified_scores[score['inference_candidate_id']],label))
    assert len(group)==int(r['n_unique_single_mutants'])
    x,y=zip(*group)
    result=rho(x,y)
    error=abs(result-float(r['spearman_rho']));max_rho_error=max(max_rho_error,error)
    assert error<1e-12
    pairs=[(math.copysign(1,a-c) if a!=c else 0,math.copysign(1,b-d)) for i,(a,b) in enumerate(group) for c,d in group[:i] if b!=d]
    correct=sum(float(u*v>0)+.5*float(u==0) for u,v in pairs)/len(pairs)
    assert abs(correct-float(r['pair_order_agreement_with_half_credit_ties']))<1e-12
    checked+=1
assert checked==39
assert len(rows(ROOT/'primary_26_single_mutants_ranked.csv'))==26
assert len(rows(ROOT/'per_site_rmsf_profiles_and_drops.csv'))==684
report={'status':'passed','implementation':'Independent Python csv/math/statistics; no pandas/scipy or analyzer helper imports',
  'checked_profiles':len(export),'checked_site_values':76*9*3,
  'checked_Spearman_and_pair_order_groups':checked,'maximum_profile_absolute_error':max_error,
  'maximum_Spearman_absolute_error':max_rho_error,
  'original_KD_pKD_parent_relative_labels_recomputed':True,
  'source_protocol_sha256':digest(ROOT/'source_only_protocol.json')}
(ROOT/'independent_numeric_verification.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(report,ensure_ascii=False,indent=2))
