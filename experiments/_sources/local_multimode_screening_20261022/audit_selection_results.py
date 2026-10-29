"""Recompute reported yields from archived sequence IDs and source labels."""
from pathlib import Path
from collections import defaultdict
from datetime import datetime,timezone
import csv,json,hashlib
import numpy as np

ROOT=Path(__file__).resolve().parent
read=lambda p:list(csv.DictReader(p.open(encoding='utf-8-sig',newline='')))
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    checks=[];total=0
    for setting,folder in [('primary','results_primary'),('top_decile','results_top_decile'),('proteinmpnn002','results_proteinmpnn002'),('three_dose_AUC','results_three_dose_AUC'),('primary','results_controls'),('primary','results_stronger_baselines')]:
        contracts=json.loads((ROOT/'features'/f'{setting}__contracts.json').read_text())['panels']
        pp={p['panel_id']:p for p in contracts};lookup={};hitlookup={};panel_ids={}
        for p in contracts:
            rr=read(Path(p['labels_csv']));values={r['candidate_id']:float(r['value_oriented_higher_better']) for r in rr}
            assert len(values)==len(rr)
            cutoff=float(p['declared_hit_threshold']) if p['hit_kind']=='fixed_threshold' else float(np.quantile(list(values.values()),.9))
            lookup.update(values);hitlookup.update({k:v>=cutoff for k,v in values.items()});panel_ids[p['panel_id']]=set(values)
        rows=read(ROOT/'evaluation'/folder/'all_candidate_budget_metrics.csv')
        initials={}
        for r in rows:
            selected=json.loads(r['selected_candidate_ids']);initial=json.loads(r['initial_candidate_ids'])
            assert len(selected)==len(set(selected))==int(r['extra_assays'])
            assert len(initial)==len(set(initial))==int(r['initial_n'])
            assert not set(selected)&set(initial) and set(selected+initial)<=panel_ids[r['panel_id']]
            key=r['panel_id'],r['mode'],r['seed']
            if key in initials:assert initial==initials[key]
            initials[key]=initial
            assert int(r['new_hits'])==sum(hitlookup[c] for c in selected)
            assert int(r['initial_hits'])==sum(hitlookup[c] for c in initial)
            assert int(r['cumulative_hits'])==sum(hitlookup[c] for c in initial+selected)
            assert float(r['cumulative_best'])==max(lookup[c] for c in initial+selected)
            if selected:assert float(r['best_extra'])==max(lookup[c] for c in selected)
        sites=read(ROOT/'evaluation'/folder/'all_site_nomination_metrics.csv')
        for r in sites:
            ids=json.loads(r['selected_candidate_ids']);positions=json.loads(r['selected_positions'])
            assert len(ids)==len(set(ids))==int(r['assay_budget'])==len(positions)*19
            assert int(r['new_hits'])==sum(hitlookup[c] for c in ids)
            assert float(r['best_effect'])==max(lookup[c] for c in ids)
            assert abs(float(r['mean_effect'])-np.mean([lookup[c] for c in ids]))<1e-12
        total+=len(rows)+len(sites)
        checks.append(dict(setting=setting,folder=folder,candidate_rows=len(rows),site_rows=len(sites),all_yields_independently_recomputed=True))
    # Threshold sensitivity must never change any selection, since no hit flags
    # are passed to policies and continuous observed responses are unchanged.
    for i in range(1,8):
        pid=f'MSKCC_TCR{i}'
        p=read(ROOT/'evaluation/results_primary'/pid/'selected_candidates.csv')
        t=read(ROOT/'evaluation/results_top_decile'/pid/'selected_candidates.csv')
        assert p==t
        p=read(ROOT/'evaluation/results_primary'/pid/'site_nomination_metrics.csv')
        t=read(ROOT/'evaluation/results_top_decile'/pid/'site_nomination_metrics.csv')
        assert len(p)==len(t)
        assert all(a['selected_candidate_ids']==b['selected_candidate_ids'] for a,b in zip(p,t))
    # Feedback batches can depend only on the prior measured set. Verify the
    # trace itself has no repeats and every candidate becomes visible later.
    for folder in ['results_primary','results_controls','results_stronger_baselines']:
        for i in range(1,8):
            rr=read(ROOT/'evaluation'/folder/f'MSKCC_TCR{i}'/'feedback_information_trace.csv')
            history={}
            for r in rr:
                key=r['strategy'],r['seed'];before=json.loads(r['observed_candidate_ids_before']);batch=json.loads(r['selected_batch'])
                assert len(batch)==5 and len(set(batch))==5 and not set(batch)&set(before)
                if key in history:assert before==history[key]
                history[key]=before+batch
    for i in range(1,8):
        pid=f'MSKCC_TCR{i}'
        original=read(ROOT/'evaluation/results_primary'/pid/'selected_candidates.csv')
        stronger=read(ROOT/'evaluation/results_stronger_baselines'/pid/'selected_candidates.csv')
        existing={(r['mode'],r['strategy'],r['seed']):r for r in stronger}
        assert all(existing[r['mode'],r['strategy'],r['seed']]==r for r in original)
    out=dict(status='passed_independent_recomputation',timestamp_utc=datetime.now(timezone.utc).isoformat(),
        checks=checks,total_metric_rows_recomputed=total,threshold_sensitivity_selection_identical=True,
        all_initial_candidate_ids_and_labels_matched=True,all_feedback_traces_no_within_batch_updates=True,
        auditor_sha256=sha(Path(__file__)),not_a_test_of_prospective_wet_lab_efficacy=True)
    (ROOT/'evaluation/independent_result_audit.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf8')
    print(json.dumps(out))

if __name__=='__main__':main()
