"""Mean-metric effect attributable to the nine retained degenerate GNM records.

Uses existing per-record diagnostics only; no model or eigenproblem is rerun.
The old/new difference elsewhere below numeric tolerance is outside this
attribution, hence it must not be called an exact full-cohort old/new rescore.
"""
import csv
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parent
BENCH=ROOT.parent


def read(path):
    with path.open(encoding='utf-8-sig',newline='') as handle:
        return list(csv.DictReader(handle))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    cohort=read(BENCH/'cohort/records.csv')
    group_by_sample={r['sample_id']:r['main_group_candidate'] for r in cohort}
    impacts=read(ROOT/'gnm_ckdtree_v2/rho_ap_metric_impact.csv')
    nine=[r for r in impacts if r['degenerate_cutoff']=='True']
    assert len(nine)==9
    assert all(int(r['changed_adjacency_edges'])==0 for r in nine)
    group_sets={g:{r['sample_id'] for r in cohort if r['main_group_candidate']==g} for g in ['single_chain','other_multichain','tcr_pmhc_itrace']}
    group_sets['all_multichain']=group_sets['other_multichain']|group_sets['tcr_pmhc_itrace']
    group_sets['all']=set(group_by_sample)
    summary=[]
    for group,ids in group_sets.items():
        changed=[r for r in nine if r['sample_id'] in ids]
        result=dict(group=group,n_records=len(ids),n_degenerate_records=len(changed),
            affected_sample_ids=';'.join(r['sample_id'] for r in changed))
        for old_key,metric in [('rho','spearman'),('ap90','top10_ap'),('standardized_mae','z_mae')]:
            contribution=sum(float(r['v2_minus_old_'+old_key]) for r in changed)
            result[metric+'_delta_sum']=contribution
            result[metric+'_equal_record_mean_delta']=contribution/len(ids)
        summary.append(result)
    output=ROOT/'gnm_ckdtree_v2'
    with (output/'degenerate_group_mean_impact.csv').open('w',encoding='utf-8',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(summary[0]));writer.writeheader();writer.writerows(summary)
    audit=dict(groups=summary,diagnosed_records=9,
        comparison='Contribution from retained nine graph-identical degenerate-cutoff records: fresh local v2 minus frozen original ProDy values, using equal record means.',
        attribution_limit='Uses only the nine diagnosed records, assigning zero contribution to other records. This is not an exact full-cohort old/new AP rescore; rounding/ties in other records are not evaluated here.',
        integration_review='Read-only review found v2 overlay routing and SHA checks correct: ANM always base v1; GNM changed-graph sample overlay; all other GNM base v1. Full original v1 remains untouched.',
        sources={str(path):sha(path) for path in [BENCH/'cohort/records.csv',output/'rho_ap_metric_impact.csv',BENCH/'evaluate_benchmark.py']},
        script_sha256=sha(Path(__file__)))
    (output/'degenerate_group_mean_impact.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
    print(json.dumps(audit,indent=2),flush=True)


if __name__=='__main__':
    main()
