"""Resolve complete v2 provenance and publish transparent portability/metric audit."""
import csv
import gzip
from itertools import groupby
import json
import os
import sys
for name in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:
    os.environ[name]='1'
sys.dont_write_bytecode=True
import numpy as np
from cohort_common import ROOT, OLD, sha, load_records, record_key, write_json


def read_csv(path):
    with path.open(encoding='utf-8',newline='') as handle:
        return list(csv.DictReader(handle))


def write_csv(path,rows):
    with path.open('w',encoding='utf-8',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)


def main():
    output=ROOT/'gnm_ckdtree_v2'
    manifest=json.loads((output/'overlay_manifest.json').read_text())
    overlays={r['sample_id']:r for r in manifest['overlays']}
    v1=ROOT/'nma_scipy/predictions'
    records,_=load_records(False)
    resolved=[]
    for record in records:
        sample=record['sample_id'];key=record_key(sample)
        base_path=v1/(key+'.npz');base_status_path=v1/(key+'.json')
        base_status=json.loads(base_status_path.read_text())
        path=output/'predictions'/(key+'.npz') if sample in overlays else base_path
        status_path=output/'predictions'/(key+'.json') if sample in overlays else base_status_path
        status=json.loads(status_path.read_text())
        assert status['input_sha256']==record['input_sha256'] and status['status']=='complete'
        if sample in overlays:
            assert status['source_v1_sha256']==base_status['output_sha256']
            assert sha(path)==status['output_sha256']
        resolved.append(dict(sample_id=sample,dataset=record['cohort'],residues=record['n_residues'],
            input_sha256=record['input_sha256'],gnm_version='ckdtree_float32_v2' if sample in overlays else 'fresh_local_v1_unchanged_graph',
            gnm_path=str(path),gnm_array_key='gnm',gnm_status_path=str(status_path),gnm_output_sha256=status['output_sha256'],
            anm_path=str(base_path),anm_array_key='anm',anm_status_path=str(base_status_path),anm_output_sha256=base_status['output_sha256']))
    write_csv(output/'resolved_prediction_index.csv',resolved)
    compare=read_csv(ROOT/'nma_scipy/full_reference_comparison.csv')
    comparison={r['sample_id']:r for r in compare}
    with gzip.open(OLD/'prody_common_evaluation/common_residue_predictions.csv.gz','rt',encoding='utf-8',newline='') as handle:
        remaining=set(overlays)
        for sample,group in groupby(csv.DictReader(handle),lambda row:row['sample_id']):
            if sample not in remaining:
                continue
            refs=list(group)
            with np.load(overlays[sample]['path'],allow_pickle=False) as saved:
                old=np.asarray([float(r['gnm']) for r in refs]);new=saved['gnm']
                assert saved['chain'].tolist()==[r['chain'] for r in refs]
                assert saved['author_residue_id'].tolist()==[r['author_residue_id'] for r in refs]
                diff=np.abs(new-old)
                comparison[sample].update(gnm_maxabs=float(diff.max()),gnm_meanabs=float(diff.mean()),
                    **{'gnm_n_above_1e-6':int(np.sum(diff>1e-6))})
            remaining.remove(sample)
            if not remaining:
                break
        assert not remaining
    for row in compare:
        row['pass_maxabs_1e-6']=float(row['anm_maxabs'])<1e-6 and float(row['gnm_maxabs'])<1e-6
    write_csv(output/'full_reference_comparison_v2.csv',compare)
    diagnoses=json.loads((output/'difference_diagnosis_and_metric_impact.json').read_text())['records']
    impact=[]
    for row in diagnoses:
        item=dict(sample_id=row['sample_id'],dataset=row['dataset'],residues=row['residues'],
            changed_adjacency_edges=row['changed_edge_count'],edges_float64=row['edges_float64'],edges_float32=row['edges_ckdtree_float32'],
            eigenvalue20=row['eigenvalue20'],eigenvalue21=row['eigenvalue21'],spectrum_gap=row['spectrum_gap'],
            zero_modes_v1=row['zero_modes_v1'],zero_modes_float32=row['zero_modes_ckdtree_float32'],degenerate_cutoff=row['degenerate_cutoff'])
        for label,metrics in [('old',row['old_metrics']),('v1',row['v1_metrics']),('v2',row['v2_metrics']),('v2_minus_old',row['v2_minus_old_metrics']),('v2_minus_v1',row['v2_minus_v1_metrics'])]:
            item.update({label+'_'+key:value for key,value in metrics.items()})
        impact.append(item)
    write_csv(output/'rho_ap_metric_impact.csv',impact)
    audit=dict(records=len(compare),residues=sum(int(r['residues']) for r in compare),
        anm_maxabs=max(float(r['anm_maxabs']) for r in compare),gnm_maxabs=max(float(r['gnm_maxabs']) for r in compare),
        anm_positions_above_1e_6=sum(int(r['anm_n_above_1e-6']) for r in compare),
        gnm_positions_above_1e_6=sum(int(r['gnm_n_above_1e-6']) for r in compare),
        corrected_adjacency_records=list(overlays),remaining_different_records=[r['sample_id'] for r in compare if not r['pass_maxabs_1e-6']],
        interpretation='ANM passes all-position tolerance. Corrected default GNM CKDTree float32 adjacency; nine remaining differences are graph-identical degenerate20-mode truncation basis choices, retained as fresh local outputs. Do not label GNM fully pointwise equivalent.',
        largest_absolute_rho_delta_among_diagnosed_records=max(abs(r['v2_minus_old_rho']) for r in impact),
        largest_absolute_ap90_delta_among_diagnosed_records=max(abs(r['v2_minus_old_ap90']) for r in impact),
        impact_scope='Metrics compare the ten diagnosed records only; the uniform benchmark independently scores all3063 records. No reference score chose a prediction.',
        index_path=str(output/'resolved_prediction_index.csv'),index_sha256=sha(output/'resolved_prediction_index.csv'),
        overlay_manifest_sha256=sha(output/'overlay_manifest.json'),script_sha256=sha(__file__))
    write_json(output/'final_reference_audit_v2.json',audit)
    print(json.dumps(audit,indent=2),flush=True)


if __name__=='__main__':
    main()
