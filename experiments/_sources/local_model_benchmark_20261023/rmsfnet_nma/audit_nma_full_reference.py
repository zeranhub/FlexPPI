"""Streaming all-position comparison; one record at a time to conserve memory."""
import csv
import gzip
from itertools import groupby
import json
import os
from pathlib import Path
import sys
import time
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:
    os.environ[key]='1'
sys.dont_write_bytecode=True
import numpy as np
from cohort_common import ROOT, OLD, sha, record_key, write_json


def main():
    directory=ROOT/'nma_scipy/predictions'
    summary=json.loads((directory/'summary.json').read_text())
    assert summary['complete']==3063 and summary['failed']==0
    reference=OLD/'prody_common_evaluation/common_residue_predictions.csv.gz'
    audit_rows=[]
    tick=time.perf_counter()
    with gzip.open(reference,'rt',encoding='utf-8',newline='') as handle:
        for sample,group in groupby(csv.DictReader(handle),lambda row:row['sample_id']):
            rows=list(group)
            key=record_key(sample)
            path=directory/(key+'.npz')
            status=json.loads((directory/(key+'.json')).read_text())
            assert status['status']=='complete' and sha(path)==status['output_sha256']
            with np.load(path,allow_pickle=False) as saved:
                assert np.array_equal(saved['residue_index'],np.asarray([int(r['residue_index']) for r in rows]))
                assert saved['chain'].tolist()==[r['chain'] for r in rows]
                assert saved['author_residue_id'].tolist()==[r['author_residue_id'] for r in rows]
                item=dict(sample_id=sample,dataset=rows[0]['dataset'],residues=len(rows))
                for method in ['anm','gnm']:
                    values=np.asarray([float(r[method]) for r in rows],dtype=np.float64)
                    difference=np.abs(saved[method]-values)
                    item[method+'_maxabs']=float(difference.max())
                    item[method+'_meanabs']=float(difference.mean())
                    item[method+'_n_above_1e-6']=int(np.sum(difference>1e-6))
            item['pass_maxabs_1e-6']=item['anm_maxabs']<1e-6 and item['gnm_maxabs']<1e-6
            audit_rows.append(item)
    assert len(audit_rows)==3063 and sum(r['residues'] for r in audit_rows)==961321
    output=ROOT/'nma_scipy'
    with (output/'full_reference_comparison.csv').open('w',encoding='utf-8',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(audit_rows[0]));writer.writeheader();writer.writerows(audit_rows)
    audit=dict(records=len(audit_rows),residues=sum(r['residues'] for r in audit_rows),
        reference_path=str(reference),reference_sha256=sha(reference),
        tolerance_maxabs=1e-6,all_pass=all(r['pass_maxabs_1e-6'] for r in audit_rows),
        anm_maxabs=max(r['anm_maxabs'] for r in audit_rows),gnm_maxabs=max(r['gnm_maxabs'] for r in audit_rows),
        anm_positions_above_tolerance=sum(r['anm_n_above_1e-6'] for r in audit_rows),
        gnm_positions_above_tolerance=sum(r['gnm_n_above_1e-6'] for r in audit_rows),
        failed_records=[r for r in audit_rows if not r['pass_maxabs_1e-6']],
        audit_seconds=time.perf_counter()-tick,script_sha256=sha(__file__),
        purpose='Independent audit only: fresh local outputs were not replaced, fit, scaled or recalibrated using the reference.')
    write_json(output/'full_reference_audit.json',audit)
    print(json.dumps(audit,indent=2),flush=True)
    assert audit['all_pass'], 'Some records differ beyond float64 numeric tolerance; inspect audit before treating implementations as equivalent'


if __name__=='__main__':
    main()
