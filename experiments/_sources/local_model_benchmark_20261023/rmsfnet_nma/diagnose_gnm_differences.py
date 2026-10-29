"""Diagnose GNM numerical portability; preserve every original/new output."""
import csv
import gzip
import json
import os
import sys
for name in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:
    os.environ[name]='2'
sys.dont_write_bytecode=True
import numpy as np
from scipy.linalg import eigvalsh
from cohort_common import ROOT, OLD, sha, record_key, write_json
from run_nma_scipy import build_matrices


def main():
    audit=json.loads((ROOT/'nma_scipy/full_reference_audit.json').read_text())
    failures={r['sample_id']:r for r in audit['failed_records']}
    refs={sample:[] for sample in failures}
    with gzip.open(OLD/'prody_common_evaluation/common_residue_predictions.csv.gz','rt',encoding='utf-8',newline='') as handle:
        for row in csv.DictReader(handle):
            if row['sample_id'] in refs:
                refs[row['sample_id']].append(row)
    result=[]
    for sample,reference_rows in refs.items():
        with np.load(ROOT/'nma_scipy/predictions'/(record_key(sample)+'.npz'),allow_pickle=False) as saved:
            coords=saved['ca_coords'].copy()
            predicted=saved['gnm'].copy()
            chain=saved['chain'].copy();author=saved['author_residue_id'].copy()
        hessian,kirchhoff=build_matrices(coords)
        values=eigvalsh(kirchhoff)
        nonzero=values[values>=1e-6]
        count20=min(20,len(coords)-1,3*len(coords)-6)
        lower,upper=float(nonzero[count20-1]),float(nonzero[count20])
        degenerate=np.isclose(nonzero,lower,atol=1e-10,rtol=1e-12)
        old=np.asarray([float(row['gnm']) for row in reference_rows])
        boundary_pairs=[]
        for i in range(len(coords)):
            deltas=coords[i+1:]-coords[i]
            dist2=(deltas**2).sum(1)
            for off in np.flatnonzero(np.abs(dist2-225.0)<1e-8):
                j=i+1+int(off)
                boundary_pairs.append(dict(i=i,j=j,chain_i=str(chain[i]),chain_j=str(chain[j]),
                    author_i=str(author[i]),author_j=str(author[j]),distance2=float(dist2[off])))
        difference=old**2-predicted**2
        item=dict(sample_id=sample,retained_nonzero_modes=count20,eigenvalue_at_cutoff=lower,next_eigenvalue=upper,
            eigenvalue_gap=upper-lower,boundary_eigenvalue_multiplicity=int(degenerate.sum()),
            gap_effectively_zero=bool(np.isclose(lower,upper,atol=1e-10,rtol=1e-12)),
            old_total_sq_fluct=float(np.sum(old**2)),new_total_sq_fluct=float(np.sum(predicted**2)),
            total_sq_fluct_difference=float(difference.sum()),near_15A_contact_pairs=boundary_pairs,
            gnm_maxabs=failures[sample]['gnm_maxabs'],gnm_positions_above_tolerance=failures[sample]['gnm_n_above_1e-6'])
        result.append(item)
    write_json(ROOT/'nma_scipy/gnm_difference_diagnosis.json',dict(records=result,script_sha256=sha(__file__),
        output_policy='Diagnostics only; no new output replaced or altered. Degenerate-mode basis is a hypothesis only unless cutoff gap and total variance support it.'))
    print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':
    main()
