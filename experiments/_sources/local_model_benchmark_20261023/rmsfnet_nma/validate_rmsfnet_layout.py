"""Compare float32 channels-last inference to the unchanged three-source baseline."""
import csv
import gzip
import json
import math
from pathlib import Path
import time
import numpy as np
from scipy.stats import spearmanr
from profile_rmsfnet import load_adapter
from cohort_common import ROOT, COHORT, sha, load_records, local_input, record_key, write_json


def average_precision(target, scores):
    order=np.argsort(-scores,kind='stable')
    truth=target[order].astype(int)
    ordered=scores[order]
    boundaries=np.r_[np.flatnonzero(ordered[1:]!=ordered[:-1]),len(scores)-1]
    tp=np.cumsum(truth)[boundaries]
    recall=tp/max(1,np.sum(truth))
    precision=tp/(boundaries+1)
    return float(np.sum(np.diff(np.r_[0.0,recall])*precision))


def metrics(target,scores):
    true_std=(target-target.mean())/target.std()
    pred_std=(scores-scores.mean())/scores.std()
    return dict(rho=float(spearmanr(target,scores).statistic),
        standardized_mae=float(np.mean(np.abs(true_std-pred_std))),
        ap90=average_precision(target>=np.quantile(target,0.9),scores))


def main():
    adapter=load_adapter()
    torch=adapter.torch
    rows,_=load_records(True)
    targets={r['sample_id']:[] for r in rows}
    with gzip.open(COHORT/'truth_residues.csv.gz','rt',encoding='utf-8',newline='') as handle:
        for row in csv.DictReader(handle):
            if row['sample_id'] in targets:
                assert row['valid_rmsf'].lower() in ('true','1','1.0')
                targets[row['sample_id']].append(float(row['target_rmsf']))
    output=ROOT/'rmsfnet/channels_last_validation'
    output.mkdir(parents=True,exist_ok=True)
    weights=Path('E:/Project/FlexPPI/tmp/third_party/RMSF-net/model_cpu/sim_model.pth')
    model=adapter.load_model(weights).to(memory_format=torch.channels_last_3d)
    audit=dict(threads=2,batch_size=2,dtype='float32',layout='channels_last_3d',
        original_model_unmodified=True,weights_sha256=sha(weights),script_sha256=sha(__file__),records=[])
    # This forward pre-hook changes memory strides only, retaining the original model and adapter mathematics.
    hook=model.register_forward_pre_hook(lambda module,args: (args[0].contiguous(memory_format=torch.channels_last_3d),))
    for record in rows:
        sample=record['sample_id']
        key=record_key(sample)
        tick=time.perf_counter()
        residues,values,timing=adapter.run_sample(model,str(local_input(record)),2)
        elapsed=time.perf_counter()-tick
        values=np.asarray(values,dtype=np.float64)
        with np.load(ROOT/'rmsfnet/smoke'/(key+'.npz'),allow_pickle=False) as original:
            old=original['score'].copy()
            chain=np.asarray([r.chain for r in residues])
            author=np.asarray([str(int(r.number))+r.insertion for r in residues])
            assert np.array_equal(chain,original['chain']) and np.array_equal(author,original['author_residue_id'])
        true=np.asarray(targets[sample],dtype=np.float64)
        assert true.shape==values.shape==old.shape
        before,after=metrics(true,old),metrics(true,values)
        item=dict(sample_id=sample,dataset=record['cohort'],residues=len(values),
            score_maxabs=float(np.max(np.abs(old-values))),score_meanabs=float(np.mean(np.abs(old-values))),
            original_metrics=before,layout_metrics=after,
            metric_difference={k:after[k]-before[k] for k in before},total_seconds=elapsed,**timing)
        np.savez_compressed(output/(key+'.npz'),score=values,chain=chain,author_residue_id=author)
        item['output_sha256']=sha(output/(key+'.npz'))
        audit['records'].append(item)
        write_json(output/'equivalence_audit.json',audit)
        print(json.dumps(item),flush=True)
    hook.remove()
    audit['all_residues_mapped']=True
    audit['largest_rho_difference']=max(abs(r['metric_difference']['rho']) for r in audit['records'])
    audit['largest_standardized_mae_difference']=max(abs(r['metric_difference']['standardized_mae']) for r in audit['records'])
    audit['largest_ap90_difference']=max(abs(r['metric_difference']['ap90']) for r in audit['records'])
    audit['float32_equivalence_pass']=bool(max(r['score_maxabs'] for r in audit['records'])<1e-4 and audit['largest_rho_difference']<1e-5)
    write_json(output/'equivalence_audit.json',audit)
    assert audit['float32_equivalence_pass']


if __name__=='__main__':
    main()
