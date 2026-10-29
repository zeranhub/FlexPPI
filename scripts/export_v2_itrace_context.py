"""Fixed-coordinate ITrace context ablations with the parent system's OOF model."""
from __future__ import annotations
import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch

from export_v2_oof_predictions import (
    ROOT, METRICS, load_fold, prepare_input, load_cache, save_cache, cache_name,
    sha256, digest_json, canonical_json, atomic_json, implementation_hash,
    collate_states, load_state_example)
from analyze_v2_oof import regression

TENSOR_KEYS=('aatype','residue_index','chain_index','backbone','backbone_mask',
             'geometry_mask','target_mask','label_values','label_variances','label_mask')


def subset_example(example, indices):
    """Keep complete chains; their within-chain positions and coordinates stay unchanged."""
    output=dict(example)
    for key in TENSOR_KEYS:
        output[key]=example[key][indices].clone()
    ordered={int(chain):i for i,chain in enumerate(dict.fromkeys(output['chain_index'].tolist()))}
    output['chain_index']=torch.tensor([ordered[int(c)] for c in output['chain_index']],dtype=torch.long)
    return output


def specs(chains, roles):
    if set(chains)!=set(roles) or set(roles.values())!={'mhc','b2m','peptide','tcr_alpha','tcr_beta'}:
        raise ValueError('Unexpected ITrace chains/roles')
    return [('pmhc',[c for c in chains if roles[c] in {'mhc','b2m','peptide'}]),
            ('tcr',[c for c in chains if roles[c] in {'tcr_alpha','tcr_beta'}])]+[
                ('isolated_'+c,[c]) for c in chains]


def context_arrays(model, example, indices, full, identity, condition):
    subset=subset_example(example,indices)
    with torch.inference_mode():
        output=model.forward_state(collate_states([subset]))
    arrays={k:v[indices].copy() for k,v in full.items()
            if k not in ['metadata_json','metric_names']}
    arrays['metric_names']=full['metric_names'].copy()
    arrays['mean']=output['mean'][0].numpy()
    arrays['model_variance']=output['log_variance'][0].exp().numpy()
    arrays['output_mask']=output['mask'][0].numpy().astype(bool)
    arrays['valid']=arrays['label_mask'] & arrays['target_mask'][:,None]
    arrays['valid'] &= arrays['geometry_mask'][:,None] & arrays['output_mask'][:,None]
    if not np.array_equal(arrays['valid'],full['valid'][indices]):
        raise ValueError('Context changed valid target coverage')
    metadata=json.loads(str(full['metadata_json'].item()))
    metadata.update(identity=identity,context=condition,residue_count=len(indices))
    metadata['protocol']={**metadata['protocol'],'context':condition,
                          'intervention':'delete entire chains at fixed coordinates; not apo/holo MD'}
    arrays['metadata_json']=np.asarray(canonical_json(metadata))
    return arrays


def summarize(frame, out):
    rows=[]
    for (sid,fold,chain,role,condition),part in frame.groupby(
            ['sample_id','fold','chain','chain_role','context'],sort=False):
        for metric in METRICS:
            s=part.loc[part['valid_'+metric]]
            if s.empty:continue
            rows.append(dict(sample_id=sid,fold=fold,chain=chain,chain_role=role,context=condition,
                             metric=metric,n_residues=len(s),**regression(s['pred_'+metric],s['target_'+metric])))
    perchain=pd.DataFrame(rows);perchain.to_csv(out/'per_chain_context_metrics.csv',index=False)
    # ITrace roles map to one chain per complex; aggregate explicitly before resampling.
    stats=['mae','rmse','bias','pearson','spearman']
    system=perchain.groupby(['sample_id','fold','chain_role','context','metric'])[stats].mean().reset_index()
    system.groupby(['fold','chain_role','context','metric'])[stats].mean().to_csv(out/'fold_context_metrics.csv')
    summaries=[];rng=np.random.default_rng(2026)
    for keys,part in system.groupby(['chain_role','context','metric']):
        for stat in stats:
            values=part[stat].dropna().to_numpy(float)
            if not len(values):continue
            boot=rng.choice(values,size=(2000,len(values)),replace=True).mean(1)
            summaries.append(dict(chain_role=keys[0],context=keys[1],metric=keys[2],statistic=stat,
                                   n_complexes=len(values),undefined_complexes=int(part[stat].isna().sum()),
                                   mean=float(values.mean()),ci_low=float(np.quantile(boot,.025)),
                                   ci_high=float(np.quantile(boot,.975))))
    pd.DataFrame(summaries).to_csv(out/'context_summary.csv',index=False)
    paired=[]
    for (role,metric),part in system.groupby(['chain_role','metric']):
        for ref in ['partner_removed','isolated']:
            for stat in stats:
                pivot=part.pivot(index='sample_id',columns='context',values=stat)
                pair=pivot[['full',ref]].dropna();delta=(pair.full-pair[ref]).to_numpy()
                if not len(delta):continue
                boot=rng.choice(delta,size=(2000,len(delta)),replace=True).mean(1)
                paired.append(dict(chain_role=role,metric=metric,statistic=stat,comparison='full_minus_'+ref,
                                    n_complexes=len(delta),mean_difference=float(delta.mean()),
                                    ci_low=float(np.quantile(boot,.025)),ci_high=float(np.quantile(boot,.975))))
    pd.DataFrame(paired).to_csv(out/'paired_context_differences.csv',index=False)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input-dir',type=Path,default=ROOT/'artifacts/evaluation_v2/oof')
    p.add_argument('--output-dir',type=Path,default=ROOT/'artifacts/evaluation_v2/itrace_context')
    p.add_argument('--limit',type=int);p.add_argument('--threads',type=int,default=6)
    args=p.parse_args()
    if args.limit is not None and args.limit<1:p.error('limit must be positive')
    torch.set_num_threads(args.threads)
    run=json.loads((args.input_dir/'run_summary.json').read_text())
    if run['status']!='complete':raise ValueError('Full OOF export must finish first')
    entries=[json.loads(x) for x in (args.input_dir/'prediction_manifest.jsonl').read_text().splitlines()]
    entries=[r for r in entries if r['dataset']=='itrace']
    if len(entries)!=244:raise ValueError('Expected 244 ITrace parents')
    roles_path=ROOT/'output/itrace_benchmark/raw/itrace_chain_roles.json'
    roles=json.loads(roles_path.read_text());roles_hash=sha256(roles_path)
    implementation=digest_json([implementation_hash(),sha256(Path(__file__))])
    selected=sorted(entries,key=lambda r:(r['fold'],r['sample_id']))
    if args.limit is not None:selected=selected[:args.limit]
    out=args.output_dir;out.mkdir(parents=True,exist_ok=True)
    status=dict(status='running',expected_parents=244,processed_parents=0,conditions=0,
                partial=args.limit is not None,full_oof_run_sha256=sha256(args.input_dir/'run_summary.json'),
                semantics='fixed-coordinate input ablation, not experimental apo/holo changes')
    atomic_json(out/'run_summary.json',status)
    frames=[];model=None;loaded_fold=None;start=time.monotonic()
    try:
        for row in selected:
            fold=row['fold'];sid=row['sample_id'];pdb=sid.removeprefix('itrace_').lower()
            if loaded_fold!=fold:
                model,config,_,checkpoint_sha=load_fold(ROOT,fold,torch.device('cpu'),True);loaded_fold=fold
            if checkpoint_sha!=row['identity']['checkpoint_sha256']:raise ValueError('Parent checkpoint changed')
            full=load_cache(Path(row['cache_path']),row['identity'])
            prepared=prepare_input(row['record'],config['data'].get('path_mappings',{}))
            if sha256(prepared[3])!=row['identity']['input_sha256']:raise ValueError('Parent input changed')
            example=load_state_example(prepared[0]);chain_roles=roles[pdb]
            conditions=[('full',list(dict.fromkeys(full['chain'].tolist())))]+specs(
                list(dict.fromkeys(full['chain'].tolist())),chain_roles)
            for condition,chains in conditions:
                indices=np.flatnonzero(np.isin(full['chain'],chains))
                if condition=='full':arrays=full
                else:
                    identity=dict(parent_identity=row['identity'],context=condition,chains=chains,
                                  selected_full_indices_sha256=digest_json(indices.tolist()),
                                  implementation_sha256=implementation,roles_sha256=roles_hash)
                    path=out/f'fold_{fold}'/cache_name(sid+'__'+condition)
                    if path.exists():arrays=load_cache(path,identity)
                    else:
                        arrays=context_arrays(model,example,indices,full,identity,condition)
                        save_cache(path,arrays,identity)
                context='full' if condition=='full' else ('isolated' if condition.startswith('isolated_') else 'partner_removed')
                frame=pd.DataFrame(dict(sample_id=sid,pdb_id=pdb,fold=fold,condition=condition,context=context,
                                        chain=arrays['chain'],author_residue_id=arrays['author_residue_id'],
                                        amino_acid=arrays['amino_acid']))
                frame['chain_role']=frame.chain.map(chain_roles)
                for i,metric in enumerate(METRICS):
                    frame['pred_'+metric]=arrays['mean'][:,i];frame['target_'+metric]=arrays['target'][:,i]
                    frame['valid_'+metric]=arrays['valid'][:,i]
                frames.append(frame);status['conditions']+=1
            status['processed_parents']+=1;status['elapsed_seconds']=time.monotonic()-start
            atomic_json(out/'run_summary.json',status)
            if status['processed_parents']%10==0:print('Context parents',status['processed_parents'],'elapsed',round(status['elapsed_seconds']),flush=True)
        data=pd.concat(frames,ignore_index=True)
        if data.duplicated(['sample_id','chain','author_residue_id','context']).any():raise ValueError('Duplicate context rows')
        if not args.limit:
            for context in ['full','partner_removed','isolated']:
                s=data.loc[data.context==context]
                if (int(s.valid_rmsf.sum()),int(s.valid_chi_entropy.sum()))!=(199668,168607):
                    raise ValueError('Context valid residue coverage changed')
        data.to_csv(out/'all_context_predictions.csv.gz',index=False,compression='gzip')
        summarize(data,out)
        status['status']='partial_smoke' if args.limit else 'complete'
        status['elapsed_seconds']=time.monotonic()-start;atomic_json(out/'run_summary.json',status)
        print(json.dumps(status,indent=2),flush=True)
    except Exception as exc:
        status.update(status='failed',error=repr(exc));atomic_json(out/'run_summary.json',status);raise


if __name__=='__main__':main()
