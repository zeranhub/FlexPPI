"""Frozen 5-checkpoint predictions for real B4.2.3 and 1G4 structural states."""
from pathlib import Path
import argparse
import hashlib
import json
import sys
import time
import gc
import numpy as np
import torch

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT/'src'), str(ROOT/'scripts')]
from flexppi.inference import _structure_example
from flexppi.data.structures import parse_pdb
from flexppi.data.dataset import collate_states
from export_atlas_hchi_v1_oof_predictions import load_fold, implementation_hash
from run_apo_holo_exploration import subset
from chunked_inference import forward_chunked

BASE = ROOT/'artifacts/tcr_real_states_20260915'
KEYS = ['aatype','residue_index','chain_index','backbone','backbone_mask','geometry_mask','target_mask',
        'label_values','label_variances','label_mask']


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_state(spec, specs):
    structure = parse_pdb(spec['path'])
    example = _structure_example(spec['path'], spec['target_chains'], 300., 100., 0.)
    # Keep the same within-chain sequence positions when testing a subset of a
    # state's observed context. Author IDs are not blindly used as model indices.
    parent = specs.get(spec.get('full_state_id'), spec)
    original = parse_pdb(parent['path'])
    lookup = {(c,r): int(i) for c,r,i in zip(original.chain_ids,original.author_residue_ids,original.residue_index)}
    indices = [lookup[(c,r)] for c,r in zip(structure.chain_ids,structure.author_residue_ids)]
    example['residue_index'] = torch.tensor(indices, dtype=torch.long)
    assert structure.num_residues >= 2
    assert np.all(structure.backbone_mask), 'Prepared residues must have complete backbone'
    return structure, example


def validate(model, example):
    # Actual deposited coordinates, mixed chains if present; no synthetic evidence.
    n = len(example['aatype'])
    ix = np.unique(np.linspace(0, n-1, min(n, 144), dtype=int))
    small = subset(example, ix)
    batch = collate_states([small])
    with torch.inference_mode():
        native = model.forward_state(batch)
        rows = []
        for size in [17,48,144]:
            trial = forward_chunked(model, batch, size)
            for key in ['mean','log_variance','representation']:
                error = float((native[key]-trial[key]).abs().max())
                torch.testing.assert_close(native[key], trial[key], rtol=2e-5, atol=2e-5)
                rows.append({'query_chunk':size,'tensor':key,'maximum_absolute_error':error})
            assert torch.equal(native['mask'],trial['mask'])
        # A frame/padding mask test uses the unchanged same physical coordinates.
        batch['geometry_mask'][0,3] = False
        batch['padding_mask'][0,-1] = False
        left=model.forward_state(batch);right=forward_chunked(model,batch,23)
        torch.testing.assert_close(left['mean'],right['mean'],rtol=2e-5,atol=2e-5)
        rows.append({'test':'invalid_frame_and_padding','maximum_absolute_error':float((left['mean']-right['mean']).abs().max())})
    return rows


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--folds',default='0,1,2,3,4')
    parser.add_argument('--state-filter',default='')
    parser.add_argument('--query-chunk',type=int,default=48)
    parser.add_argument('--validate-only',action='store_true')
    args=parser.parse_args()
    torch.set_num_threads(4)
    specs_list=json.loads((BASE/'states.json').read_text(encoding='utf-8'))
    specs={s['id']:s for s in specs_list}
    selected=[s for s in specs_list if args.state_filter in s['id']]
    if not selected:
        raise ValueError('No states selected')
    cache=BASE/'prediction_cache';cache.mkdir(parents=True,exist_ok=True)
    native_hash=implementation_hash()
    wrapper_hash=digest(Path(__file__).with_name('chunked_inference.py'))
    runner_hash=digest(__file__)
    provenance, checks=[] , []
    for fold in map(int,args.folds.split(',')):
        model,config,path,sha=load_fold(ROOT,fold,torch.device('cpu'),True)
        st,example=load_state(selected[0],specs)
        check=validate(model,example)
        checks.append({'fold':fold,'tests':check})
        (BASE/f'chunk_validation_fold{fold}.json').write_text(json.dumps(check,indent=2),encoding='utf-8')
        print('Validated chunked/native equivalence, fold',fold,flush=True)
        if args.validate_only:
            continue
        for s in selected:
            st,example=load_state(s,specs)
            expected=json.dumps({'input_sha256':digest(s['path']), 'checkpoint_sha256':sha,
                'native_implementation':native_hash,'chunked_implementation':wrapper_hash,
                'runner':runner_hash,'index_mode':'parent_full_observed_index',
                'index_sha256':hashlib.sha256(example['residue_index'].numpy().tobytes()).hexdigest()},sort_keys=True)
            target=cache/(s['id']+f'_fold{fold}.npz')
            if target.exists():
                with np.load(target,allow_pickle=False) as z:
                    if str(z['provenance'])==expected:
                        print('Cached',fold,s['id'],flush=True);continue
            start=time.time()
            prediction=forward_chunked(model,collate_states([example]),args.query_chunk)
            values=prediction['mean'][0].numpy().copy()
            mask=prediction['mask'][0].numpy().copy()
            assert mask.all() and np.isfinite(values).all()
            np.savez_compressed(target,prediction=values,mask=mask,
                chain=np.array(st.chain_ids),author=np.array(st.author_residue_ids),
                aa=np.array(list(st.sequence)),residue_index=example['residue_index'].numpy(),provenance=expected)
            provenance.append({'fold':fold,'state':s['id'],'n_residues':st.num_residues,
                'above_training_max_residues':st.num_residues>config['data']['max_residues'],
                'seconds':time.time()-start,'checkpoint':str(path),'checkpoint_sha256':sha})
            print('Predicted',fold,s['id'],st.num_residues,round(time.time()-start,2),'s',flush=True)
            del prediction;gc.collect()
        del model;gc.collect()
    (BASE/'latest_prediction_run.json').write_text(json.dumps({'provenance':provenance,'validation':checks,
        'arguments':vars(args),'no_truncation':True,'all_keys_used':True,
        'caveat':'Five checkpoint spread is not biological replication; inputs longer than training maximum are length extrapolation.'},indent=2),encoding='utf-8')


if __name__=='__main__':
    main()
