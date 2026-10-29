"""Frozen five-checkpoint inference, exact residue matching, no fitting or MD.

Caches are tied to checkpoint, implementation, and input PDB SHA256 values.
Quantitative source data retain every fold and every deposited NMR conformer.
"""
from pathlib import Path
import json
import re
import time
import numpy as np
import pandas as pd
import torch
from scipy.spatial import cKDTree
from run_apo_holo_exploration import state, infer, heavy_atoms, digest, write_json
from export_atlas_hchi_v1_oof_predictions import load_fold, implementation_hash
from flexppi.data.pairing import global_sequence_alignment

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'artifacts/mechanistic_case_studies_20260913'
METRICS = ['rmsf','chi_entropy']

def match(before, after, chainmap):
    ib, ia, audit = [], [], []
    for bc, ac in chainmap.items():
        bi=np.flatnonzero(before['chain']==bc); ai=np.flatnonzero(after['chain']==ac)
        aln=global_sequence_alignment(''.join(before['aa'][bi]), ''.join(after['aa'][ai]))
        assert aln.identity > .9, (bc,ac,aln.identity)
        b=bi[aln.indices_a]; a=ai[aln.indices_b]
        valid=before['example']['geometry_mask'][b].numpy() & after['example']['geometry_mask'][a].numpy()
        ib.extend(b[valid]);ia.extend(a[valid])
        audit.append(dict(before_chain=bc,after_chain=ac,exact_identity=aln.identity,
                          before_n=len(bi),after_n=len(ai),valid_exact_n=int(valid.sum())))
    ib,ia=np.asarray(ib),np.asarray(ia)
    assert np.array_equal(before['aa'][ib],after['aa'][ia])
    return ib,ia,audit

def distances(st, targets):
    xyz,ri,_=heavy_atoms(st)
    is_target=np.isin(st['chain'][ri],targets)
    out=np.full(len(st['chain']),np.inf)
    if (~is_target).any():
        ds=cKDTree(xyz[~is_target]).query(xyz[is_target])[0]
        np.minimum.at(out,ri[is_target],ds)
    return out

def main():
    torch.set_num_threads(6)
    specs=json.loads((BASE/'states.json').read_text())
    cons=json.loads((BASE/'contrasts.json').read_text())
    states={s['id']:state(s) for s in specs}
    preds={s['id']:[] for s in specs}
    cache=BASE/'prediction_cache';cache.mkdir(exist_ok=True)
    impl=implementation_hash()
    provenance=[]
    for fold in range(5):
        model,config,cp,sha=load_fold(ROOT,fold,torch.device('cpu'),True)
        provenance.append(dict(fold=fold,checkpoint=str(cp),sha256=sha,implementation_sha256=impl))
        started=time.time()
        for i,s in enumerate(specs):
            st=states[s['id']];path=cache/f"{s['id']}_fold{fold}.npz"
            expected=json.dumps(dict(input=digest(s['path']),checkpoint=sha,implementation=impl,protocol='mechanistic_inputs_v2_60I_UCH'),sort_keys=True)
            p=None
            if path.exists():
                with np.load(path,allow_pickle=False) as z:
                    if str(z['provenance'])==expected:
                        p=z['prediction'].copy()
            if p is None:
                p=infer(model,st)
                np.savez_compressed(path,prediction=p,provenance=expected,
                                    chain=st['chain'],author=st['author'],aa=st['aa'])
            assert p.shape==(len(st['chain']),2)
            preds[s['id']].append(p)
            if i%10==0:
                print(f"Fold {fold+1}/5: {i+1}/{len(specs)} states; {time.time()-started:.1f}s",flush=True)
        del model
    write_json(BASE/'checkpoint_provenance.json',provenance)
    allrows,regions,audits=[],[],[]
    for c in cons:
        b,a=states[c['before']],states[c['after']]
        ib,ia,audit=match(b,a,c['chainmap']);audits.append(dict(contrast=c['id'],mapping=audit))
        distance=distances(a,list(c['chainmap'].values()))[ia]
        num=np.array([int(re.match(r'-?\d+',str(x)).group()) for x in a['author'][ia]])
        pb,pa=np.stack(preds[c['before']])[:,ib],np.stack(preds[c['after']])[:,ia]
        metadata=dict(contrast=c['id'],case=c['case'],chain=a['chain'][ia],author_id=a['author'][ia],
                      author_number=num,aa=a['aa'][ia],before_chain=b['chain'][ib],before_author=b['author'][ib],
                      partner_distance=distance,source_model=a['spec']['source_model'])
        rows=[]
        for f in range(5):
            frame=pd.DataFrame(metadata);frame['fold']=f
            for j,m in enumerate(METRICS):
                frame['before_'+m]=pb[f,:,j];frame['after_'+m]=pa[f,:,j]
                frame['delta_'+m]=pa[f,:,j]-pb[f,:,j]
                # A/G have no chi torsions; P is retained in raw data but
                # excluded from primary regional Hchi summary (near-constant labels).
            rows.append(frame)
        frame=pd.concat(rows,ignore_index=True);allrows.append(frame)
        for reg in c['regions']+[dict(name='all shared',chain=None,low=-9999,high=9999)]:
            sel=(frame.author_number>=reg['low'])&(frame.author_number<=reg['high'])
            if reg['chain'] is not None:
                sel &= frame.chain==reg['chain']
            for m in METRICS:
                selected=frame[sel & (~frame.aa.isin(['A','G','P']) if m=='chi_entropy' else True)]
                for f,g in selected.groupby('fold'):
                    values=g['delta_'+m].dropna()
                    regions.append(dict(contrast=c['id'],case=c['case'],region=reg['name'],metric=m,
                                        fold=int(f),source_model=a['spec']['source_model'],n=len(values),
                                        delta=float(values.mean()),median=float(values.median()),
                                        mean_partner_distance=float(g.partner_distance.replace(np.inf,np.nan).mean())))
    data=pd.concat(allrows,ignore_index=True)
    data.to_csv(BASE/'residue_predictions.csv',index=False)
    pd.DataFrame(regions).to_csv(BASE/'region_predictions.csv',index=False)
    write_json(BASE/'residue_mapping_audit.json',audits)
    write_json(BASE/'analysis_contract.json',dict(status='completed',n_states=len(specs),n_contrasts=len(cons),
       n_checkpoints=5,kix_conformers=20,model_temperature_K=300,model_duration_ns=100,
       conditions_used_by_model=False,reference='after minus before',
       uncertainty='Paired checkpoint SD is model spread, not biological replication or a confidence interval.',
       hchi_region_exclusions=['A','G','P'],no_test_labels=True,
       limitations=['Fixed-coordinate contexts are not experimental apo/holo ensembles.',
                    'The model sees amino-acid identities and N/CA/C/O; no ligands, covalent chemistry or side-chain atom positions.',
                    'Hchi is a marginal chi-dihedral entropy proxy, not S2, HDX, thermodynamic entropy or correlated motion.',
                    'Per-chain RMSF is not automatically comparable to whole-complex MD alignment.',
                    'TCR FG-loop boundary provisional pending primary-paper numbering confirmation.']))
    print('COMPLETE',data.shape,flush=True)

if __name__=='__main__':
    main()
