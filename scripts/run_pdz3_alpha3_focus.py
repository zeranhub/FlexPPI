"""Frozen five-checkpoint prediction and exact-ID no-ligand residue contrasts."""
from pathlib import Path
import json
import time
import numpy as np
import pandas as pd
import torch
from run_apo_holo_exploration import state, infer, digest, write_json
from export_atlas_hchi_v1_oof_predictions import load_fold, implementation_hash

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'artifacts/pdz3_alpha3_focus_20260914'
METRICS=['rmsf','chi_entropy']

def main():
    torch.set_num_threads(6)
    specs=json.loads((BASE/'states.json').read_text(encoding='utf-8'))
    states={s['id']:state(s) for s in specs}
    cache=BASE/'prediction_cache';cache.mkdir(exist_ok=True)
    impl=implementation_hash();runner_sha=digest(__file__)
    provenance=[];rows=[]
    for f in range(5):
        model,config,cp,sha=load_fold(ROOT,f,torch.device('cpu'),True)
        provenance.append(dict(fold=f,checkpoint=str(cp),sha256=sha,implementation_sha256=impl,runner_sha256=runner_sha))
        t=time.time()
        for i,s in enumerate(specs):
            st=states[s['id']]
            expected=json.dumps(dict(input=digest(s['path']),checkpoint=sha,implementation=impl,
                                     protocol='pdz3_alpha3_no_ligand_v1'),sort_keys=True)
            path=cache/f"{s['id']}_fold{f}.npz";p=None
            if path.exists():
                with np.load(path,allow_pickle=False) as z:
                    if str(z['provenance'])==expected:
                        p=z['prediction'].copy()
            if p is None:
                p=infer(model,st)
                np.savez_compressed(path,prediction=p,provenance=expected,chain=st['chain'],author=st['author'],aa=st['aa'])
            assert p.shape==(s['n_residues'],2) and np.isfinite(p).all(),s['id']
            d=pd.DataFrame(dict(state=s['id'],background=s['background'],context=s['context'],mutant=s['mutant'],
                 fold=f,author_number=st['author'].astype(int),aa=st['aa'],rmsf=p[:,0],chi_entropy=p[:,1]))
            rows.append(d)
            if i%50==0:
                print(f'Fold {f+1}/5: {i+1}/{len(specs)} states; {time.time()-t:.1f}s',flush=True)
        del model
    write_json(BASE/'checkpoint_provenance.json',provenance)
    data=pd.concat(rows,ignore_index=True)
    assert not data.duplicated(['state','fold','author_number']).any()
    data.to_csv(BASE/'absolute_predictions.csv',index=False)
    grouped={s:g for s,g in data.groupby('state',sort=False)}
    specmap={s['id']:s for s in specs}
    contrasts=[]
    def diff(before,after,kind):
        b=grouped[before];a=grouped[after]
        m=b.merge(a,on=['fold','author_number','aa'],suffixes=('_before','_after'),validate='one_to_one')
        changed=set(int(x[1:-1]) for x in specmap[before]['mutations']+specmap[after]['mutations'])
        m=m[~m.author_number.isin(changed)].copy()
        out=m[['fold','author_number','aa']].copy()
        out['before']=before;out['after']=after;out['kind']=kind
        out['background']=specmap[after]['background'];out['context']=specmap[after]['context'];out['mutant']=specmap[after]['mutant']
        out['hchi_valid']=~out.aa.isin(['A','G','P'])
        out['boundary_near']=out.author_number>=391
        for metric in METRICS:
            out['before_'+metric]=m[metric+'_before'];out['after_'+metric]=m[metric+'_after']
            out['delta_'+metric]=m[metric+'_after']-m[metric+'_before']
        contrasts.append(out)
    for s in specs:
        bg,ctx,mut=s['background'],s['context'],s['mutant']
        if mut!='WT':
            diff(f'{bg}_{ctx}_WT',s['id'],'mutation')
        if ctx in ('d7','d8','dN','dNd8') and f'{bg}_full_{mut}' in specmap:
            diff(f'{bg}_full_{mut}',s['id'],'deletion')
    delta=pd.concat(contrasts,ignore_index=True)
    delta.to_csv(BASE/'residue_deltas.csv',index=False)
    inter=[]
    for bg in ('apo1bfe','strip5heb'):
        for ctx in ('d7','d8'):
            for mut in sorted(data[(data.background==bg)&(data.context==ctx)].mutant.unique()):
                if mut=='WT':continue
                a=delta[(delta.background==bg)&(delta.context=='full')&(delta.kind=='mutation')&(delta.mutant==mut)]
                b=delta[(delta.background==bg)&(delta.context==ctx)&(delta.kind=='mutation')&(delta.mutant==mut)]
                m=a.merge(b,on=['author_number','fold','aa'],suffixes=('_full','_deleted'),validate='one_to_one')
                if m.empty:continue
                out=m[['author_number','fold','aa']].copy();out['background']=bg;out['context']=ctx;out['mutant']=mut
                out['hchi_valid']=~out.aa.isin(['A','G','P'])
                for metric in METRICS:
                    out['interaction_'+metric]=m['delta_'+metric+'_deleted']-m['delta_'+metric+'_full']
                inter.append(out)
    pd.concat(inter,ignore_index=True).to_csv(BASE/'mutation_alpha3_interactions.csv',index=False)
    write_json(BASE/'run_status.json',dict(status='complete',states=len(specs),evaluations=len(specs)*5,
        absolute_rows=len(data),delta_rows=len(delta),contrast_pairs=delta[['before','after']].drop_duplicates().shape[0],
        units={'rmsf':'Angstrom','chi_entropy':'nats'},interaction='Mutation effect without alpha3 minus mutation effect with alpha3; output interaction, not energy.'))
    print('COMPLETE',len(data),len(delta),flush=True)

if __name__=='__main__':main()
