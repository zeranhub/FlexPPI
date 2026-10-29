"""Model-first new inference, all 24 available Zenodo states and all 5 folds.

Legacy MD labels are provisional screening only until raw-trajectory chemical
and population checks pass for any selected case. No fitting or model choice.
"""
import os
for k in ['OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS']: os.environ[k]='2'
os.environ['CUDA_VISIBLE_DEVICES']=''
import sys,json,csv,hashlib,time
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import torch
ROOT=Path(r'E:\Project\FlexPPI\experiments\core_corrected_pipeline')
DATA=Path(r'E:\Project\FlexPPI\artifacts\zenodo_6638504\data')
OUT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'training_payload/src'))
from flexppi.config import ModelConfig
from flexppi.constants import METRIC_NAMES
from flexppi.data.structures import parse_pdb,target_mask
from flexppi.data.dataset import collate_states
from flexppi.model.information_ablation import build_information_ablation
torch.set_num_threads(2);torch.set_num_interop_threads(1)
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p,x): p.write_text(json.dumps(x,indent=2,ensure_ascii=False),encoding='utf-8')
def rho(a,b): return float(spearmanr(a,b).statistic) if len(a)>3 and np.ptp(a)>0 and np.ptp(b)>0 else float('nan')
def main():
    start=time.time();(OUT/'predictions').mkdir(exist_ok=True)
    records=[json.loads(x) for x in (DATA/'zenodo_6638504_states.jsonl').read_text().splitlines()]
    prepared={};input_audit=[];weights=[]
    for r in records:
        s=parse_pdb(r['structure_path']);n=s.num_residues
        e=dict(sample_id=r['sample_id'],dataset='zenodo',source_group='new_inference',metadata={},
            aatype=torch.from_numpy(s.aatype),residue_index=torch.from_numpy(s.residue_index),chain_index=torch.from_numpy(s.chain_index),
            backbone=torch.from_numpy(s.backbone),backbone_mask=torch.from_numpy(s.backbone_mask),geometry_mask=torch.from_numpy(s.frame_mask),
            target_mask=torch.from_numpy(target_mask(s,[])),condition=torch.zeros(3),
            label_values=torch.full((n,3),float('nan')),label_variances=torch.zeros(n,3),label_mask=torch.zeros(n,3,dtype=torch.bool))
        prepared[r['sample_id']]=(s,collate_states([e]),r)
        input_audit.append(dict(sample_id=r['sample_id'],n=n,path=r['structure_path'],sha256=sha(r['structure_path']),sequence=s.sequence))
    dump(OUT/'input_audit.json',input_audit)
    for variant in ['full','intrachain']:
        for fold in range(5):
            cp=ROOT/'trained_checkpoints'/variant/f'fold_{fold}'/'epoch_040.pt'
            c=torch.load(cp,map_location='cpu',weights_only=False);conf=c['config']['model']
            assert conf['information_ablation']==variant and conf['metrics']==['rmsf','chi_entropy']
            model=build_information_ablation(ModelConfig(**conf),variant).eval()
            model.load_state_dict(c['model'],strict=True);model.requires_grad_(False)
            weights.append(dict(variant=variant,fold=fold,path=str(cp),sha256=sha(cp)))
            for sid,(s,b,r) in prepared.items():
                with torch.inference_mode(): z=model.forward_state(b)
                mean=z['mean'][0].numpy();mask=z['mask'][0].numpy()
                assert mean.shape==(s.num_residues,2) and np.isfinite(mean).all() and np.array_equal(mask,s.frame_mask)
                np.savez_compressed(OUT/'predictions'/f'{variant}_f{fold}_{sid}.npz',mean=mean,mask=mask,
                    chain=np.array(s.chain_ids),author=np.array(s.author_residue_ids),aa=np.array(list(s.sequence)))
            print('INFERRED',variant,fold,flush=True)
            del model,c
    dump(OUT/'checkpoint_manifest.json',weights)
    states={};state_metrics=[];rows=[]
    for sid,(s,b,r) in prepared.items():
        with np.load(r['label_path']) as z: y=z['values'][:,:2];valid=z['mask'][:,:2]
        if len(y)!=s.num_residues: raise ValueError(('Legacy label population mismatch',sid,len(y),s.num_residues))
        means={v:np.mean([np.load(OUT/'predictions'/f'{v}_f{f}_{sid}.npz')['mean'] for f in range(5)],axis=0) for v in ['full','intrachain']}
        states[sid]=dict(s=s,y=y,valid=valid,means=means)
        for chain in s.unique_chains():
            for m,metric in enumerate(['rmsf','hchi']):
                ix=np.array(s.chain_ids)==chain;ix&=valid[:,m]&s.frame_mask
                if m==1: ix&=~np.isin(list(s.sequence),list('AGP'))
                a=dict(sample_id=sid,chain=chain,metric=metric,n=int(ix.sum()),legacy_label_status='provisional_screen')
                for v in means: a.update({v+'_mae':float(np.abs(means[v][ix,m]-y[ix,m]).mean()),v+'_rho':rho(means[v][ix,m],y[ix,m])})
                state_metrics.append(a)
        for i,(chain,author,aa) in enumerate(zip(s.chain_ids,s.author_residue_ids,s.sequence)):
            row=dict(sample_id=sid,index=i,chain=chain,author=author,aa=aa)
            for m,metric in enumerate(['rmsf','hchi']):
                row['md_'+metric]=float(y[i,m]);row['valid_'+metric]=bool(valid[i,m] and s.frame_mask[i])
                for v in means:row['pred_'+v+'_'+metric]=float(means[v][i,m])
            rows.append(row)
    pd.DataFrame(rows).to_csv(OUT/'all_state_residues_provisional.csv',index=False)
    pd.DataFrame(state_metrics).to_csv(OUT/'chain_state_metrics_provisional.csv',index=False)
    pairs=[json.loads(x) for x in (DATA/'zenodo_6638504_pairs.jsonl').read_text().splitlines()]
    result=[];deltas=[]
    for p in pairs:
        a=states[p['apo']['sample_id']];h=states[p['holo']['sample_id']]
        ai=np.array(p['apo_indices']);hi=np.array(p['holo_indices'])
        assert ''.join(a['s'].sequence[i] for i in ai)==''.join(h['s'].sequence[i] for i in hi)
        chain=p['metadata']['holo_target_chain'];s=h['s']
        ci=np.where(np.array(s.chain_ids)==chain)[0]; oi=np.where(np.array(s.chain_ids)!=chain)[0]
        xyz=s.backbone[:,1];dist=np.sqrt(((xyz[ci,None,:]-xyz[None,oi,:])**2).sum(-1)).min(1)
        interfaces=dict(zip(ci,(dist<=8).tolist()))
        for j,(i,k) in enumerate(zip(ai,hi)):
            row=dict(pair=p['pair_id'],apo_index=int(i),holo_index=int(k),chain=chain,author=s.author_residue_ids[k],aa=s.sequence[k],interface_ca8=interfaces.get(k,False))
            for m,metric in enumerate(['rmsf','hchi']):
                row['md_delta_'+metric]=float(h['y'][k,m]-a['y'][i,m]);row['valid_'+metric]=bool(a['valid'][i,m] and h['valid'][k,m] and a['s'].frame_mask[i] and s.frame_mask[k])
                for v in ['full','intrachain']: row[v+'_delta_'+metric]=float(h['means'][v][k,m]-a['means'][v][i,m])
            deltas.append(row)
        for m,metric in enumerate(['rmsf','hchi']):
            good=a['valid'][ai,m]&h['valid'][hi,m]&a['s'].frame_mask[ai]&s.frame_mask[hi]
            if m==1:good&=~np.isin([s.sequence[k] for k in hi],list('AGP'))
            for region in ['all','interface','noninterface']:
                g=good.copy()
                if region!='all': g&=np.array([interfaces[k]==(region=='interface') for k in hi])
                if g.sum()<4:continue
                y=h['y'][hi[g],m]-a['y'][ai[g],m]
                row=dict(pair=p['pair_id'],metric=metric,region=region,n=int(g.sum()),zero_mae=float(np.abs(y).mean()),md_mean_delta=float(y.mean()))
                for v in ['full','intrachain']:
                    z=h['means'][v][hi[g],m]-a['means'][v][ai[g],m]
                    row[v+'_mae']=float(np.abs(z-y).mean());row[v+'_rho']=rho(y,z);row[v+'_pred_mean_delta']=float(z.mean())
                    row[v+'_direction_all']=float((np.sign(y)==np.sign(z)).mean());row[v+'_gain_zero']=row['zero_mae']-row[v+'_mae']
                result.append(row)
    df=pd.DataFrame(result);df.to_csv(OUT/'all_pair_metrics_provisional.csv',index=False)
    pd.DataFrame(deltas).to_csv(OUT/'all_pair_residues_provisional.csv',index=False)
    print(df[(df.region=='all')].sort_values('full_gain_zero',ascending=False).to_string(index=False),flush=True)
    dump(OUT/'completion.json',dict(status='all_inference_complete',states=24,neural_forwards=240,elapsed_seconds=time.time()-start,
        screening_labels='Legacy labels are provisional, selected cases require raw canonical recomputation before claims.',script_sha256=sha(__file__)))
if __name__=='__main__':main()
