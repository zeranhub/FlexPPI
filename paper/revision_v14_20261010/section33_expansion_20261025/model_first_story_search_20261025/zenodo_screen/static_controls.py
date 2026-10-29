"""Apply existing fivefold SASA-aware regressors to the selected RNase case."""
import os
for k in ['OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS']:os.environ[k]='2'
from pathlib import Path
import importlib.util,json,sys
import numpy as np
import pandas as pd
import joblib
import mdtraj as md
from threadpoolctl import threadpool_limits
OUT=Path(__file__).resolve().parent;R=OUT.parent.parent;CORE=Path(r'E:\Project\FlexPPI\experiments\core_corrected_pipeline')
sys.path.insert(0,str(CORE/'training_payload/src'))
from flexppi.data.structures import parse_pdb
def loadmod(name,path):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
feat=loadmod('area_features',R/'feature_baseline/run_feature_baselines.py')
static=loadmod('static_features',CORE/'static_feature_definitions.py')
def main():
    contract=json.loads((R/'feature_baseline/baseline_contract.json').read_text())
    names=['1AY7_apo_A','1AY7_apo_B','1AY7_holo_complex'];preds={};summary=[]
    for sid in names:
        p=OUT/'raw_recomputed'/sid/'canonical_topology.pdb';s=parse_pdb(str(p));x,n,_,_=static.features(s);tab=pd.DataFrame(x,columns=n)
        traj,keys,rids,sc,quality=feat.heavy_structure(p);count=len(keys);assert count==s.num_residues
        area=md.shrake_rupley(traj,probe_radius=.14,n_sphere_points=960,mode='atom')[0]*100
        cw=np.bincount(rids,weights=area,minlength=count);cs=np.bincount(rids[sc],weights=area[sc],minlength=count);iw=cw.copy();isc=cs.copy()
        for chain in traj.topology.chains:
            atoms=np.array([a.index for a in chain.atoms],int);aa=md.shrake_rupley(traj.atom_slice(atoms),probe_radius=.14,n_sphere_points=960,mode='atom')[0]*100
            rr=rids[atoms];scm=sc[atoms];unique=np.unique(rr)
            iw[unique]=np.bincount(rr,weights=aa,minlength=count)[unique];isc[unique]=np.bincount(rr[scm],weights=aa[scm],minlength=count)[unique]
        for col,val in zip(feat.SASA_COLUMNS,[cs,cw,isc,iw,isc-cs,iw-cw]):tab[col]=val
        tab.to_csv(OUT/(sid+'_static_features.csv'),index=False)
        with threadpool_limits(limits=2):
            for model in ['aa_sasa_contacts','static47_sasa']:
                cols=contract['feature_sets'][model]
                pred=[]
                for metric in ['rmsf','hchi']:
                    pred.append(np.mean([joblib.load(R/'feature_baseline/models'/f'{model}_fold{f}_{metric}.joblib').predict(tab[cols].to_numpy(np.float32)) for f in range(5)],axis=0))
                preds[sid,model]=np.array(pred).T
    pairs=[json.loads(x) for x in (Path(r'E:\Project\FlexPPI\artifacts\zenodo_6638504\data')/'zenodo_6638504_pairs.jsonl').read_text().splitlines()]
    original=pd.read_csv(OUT/'verified_all_residues.csv');rows=[]
    for p in pairs:
        if p['metadata']['pdb_id']!='1AY7':continue
        part=original[original.pair==p['pair_id']].copy();ai=np.array(p['apo_indices']);hi=np.array(p['holo_indices'])
        for model in ['aa_sasa_contacts','static47_sasa']:
            for j,metric in enumerate(['rmsf','hchi']):part[model+'_delta_'+metric]=preds[p['holo']['sample_id'],model][hi,j]-preds[p['apo']['sample_id'],model][ai,j]
        rows.append(part)
        for metric in ['rmsf','hchi']:
            good=part['valid_'+metric].to_numpy(bool,copy=True)
            if metric=='hchi':good&=~part.aa.isin(list('AGP')).to_numpy()
            for region in ['all','interface','noninterface']:
                keep=good.copy()
                if region!='all':keep&=part.interface_ca8.to_numpy(bool)==(region=='interface')
                g=part[keep];y=g['md_delta_'+metric]
                for model in ['full','intrachain','aa_sasa_contacts','static47_sasa']:
                    z=g[model+'_delta_'+metric]
                    summary.append(dict(pair=p['pair_id'],metric=metric,region=region,model=model,n=len(g),mae=float((z-y).abs().mean()),zero_mae=float(y.abs().mean())))
    pd.concat(rows).to_csv(OUT/'verified_all_residues_with_static_controls.csv',index=False)
    q=pd.DataFrame(summary);q.to_csv(OUT/'verified_all_model_comparison.csv',index=False);print(q.to_string(index=False))
if __name__=='__main__':main()
