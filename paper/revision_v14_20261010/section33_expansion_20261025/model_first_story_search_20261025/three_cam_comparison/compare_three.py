"""Fair three-state CaM comparison with one frozen fold0 model per arm."""
import os
for k in ['OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS']:os.environ[k]='2'
from pathlib import Path
import json,itertools
import numpy as np
import pandas as pd
from Bio.Align import PairwiseAligner
from Bio.SeqUtils import seq1
from scipy.stats import spearmanr
import joblib
from threadpoolctl import threadpool_limits
OUT=Path(__file__).resolve().parent;P=OUT.parent/'partner_screen';R=OUT.parent.parent
CORE=Path(r'E:\Project\FlexPPI\experiments\core_corrected_pipeline')
SIDS=['A001K','A004C','A007U'];LABELS=dict(zip(SIDS,['CaMKI_1MXE','NHE1_2YGG','IQCG_4LZX']))
REF='ADQLTEEQIAEFKEAFSLFDKDGDGTITTKELGTVMRSLGQNPTEAELQDMINEVDADGNGTIDFPEFLTMMARKMKDTDSEEEIREAFRVFDKDGNGYISAAELRHVMTNLGEKLTDEEVDEMIREADIDGDGQVNYEEFVQMMTAK'
def score(y,z):
    y=np.asarray(y);z=np.asarray(z);large=np.abs(y)>=.1
    return dict(n=len(y),mae=float(np.mean(np.abs(y-z))),rho=float(spearmanr(y,z).statistic),
        zero_mae=float(np.abs(y).mean()),MD_mean=float(y.mean()),pred_mean=float(z.mean()),
        large_n=int(large.sum()),large_correct=int(((np.sign(y)==np.sign(z))&large).sum()))
def main():
    a=PairwiseAligner(mode='global',match_score=3,mismatch_score=-5,open_gap_score=-8,extend_gap_score=-1)
    states={};mapping=[];audit=[]
    for sid in SIDS:
        z=np.load(CORE/'neural_exports/evaluation/full'/(sid+'.npz'));chain=max(set(z['chain']),key=lambda c:sum(z['chain']==c))
        ix=np.flatnonzero(z['chain']==chain);seq=''.join(z['amino_acid'][ix]);al=a.align(REF,seq)[0];mm={}
        mismatch=[];missing=[]
        pairs=list(al.indices.T)
        if sid=='A001K':
            # The truncated C-terminal sequence permits an equally scoring but
            # biologically wrong alignment. Native PDB numbering resolves it:
            # REMARK 465 explicitly identifies Lys148 as absent, Ser147 present.
            raw=Path(r'E:\Project\FlexPPI\artifacts\apo_holo_exploration_20260913\raw\1MXE.pdb')
            native={}
            for line in raw.read_text().splitlines():
                if line[:6].strip()=='ATOM' and line[21:22]=='A' and line[12:16].strip()=='CA' and line[16:17] in [' ','A']:
                    native.setdefault(int(line[22:26]),seq1(line[17:20]))
            assert ''.join(native.values())==seq
            pairs=[(pos-1,j) for j,pos in enumerate(native)]
            missing=[p for p in range(1,149) if p not in native]
        for rr,ss in pairs:
            if rr<0:continue
            if ss<0:missing.append(int(rr+1));continue
            index=int(ix[ss]);mm[int(rr+1)]=index
            if REF[rr]!=seq[ss]:mismatch.append(dict(mature_position=int(rr+1),reference=REF[rr],observed=seq[ss]))
            mapping.append(dict(sample_id=sid,label=LABELS[sid],chain=chain,residue_index=index,mature_position=int(rr+1),amino_acid=seq[ss],reference_amino_acid=REF[rr],exact=REF[rr]==seq[ss]))
        preds={v:np.load(P/'training_member_predictions'/(v+'_'+sid+'.npz'))['mean'] for v in ['full','intrachain']}
        states[sid]=dict(cache=z,map=mm,preds=preds,chain=chain)
        audit.append(dict(sample_id=sid,label=LABELS[sid],cam_chain=chain,observed_length=len(seq),mature_missing=missing,mature_substitutions=mismatch,sequence=seq,
            numbering_method='Native 1MXE chain A numbering, sequence exactly equal to prepared CaM; missing Lys148 established by REMARK465.' if sid=='A001K' else 'Exact mature-CaM sequence alignment; extra expression residue excluded.'))
    pd.DataFrame(mapping).to_csv(OUT/'all_three_mature_position_mapping.csv',index=False)
    (OUT/'sequence_audit.json').write_text(json.dumps(audit,indent=2))
    feat=pd.read_parquet(R/'feature_baseline/features_all.parquet',filters=[('sample_id','in',SIDS)]);contract=json.loads((R/'feature_baseline/baseline_contract.json').read_text())
    with threadpool_limits(limits=2):
        for sid,state in states.items():
            f=feat[feat.sample_id==sid].sort_values('residue_index');assert np.array_equal(f.residue_index,np.arange(len(f)))
            for model in ['aa_sasa_contacts','static47_sasa']:
                pred=np.column_stack([joblib.load(R/'feature_baseline/models'/f'{model}_fold0_{m}.joblib').predict(f[contract['feature_sets'][model]].to_numpy(np.float32)) for m in ['rmsf','hchi']])
                state['preds'][model]=pred
    common=sorted(set.intersection(*(set(s['map']) for s in states.values())))
    exact=[p for p in common if all(s['cache']['amino_acid'][s['map'][p]]==REF[p-1] for s in states.values())]
    print('EXACT COMMON',len(exact),'excluded mismatches',[p for p in common if p not in exact],flush=True)
    absolute=[];pairs=[];long=[]
    for metric,m in [('rmsf',0),('hchi',1)]:
        positions=[p for p in exact if all(s['cache']['valid'][s['map'][p],m] for s in states.values()) and (m==0 or REF[p-1] not in 'AGP')]
        for sid,s in states.items():
            for population,pp in [('all_own_CaM',[p for p,i in s['map'].items() if s['cache']['valid'][i,m] and (m==0 or s['cache']['amino_acid'][i] not in 'AGP')]),('three_state_common_exact',positions)]:
                ix=np.array([s['map'][p] for p in pp]);y=s['cache']['target'][ix,m]
                for model,pr in s['preds'].items():absolute.append(dict(state=LABELS[sid],population=population,metric=metric,model=model,**score(y,pr[ix,m])))
        for sid1,sid2 in itertools.combinations(SIDS,2):
            s1,s2=states[sid1],states[sid2]
            for population,pp in [('pair_common_exact',[p for p in sorted(set(s1['map'])&set(s2['map'])) if s1['cache']['amino_acid'][s1['map'][p]]==s2['cache']['amino_acid'][s2['map'][p]] and s1['cache']['valid'][s1['map'][p],m] and s2['cache']['valid'][s2['map'][p],m] and (m==0 or s1['cache']['amino_acid'][s1['map'][p]] not in 'AGP')]),('three_state_common_exact',positions)]:
                i1=np.array([s1['map'][p] for p in pp]);i2=np.array([s2['map'][p] for p in pp]);y=s2['cache']['target'][i2,m].astype(float)-s1['cache']['target'][i1,m].astype(float)
                for model in s1['preds']:
                    pr=s2['preds'][model][i2,m].astype(float)-s1['preds'][model][i1,m].astype(float)
                    pairs.append(dict(state_a=LABELS[sid1],state_b=LABELS[sid2],population=population,metric=metric,model=model,**score(y,pr)))
                for p,i,j in zip(pp,i1,i2):
                    row=dict(state_a=LABELS[sid1],state_b=LABELS[sid2],population=population,metric=metric,mature_position=p,aa=REF[p-1],
                        index_a=int(i),index_b=int(j),md_a=float(s1['cache']['target'][i,m]),md_b=float(s2['cache']['target'][j,m]))
                    row['md_delta']=row['md_b']-row['md_a']
                    for model in s1['preds']:row[model+'_delta']=float(s2['preds'][model][j,m])-float(s1['preds'][model][i,m])
                    long.append(row)
    pd.DataFrame(absolute).to_csv(OUT/'three_state_absolute_metrics.csv',index=False)
    pd.DataFrame(pairs).to_csv(OUT/'three_state_pairwise_metrics.csv',index=False)
    pd.DataFrame(long).to_csv(OUT/'three_state_pairwise_all_residues.csv',index=False)
    print('\nABSOLUTE COMMON\n',pd.DataFrame(absolute).query("population=='three_state_common_exact'").to_string(index=False))
    print('\nPAIR COMMON\n',pd.DataFrame(pairs).query("population=='three_state_common_exact'").to_string(index=False))
    (OUT/'completion.json').write_text(json.dumps(dict(status='complete',fold=0,training_members=SIDS,models=['full','intrachain','aa_sasa_contacts','static47_sasa'],common_exact_positions=len(exact),
        rationale='Same frozen model and same population for all pair comparisons. No new fitting; source data values retained.'),indent=2))
if __name__=='__main__':main()
