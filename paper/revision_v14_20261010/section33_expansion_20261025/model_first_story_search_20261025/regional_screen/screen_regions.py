"""Inspect existing current core OOF predictions before biological interpretation."""
import os
os.environ['OMP_NUM_THREADS']='1'
os.environ['OPENBLAS_NUM_THREADS']='1'
from pathlib import Path
import json
from difflib import SequenceMatcher
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy.stats import spearmanr

OUT=Path(__file__).resolve().parent
ROOT=Path('E:/Project/FlexPPI')
BASE=ROOT/'paper/revision_v14_20261010/section33_expansion_20261025/feature_baseline'
IDS=['A004F','A007U','A006F','A00AQ','A000Y','A00BG','A001D']
METHODS=['full','intra','aa_sasa_contacts','static47_sasa']

def rho(a,b):
    return float(spearmanr(a,b).statistic) if len(a)>3 and np.ptp(a)>0 and np.ptp(b)>0 else np.nan

def summarize(d, sid, chain, region):
    row={'sample_id':sid,'chain':chain,'region':region,'residues':len(d)}
    for metric in ['rmsf','hchi']:
        q=d[d['valid_'+metric] & (d.amino_acid!='P') & d.sidechain_complete] if metric=='hchi' else d[d.valid_rmsf]
        y=q['md_'+metric].to_numpy()
        row[metric+'_n']=len(q)
        row[metric+'_md_mean']=float(y.mean()) if len(q) else np.nan
        for method in METHODS:
            p=q['pred_'+method+'_'+metric].to_numpy()
            row[f'{metric}_{method}_mae']=float(abs(y-p).mean()) if len(q) else np.nan
            row[f'{metric}_{method}_rho']=rho(y,p)
    return row

def main():
    f=pq.read_table(BASE/'features_all.parquet',filters=[('sample_id','in',IDS)]).to_pandas()
    predcols=['sample_id','residue_index']+['pred_'+m+'_'+s for m in ['aa_sasa_contacts','static47_sasa'] for s in ['rmsf','hchi']]
    p=pq.read_table(BASE/'oof_predictions_all.parquet',columns=predcols,filters=[('sample_id','in',IDS)]).to_pandas()
    d=f.merge(p,on=['sample_id','residue_index'],validate='one_to_one')
    metadata=pd.read_csv(ROOT/'experiments/non_tcr_case_screen_20261025/all_chain_metrics.csv')[['sample_id','chain','pdb_id','description','identity']]
    d=d.merge(metadata,on=['sample_id','chain'],how='left',validate='many_to_one')
    d.to_csv(OUT/'seven_candidate_residue_predictions.csv',index=False)
    summaries=[]
    for (sid,chain), q in d.groupby(['sample_id','chain']):
        for region,s in [('all',q),('interface_ca8',q[q.interface_ca8]),('outside_interface',q[~q.interface_ca8]),('interior_minus_5',q[q.distance_to_chain_end>=5])]:
            z=summarize(s,sid,chain,region)
            z.update(pdb_id=q.pdb_id.iloc[0],description=q.description.iloc[0])
            summaries.append(z)
    pd.DataFrame(summaries).to_csv(OUT/'seven_candidate_chain_and_region_metrics.csv',index=False)
    g=d[d.sample_id=='A004F']
    ga=g[g.chain=='A'].reset_index(drop=True).copy()
    gb=g[g.chain=='B'].reset_index(drop=True).copy()
    aseq=''.join(ga.amino_acid);bseq=''.join(gb.amino_acid)
    matches=SequenceMatcher(None,aseq,bseq,autojunk=False).get_matching_blocks()
    ia=[a for z in matches for a in range(z.a,z.a+z.size)]
    ib=[b for z in matches for b in range(z.b,z.b+z.size)]
    aa=pd.concat([ga.iloc[ia].reset_index(drop=True).add_suffix('_A'),gb.iloc[ib].reset_index(drop=True).add_suffix('_B')],axis=1)
    assert (aa.amino_acid_A==aa.amino_acid_B).all()
    aa['amino_acid']=aa.amino_acid_A
    (OUT/'geminin_chain_sequence_alignment.json').write_text(json.dumps({'chain_A_sequence':aseq,'chain_B_sequence':bseq,'exact_matching_blocks':[list(z) for z in matches],'note':'Input author_residue_id is globally renumbered and cannot be merged across chains; use exact sequence matching.'},indent=2))
    dr=[]
    for met in ['rmsf','hchi']:
        v=aa['valid_'+met+'_A']&aa['valid_'+met+'_B']
        if met=='hchi':v&=aa.sidechain_complete_A & aa.sidechain_complete_B & (aa.amino_acid!='P')
        t=aa.loc[v].copy()
        y=t['md_'+met+'_B'].to_numpy()-t['md_'+met+'_A'].to_numpy()
        for method in METHODS:
            z=t[f'pred_{method}_{met}_B'].to_numpy()-t[f'pred_{method}_{met}_A'].to_numpy()
            dr.append({'metric':met,'method':method,'n':len(t),'md_B_minus_A_mean':y.mean(),'pred_B_minus_A_mean':z.mean(),'delta_mae':abs(y-z).mean(),'zero_delta_mae':abs(y).mean(),'delta_rho':rho(y,z),'direction_accuracy_all':np.mean(np.sign(y)==np.sign(z))})
        aa.loc[v,'md_delta_'+met]=y
        aa.loc[v,'full_delta_'+met]=t[f'pred_full_{met}_B'].to_numpy()-t[f'pred_full_{met}_A'].to_numpy()
    pd.DataFrame(dr).to_csv(OUT/'geminin_same_position_chain_asymmetry_metrics.csv',index=False)
    aa.to_csv(OUT/'geminin_same_position_chain_asymmetry_residues.csv',index=False)
    pd.DataFrame(summaries).query("region=='all'").to_csv(OUT/'top_candidate_chain_summary.csv',index=False)
    print(pd.DataFrame(summaries).query("region=='all'")[["sample_id","pdb_id","chain","description","residues","rmsf_full_mae","rmsf_intra_mae","rmsf_full_rho","hchi_full_mae","hchi_intra_mae","hchi_full_rho","hchi_aa_sasa_contacts_mae","hchi_static47_sasa_mae"]].to_string(index=False))
    print(pd.DataFrame(dr).to_string(index=False))

if __name__=='__main__':main()
