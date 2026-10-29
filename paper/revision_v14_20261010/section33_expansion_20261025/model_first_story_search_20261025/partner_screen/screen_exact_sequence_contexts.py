"""Model-first exploratory screen of exact-sequence chains in existing exports.

No training; primary screen subtracts predictions from the SAME OOF fold.
All exact-sequence residues, unselected by model agreement, are retained.
"""
from pathlib import Path
import json, itertools, hashlib
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

R=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010\section33_expansion_20261025')
OUT=Path(__file__).resolve().parent
MODELS=['full','intra','aa_sasa_contacts','static47_sasa']

def score(y,p):
    ok=np.isfinite(y)&np.isfinite(p); y,p=y[ok],p[ok]
    s=np.abs(y)>=.1
    return dict(n=len(y),mae=np.mean(np.abs(y-p)),zero_mae=np.mean(np.abs(y)),
                rho=spearmanr(y,p).statistic if len(y)>3 and np.std(y)>1e-9 and np.std(p)>1e-9 else np.nan,
                direction=np.mean(np.sign(y[s])==np.sign(p[s])) if s.any() else np.nan,
                strong_n=int(s.sum()),mean_md=y.mean(),mean_pred=p.mean(),
                up_recall=np.mean(p[y>=.1]>0) if (y>=.1).any() else np.nan,
                down_recall=np.mean(p[y<=-.1]<0) if (y<=-.1).any() else np.nan)

def main():
    cols=['sample_id','residue_index','chain','author_residue_id','amino_acid',
          'sidechain_sasa_complex_A2','sidechain_sasa_isolated_A2','partner_ca_contacts_8A','interface_ca8']
    f=pd.read_parquet(R/'feature_baseline/features_all.parquet',columns=cols)
    p=pd.read_parquet(R/'feature_baseline/oof_predictions_all.parquet')
    f=f.merge(p,on=['sample_id','residue_index'],validate='one_to_one',suffixes=('','_p'))
    rec=pd.read_csv(R/'feature_baseline/cohort_records.csv').set_index('sample_id')
    groups={}; chain_meta=[]
    for (sid,ch),g in f.groupby(['sample_id','chain'],sort=False):
        g=g.sort_values('residue_index').reset_index(drop=True)
        seq=''.join(g.amino_acid)
        if len(seq)<40:continue
        key=(seq,int(g.fold.iloc[0]))
        groups.setdefault(key,[]).append((sid,ch,g))
        chain_meta.append(dict(sample_id=sid,chain=ch,n=len(seq),fold=key[1],sequence=seq,
                               pdb_id=rec.loc[sid,'pdb_id'],source_group=rec.loc[sid,'source_group'],
                               interface_n=int(g.interface_ca8.sum()),contact_total=int(g.partner_ca_contacts_8A.sum())))
    pd.DataFrame(chain_meta).to_csv(OUT/'all_chain_inventory.csv',index=False)
    results=[]; pairs=[]
    for (seq,fold),members in groups.items():
        if len(members)<2:continue
        for (s1,c1,a),(s2,c2,b) in itertools.combinations(members,2):
            # Same simulation copies are equivalent-state contexts, not state changes.
            if s1==s2:continue
            interface=a.interface_ca8.to_numpy()|b.interface_ca8.to_numpy()
            dc=b.partner_ca_contacts_8A.to_numpy()-a.partner_ca_contacts_8A.to_numpy()
            if not interface.any(): continue
            # Skip nominal replicas with effectively identical PPI geometry.
            if np.abs(dc).sum()<3:continue
            pid=f'{s1}_{c1}__{s2}_{c2}'
            n1=int(a.interface_ca8.sum());n2=int(b.interface_ca8.sum())
            pair=dict(pair_id=pid,sample_a=s1,chain_a=c1,sample_b=s2,chain_b=c2,
                      pdb_a=rec.loc[s1,'pdb_id'],pdb_b=rec.loc[s2,'pdb_id'],
                      source_a=rec.loc[s1,'source_group'],source_b=rec.loc[s2,'source_group'],
                      fold=fold,chain_n=len(seq),interface_a=n1,interface_b=n2,
                      contact_delta_sum=int(dc.sum()),contact_change_abs=int(np.abs(dc).sum()),
                      sequence_sha256=hashlib.sha256(seq.encode()).hexdigest())
            pairs.append(pair)
            for m in ['rmsf','hchi']:
                valid=a['valid_'+m].to_numpy()&b['valid_'+m].to_numpy()
                if m=='hchi':valid &= ~a.amino_acid.isin(['A','G','P']).to_numpy()
                y=b['md_'+m].to_numpy()-a['md_'+m].to_numpy()
                for region,mask in [('all',valid),('interface_union',valid&interface),('noninterface',valid&~interface)]:
                    if mask.sum()<10:continue
                    for model in MODELS:
                        pred=b['pred_'+model+'_'+m].to_numpy()-a['pred_'+model+'_'+m].to_numpy()
                        results.append(dict(**pair,metric=m,region=region,model=model,**score(y[mask],pred[mask])))
    results=pd.DataFrame(results);results['gain_vs_zero']=results.zero_mae-results.mae
    results.to_csv(OUT/'all_samefold_exact_sequence_context_metrics.csv',index=False)
    pd.DataFrame(pairs).to_csv(OUT/'all_samefold_exact_sequence_context_pairs.csv',index=False)
    q=results.query("region=='all'")
    wide=q.pivot(index=['pair_id','metric'],columns='model',values=['mae','rho','gain_vs_zero']).reset_index()
    wide.columns=['_'.join([v for v in x if v]) if isinstance(x,tuple) else x for x in wide.columns]
    wide['gain_vs_intra']=wide.mae_intra-wide.mae_full
    wide['gain_vs_strong']=wide[['mae_aa_sasa_contacts','mae_static47_sasa']].min(axis=1)-wide.mae_full
    wide=wide.merge(pd.DataFrame(pairs),on='pair_id',validate='many_to_one')
    wide.to_csv(OUT/'all_samefold_context_ranking.csv',index=False)
    print('Exact chain contexts',len(pairs))
    for m in ['rmsf','hchi']:
        print('\n',m)
        z=wide[wide.metric.eq(m)].sort_values(['gain_vs_strong','rho_full'],ascending=False).head(15)
        print(z[['pair_id','pdb_a','pdb_b','chain_n','interface_a','interface_b','mae_full','gain_vs_zero_full','gain_vs_intra','gain_vs_strong','rho_full']].to_string(index=False))
    # Save rows of strongest broad candidates for concrete, reproducible follow-up.
    top=set(wide.query('gain_vs_strong>0 & gain_vs_zero_full>0 & rho_full>.25').pair_id)
    retained=[]
    bychain={(sid,ch):g for members in groups.values() for sid,ch,g in members}
    for pair in pairs:
        if pair['pair_id'] not in top:continue
        a=bychain[pair['sample_a'],pair['chain_a']];b=bychain[pair['sample_b'],pair['chain_b']]
        for i in range(len(a)):
            row=dict(pair_id=pair['pair_id'],position=i,amino_acid=a.amino_acid.iloc[i],author_a=a.author_residue_id.iloc[i],author_b=b.author_residue_id.iloc[i])
            for prefix,g in [('a',a),('b',b)]:
                for col in ['residue_index','interface_ca8','partner_ca_contacts_8A','sidechain_sasa_complex_A2','sidechain_sasa_isolated_A2','valid_rmsf','valid_hchi','md_rmsf','md_hchi']+[f'pred_{model}_{m}' for model in MODELS for m in ['rmsf','hchi']]:
                    row[prefix+'_'+col]=g[col].iloc[i]
            retained.append(row)
    pd.DataFrame(retained).to_csv(OUT/'positive_candidate_all_residues.csv.gz',index=False,compression='gzip')

if __name__=='__main__':main()
