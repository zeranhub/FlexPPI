"""Exploratory case discovery after fixed OOF baseline training; no manuscript edits.

All 616 previously identity-audited non-TCR heteromer records are screened.
Selection is explicitly exploratory; no case-selection p-values are produced.
"""
import os
os.environ['OMP_NUM_THREADS']='1'
os.environ['MKL_NUM_THREADS']='1'
os.environ['OPENBLAS_NUM_THREADS']='1'
from pathlib import Path
import json
from collections import Counter
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy.stats import spearmanr

ROOT=Path('E:/Project/FlexPPI')
BASE=Path(__file__).resolve().parent
OUT=BASE/'case_candidates'
METHODS=['full','intra','static','aa_only','aa_sasa','aa_sasa_contacts','static47_sasa']
STRONG=['aa_sasa_contacts','static47_sasa']


def rho(y,p):
    return float(spearmanr(y,p).statistic) if len(y)>=4 and np.ptp(y)>0 and np.ptp(p)>0 else np.nan


def centered_rho(part,method):
    counts=Counter(part.amino_acid)
    keep=np.array([counts[a]>=4 and a!='P' for a in part.amino_acid])
    q=part.loc[keep]
    if len(q)<8:
        return np.nan,len(q)
    aa=q.amino_acid.to_numpy()
    y=q.md_hchi.to_numpy(dtype=float).copy()
    p=q['pred_'+method+'_hchi'].to_numpy(dtype=float).copy()
    for typ in set(aa):
        ix=aa==typ
        y[ix]=0 if np.ptp(y[ix])==0 else y[ix]-y[ix].mean()
        p[ix]=0 if np.ptp(p[ix])==0 else p[ix]-p[ix].mean()
    return rho(y,p),len(q)


def main():
    OUT.mkdir(exist_ok=True)
    cohort=pd.read_csv(BASE/'cohort_records.csv')
    previous=pd.read_csv(ROOT/'experiments/non_tcr_case_screen_20261025/all_616_record_screen.csv')
    chainlabels=pd.read_csv(ROOT/'experiments/non_tcr_case_screen_20261025/all_chain_metrics.csv')[['sample_id','chain','description','identity']].drop_duplicates(['sample_id','chain'])
    assert len(previous)==616
    allowed=set(previous.sample_id)
    columns=['sample_id','record_index','residue_index','dataset','fold','chain','author_residue_id','amino_acid','sidechain_complete',
             'valid_hchi','md_hchi','md_rmsf','valid_rmsf','interface_ca8','sidechain_sasa_complex_A2','sidechain_sasa_isolated_A2',
             'delta_sidechain_sasa_A2','intra_ca_contacts_8A','partner_ca_contacts_8A','min_partner_ca_distance_A',
             'chain_length','distance_to_chain_end']
    frame=pq.read_table(BASE/'features_all.parquet',columns=columns,
                        filters=[('sample_id','in',list(allowed))]).to_pandas()
    predcols=['sample_id','residue_index']+['pred_'+method+'_hchi' for method in METHODS]
    predictions=pq.read_table(BASE/'oof_predictions_all.parquet',columns=predcols,
                              filters=[('sample_id','in',list(allowed))]).to_pandas()
    frame=frame.merge(predictions,on=['sample_id','residue_index'],validate='one_to_one')
    inventory=pd.read_csv(BASE.parent/'paired_and_experimental/paired_inventory.csv')
    paired=set(inventory.loc[inventory.status=='evaluated','apo_sample'])|set(inventory.loc[inventory.status=='evaluated','holo_sample'])
    family=pd.read_csv(ROOT/'experiments/_sources/family_folds/record_family_folds.csv')
    clean=set(family.groupby('component_id').original_fold.nunique().loc[lambda x:x==1].index)
    family['family_clean_original_fold']=family.component_id.isin(clean)
    records=cohort.merge(previous[['sample_id','description','n_protein_identities']],on='sample_id',how='inner',validate='one_to_one')
    records=records.merge(family[['sample_id','component_id','family_clean_original_fold']],on='sample_id',validate='one_to_one')
    records['true_paired_available']=records.sample_id.isin(paired)
    metadata=records.set_index('sample_id').to_dict('index')
    metric_rows=[];pair_rows=[]
    for (sid,chain),chainpart in frame.groupby(['sample_id','chain'],sort=False):
        part=chainpart[chainpart.valid_hchi&chainpart.sidechain_complete&(chainpart.amino_acid!='P')].copy()
        for subset,subpart in [('chain',part),('interface_ca8',part[part.interface_ca8])]:
            if len(subpart)<8:
                continue
            result=dict(sample_id=sid,chain=chain,subset=subset,n=len(subpart),
                        true_hchi_sd=float(subpart.md_hchi.std()),true_hchi_span=float(subpart.md_hchi.max()-subpart.md_hchi.min()),
                        true_rmsf_mean=float(subpart.md_rmsf.mean()))
            for method in METHODS:
                y=subpart.md_hchi.to_numpy();p=subpart['pred_'+method+'_hchi'].to_numpy()
                cr,n=centered_rho(subpart,method)
                result[method+'_mae']=float(np.abs(y-p).mean())
                result[method+'_rho']=rho(y,p)
                result[method+'_typed_rho']=cr
                result['typed_n']=n
            for method in ['intra']+STRONG:
                result['mae_gain_vs_'+method]=result[method+'_mae']-result['full_mae']
                result['typed_rho_gain_vs_'+method]=result['full_typed_rho']-result[method+'_typed_rho']
            metric_rows.append(result)
        for aa,group in part.groupby('amino_acid',sort=False):
            if len(group)<4:
                continue
            y=group.md_hchi.to_numpy();rmsf=group.md_rmsf.to_numpy();sasa=group.sidechain_sasa_complex_A2.to_numpy()
            sasa_rank=pd.Series(sasa).rank(method='average',pct=True).to_numpy()
            rows=group.to_dict('records')
            for i in range(len(group)):
                for j in range(i+1,len(group)):
                    if abs(sasa[i]-sasa[j])>10 or abs(sasa_rank[i]-sasa_rank[j])>.25 or abs(rmsf[i]-rmsf[j])>.10:
                        continue
                    if min(rows[i]['distance_to_chain_end'],rows[j]['distance_to_chain_end'])<=3:
                        continue
                    if not(rows[i]['interface_ca8'] or rows[j]['interface_ca8']):
                        continue
                    # Always order by true entropy, never by model prediction.
                    lo,hi=(i,j) if y[i]<y[j] else (j,i)
                    low,high=rows[lo],rows[hi]
                    truth=float(y[hi]-y[lo])
                    result=dict(sample_id=sid,chain=chain,amino_acid=aa,aa_count_in_chain=len(group),
                                low_residue_index=low['residue_index'],high_residue_index=high['residue_index'],
                                low_residue=aa+str(low['author_residue_id']),high_residue=aa+str(high['author_residue_id']),
                                low_sasa=float(sasa[lo]),high_sasa=float(sasa[hi]),sasa_abs_gap=float(abs(sasa[hi]-sasa[lo])),
                                sasa_rank_gap=float(abs(sasa_rank[hi]-sasa_rank[lo])),
                                low_rmsf=float(rmsf[lo]),high_rmsf=float(rmsf[hi]),rmsf_abs_gap=float(abs(rmsf[hi]-rmsf[lo])),
                                low_true_hchi=float(y[lo]),high_true_hchi=float(y[hi]),true_gap=truth,
                                low_interface=low['interface_ca8'],high_interface=high['interface_ca8'],
                                low_intra_contacts=low['intra_ca_contacts_8A'],high_intra_contacts=high['intra_ca_contacts_8A'],
                                low_partner_contacts=low['partner_ca_contacts_8A'],high_partner_contacts=high['partner_ca_contacts_8A'],
                                low_partner_distance=low['min_partner_ca_distance_A'],high_partner_distance=high['min_partner_ca_distance_A'])
                    for method in METHODS:
                        a=low['pred_'+method+'_hchi'];b=high['pred_'+method+'_hchi']
                        result[method+'_low']=a;result[method+'_high']=b
                        result[method+'_gap']=b-a
                        result[method+'_gap_error']=float(abs((b-a)-truth))
                        result[method+'_site_mae']=float((abs(a-y[lo])+abs(b-y[hi]))/2)
                    result['dynamic_gap']=truth>=.40
                    result['full_restores_gap']=(truth>=.40 and result['full_gap']>=.25 and result['full_gap_error']<=.20 and result['full_site_mae']<=.20)
                    result['stronger_than_both_baselines']=(result['full_restores_gap'] and
                        min(result[m+'_gap_error'] for m in STRONG)-result['full_gap_error']>=.10 and
                        min(result[m+'_site_mae'] for m in STRONG)-result['full_site_mae']>=.05)
                    pair_rows.append(result)
    metrics=pd.DataFrame(metric_rows)
    pairs=pd.DataFrame(pair_rows)
    metrics=metrics.merge(chainlabels,on=['sample_id','chain'],how='left',validate='many_to_one')
    metrics=metrics.merge(records[['sample_id','pdb_id','true_paired_available','family_clean_original_fold','component_id']],on='sample_id',validate='many_to_one')
    metrics.to_csv(OUT/'all_chain_and_interface_metrics.csv',index=False)
    pairs=pairs.merge(records[['sample_id','pdb_id','description','true_paired_available','family_clean_original_fold']],on='sample_id',validate='many_to_one')
    pairs.to_csv(OUT/'all_prespecified_similar_environment_pairs.csv.gz',index=False,compression='gzip')
    selected=pairs[pairs.stronger_than_both_baselines].copy()
    selected.to_csv(OUT/'successful_pairs_against_both_strong_baselines.csv',index=False)
    summaries=[]
    for (sid,chain),part in pairs.groupby(['sample_id','chain'],sort=False):
        success=part[part.stronger_than_both_baselines].sort_values('true_gap',ascending=False)
        used=set();disjoint=[]
        for _,row in success.iterrows():
            ids={row.low_residue_index,row.high_residue_index}
            if not used&ids:
                disjoint.append(row.to_dict());used|=ids
        meta=metadata[sid]
        summaries.append(dict(sample_id=sid,chain=chain,pdb_id=meta['pdb_id'],description=meta['description'],
                              component_id=meta['component_id'],family_clean_original_fold=meta['family_clean_original_fold'],
                              true_paired_available=meta['true_paired_available'],n_matched_pairs=len(part),
                              n_dynamic_pairs=int(part.dynamic_gap.sum()),n_full_restored_pairs=int(part.full_restores_gap.sum()),
                              n_pairs_superior_to_both=len(success),n_disjoint_superior_pairs=len(disjoint),
                              n_aa_in_superior_pairs=success.amino_acid.nunique(),
                              selected_residue_pairs='; '.join(r['low_residue']+'/'+r['high_residue'] for r in disjoint[:8])))
    ranking=pd.DataFrame(summaries)
    chainm=metrics[metrics.subset=='chain'].drop(columns=['description','pdb_id','true_paired_available','family_clean_original_fold','component_id','subset','identity'])
    ranking=ranking.merge(chainm,on=['sample_id','chain'],how='left',validate='one_to_one')
    inter=metrics[metrics.subset=='interface_ca8'][['sample_id','chain','n','full_mae','intra_mae','full_typed_rho','intra_typed_rho','mae_gain_vs_intra']].rename(columns={c:'interface_'+c for c in ['n','full_mae','intra_mae','full_typed_rho','intra_typed_rho','mae_gain_vs_intra']})
    ranking=ranking.merge(inter,on=['sample_id','chain'],how='left',validate='one_to_one')
    ranking['passes_chain_gains']=(ranking.mae_gain_vs_aa_sasa_contacts>.02)&(ranking.mae_gain_vs_static47_sasa>.02)&(ranking.typed_rho_gain_vs_aa_sasa_contacts>.10)&(ranking.typed_rho_gain_vs_static47_sasa>.10)
    ranking['passes_interface_gain']=(ranking.interface_mae_gain_vs_intra>.02)&(ranking.interface_n>=15)
    ranking['strict_case_candidate']=ranking.passes_chain_gains&ranking.passes_interface_gain&(ranking.n_disjoint_superior_pairs>=3)&(ranking.typed_n>=30)
    ranking['case_score']=ranking.n_disjoint_superior_pairs.clip(upper=8)+ranking.n_aa_in_superior_pairs.clip(upper=4)+3*ranking.passes_chain_gains+2*ranking.passes_interface_gain+ranking.true_paired_available
    ranking=ranking.sort_values(['strict_case_candidate','case_score','n_disjoint_superior_pairs','mae_gain_vs_static47_sasa'],ascending=False)
    ranking.to_csv(OUT/'candidate_chain_ranking.csv',index=False)
    ranking[ranking.true_paired_available].to_csv(OUT/'true_paired_candidate_chain_ranking.csv',index=False)
    top=ranking.head(15)
    topids=set(top.sample_id)
    frame[frame.sample_id.isin(topids)].to_parquet(OUT/'top15_candidate_residues.parquet',index=False,compression='zstd')
    selected[selected.sample_id.isin(topids)].to_csv(OUT/'top15_candidate_selected_pairs.csv',index=False)
    write=dict(status='complete',screened_records=len(allowed),screened_chains=int(metrics[metrics.subset=='chain'].shape[0]),
               source_scope='616 previously protein-identity-audited non-TCR heteromer records; all retained',
               environment_match='Same record, same chain, same AA group n>=4, nonPro; absolute sidechain SASA difference <=10A2, within-AA SASA rank difference <=0.25, true C-alpha RMSF difference <=0.10A; at least one partner CA8 contact; >3 residues from either chain end',
               dynamic_pair='True Hchi gap >=0.40 nats',
               restored='Full predicted gap >=0.25, gap error<=0.20 and two-site MAE<=0.20 nats',
               superior='Gap error improves by >=0.10 and two-site MAE by >=0.05 against each AA+SASA+contacts and Static47+SASA',
               chain_gains='Full MAE gain>0.02 and type-centered rho gain>0.10 over both strong baselines',
               strict_case='Chain gains, interface Full vs Intra MAE gain>0.02 with n>=15, >=3 mutually disjoint successful pairs, >=30 within-AA eligible chain positions',
               all_prespecified_pairs=len(pairs),dynamic_pairs=int(pairs.dynamic_gap.sum()),superior_pairs=len(selected),strict_candidates=int(ranking.strict_case_candidate.sum()),
               limitations='Exploratory post-training case discovery; thresholds are descriptive screening criteria, not heldout significance tests; no p-values. Similar SASA/RMSF are label/structure matching, not mechanistic causality.')
    (OUT/'screening_contract.json').write_text(json.dumps(write,indent=2),encoding='utf-8')
    print(json.dumps(write,indent=2))
    print(top[['sample_id','pdb_id','chain','strict_case_candidate','n_disjoint_superior_pairs','n_aa_in_superior_pairs','full_mae','aa_sasa_contacts_mae','static47_sasa_mae','full_typed_rho','aa_sasa_contacts_typed_rho','static47_sasa_typed_rho','interface_mae_gain_vs_intra','true_paired_available','selected_residue_pairs']].to_string(index=False))


if __name__=='__main__':main()
