"""Paired Full-minus-control comparisons; whole-family cluster bootstrap."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd

OUT=Path(__file__).resolve().parent
CORE=Path('E:/Project/FlexPPI/experiments/core_corrected_pipeline')
FAMILY=Path('E:/Project/FlexPPI/experiments/_sources/family_folds/record_family_folds.csv')
BOOTSTRAPS=2000


def interval(frame,value,key):
    data=frame[['component_id',value]].dropna()
    n=len(data)
    if not n:
        return dict(n_records=0,n_families=0,mean_gain=np.nan,ci_low=np.nan,ci_high=np.nan)
    grouped=data.groupby('component_id')[value].agg(['sum','count'])
    g=len(grouped)
    result=dict(n_records=n,n_families=g,mean_gain=float(data[value].mean()),ci_low=np.nan,ci_high=np.nan,
                fraction_full_better=float((data[value]>0).mean()))
    if g<2:
        result['ci_note']='One family: no nondegenerate family confidence interval is identifiable.'
        return result
    seed=20261025+int(hashlib.sha256(key.encode()).hexdigest()[:7],16)
    rng=np.random.default_rng(seed)
    totals=grouped['sum'].to_numpy();counts=grouped['count'].to_numpy()
    # Batch bootstrap draws to keep memory bounded on ordinary desktop hardware.
    draws=[]
    for start in range(0,BOOTSTRAPS,100):
        chosen=rng.integers(0,g,size=(min(100,BOOTSTRAPS-start),g))
        draws.extend((totals[chosen].sum(axis=1)/counts[chosen].sum(axis=1)).tolist())
    result['ci_low'],result['ci_high']=np.quantile(draws,[.025,.975]).tolist()
    result['ci_note']='Whole-family resampling; record-weighted point estimate and pooled-record bootstrap statistic.'
    return result


def main():
    family=pd.read_csv(FAMILY)
    assert not family.sample_id.duplicated().any()
    cohort=pd.read_csv(OUT/'cohort_records.csv')
    assert len(family)==len(cohort)==3063 and set(family.sample_id)==set(cohort.sample_id)
    clean_components=set(family.groupby('component_id').original_fold.nunique().loc[lambda v:v==1].index)
    family['family_clean_original_fold']=family.component_id.isin(clean_components)
    metrics=pd.read_csv(OUT/'per_record_metrics.csv')
    metric_rows=[]
    for metric in ['rmsf','hchi']:
        for subset in ['all','complete_sidechains','interface_ca8','noninterface_ca8']:
            part=metrics[(metrics.metric==metric)&(metrics.subset==subset)]
            full=part[part.variant=='full'][['sample_id','mae','rho']].rename(columns={'mae':'full_mae','rho':'full_rho'})
            joined=part.merge(full,on='sample_id',validate='many_to_one').merge(family[['sample_id','component_id','family_clean_original_fold']],on='sample_id',validate='many_to_one')
            joined['mae_gain']=joined.mae-joined.full_mae
            joined['rho_gain']=joined.full_rho-joined.rho
            metric_rows.append(joined)
    pairs=pd.concat(metric_rows,ignore_index=True)
    pairs.to_csv(OUT/'paired_full_vs_baseline_per_record.csv',index=False)
    same=pd.read_csv(OUT/'same_type_per_record_metrics.csv')
    full=same[same.variant=='full'][['sample_id','metric','subset','rho_train_aa_centered','rho_within_record_aa_centered']].rename(columns={
        'rho_train_aa_centered':'full_rho_train_aa_centered','rho_within_record_aa_centered':'full_rho_within_record_aa_centered'})
    typed=same.merge(full,on=['sample_id','metric','subset'],validate='many_to_one').merge(family[['sample_id','component_id','family_clean_original_fold']],on='sample_id',validate='many_to_one')
    for statistic in ['rho_train_aa_centered','rho_within_record_aa_centered']:
        typed[statistic+'_gain']=typed['full_'+statistic]-typed[statistic]
    typed.to_csv(OUT/'paired_full_vs_baseline_typed_per_record.csv',index=False)
    results=[]
    for source in ['all','atlas','dynarepo','itrace']:
        for population in ['all_records','family_clean_original_fold']:
            p=pairs if source=='all' else pairs[pairs.dataset==source]
            t=typed if source=='all' else typed[typed.dataset==source]
            if population!='all_records':
                p=p[p.family_clean_original_fold];t=t[t.family_clean_original_fold]
            for (variant,metric,subset),block in p.groupby(['variant','metric','subset'],sort=False):
                if variant=='full':continue
                for value in ['mae_gain','rho_gain']:
                    label=f'{source}|{population}|{variant}|{metric}|{subset}|{value}'
                    results.append(dict(source=source,population=population,baseline=variant,metric=metric,subset=subset,statistic=value,
                                        **interval(block,value,label)))
            for (variant,metric,subset),block in t.groupby(['variant','metric','subset'],sort=False):
                if variant=='full':continue
                for value in ['rho_train_aa_centered_gain','rho_within_record_aa_centered_gain']:
                    label=f'{source}|{population}|{variant}|{metric}|{subset}|{value}'
                    results.append(dict(source=source,population=population,baseline=variant,metric=metric,subset=subset,statistic=value,
                                        **interval(block,value,label)))
    pd.DataFrame(results).to_csv(OUT/'paired_full_vs_baseline_family_ci.csv',index=False)
    # Reproduce original AA/Static/Full/Intra per-record values rather than compare rounded narrative numbers.
    checks=[]
    for origin,name in [(CORE/'static_corrected/oof_per_record_metrics.csv','original_static'),
                        (CORE/'neural_exports/evaluation/per_record_metrics.csv','original_neural')]:
        old=pd.read_csv(origin).rename(columns={'spearman':'rho'})
        old['metric']=old.metric.replace({'chi_entropy':'hchi'})
        old['variant']=old.variant.replace({'static_sequence_complex':'static','intrachain':'intra'})
        new=metrics[metrics.subset=='all']
        joined=old.merge(new,on=['sample_id','variant','metric'],suffixes=('_old','_new'),validate='one_to_one')
        for variant,part in joined.groupby('variant'):
            checks.append(dict(origin=name,variant=variant,matched_rows=len(part),
                               max_mae_abs_difference=float(np.nanmax(abs(part.mae_old-part.mae_new))),
                               max_rho_abs_difference=float(np.nanmax(abs(part.rho_old-part.rho_new)))))
    pd.DataFrame(checks).to_csv(OUT/'original_oof_reproduction_audit.csv',index=False)
    assert all(check['max_mae_abs_difference']<1e-6 for check in checks)
    assert all(check['max_rho_abs_difference']<1e-10 for check in checks)
    report=dict(status='complete',bootstrap_replicates=BOOTSTRAPS,seed_rule='20261025 plus SHA256 comparison label',
                cluster='component_id from frozen record_family_folds.csv',
                family_clean_definition='Every record in component has the same original fold; do not call this a newly trained family split',
                n_clean_records=int(family.family_clean_original_fold.sum()),n_clean_families=len(clean_components),
                positive_mae_gain='baseline MAE minus Full MAE',positive_rho_gain='Full rho minus baseline rho',
                limitations='Family bootstrap represents record dependence, not family-isolated retraining. ITrace has one large family; source-specific CI withheld when fewer than two families.',
                original_reproduction_checks=checks)
    (OUT/'cluster_ci_contract.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
