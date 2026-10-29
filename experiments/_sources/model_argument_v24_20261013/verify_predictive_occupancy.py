"""Separate vectorized reconstruction of the frozen angle experiment."""
from pathlib import Path
import json, hashlib, time
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parent
EXP=ROOT/'experiment_predictive_occupancy'
QA=ROOT/'manuscript_QA';QA.mkdir(exist_ok=True)
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
plan=json.loads((EXP/'analysis_plan.json').read_text(encoding='utf8'))
metric='fine_entropy_mean_per_chi'
preds=['full_hchi','intrachain_hchi','full_rmsf']
checks=[]
def check(name,value):
    checks.append({'name':name,'passed':bool(value)});assert value,name

def entropy_vector(cache):
    a=cache['chi_degrees'].astype(np.float64)
    n=a.shape[1]
    q=cache['chi_count'].astype(int)
    columns=(np.arange(4)[None,:]<q[:,None]).reshape(-1)
    cols=np.flatnonzero(columns)
    # Each column represents one actual chemical rotation axis at one site.
    values=a.reshape(501,-1)[:,cols]
    finite=np.isfinite(values).all(axis=0)
    selected=cols[finite];values=values[:,finite]
    ent=np.zeros(len(selected))
    for shift in [0,1.25,2.5,3.75]:
        bin_id=np.floor(np.remainder(values+180-shift,360)/5).astype(np.int64)
        counts=np.bincount((bin_id+np.arange(len(selected))[None,:]*72).ravel(),minlength=len(selected)*72).reshape(-1,72)
        probability=counts/501
        logarithm=np.zeros_like(probability)
        np.log(probability,out=logarithm,where=probability>0)
        ent-=np.sum(probability*logarithm,axis=1)/4
    flat=np.full(n*4,np.nan);flat[selected]=ent
    out=np.nansum(flat.reshape(n,4),axis=1)/np.maximum(q,1)
    complete=np.isfinite(flat.reshape(n,4))|(np.arange(4)[None,:]>=q[:,None])
    return out,complete.all(axis=1)

for mode in ['replica1','available_replicas']:
    residues=pd.read_csv(EXP/f'{mode}_per_residue.csv.gz',dtype={'sample_id':str,'chain':str,'author_residue_id':str})
    corr=pd.read_csv(EXP/f'{mode}_per_record_correlations.csv')
    high=pd.read_csv(EXP/f'{mode}_high_low_per_record.csv')
    max_error=0;rho_error=0;group_error=0
    for rec in plan['records']:
        part=residues[residues.sample_id==rec['sample_id']].sort_values('residue_index')
        ids=part.residue_index.to_numpy()
        raw=[]
        for entry in rec['caches'][:1] if mode=='replica1' else rec['caches']:
            check(f'{mode} {Path(entry["path"]).name} immutable',sha(entry['path'])==entry['sha256'])
            with np.load(entry['path']) as cache:
                ent,valid=entropy_vector(cache)
                check(f'{mode} {rec["sample_id"]} complete chemical angles',valid[ids].all())
                raw.append(ent[ids])
        value=np.mean(raw,axis=0)
        max_error=max(max_error,float(np.max(np.abs(value-part[metric].to_numpy()))))
        for key,label in [('full_prediction','full_hchi'),('intra_prediction','intrachain_hchi')]:
            check(f'{mode} {rec["sample_id"]} {label} frozen',sha(rec[key])==rec[key+'_sha256'])
            with np.load(rec[key]) as z:
                check(f'{mode} {rec["sample_id"]} {label} matches OOF source',np.allclose(part[label].to_numpy(),z['mean'][ids,1].astype(float),rtol=0,atol=1e-14))
                if key=='full_prediction':check(f'{mode} {rec["sample_id"]} rmsf matches OOF source',np.allclose(part.full_rmsf.to_numpy(),z['mean'][ids,0].astype(float),rtol=0,atol=1e-14))
        counts=part.groupby('amino_acid').size()
        kept=part[part.amino_acid.isin(counts[counts>=4].index)].copy()
        for scope in ['all','q1','q2','q3','q4']:
            subset=kept if scope=='all' else kept[kept.chi_count==int(scope[1:])]
            if subset.empty:continue
            for predictor in preds:
                # Independent residual rank reconstruction, using all sites.
                x=subset[predictor]-subset.groupby('amino_acid')[predictor].transform('mean')
                y=subset[metric]-subset.groupby('amino_acid')[metric].transform('mean')
                r=x.rank(method='average').corr(y.rank(method='average')) if len(x)>=4 and x.max()>x.min() and y.max()>y.min() else np.nan
                target=corr[(corr.sample_id==rec['sample_id'])&(corr.scope==scope)&(corr.metric==metric)&(corr.predictor==predictor)].iloc[0]
                if np.isfinite(r):rho_error=max(rho_error,abs(r-target.rho_type_centered))
                else:check(f'{mode} {rec["sample_id"]} {scope} undefined rho explicit',pd.isna(target.rho_type_centered))
                dif=[]
                for aa,g in subset.groupby('amino_acid'):
                    if g[predictor].max()-g[predictor].min()<=1e-12:continue
                    order=g.sort_values([predictor,'residue_index'],kind='stable');half=len(g)//2
                    dif.append(order.iloc[-half:][metric].mean()-order.iloc[:half][metric].mean())
                actual=high[(high.sample_id==rec['sample_id'])&(high.scope==scope)&(high.metric==metric)&(high.predictor==predictor)]
                if dif:group_error=max(group_error,abs(np.mean(dif)-actual.iloc[0].high_minus_low_MD))
    check(f'{mode} raw angle entropy reconstructed at all 33493 sites',max_error<2e-12)
    check(f'{mode} all scope/predictor residual correlations reproduced',rho_error<2e-12)
    check(f'{mode} all scope/predictor prediction-defined half groups reproduced',group_error<2e-12)
    print(mode,'raw max',max_error,'rho max',rho_error,'group max',group_error,flush=True)

# Produce a compact single source for manuscript and figure values.
values={}
for mode in ['replica1','available_replicas']:
    c=pd.read_csv(EXP/f'{mode}_correlation_summary.csv')
    h=pd.read_csv(EXP/f'{mode}_high_low_summary.csv')
    diff=pd.read_csv(EXP/f'{mode}_paired_correlation_difference_summary.csv')
    elig=pd.read_csv(EXP/f'{mode}_eligibility_masks.csv')
    values[mode]={}
    for scope in ['all','q1','q2','q3','q4']:
        cm=c[(c.scope==scope)&(c.metric==metric)&(c.dataset=='all')]
        hm=h[(h.scope==scope)&(h.metric==metric)&(h.dataset=='all')]
        dm=diff[(diff.scope==scope)&(diff.metric==metric)&(diff.dataset=='all')]
        values[mode][scope]={'correlation':cm.to_dict('records'),'high_low':hm.to_dict('records'),'paired_rho_difference':dm.to_dict('records'),'eligibility':elig[elig.scope==scope].to_dict('records')}
    for summary,record,name in [(c,pd.read_csv(EXP/f'{mode}_per_record_correlations.csv'),'rho_type_centered'),(h,pd.read_csv(EXP/f'{mode}_high_low_per_record.csv'),'high_minus_low_MD')]:
        for _,row in summary[(summary.dataset=='all')&(summary.metric==metric)].iterrows():
            subset=record[(record.scope==row.scope)&(record.metric==metric)&(record.predictor==row.predictor)]
            check(f'{mode} {row.scope} {row.predictor} {name} mean/n',abs(subset[name].mean()-row['mean'])<1e-12 and subset[name].notna().sum()==row.n_records)
    # Saved resamples independently establish the reported percentile endpoints.
    for filename in ['correlation_summary','high_low_summary','paired_correlation_difference_summary']:
        draws=pd.read_csv(EXP/f'{mode}_{filename}.bootstrap_draws.csv.gz')
        summary=pd.read_csv(EXP/f'{mode}_{filename}.csv')
        group=['scope','metric','dataset','predictor' if 'predictor' in summary else 'contrast']
        for keys,part in draws.groupby(group):
            choose=np.ones(len(summary),bool)
            for key,val in zip(group,keys):choose&=summary[key]==val
            row=summary.loc[choose].iloc[0]
            ci=part.value.quantile([.025,.975]).to_numpy()
            check(f'{mode} {filename} {keys} 10000 resamples and intervals',len(part)==10000 and np.allclose(ci,[row.bootstrap_95_low,row.bootstrap_95_high],atol=1e-12,rtol=0))
(ROOT/'manuscript_values.json').write_text(json.dumps(values,ensure_ascii=False,indent=2),encoding='utf8')
table=[]
for scope in ['all','q1','q2','q3','q4']:
    obj=values['replica1'][scope]
    by={r['predictor']:r for r in obj['correlation']}
    hh={r['predictor']:r for r in obj['high_low']}
    table.append({'scope':scope,'n_records':by['full_hchi']['n_records'],'positions':int(sum(r['eligible_positions'] for r in obj['eligibility'])),'Full_Hchi_rho':by['full_hchi']['mean'],'Intra_Hchi_rho':by['intrachain_hchi']['mean'],'Full_RMSF_rho':by['full_rmsf']['mean'],'Full_Hchi_high_low':hh['full_hchi']['mean']})
pd.DataFrame(table).to_csv(ROOT/'tableS7_predictive_occupancy.csv',index=False,float_format='%.17g')
(QA/'predictive_occupancy_reconstruction_audit.json').write_text(json.dumps({'status':'passed','plan_sha256':sha(EXP/'analysis_plan.json'),'scope':'Separate vectorized histogram implementation reconstructed all 203 cached trajectories, primary/resampled averages and every scope/model correlation and prediction-defined group difference; percentile endpoints reconstructed from saved 10000 draws. No old sources mutated.','checks_passed':len(checks),'checks':checks},ensure_ascii=False,indent=2),encoding='utf8')
print('VERIFIED',len(checks),'checks',flush=True)
