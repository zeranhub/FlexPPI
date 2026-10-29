from pathlib import Path
import sys,json,ast
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'completed_results_20261012/plot_runtime'))
import pandas as pd,numpy as np
SOURCE=ROOT/'tcr_case_v40_20261017'
OUT=Path(__file__).resolve().parent
meta=pd.read_csv(SOURCE/'local_pair_screen/biological_candidate_pairs.csv').set_index(['a','b'])
catalog={r['pdb_id']:r for r in json.loads((SOURCE/'biology_sources/complexes_catalog_summary.json').read_text(encoding='utf8'))}

def corr(x,y):
    return np.corrcoef(pd.Series(x).rank(),pd.Series(y).rank())[0,1] if len(x)>=4 and np.ptp(x)>0 and np.ptp(y)>0 else np.nan

def summarize(f,scope,mask):
    out=[]
    for metric in ['rmsf','hchi']:
        g=f[mask&f[metric+'_valid']].copy()
        x=g[metric+'_md_delta'].to_numpy();y=g[metric+'_full_delta'].to_numpy();sem=g[metric+'_delta_sem'].to_numpy()
        if not len(g):continue
        robust=abs(x)>2*sem; nz=abs(x)>(.1 if metric=='rmsf' else .2)
        z={'a':f.a.iloc[0],'b':f.b.iloc[0],'category':f.category.iloc[0],'scope':scope,'metric':metric,'n':len(g),'rho_delta':corr(x,y),'mae_delta':np.mean(abs(x-y)),'zero_mae':np.mean(abs(x)),'gain_over_zero':1-np.mean(abs(x-y))/np.mean(abs(x)) if np.mean(abs(x))>0 else np.nan,'mean_md_delta':np.mean(x),'mean_model_delta':np.mean(y),'robust_n':int(robust.sum()),'robust_direction_fraction':np.mean(x[robust]*y[robust]>0) if robust.any() else np.nan,'robust_mae':np.mean(abs(x[robust]-y[robust])) if robust.any() else np.nan,'robust_zero_mae':np.mean(abs(x[robust])) if robust.any() else np.nan,'substantial_n':int(nz.sum()),'substantial_direction_fraction':np.mean(x[nz]*y[nz]>0) if nz.any() else np.nan}
        z['robust_gain_over_zero']=1-z['robust_mae']/z['robust_zero_mae'] if robust.any() else np.nan
        # Mean changes can look good when all sites shift similarly. This separates spatial pattern from a common shift.
        z['centered_mae']=np.mean(abs((x-x.mean())-(y-y.mean())))
        z['constant_shift_mae']=np.mean(abs(x-x.mean()))
        z['gain_over_observed_constant_shift']=1-z['centered_mae']/z['constant_shift_mae'] if z['constant_shift_mae']>0 else np.nan
        z['md_delta_sd']=np.std(x);z['model_delta_sd']=np.std(y)
        out.append(z)
    return out

rows=[]
def process(f):
    same=f.same_aa
    scopes={'all_unchanged':same,'peptide_unchanged':same&(f.role=='peptide'),'tcr_interface_unchanged':same&f.interface&f.role.isin(['tcr_alpha','tcr_beta']),'cdr3_unchanged':same&f.region_a.str.contains('CDR3')&f.region_b.str.contains('CDR3'),'mhc_groove_unchanged':same&(f.role=='mhc')&(f.position_a<=180)&(f.position_b<=180),'all_interface_unchanged':same&f.interface}
    for region in ['CDR1α','CDR1β','CDR2α','CDR2β','CDR3α','CDR3β','α1 helix','α2 helix','MHC floor']:
        scopes[region+'_unchanged']=same&(f.region_a==region)&(f.region_b==region)
        if region=='MHC floor':
            scopes[region+'_unchanged'] &= (f.role=='mhc')&(f.position_a<=180)&(f.position_b<=180)
    for scope,mask in scopes.items():rows.extend(summarize(f,scope,mask))

cols=['a','b','category','role','position_a','position_b','same_aa','region_a','region_b','interface']
cols += [metric+suffix for metric in ['rmsf','hchi'] for suffix in ['_valid','_md_delta','_full_delta','_delta_sem']]
pending=pd.DataFrame()
for chunk in pd.read_csv(SOURCE/'change_screen/all_candidate_changes.csv.gz',usecols=cols,chunksize=8000):
    chunk=pd.concat([pending,chunk],ignore_index=True)
    last=(chunk.a.iloc[-1],chunk.b.iloc[-1]); pending=chunk[(chunk.a==last[0])&(chunk.b==last[1])]
    for _,frame in chunk[~((chunk.a==last[0])&(chunk.b==last[1]))].groupby(['a','b'],sort=False):process(frame)
if len(pending):process(pending)
res=pd.DataFrame(rows)
res.to_csv(OUT/'all_predefined_region_metrics.csv',index=False,encoding='utf-8-sig')
pivot=res.pivot(index=['a','b','category','scope'],columns='metric')
pivot.columns=['_'.join(c) for c in pivot.columns]
pivot=pivot.reset_index()
pivot['min_gain']=pivot[['gain_over_zero_rmsf','gain_over_zero_hchi']].min(axis=1)
pivot['mean_gain']=pivot[['gain_over_zero_rmsf','gain_over_zero_hchi']].mean(axis=1)
pivot['min_rho']=pivot[['rho_delta_rmsf','rho_delta_hchi']].min(axis=1)
pivot['mean_rho']=pivot[['rho_delta_rmsf','rho_delta_hchi']].mean(axis=1)
pivot['min_robust_n']=pivot[['robust_n_rmsf','robust_n_hchi']].min(axis=1)
pivot['min_n']=pivot[['n_rmsf','n_hchi']].min(axis=1)
pivot['dual_gain_positive']=(pivot.gain_over_zero_rmsf>0)&(pivot.gain_over_zero_hchi>0)
pivot['dual_robust_two']=(pivot.robust_n_rmsf>=2)&(pivot.robust_n_hchi>=2)
pivot['rank_score']=pivot.min_gain+.25*pivot.mean_rho+.05*np.log1p(pivot.min_robust_n)
def metadata(row):
    vals={}
    for side in ['a','b']:
        c=catalog.get(row[side],{})
        vals['title_'+side]=c.get('title',''); vals['peptide_'+side]=c.get('peptide_seq',''); vals['hla_'+side]=c.get('hla_allele','')
        citation=ast.literal_eval(c.get('citation','{}'));vals['doi_'+side]=citation.get('doi','');vals['citation_title_'+side]=citation.get('title','')
    if (row.a,row.b) in meta.index:
        for k in ['peptide_diff','mhc_diff','tcr_alpha_diff','tcr_beta_diff']: vals[k]=meta.loc[(row.a,row.b),k]
    return pd.Series(vals)
pivot=pd.concat([pivot,pivot.apply(metadata,axis=1)],axis=1)
pivot.sort_values('rank_score',ascending=False).to_csv(OUT/'joint_region_candidate_rank.csv',index=False,encoding='utf-8-sig')
qualified=pivot[pivot.dual_gain_positive&(pivot.min_n>=5)&(pivot.rho_delta_rmsf>=.3)&(pivot.rho_delta_hchi>=.3)&(pivot.zero_mae_rmsf>=.08)&(pivot.zero_mae_hchi>=.12)].sort_values(['dual_robust_two','rank_score'],ascending=False)
qualified.to_csv(OUT/'joint_region_shortlist.csv',index=False,encoding='utf-8-sig')
print(qualified[['a','b','scope','n_rmsf','n_hchi','rho_delta_rmsf','rho_delta_hchi','gain_over_zero_rmsf','gain_over_zero_hchi','robust_n_rmsf','robust_n_hchi','rank_score']].head(35).to_string(index=False),flush=True)
(OUT/'rank_audit.json').write_text(json.dumps({'status':'complete','pair_n':res[['a','b']].drop_duplicates().shape[0],'predefined_scopes':res.scope.unique().tolist(),'source':'existing 251 pair screen, current Full family-disjoint model','change_definition':'b-a; only identical-amino-acid aligned positions enter primary ranking','zero_baseline':'All differences predicted zero; gain=1-MAE_model/mean_abs_MD_change','robust_screen':'abs(mean_MD_change)>2*sqrt((SD_a^2+SD_b^2)/3); descriptive screening, not formal significance','substantial_screen':'abs(mean_MD_change)>.1 Angstrom RMSF or>.2 nats Hchi','selection_caveat':'This is retrospective candidate discovery. A candidate passing one local region is not evidence of accuracy across the whole complex.','qualified_region_count':len(qualified)},ensure_ascii=False,indent=2),encoding='utf8')
