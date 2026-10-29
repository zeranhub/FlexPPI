from pathlib import Path
import os, sys, json, hashlib, warnings, time
os.environ['OMP_NUM_THREADS']='2'
os.environ['OPENBLAS_NUM_THREADS']='2'
BASE=Path('E:/Project/FlexPPI/experiments/_sources')
sys.path.insert(0,str(BASE/'completed_results_20261012/plot_runtime'))
sys.path.append(str(BASE/'label_sensitivity/runtime'))
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

CORE=Path('E:/Project/FlexPPI/experiments/core_corrected_pipeline')
OUT=Path(__file__).resolve().parent
METHODS=['Full','Intra','Static','AA']
TARGETS=['rmsf','chi_entropy']
GRID=np.linspace(0,1,101)
SEED=20261018

def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()

def curves(y,score):
 # Whole score-tie groups in descending score order, never arbitrary tie ranks.
 order=np.argsort(-score,kind='stable')
 s=score[order]; z=y[order].astype(np.int64)
 end=np.r_[np.flatnonzero(np.diff(s)!=0),len(s)-1]
 tp=np.cumsum(z,dtype=np.int64)[end]; fp=(end+1)-tp
 p=int(y.sum()); n=len(y)-p
 recall=tp/p; precision=tp/(end+1); fpr=fp/n
 ap=float(np.sum(np.diff(np.r_[0.,recall])*precision))
 auc=float(np.trapezoid(np.r_[0.,recall],np.r_[0.,fpr]))
 rr=np.r_[0.,recall]; pp=np.r_[1.,precision]
 sampled=pp[np.searchsorted(rr,GRID,side='left')]
 sampled[0]=1.
 return ap,auc,sampled,dict(threshold=np.r_[np.inf,s[end]],precision=pp,recall=rr,fpr=np.r_[0.,fpr],tpr=rr)

def stratified_draws(sources,n_draw=10000):
 rng=np.random.default_rng(SEED)
 groups=[np.flatnonzero(sources==s) for s in sorted(set(sources))]
 return np.concatenate([rng.choice(ix,size=(n_draw,len(ix)),replace=True) for ix in groups],axis=1).astype(np.int32)

def cis(values,draws):
 # Batch indexing bounds peak RAM and preserves identical paired draws.
 results=[]
 for begin in range(0,len(draws),250):
  results.append(np.nanmean(values[draws[begin:begin+250]],axis=1))
 b=np.concatenate(results)
 return np.quantile(b,[.025,.975],axis=0)

def sym(x,t=1.):
 a=np.abs(x);return np.sign(x)*np.where(a<=t,a/0.9,t*(1/0.9+np.log10(np.maximum(a/t,1))))
def invsym(x,t=1.):
 a=np.abs(x);return np.sign(x)*np.where(a<=t/0.9,a*.9,t*np.power(10,a/t-1/.9))

def main():
 started=time.time()
 assert (OUT/'PRESPECIFIED_PROTOCOL.md').exists()
 inputs={}
 def hashit(p):inputs[str(p)]=sha(p)
 schema_path=CORE/'static_corrected/training_schema.json'
 truth_path=BASE/'baselines/cohort/truth_residues.csv.gz'
 schema=json.loads(schema_path.read_text(encoding='utf8'))['records']
 truth=pd.read_csv(truth_path,dtype={'sample_id':str,'dataset':str,'chain':str,'author_residue_id':str,'amino_acid':str})
 hashit(schema_path);hashit(truth_path)
 assert len(truth)==961321 and len(schema)==3063
 assert not truth[['sample_id','residue_index']].duplicated().any()
 static={}
 for method,fn in [('Static','static_sequence_complex_oof.npz'),('AA','aa_only_oof.npz')]:
  p=CORE/'static_corrected'/fn;hashit(p)
  static[method]=np.load(p)['prediction']
  assert static[method].shape==(len(truth),2)
 feature_path=CORE/'static_corrected/training_features.npz';hashit(feature_path)
 with np.load(feature_path) as z:
  print('Static training arrays',[(k,z[k].shape) for k in z.files],flush=True)
  assert np.array_equal(z['target'][:,:2],truth[['target_rmsf','target_chi_entropy']].to_numpy(dtype=np.float32),equal_nan=True)
  assert np.array_equal(z['valid'][:,:2],truth[['valid_rmsf','valid_chi_entropy']].to_numpy(dtype=bool))
  static_fold=z['fold'].copy();static_record=z['record_index'].copy()
 manifest_path=CORE/'neural_exports/evaluation/prediction_manifest.jsonl';hashit(manifest_path)
 manifest=[json.loads(x) for x in manifest_path.read_text(encoding='utf8').splitlines()]
 manifest={(r['variant'],r['sample_id']):r for r in manifest}
 old_path=CORE/'neural_exports/evaluation/per_record_metrics.csv';hashit(old_path)
 old=pd.read_csv(old_path).set_index(['sample_id','variant','metric'])
 oldstatic_path=CORE/'static_corrected/oof_per_record_metrics.csv';hashit(oldstatic_path)
 oldstatic=pd.read_csv(oldstatic_path).set_index(['sample_id','variant','metric'])
 rows=[];typed=[];records=[];audit=[]
 raw={(t,m):{k:[] for k in ['record_code','threshold','precision','recall','fpr','tpr']} for t in TARGETS for m in METHODS}
 curve_grid={(t,m):[] for t in TARGETS for m in METHODS}
 points={t:{k:[] for k in ['reference','prediction','record_code','amino_acid','fold','residue_index']} for t in TARGETS}
 exclusions=[]
 for code,r in enumerate(schema):
  sid=r['sample_id'];start,stop=r['start'],r['stop'];g=truth.iloc[start:stop];n=stop-start
  assert np.all(g.sample_id==sid) and np.all(g.dataset==r['dataset']) and np.all(g.fold==r['fold'])
  assert np.array_equal(g.residue_index,np.arange(n))
  assert np.all(static_fold[start:stop]==r['fold']) and np.all(static_record[start:stop]==r['record_index'])
  with np.load(CORE/'neural_exports/evaluation/full'/(sid+'.npz')) as z:
   full={k:z[k].copy() for k in z.files}
  with np.load(CORE/'neural_exports/evaluation/intrachain'/(sid+'.npz')) as z:
   intra={k:z[k].copy() for k in z.files}
  for variant,cache in [('full',full),('intrachain',intra)]:
   p=CORE/'neural_exports/evaluation'/variant/(sid+'.npz');hashit(p)
   assert cache['mean'].shape==(n,2) and cache['target'].shape==(n,3)
   assert np.array_equal(cache['amino_acid'].astype(str),g.amino_acid.astype(str))
   assert np.array_equal(cache['chain'].astype(str),g.chain.astype(str))
   assert np.array_equal(cache['author_residue_id'].astype(str),g.author_residue_id.astype(str))
   assert np.array_equal(cache['target'][:,:2],g[['target_rmsf','target_chi_entropy']].to_numpy(dtype=np.float32),equal_nan=True)
   assert np.array_equal(cache['valid'][:,:2],g[['valid_rmsf','valid_chi_entropy']].to_numpy(dtype=bool))
   mr=manifest[(variant,sid)]
   assert mr['dataset']==r['dataset'] and mr['fold']==r['fold'] and mr['label_sha256']==r['label_sha256']
  predictions={'Full':full['mean'],'Intra':intra['mean'],'Static':static['Static'][start:stop],'AA':static['AA'][start:stop]}
  records.append(dict(record_code=code,sample_id=sid,dataset=r['dataset'],fold=r['fold'],n_residues=n))
  for j,t in enumerate(TARGETS):
   valid=full['valid'][:,j]&intra['valid'][:,j]&np.isfinite(full['target'][:,j])
   original_n=int(valid.sum())
   for pred in predictions.values():valid&=np.isfinite(pred[:,j])
   assert valid.any()
   if int(valid.sum())!=original_n:exclusions.append(dict(sample_id=sid,target=t,excluded=original_n-int(valid.sum()),reason='nonfinite_method_prediction'))
   y=full['target'][valid,j].astype(np.float64);aa=full['amino_acid'][valid].astype('U1');nr=len(y)
   q=float(np.quantile(y,.9,method='linear'));binary=y>=q;p=int(binary.sum());neg=nr-p
   reason='' if p and neg else 'one_class_after_reference_q90_ties'
   points[t]['reference'].append(y.astype(np.float32));points[t]['prediction'].append(full['mean'][valid,j])
   points[t]['record_code'].append(np.full(nr,code,dtype=np.int32));points[t]['amino_acid'].append(aa)
   points[t]['fold'].append(np.full(nr,r['fold'],dtype=np.int8));points[t]['residue_index'].append(np.flatnonzero(valid).astype(np.int32))
   for method in METHODS:
    pred=predictions[method][valid,j].astype(np.float64)
    mae=float(np.mean(np.abs(pred-y)))
    rho=float(spearmanr(y,pred).statistic) if nr>=4 and np.ptp(y)>0 and np.ptp(pred)>0 else np.nan
    if method in ['Full','Intra']:
     previous=old.loc[(sid,{'Full':'full','Intra':'intrachain'}[method],t)]
    else:previous=oldstatic.loc[(sid,{'Static':'static_sequence_complex','AA':'aa_only'}[method],t)]
    assert int(previous['count'])==nr and abs(mae-float(previous['mae']))<2e-6,(sid,t,method,mae,previous['mae'])
    if np.isfinite(rho):assert abs(rho-float(previous['spearman']))<1e-9,(sid,t,method,'rho')
    if not reason:
     ap,auc,grid,curve=curves(binary,pred)
     curve_grid[(t,method)].append((code,grid))
     for k,a in curve.items():raw[(t,method)][k].append(a)
     raw[(t,method)]['record_code'].append(np.full(len(curve['recall']),code,dtype=np.int32))
    else:ap=auc=np.nan
    rows.append(dict(record_code=code,sample_id=sid,dataset=r['dataset'],fold=r['fold'],target=t,method=method,n_positions=nr,threshold_q90=q,n_positive=p,n_negative=neg,prevalence=p/nr,eligible=not bool(reason),undefined_reason=reason,ap=ap,auroc=auc,mae=mae,spearman=rho))
   if t=='chi_entropy':
    for a in 'CDEFHIKLMNPQRSTVWY':
     keep=aa==a;ny=int(keep.sum());yy=y[keep]
     tq=float(np.quantile(yy,.9,method='linear')) if ny else np.nan;by=yy>=tq;np_=int(by.sum())
     why='fewer_than_four_positions' if ny<4 else ('one_class_after_reference_q90_ties' if np_==ny or np_==0 else '')
     for method in METHODS:
      sc=predictions[method][valid,j][keep].astype(np.float64)
      if not why:
       aap,aauc,_,_=curves(by,sc)
       if method=='AA':
        assert np.ptp(sc)==0 and abs(aap-np_/ny)<1e-12 and abs(aauc-.5)<1e-12
      else:aap=aauc=np.nan
      typed.append(dict(record_code=code,sample_id=sid,dataset=r['dataset'],fold=r['fold'],amino_acid=a,n_positions=ny,threshold_q90=tq,n_positive=np_,prevalence=np_/ny if ny else np.nan,method=method,eligible=not bool(why),undefined_reason=why,ap=aap,auroc=aauc))
  if (code+1)%250==0:print('Identity, mask and per-record metrics verified',code+1,'/',len(schema),'seconds',round(time.time()-started,1),flush=True)
 pd.DataFrame(records).to_csv(OUT/'record_codes.csv',index=False)
 df=pd.DataFrame(rows);df.to_csv(OUT/'per_record_metrics.csv',index=False)
 typ=pd.DataFrame(typed);typ.to_csv(OUT/'within_type_per_record_metrics.csv',index=False)
 print('Saving threshold curves',flush=True)
 for (t,m),data in raw.items():
  np.savez_compressed(OUT/f'raw_curves_{t}_{m}.npz',**{k:np.concatenate(v) for k,v in data.items()})
 raw.clear()
 sources=np.array([r['dataset'] for r in records]);draws=stratified_draws(sources)
 curve_weights=np.zeros((1000,len(schema)),dtype=np.float64)
 for b,ix in enumerate(draws[:1000]):curve_weights[b]=np.bincount(ix,minlength=len(schema))
 summary={'status':'complete','protocol_sha256':sha(OUT/'PRESPECIFIED_PROTOCOL.md'),'code_sha256':sha(__file__),'seed':SEED,'bootstrap_replicates':10000,'curve_bootstrap_replicates':1000,'statistical_unit':'record','bootstrap_strata':sorted(set(sources)),'metrics':{},'within_type':{},'coverage':{}}
 macro_rows=[]
 for t in TARGETS:
  summary['metrics'][t]={}
  for m in METHODS:
   g=df[(df.target==t)&(df.method==m)].sort_values('record_code')
   eligible=g.eligible.to_numpy();vals=g[['ap','auroc']].to_numpy(dtype=float)
   ci=cis(vals,draws)
   codes=np.array([x[0] for x in curve_grid[(t,m)]]);grids=np.stack([x[1] for x in curve_grid[(t,m)]])
   w=curve_weights[:,codes];macro=grids.mean(axis=0);boot=(w@grids)/w.sum(axis=1)[:,None];cic=np.quantile(boot,[.025,.975],axis=0)
   prevalence=float(g.loc[g.eligible,'prevalence'].mean())
   entry=dict(mean_ap=float(np.nanmean(vals[:,0])),ap_ci=ci[:,0].tolist(),mean_auroc=float(np.nanmean(vals[:,1])),auroc_ci=ci[:,1].tolist(),macro_mae=float(g.mae.mean()),mean_spearman=float(g.spearman.mean()),n_records=int(eligible.sum()),n_total_records=len(g),n_positions=int(g.n_positions.sum()),prevalence_mean=prevalence,prevalence_min=float(g.loc[g.eligible,'prevalence'].min()),prevalence_max=float(g.loc[g.eligible,'prevalence'].max()),undefined_records=int((~eligible).sum()))
   summary['metrics'][t][m]=entry
   for ri,re in enumerate(GRID):macro_rows.append(dict(target=t,method=m,recall=float(re),precision=float(macro[ri]),precision_ci_low=float(cic[0,ri]),precision_ci_high=float(cic[1,ri]),n_records=int(eligible.sum()),no_skill_precision=prevalence))
   print(t,m,entry,flush=True)
  gg=df[(df.target==t)&(df.method=='Full')]
  summary['coverage'][t]={'n_total_records':len(gg),'n_eligible_records':int(gg.eligible.sum()),'n_valid_positions':int(gg.n_positions.sum()),'undefined_reasons':gg.loc[~gg.eligible,'undefined_reason'].value_counts().to_dict()}
 pd.DataFrame(macro_rows).to_csv(OUT/'macro_pr_curve.csv',index=False)
 typed_macro=typ[typ.eligible].groupby(['record_code','sample_id','dataset','fold','method'],as_index=False)[['ap','auroc','prevalence']].mean()
 typed_macro.to_csv(OUT/'within_type_record_macro_metrics.csv',index=False)
 for m in METHODS:
  g=typed_macro[typed_macro.method==m].set_index('record_code').reindex(range(len(schema)))
  vals=g[['ap','auroc']].to_numpy(dtype=float);ci=cis(vals,draws)
  summary['within_type'][m]={'mean_ap':float(np.nanmean(vals[:,0])),'ap_ci':ci[:,0].tolist(),'mean_auroc':float(np.nanmean(vals[:,1])),'auroc_ci':ci[:,1].tolist(),'n_records':int(np.isfinite(vals[:,0]).sum()),'n_eligible_record_types':int(((typ.method==m)&typ.eligible).sum()),'undefined_record_types':typ.loc[(typ.method==m)&~typ.eligible,'undefined_reason'].value_counts().to_dict(),'mean_record_type_prevalence':float(g.prevalence.mean())}
 print('Saving record-equal full density and conditional trends',flush=True)
 trend_rows=[]
 for t in TARGETS:
  data={k:np.concatenate(v) for k,v in points[t].items()}
  ref=data['reference'].astype(np.float64);pred=data['prediction'].astype(np.float64);codes=data['record_code']
  counts=np.bincount(codes,minlength=len(schema));weights=1./(len(schema)*counts[codes])
  lo=float(min(0.,ref.min(),pred.min()));hi=float(max(ref.max(),pred.max()))
  if t=='rmsf':
   edges=invsym(np.linspace(sym(np.array(lo)),sym(np.array(hi)),91));scale='symlog';linthresh=1.
  else:edges=np.linspace(lo,hi,91);scale='linear';linthresh=0.
  edges[0]=np.nextafter(lo,-np.inf);edges[-1]=np.nextafter(hi,np.inf)
  hist,xe,ye=np.histogram2d(ref,pred,bins=(edges,edges),weights=weights)
  assert abs(hist.sum()-1)<1e-10 and len(ref)==summary['metrics'][t]['Full']['n_positions']
  # Reference bins for conditional means: equal spacing in the same displayed axis transform.
  trend_edges=invsym(np.linspace(sym(np.array(lo)),sym(np.array(hi)),21)) if t=='rmsf' else np.linspace(lo,hi,21)
  trend_edges[0]=np.nextafter(lo,-np.inf);trend_edges[-1]=np.nextafter(hi,np.inf)
  bins=np.searchsorted(trend_edges,ref,side='right')-1
  den=np.zeros((len(schema),20));num=den.copy();xnum=den.copy()
  np.add.at(den,(codes,bins),1./counts[codes]);np.add.at(num,(codes,bins),pred/counts[codes]);np.add.at(xnum,(codes,bins),ref/counts[codes])
  db=curve_weights@den;nb=curve_weights@num
  with warnings.catch_warnings():
   warnings.simplefilter('ignore',RuntimeWarning)
   means=num.sum(axis=0)/den.sum(axis=0);xm=xnum.sum(axis=0)/den.sum(axis=0)
   trci=np.nanquantile(nb/db,[.025,.975],axis=0)
  for b in range(20):
   trend_rows.append(dict(target=t,bin_index=b,bin_left=float(trend_edges[b]),bin_right=float(trend_edges[b+1]),x_mean=float(xm[b]),prediction_mean=float(means[b]),ci_low=float(trci[0,b]),ci_high=float(trci[1,b]),n_positions=int((bins==b).sum()),n_records=int((den[:,b]>0).sum()),record_equal_mass=float(den[:,b].sum()/len(schema))))
  np.savez_compressed(OUT/f'prediction_reference_{t}.npz',**data,histogram_mass=hist,x_edges=xe,y_edges=ye,density_axis_scale=np.array(scale),linthresh=np.array(linthresh))
  summary['coverage'][t].update(reference_min=float(ref.min()),reference_max=float(ref.max()),prediction_min=float(pred.min()),prediction_max=float(pred.max()),density_weight_sum=float(hist.sum()),density_axis_scale=scale,density_linthresh=linthresh)
 pd.DataFrame(trend_rows).to_csv(OUT/'reference_binned_trend.csv',index=False)
 (OUT/'summary.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False,allow_nan=False),encoding='utf8')
 audit_data={'status':'passed','records':len(schema),'rows':len(truth),'all_identity_target_mask_assertions_passed':True,'static_schema_and_features_order_verified':True,'all_4_methods_original_per_record_MAE_and_spearman_reproduced':True,'exclusions':exclusions,'within_type_AA_constant_score_AP_prevalence_AUC_half_verified':True,'inputs_sha256':inputs,'code_sha256':sha(__file__),'protocol_sha256':sha(OUT/'PRESPECIFIED_PROTOCOL.md'),'runtime_seconds':time.time()-started}
 (OUT/'input_alignment_audit.json').write_text(json.dumps(audit_data,indent=2,ensure_ascii=False),encoding='utf8')
 print('COMPLETE',round(time.time()-started,1),'seconds',flush=True)

if __name__=='__main__':main()
