"""Fixed, group-held-out downstream comparison. No FlexPPI retraining."""
from pathlib import Path
import argparse,hashlib,json,sys,time,warnings
import numpy as np,pandas as pd,scipy,sklearn
from scipy.stats import spearmanr
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score,roc_auc_score
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
warnings.filterwarnings('ignore',category=FutureWarning,module='sklearn')

R=Path(r'E:\Project\FlexPPI\experiments\_sources\antigen_incremental_value_20261028')
p=argparse.ArgumentParser();p.add_argument('--variant',choices=['full','intrachain'],required=True);p.add_argument('--template-version',choices=['primary','secondary'],required=True);args=p.parse_args()
V=args.variant;TV=args.template_version
O=R/'evaluation'/f'{V}_{TV}';O.mkdir(parents=True,exist_ok=True)
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
states=json.loads((R/'states.json').read_text());states=[s for s in states if s['template_version']==TV]
assert len(states)==282
missing=[s['id'] for s in states if not (R/'predictions'/V/(s['id']+'.npz')).exists()]
if missing:raise RuntimeError(f'{len(missing)} required predictions absent; first={missing[:3]}')
number=pd.read_csv(R/'tcr_numbering.csv',dtype={'author':str})
absolute={};absolute_records=[]
for s in states:
 z=np.load(R/'predictions'/V/(s['id']+'.npz'))
 assert str(z['input_sha256'])==s['sha256']
 lookup={(ch,a):i for i,(ch,a) in enumerate(zip(z['chain'],z['author']))}
 features={}
 for role,ch in [('alpha','D'),('beta','E')]:
  for name,regions in [('CDR3',['CDR3']),('CDR12',['CDR1','CDR2'])]:
   n=number[(number.template==s['template'])&(number.chain==ch)&number.region.isin(regions)]
   ix=[lookup[(r.chain,r.author)] for r in n.itertuples()]
   assert len(ix)>0
   for j,metric in enumerate(['rmsf','hchi']):
    good=z['model_mask']&z['geometry_mask'] if metric=='rmsf' else z['hchi_valid']
    sel=[i for i in ix if good[i]];assert sel
    key=f'{role}_{name}_{metric}';features[key]=float(np.mean(z['mean'][sel,j]))
    absolute_records.append({'state_id':s['id'],'candidate_id':s['candidate_id'],'system':s['system'],'template':s['template'],'feature':key,'mean':features[key],'n':len(sel)})
 absolute[s['id']]=features
dynamic=[]
for s in states:
 if s['is_reference']:continue
 dyn={'candidate_id':s['candidate_id']}
 for k,a in absolute[s['id']].items():dyn[k]=a-absolute[s['reference_state']][k]
 dynamic.append(dyn)
pd.DataFrame(absolute_records).to_csv(O/'absolute_regional_checks.csv',index=False)
d=pd.read_csv(R/'functional_candidates.csv');assert len(d)==279 and d.candidate_id.is_unique
n=pd.read_csv(R/'baselines/nettcr/scores.csv')
e=pd.read_csv(R/'baselines/ergo/all_ergo_scores.csv')
e=e.pivot(index='candidate_id',columns='model',values='binding_score').reset_index()
d=d.merge(n[['candidate_id','NetTCR22_selected']],on='candidate_id',validate='one_to_one').merge(e,on='candidate_id',validate='one_to_one').merge(pd.DataFrame(dynamic),on='candidate_id',validate='one_to_one')
assert len(d)==279
AA=list('ARNDCQEGHILKMFPSTWYV')
extra={a:None for a in []}
for a in AA:extra['frac_'+a]=d.sequence.map(lambda x:x.count(a)/len(x))
position_cols=[]
for i in range(12):
 for a in AA:
  k=f'pos_{i+1}_{a}';position_cols.append(k);extra[k]=d.sequence.map(lambda x:int(i<len(x) and x[i]==a))
d=pd.concat([d,pd.DataFrame(extra,index=d.index)],axis=1)
BASE=['NetTCR22_selected','ERGO-II_mcpas','ERGO-II_vdjdb']
STATIC=['tcrdock_pae','af3_interface_pae','mhc_el_rank_binder_0p5','mhc_el_rank_binder_2p0','min_hamming_distance_from_wildtypes']
DYNAMIC=list(dynamic[0])[1:]
FRACTION=['frac_'+a for a in AA]
MODELS={'B0_binding':BASE,'B1_binding_static':BASE+STATIC,'B2_binding_static_dynamics':BASE+STATIC+DYNAMIC,
        'B1_AA':BASE+STATIC+FRACTION,'B2_AA':BASE+STATIC+FRACTION+DYNAMIC,
        'B1_position_AA':BASE+STATIC+position_cols,'B2_position_AA':BASE+STATIC+position_cols+DYNAMIC}
features=d[['candidate_id','system','sequence','functional_positive','functional_response']+list(dict.fromkeys(BASE+STATIC+DYNAMIC+FRACTION+position_cols))]
features.to_csv(O/'features_and_labels.csv',index=False)
y=d.functional_positive.to_numpy(dtype=int)
def model():
 return Pipeline([('impute',SimpleImputer(strategy='median',add_indicator=True)),('scale',StandardScaler()),('logit',LogisticRegression(C=1,penalty='l2',solver='lbfgs',class_weight='balanced',max_iter=2000,random_state=0))])
def score_metrics(ix,scores,model_name,split,fold):
 f=d.iloc[ix].copy();f['score']=scores
 ranked=f.sort_values(['score','candidate_id'],ascending=[False,True],kind='stable');top=ranked.head(10)
 return {'split':split,'fold':str(fold),'model':model_name,'n':len(f),'n_positive':int(f.functional_positive.sum()),
    'top10_hits':int(top.functional_positive.sum()),'top10_tested':len(top),'top10_precision':float(top.functional_positive.mean()),
    'average_precision':float(average_precision_score(y[ix],scores)),
    'auroc':float(roc_auc_score(y[ix],scores)) if len(set(y[ix]))==2 else None,
    'response_spearman':float(spearmanr(f.functional_response,scores).statistic),
    'top10_candidates':'|'.join(top.candidate_id)}
preds=[];metrics=[];fit_log=[];coefficients=[]
def run_fold(train,test,split,fold):
 assert not set(d.iloc[train].candidate_id)&set(d.iloc[test].candidate_id)
 for name,cols in MODELS.items():
  m=model();m.fit(d.iloc[train][cols].to_numpy(float),y[train]);scores=m.predict_proba(d.iloc[test][cols].to_numpy(float))[:,1]
  fit_log.append({'split':split,'fold':str(fold),'model':name,'training_ids':d.iloc[train].candidate_id.tolist(),'test_ids':d.iloc[test].candidate_id.tolist(),'features':cols,'n_iter':m.named_steps['logit'].n_iter_.tolist()})
  metrics.append(score_metrics(test,scores,name,split,fold))
  for i,s in zip(test,scores):preds.append({'split':split,'fold':str(fold),'model':name,'candidate_id':d.iloc[i].candidate_id,'system':d.iloc[i].system,'functional_positive':int(y[i]),'functional_response':float(d.iloc[i].functional_response),'score':float(s)})
  if split=='leave_system_out':
   names=m.named_steps['impute'].get_feature_names_out(cols)
   coefficients.extend({'fold':str(fold),'model':name,'feature':f,'scaled_coefficient':float(c)} for f,c in zip(names,m.named_steps['logit'].coef_[0]))

for system in sorted(d.system.unique()):
 test=np.flatnonzero(d.system.to_numpy()==system);train=np.flatnonzero(d.system.to_numpy()!=system)
 run_fold(train,test,'leave_system_out',system)
 for name in BASE:
  scores=d.iloc[test][name].to_numpy(float);metrics.append(score_metrics(test,scores,'raw_'+name,'leave_system_out',system))
  for i,s in zip(test,scores):preds.append({'split':'leave_system_out','fold':system,'model':'raw_'+name,'candidate_id':d.iloc[i].candidate_id,'system':system,'functional_positive':int(y[i]),'functional_response':float(d.iloc[i].functional_response),'score':float(s)})

# Secondary: connected components of peptides within Hamming distance two.
cluster_info=[];skipped=[]
for system in sorted(d.system.unique()):
 ix=np.flatnonzero(d.system.to_numpy()==system);seq=d.iloc[ix].sequence.tolist();parent=list(range(len(ix)))
 def root(i):
  while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
  return i
 for i in range(len(ix)):
  for j in range(i):
   if len(seq[i])==len(seq[j]) and sum(a!=b for a,b in zip(seq[i],seq[j]))<=2:
    a,b=root(i),root(j);parent[a]=b
 groups=np.array([root(i) for i in range(len(ix))]);ng=len(set(groups))
 for i,g in zip(ix,groups):cluster_info.append({'candidate_id':d.iloc[i].candidate_id,'system':system,'cluster':system+'_'+str(g),'cluster_size':int((groups==g).sum())})
 if ng<2:skipped.append({'system':system,'reason':'All candidates in one Hamming<=2 connected component.'});continue
 cv=list(GroupKFold(n_splits=min(5,ng)).split(ix,y[ix],groups))
 if any(len(set(y[ix[tr]]))<2 for tr,te in cv):skipped.append({'system':system,'reason':'At least one training group split has a single label.'});continue
 local_start=len(preds)
 for fold,(tr,te) in enumerate(cv):run_fold(ix[tr],ix[te],'within_system_cluster',system+'_'+str(fold))
 local=pd.DataFrame(preds[local_start:])
 for name,g in local.groupby('model'):
  score_by=dict(zip(g.candidate_id,g.score));scores=np.array([score_by[c] for c in d.iloc[ix].candidate_id])
  metrics.append(score_metrics(ix,scores,name,'within_system_cluster_aggregated',system))

pd.DataFrame(preds).to_csv(O/'all_heldout_predictions.csv',index=False)
pd.DataFrame(metrics).to_csv(O/'metrics_by_system_and_fold.csv',index=False)
pd.DataFrame(coefficients).to_csv(O/'training_only_scaled_coefficients.csv',index=False)
pd.DataFrame(cluster_info).to_csv(O/'sequence_components.csv',index=False)
(O/'skipped_analyses.json').write_text(json.dumps(skipped,indent=2),encoding='utf-8')
(O/'fit_membership.json').write_text(json.dumps(fit_log,indent=2),encoding='utf-8')
md=pd.DataFrame(metrics)
macro=md[md.split.isin(['leave_system_out','within_system_cluster_aggregated'])].groupby(['split','model'])[['top10_hits','top10_precision','average_precision','auroc','response_spearman']].mean().reset_index()
macro.to_csv(O/'macro_metrics.csv',index=False)

# Same dimensionality: permute complete dynamic vectors within each system.
null=[]
for seed in range(200):
 rng=np.random.default_rng(seed);permuted=d[DYNAMIC].to_numpy(float).copy()
 for system in sorted(d.system.unique()):
  ix=np.flatnonzero(d.system.to_numpy()==system);permuted[ix]=permuted[rng.permutation(ix)]
 for system in sorted(d.system.unique()):
  test=np.flatnonzero(d.system.to_numpy()==system);train=np.flatnonzero(d.system.to_numpy()!=system)
  X=np.column_stack([d[BASE+STATIC].to_numpy(float),permuted]);m=model();m.fit(X[train],y[train]);scores=m.predict_proba(X[test])[:,1]
  z=score_metrics(test,scores,'B1_plus_shuffled_dynamic','leave_system_out',system);z['seed']=seed;null.append(z)
pd.DataFrame(null).to_csv(O/'shuffled_dynamic_controls_all_seeds.csv',index=False)
audit={'status':'complete','variant':V,'template_version':TV,'candidate_rows':len(d),'independent_systems':3,
 'primary_split':'Leave entire system out; same preprocessing, classifier and fixed C=1 for all feature sets.',
 'base_model_retraining':False,'downstream_model_trained':True,'same_dimension_permutation_seeds':200,
 'runtime':{'python':sys.version,'numpy':np.__version__,'pandas':pd.__version__,'scipy':scipy.__version__,'sklearn':sklearn.__version__},
 'source_sha256':{str(p):sha(p) for p in [R/'protocol_frozen_before_new_predictions.json',R/'protocol_additional_controls_before_evaluation.json',R/'functional_candidates.csv',R/'baseline_inputs.csv',R/'states.json',Path(__file__)]},
 'no_claim_of_binding_accuracy':'Labels measure T-cell response. Does not directly validate affinity, antigen processing, clinical immunogenicity or causal dynamics.'}
(O/'evaluation_audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
print(macro.to_string(index=False))
