"""Re-score frozen baseline means on the identical alternate-label intersection."""
from pathlib import Path
import sys,json
root=Path(__file__).resolve().parent
sys.path.insert(0,str(root/'runtime'))
sys.path.insert(0,str(root.parent/'single_target'))
from completion_common import sha256,statistics
import numpy as np,pandas as pd
base=Path(r'E:\Project\FlexPPI\experiments\core_corrected_pipeline')
plan=json.loads((root/'sensitivity_plan.json').read_text(encoding='utf8'))
schema_path=base/'static_corrected/training_schema.json'
schema={r['sample_id']:r for r in json.loads(schema_path.read_text(encoding='utf8'))['records']}
cache_paths={name:base/'static_corrected'/file for name,file in [('static','static_sequence_complex_oof.npz'),('aa','aa_only_oof.npz')]}
cache_paths['features']=base/'static_corrected/training_features.npz'
cache_paths['backflip2']=base/'backflip2_corrected_scoring/training_aligned.npz'
arrays={}
for name,path in cache_paths.items():
    with np.load(path,allow_pickle=False) as data:arrays[name]={k:data[k] for k in data if name!='features' or k in ['target','valid']}
rows=[];guards=[]
for record in plan['records']:
    sid=record['sample_id'];s=schema[sid];sl=slice(s['start'],s['stop'])
    assert s['input_sha256']==record['structure_sha256'] and s['label_sha256']==record['old_label_sha256']
    with np.load(record['old_label'],allow_pickle=False) as data:old=data['values'][:,:2].copy()
    pred={}
    for name,directory in [('full','full'),('intra','intrachain')]:
        p=base/'neural_exports/evaluation'/directory/(sid+'.npz')
        with np.load(p,allow_pickle=False) as data:
            pred[name]=data['mean'][:,:2].copy()
            assert np.array_equal(data['target'][:,:2],old,equal_nan=True)
            if name=='full':valid=data['valid'][:,:2].copy()
            else:assert np.array_equal(data['valid'][:,:2],valid)
        guards.append(dict(sample_id=sid,method=name,path=str(p),sha256=sha256(p)))
    assert arrays['features']['target'][sl].shape==old.shape
    assert np.array_equal(arrays['features']['target'][sl],old,equal_nan=True)
    assert np.array_equal(arrays['features']['valid'][sl],valid)
    for name in ['static','aa']:pred[name]=arrays[name]['prediction'][sl]
    bf=np.full(old.shape,np.nan);bf[:,0]=arrays['backflip2']['prediction'][sl];pred['backflip2']=bf
    assert np.array_equal(arrays['backflip2']['target'][sl],old[:,0],equal_nan=True)
    for variant in ['uniform_rep1','uniform_available_replicas']:
        path=root/'labels'/variant/(sid+'.npz')
        with np.load(path,allow_pickle=False) as data:new=data['values'][:,:2].copy();newmask=data['mask'][:,:2].copy()
        for m,metric in enumerate(['rmsf','chi_entropy']):
            methods=['full','intra','static','aa']+(['backflip2'] if m==0 else [])
            common=valid[:,m]&newmask[:,m]&np.isfinite(old[:,m])&np.isfinite(new[:,m])
            for name in methods:common&=np.isfinite(pred[name][:,m])
            assert common.sum()==valid[:,m].sum()
            for name in methods:
                for label,target in [('original',old),('common100ns',new)]:
                    stats=statistics(target[common,m],pred[name][common,m])
                    rows.append(dict(sample_id=sid,dataset=record['dataset'],variant=variant,metric=metric,method=name,label=label,**stats))
frame=pd.DataFrame(rows)
frame.to_csv(root/'controls_per_record_metrics.csv',index=False)
numeric=['mae','rmse','bias','spearman','pearson']
summary=frame.groupby(['variant','dataset','metric','method','label'])[numeric].mean().reset_index()
summary.to_csv(root/'controls_record_macro.csv',index=False)
paired=[]
for keys,group in frame.groupby(['variant','dataset','metric','sample_id','label']):
    full=group[group.method=='full'].iloc[0]
    for _,row in group[group.method!='full'].iterrows():
        paired.append(dict(zip(['variant','dataset','metric','sample_id','label'],keys),control=row.method,
            full_minus_control_mae=float(full.mae-row.mae),full_minus_control_spearman=float(full.spearman-row.spearman)))
paired=pd.DataFrame(paired)
paired.to_csv(root/'controls_paired_differences.csv',index=False)
paired.groupby(['variant','dataset','metric','label','control'])[['full_minus_control_mae','full_minus_control_spearman']].mean().to_csv(root/'controls_paired_macro.csv')
audit=dict(status='complete',selected_records=90,eligible_records=83,compared_records=83,
    input_caches={name:dict(path=str(path),sha256=sha256(path)) for name,path in cache_paths.items()},
    input_schema_sha256=sha256(schema_path),neural_guards=guards,
    plan_sha256=sha256(root/'sensitivity_plan.json'),code_sha256=sha256(__file__),
    units='RMSF Angstrom; Hchi nats',scope='Frozen means; no re-training. Source record macro, identical original/new valid residues, no outcome-based exclusions.')
(root/'controls_analysis_audit.json').write_text(json.dumps(audit,indent=2),encoding='utf8')
print(summary[summary.variant=='uniform_rep1'][['dataset','metric','method','label','mae','spearman']].to_string(index=False))
