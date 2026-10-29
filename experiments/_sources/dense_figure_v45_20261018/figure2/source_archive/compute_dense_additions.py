from pathlib import Path
import os, sys, json, hashlib, time, warnings, importlib.util
os.environ['OMP_NUM_THREADS']='2'
os.environ['OPENBLAS_NUM_THREADS']='2'
BASE=Path('E:/Project/FlexPPI/experiments/_sources')
sys.path.insert(0,str(BASE/'completed_results_20261012/plot_runtime'))
sys.path.append(str(BASE/'label_sensitivity/runtime'))
import numpy as np
import pandas as pd

OUT=Path(__file__).resolve().parent
SRC=BASE/'predictive_figure_v44_20261018/analysis'
CORE=Path('E:/Project/FlexPPI/experiments/core_corrected_pipeline')
spec=importlib.util.spec_from_file_location('v44_metrics',SRC/'compute_metrics.py')
old=importlib.util.module_from_spec(spec); spec.loader.exec_module(old)
METHODS=['Full','Intra','Static','AA']
TARGETS=['rmsf','chi_entropy']
GRID=old.GRID

def main():
    started=time.time(); inputs={}
    def register(path):
        inputs[str(path)]=old.sha(path)
    for fn in ['compute_metrics.py','PRESPECIFIED_PROTOCOL.md','summary.json','input_alignment_audit.json','reference_binned_trend.csv','record_codes.csv','per_record_metrics.csv','within_type_per_record_metrics.csv','within_type_record_macro_metrics.csv','macro_pr_curve.csv','all_selected_endpoint.csv','prediction_reference_rmsf.npz','prediction_reference_chi_entropy.npz']:
        register(SRC/fn)
    assert json.loads((SRC/'input_alignment_audit.json').read_text(encoding='utf8'))['status']=='passed'
    expected=json.loads((SRC/'summary.json').read_text(encoding='utf8'))
    protocol=OUT/'PRESPECIFIED_ADDITIONS.md'; assert protocol.exists()
    schema_path=CORE/'static_corrected/training_schema.json';register(schema_path)
    schema=json.loads(schema_path.read_text(encoding='utf8'))['records'];assert len(schema)==3063
    truth_path=BASE/'baselines/cohort/truth_residues.csv.gz';register(truth_path)
    truth=pd.read_csv(truth_path,dtype={'sample_id':str,'dataset':str,'chain':str,'author_residue_id':str,'amino_acid':str})
    assert len(truth)==961321
    static={}
    for method,fn in [('Static','static_sequence_complex_oof.npz'),('AA','aa_only_oof.npz')]:
        p=CORE/'static_corrected'/fn;register(p)
        static[method]=np.load(p)['prediction'];assert static[method].shape==(len(truth),2)
    trend_old=pd.read_csv(SRC/'reference_binned_trend.csv',float_precision='round_trip')
    old_typ=pd.read_csv(SRC/'within_type_per_record_metrics.csv',float_precision='round_trip')
    old_typ=old_typ.set_index(['record_code','amino_acid','method'])
    records=pd.read_csv(SRC/'record_codes.csv');assert len(records)==3063
    sources=records.dataset.to_numpy()
    draws=old.stratified_draws(sources,10000)
    weights=np.zeros((1000,len(schema)),dtype=np.float64)
    for b,ix in enumerate(draws[:1000]):weights[b]=np.bincount(ix,minlength=len(schema))
    del draws
    binspec={};aggregates={}
    for target in TARGETS:
        g=trend_old[trend_old.target==target].sort_values('bin_index')
        edges=np.r_[g.bin_left.to_numpy(),g.bin_right.to_numpy()[-1]]
        assert len(edges)==21 and np.array_equal(edges[1:-1],g.bin_right.to_numpy()[:-1])
        binspec[target]=edges
        aggregates[target]={'den':np.zeros((3063,20)),'xnum':np.zeros((3063,20)),'num':{m:np.zeros((3063,20)) for m in METHODS},'position_counts':np.zeros(20,dtype=np.int64),'n_positions':0}
    raw={m:{k:[] for k in ['type_code','threshold','precision','recall','fpr','tpr']} for m in METHODS}
    record_grids={m:np.full((3063,101),np.nan) for m in METHODS}
    type_map=[];metrics=[];type_code=0
    for code,r in enumerate(schema):
        sid=r['sample_id'];start,stop=r['start'],r['stop'];g=truth.iloc[start:stop];n=stop-start
        assert np.all(g.sample_id==sid) and np.all(g.dataset==r['dataset']) and np.all(g.fold==r['fold'])
        caches={}
        for variant in ['full','intrachain']:
            p=CORE/'neural_exports/evaluation'/variant/(sid+'.npz');register(p)
            with np.load(p) as z:caches[variant]={k:z[k].copy() for k in z.files}
            c=caches[variant]
            assert c['mean'].shape==(n,2)
            assert np.array_equal(c['target'][:,:2],g[['target_rmsf','target_chi_entropy']].to_numpy(dtype=np.float32),equal_nan=True)
            assert np.array_equal(c['valid'][:,:2],g[['valid_rmsf','valid_chi_entropy']].to_numpy(dtype=bool))
            assert np.array_equal(c['amino_acid'].astype(str),g.amino_acid.astype(str))
            assert np.array_equal(c['chain'].astype(str),g.chain.astype(str))
            assert np.array_equal(c['author_residue_id'].astype(str),g.author_residue_id.astype(str))
        full,intra=caches['full'],caches['intrachain']
        predictions={'Full':full['mean'],'Intra':intra['mean'],'Static':static['Static'][start:stop],'AA':static['AA'][start:stop]}
        for j,target in enumerate(TARGETS):
            valid=full['valid'][:,j]&intra['valid'][:,j]&np.isfinite(full['target'][:,j])
            n_before=int(valid.sum())
            for pred in predictions.values():valid&=np.isfinite(pred[:,j])
            assert int(valid.sum())==n_before
            y=full['target'][valid,j].astype(np.float64); aa=full['amino_acid'][valid].astype('U1'); nr=len(y)
            bins=np.searchsorted(binspec[target],y,side='right')-1;assert bins.min()>=0 and bins.max()<20
            ag=aggregates[target]
            ag['den'][code]=np.bincount(bins,weights=np.full(nr,1./nr),minlength=20)
            ag['xnum'][code]=np.bincount(bins,weights=y/nr,minlength=20)
            ag['position_counts']+=np.bincount(bins,minlength=20)
            ag['n_positions']+=nr
            for method in METHODS:ag['num'][method][code]=np.bincount(bins,weights=predictions[method][valid,j].astype(np.float64)/nr,minlength=20)
            if target!='chi_entropy':continue
            grids={m:[] for m in METHODS}
            for amino in 'CDEFHIKLMNPQRSTVWY':
                keep=aa==amino; yy=y[keep];ny=len(yy)
                if ny<4:continue
                threshold=float(np.quantile(yy,.9,method='linear')); binary=yy>=threshold;p=int(binary.sum())
                if p==0 or p==ny:continue
                type_map.append({'type_code':type_code,'record_code':code,'sample_id':sid,'dataset':r['dataset'],'fold':r['fold'],'amino_acid':amino,'n_positions':ny,'threshold_q90':threshold,'n_positive':p,'prevalence':p/ny})
                for method in METHODS:
                    score=predictions[method][valid,j][keep].astype(np.float64)
                    ap,auc,grid,curve=old.curves(binary,score)
                    prev=old_typ.loc[(code,amino,method)]
                    assert bool(prev.eligible) and int(prev.n_positions)==ny and int(prev.n_positive)==p
                    assert abs(ap-float(prev.ap))<2e-14 and abs(auc-float(prev.auroc))<2e-14
                    grids[method].append(grid)
                    metrics.append({'type_code':type_code,'record_code':code,'sample_id':sid,'dataset':r['dataset'],'fold':r['fold'],'amino_acid':amino,'method':method,'ap':ap,'auroc':auc,'prevalence':p/ny})
                    for key,values in curve.items():raw[method][key].append(values)
                    raw[method]['type_code'].append(np.full(len(curve['recall']),type_code,dtype=np.int32))
                type_code+=1
            assert len(grids['Full'])>0
            for method in METHODS:record_grids[method][code]=np.mean(grids[method],axis=0)
        if (code+1)%500==0:print('Verified and computed records',code+1,'/3063; seconds',round(time.time()-started,1),flush=True)
    assert type_code==42108
    trend_rows=[]
    for target in TARGETS:
        ag=aggregates[target];den=ag['den'];xnum=ag['xnum'];g=trend_old[trend_old.target==target].sort_values('bin_index')
        assert ag['n_positions']=={'rmsf':961321,'chi_entropy':824988}[target]
        assert np.array_equal(ag['position_counts'],g.n_positions.to_numpy())
        assert np.array_equal((den>0).sum(axis=0),g.n_records.to_numpy())
        db=weights@den
        with warnings.catch_warnings():
            warnings.simplefilter('ignore',RuntimeWarning)
            xmean=xnum.sum(axis=0)/den.sum(axis=0)
            for method in METHODS:
                num=ag['num'][method];mean=num.sum(axis=0)/den.sum(axis=0)
                if method=='Full':
                    ci=np.nanquantile((weights@num)/db,[.025,.975],axis=0)
                    for col,values in [('x_mean',xmean),('prediction_mean',mean),('ci_low',ci[0]),('ci_high',ci[1])]:
                        assert np.allclose(values,g[col].to_numpy(),rtol=0,atol=2e-12,equal_nan=True),(target,col,np.nanmax(abs(values-g[col].to_numpy())))
                else:ci=np.full((2,20),np.nan)
                for b in range(20):
                    trend_rows.append({'target':target,'method':method,'bin_index':b,'bin_left':float(binspec[target][b]),'bin_right':float(binspec[target][b+1]),'x_mean':float(xmean[b]),'prediction_mean':float(mean[b]),'ci_low':float(ci[0,b]),'ci_high':float(ci[1,b]),'n_positions':int(ag['position_counts'][b]),'n_records':int((den[:,b]>0).sum()),'record_equal_mass':float(den[:,b].sum()/3063)})
    pd.DataFrame(trend_rows).to_csv(OUT/'four_method_binned_trend.csv',index=False)
    print('Trend means and Full bootstrap intervals reproduce v44',flush=True)
    pd.DataFrame(type_map).to_csv(OUT/'within_type_codes.csv',index=False)
    metrics_df=pd.DataFrame(metrics);metrics_df.to_csv(OUT/'within_type_complete_curve_metrics.csv',index=False)
    macro_metrics=metrics_df.groupby(['record_code','sample_id','dataset','fold','method'],as_index=False)[['ap','auroc','prevalence']].mean()
    prev_macro=pd.read_csv(SRC/'within_type_record_macro_metrics.csv',float_precision='round_trip')
    merged=macro_metrics.merge(prev_macro,on=['record_code','sample_id','dataset','fold','method'],validate='one_to_one',suffixes=('_new','_old'))
    assert len(merged)==12252
    for metric in ['ap','auroc','prevalence']:assert np.allclose(merged[metric+'_new'],merged[metric+'_old'],rtol=0,atol=2e-14)
    macro_metrics.to_csv(OUT/'within_type_record_macro_metrics.csv',index=False)
    macro_rows=[];endpoints=[];within_summary={}
    for method in METHODS:
        grids=record_grids[method]; assert np.isfinite(grids).all()
        mean=grids.mean(axis=0); mm=macro_metrics[macro_metrics.method==method]
        prev=expected['within_type'][method]
        assert abs(float(mm.ap.mean())-prev['mean_ap'])<2e-14 and abs(float(mm.auroc.mean())-prev['mean_auroc'])<2e-14
        prevalence=float(mm.prevalence.mean());assert abs(prevalence-prev['mean_record_type_prevalence'])<2e-14
        for i,re in enumerate(GRID):macro_rows.append({'target':'within_type_hchi','method':method,'recall':float(re),'precision':float(mean[i]),'n_records':3063,'n_eligible_record_types':42108,'no_skill_precision':prevalence})
        endpoints.append({'target':'within_type_hchi','method':method,'recall':1.,'precision':prevalence,'n_records':3063,'n_eligible_record_types':42108,'endpoint':'all_selected_original_last_threshold'})
        within_summary[method]={'mean_ap':float(mm.ap.mean()),'mean_auroc':float(mm.auroc.mean()),'mean_record_type_prevalence':prevalence,'n_records':3063,'n_eligible_record_types':42108,'macro_curve_recall1_first_threshold_precision':float(mean[-1]),'all_selected_precision':prevalence}
        np.savez_compressed(OUT/f'raw_within_type_curves_{method}.npz',**{k:np.concatenate(v) for k,v in raw[method].items()})
        np.savez_compressed(OUT/f'within_type_record_precision_grid_{method}.npz',record_code=np.arange(3063,dtype=np.int32),recall=GRID,precision=grids)
    pd.DataFrame(macro_rows).to_csv(OUT/'within_type_macro_pr.csv',index=False)
    pd.DataFrame(endpoints).to_csv(OUT/'within_type_all_selected_endpoint.csv',index=False)
    print('All 42,108 type curves and 3,063 record AP/AUROC reproduced',flush=True)
    per_record=pd.read_csv(SRC/'per_record_metrics.csv',float_precision='round_trip')
    sources_rows=[]
    tasks={target:per_record[per_record.target==target] for target in TARGETS};tasks['within_type_hchi']=macro_metrics
    for task,df in tasks.items():
        table=df.groupby(['dataset','method'],as_index=False).agg(n_records=('record_code','nunique'),mean_ap=('ap','mean'))
        for dataset in sorted(table.dataset.unique()):
            g=table[table.dataset==dataset].set_index('method');full_ap=float(g.loc['Full','mean_ap'])
            for method in METHODS:
                sources_rows.append({'task':task,'dataset':dataset,'method':method,'n_records':int(g.loc[method,'n_records']),'mean_ap':float(g.loc[method,'mean_ap']),'full_minus_method_ap':full_ap-float(g.loc[method,'mean_ap'])})
        for method in METHODS:
            g=table[table.method==method];overall=float(np.average(g.mean_ap,weights=g.n_records))
            previous=expected['within_type'][method]['mean_ap'] if task=='within_type_hchi' else expected['metrics'][task][method]['mean_ap']
            assert int(g.n_records.sum())==3063 and abs(overall-previous)<2e-14
    pd.DataFrame(sources_rows).to_csv(OUT/'source_ap_matrix.csv',index=False)
    # Copy unchanged curves and endpoint as explicit, hash-bound sources for figure consumers.
    import shutil
    for fn in ['macro_pr_curve.csv','all_selected_endpoint.csv','summary.json']:
        shutil.copyfile(SRC/fn,OUT/('v44_'+fn))
    outputs={str(p.name):old.sha(p) for p in sorted(OUT.iterdir()) if p.is_file() and p.name not in ['addition_audit.json','inputs_sha256.json','compute_log.txt']}
    (OUT/'inputs_sha256.json').write_text(json.dumps(inputs,indent=2,ensure_ascii=False),encoding='utf8')
    audit={'status':'passed','code_sha256':old.sha(__file__),'protocol_sha256':old.sha(protocol),'input_manifest_sha256':old.sha(OUT/'inputs_sha256.json'),'n_records':3063,'n_positions':{'rmsf':961321,'chi_entropy':824988},'n_eligible_record_types':42108,'all_four_methods_original_positions_masks_identity_verified':True,'no_downsampling_no_clipping_no_calibration':True,'full_20bin_trend_and_ci_reproduced':True,'full_ci_max_tolerance':2e-12,'same_v44_stratified_bootstrap_first_1000_of_10000':True,'all_type_AP_AUROC_and_record_macro_values_reproduced':True,'source_36cells_weighted_means_reproduce_v44':True,'within_type_summary':within_summary,'outputs_sha256':outputs,'runtime_seconds':time.time()-started}
    (OUT/'addition_audit.json').write_text(json.dumps(audit,indent=2,ensure_ascii=False,allow_nan=False),encoding='utf8')
    print('PASSED',round(time.time()-started,1),'seconds',flush=True)

if __name__=='__main__':main()
