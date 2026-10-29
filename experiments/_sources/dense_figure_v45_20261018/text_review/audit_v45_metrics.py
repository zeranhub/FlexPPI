"""Independent reconstruction of added dense-Figure-2 statistics (no training)."""
from pathlib import Path
import os, sys, json, hashlib, warnings
os.environ['OMP_NUM_THREADS']='2'
os.environ['OPENBLAS_NUM_THREADS']='2'
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT.parent
OLD=BASE/'predictive_figure_v44_20261018/analysis'
NEW=ROOT/'analysis'
CORE=Path('E:/Project/FlexPPI/experiments/core_corrected_pipeline')
sys.path.insert(0,str(BASE/'completed_results_20261012/plot_runtime'))
import numpy as np
import pandas as pd

METHODS=['Full','Intra','Static','AA']
TASKS=['rmsf','chi_entropy']
SOURCES=['atlas','dynarepo','itrace']
GRID=np.linspace(0,1,101)
checks={}
details={}
inputs={}
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def bind(path):inputs[str(path)]=sha(path)
def compare(name,a,b,tol=1e-11):
    a=np.asarray(a,dtype=float);b=np.asarray(b,dtype=float)
    ok=a.shape==b.shape and np.array_equal(np.isnan(a),np.isnan(b)) and np.allclose(a,b,atol=tol,rtol=tol,equal_nan=True)
    checks[name]=bool(ok)
    details[name]={'shape_a':list(a.shape),'shape_b':list(b.shape),'max_abs_error':float(np.nanmax(np.abs(a-b))) if a.shape==b.shape and np.isfinite(a-b).any() else None}

def source_audit():
    # Source stratification is a new aggregation of already audited per-record AP.
    p=OLD/'per_record_metrics.csv';typ=OLD/'within_type_record_macro_metrics.csv'
    bind(p);bind(typ)
    main=pd.read_csv(p);within=pd.read_csv(typ)
    records=pd.read_csv(OLD/'record_codes.csv')
    checks['all_3063_records_unchanged']=len(records)==3063 and records.record_code.nunique()==3063
    checks['source_record_counts_1929_890_244']=records.groupby('dataset').size().to_dict()=={'atlas':1929,'dynarepo':890,'itrace':244}
    rows=[]
    for task in ['rmsf','chi_entropy','within_type_chi_entropy']:
        frame=main[main.target==task] if task!='within_type_chi_entropy' else within
        for source in SOURCES:
            for method in METHODS:
                sub=frame[(frame.dataset==source)&(frame.method==method)]
                rows.append({'task':task,'source':source,'method':method,'ap':float(sub.ap.mean()),'n_records':int(sub.record_code.nunique())})
    expected=pd.DataFrame(rows)
    actual=pd.read_csv(NEW/'source_ap_matrix.csv')
    # Field normalization is limited to explicit semantic synonyms.
    actual=actual.rename(columns={'target':'task','dataset':'source','mean_ap':'ap'})
    aliases={'within_type':'within_type_chi_entropy','chi_entropy_within_type':'within_type_chi_entropy','within_type_hchi':'within_type_chi_entropy','hchi_within_type':'within_type_chi_entropy'}
    actual['task']=actual.task.replace(aliases)
    index=['task','source','method']
    exp=expected.set_index(index).sort_index();act=actual.set_index(index).sort_index()
    checks['all_36_source_AP_cells_identity']=len(actual)==36 and list(exp.index)==list(act.index)
    if checks['all_36_source_AP_cells_identity']:
        compare('all_36_source_AP_values_exact',exp.ap.to_numpy(),act.ap.to_numpy())
        if 'n_records' in act:compare('all_source_AP_record_counts',exp.n_records.to_numpy(),act.n_records.to_numpy(),tol=0)
    gain=[]
    for target in TASKS:
        for source in SOURCES:
            sub=expected[(expected.task==target)&(expected.source==source)].set_index('method')
            gain.append(all(sub.loc['Full','ap']>sub.loc[m,'ap'] for m in METHODS[1:]))
    checks['all_6_source_target_Full_AP_higher_than_controls']=all(gain)
    expected.to_csv(ROOT/'text_review/independent_source_AP.csv',index=False)

def trend_audit():
    records=pd.read_csv(OLD/'record_codes.csv').sort_values('record_code')
    schema_path=CORE/'static_corrected/training_schema.json';bind(schema_path)
    schema=json.loads(schema_path.read_text(encoding='utf8'))['records']
    sources=records.dataset.to_numpy()
    rng=np.random.default_rng(20261018)
    # Reproduce first 1,000 of the original 10,000 paired record draws.
    # RNG consumption must still include all 10,000 draws in each stratum.
    groups=[np.flatnonzero(sources==s) for s in sorted(set(sources))]
    draws=np.concatenate([rng.choice(ix,size=(10000,len(ix)),replace=True)[:1000] for ix in groups],axis=1)
    bw=np.zeros((1000,len(records)),dtype=float)
    for i,row in enumerate(draws):bw[i]=np.bincount(row,minlength=len(records))
    actual=pd.read_csv(NEW/'four_method_binned_trend.csv').set_index(['target','method','bin_index']).sort_index()
    oldtrend=pd.read_csv(OLD/'reference_binned_trend.csv')
    static={}
    for method,name in [('Static','static_sequence_complex_oof.npz'),('AA','aa_only_oof.npz')]:
        path=CORE/'static_corrected'/name;bind(path);static[method]=np.load(path)['prediction']
    all_expected=[]
    for ti,target in enumerate(TASKS):
        path=OLD/f'prediction_reference_{target}.npz';bind(path)
        data=np.load(path)
        y=data['reference'].astype(float);code=data['record_code'];ri=data['residue_index']
        counts=np.bincount(code,minlength=len(records))
        offsets=np.array([s['start'] for s in schema])
        positions=offsets[code]+ri
        preds={'Full':data['prediction'].astype(float),'Static':static['Static'][positions,ti].astype(float),'AA':static['AA'][positions,ti].astype(float)}
        intra=np.empty_like(y)
        boundaries=np.r_[0,np.flatnonzero(np.diff(code))+1,len(code)]
        for left,right in zip(boundaries[:-1],boundaries[1:]):
            rid=int(code[left]);sid=records.iloc[rid].sample_id
            p=CORE/'neural_exports/evaluation/intrachain'/(sid+'.npz')
            # Bind these original source arrays in addition to derived outputs.
            bind(p)
            with np.load(p) as z:intra[left:right]=z['mean'][ri[left:right],ti]
        preds['Intra']=intra
        trend=oldtrend[oldtrend.target==target].sort_values('bin_index')
        edges=np.r_[trend.bin_left.iloc[0],trend.bin_right.to_numpy()]
        binid=np.searchsorted(edges,y,side='right')-1
        checks[target+'_all_points_inside_original_reference_bins']=bool(np.all((binid>=0)&(binid<20)))
        den=np.zeros((len(records),20));xn=np.zeros_like(den)
        np.add.at(den,(code,binid),1/counts[code]);np.add.at(xn,(code,binid),y/counts[code])
        with np.errstate(invalid='ignore',divide='ignore'):xm=xn.sum(0)/den.sum(0)
        db=bw@den
        for method,pred in preds.items():
            num=np.zeros_like(den);np.add.at(num,(code,binid),pred/counts[code])
            with warnings.catch_warnings():
                warnings.simplefilter('ignore',RuntimeWarning)
                means=num.sum(0)/den.sum(0);boot=(bw@num)/db;ci=np.nanquantile(boot,[.025,.975],axis=0)
            expected=pd.DataFrame({'bin_index':range(20),'x_mean':xm,'prediction_mean':means,'ci_low':ci[0],'ci_high':ci[1], 'n_positions':np.bincount(binid,minlength=20),'n_records':(den>0).sum(0),'record_equal_mass':den.sum(0)/len(records)})
            observed=actual.loc[(target,method)].sort_index()
            for col in expected.columns[1:]:
                if col in ['ci_low','ci_high'] and method!='Full':
                    checks[target+'_'+method+'_'+col+'_not_claimed_or_plotted']=bool(observed[col].isna().all())
                elif col in observed:compare(target+'_'+method+'_'+col,expected[col],observed[col])
            if method=='Full':
                for col in ['x_mean','prediction_mean','ci_low','ci_high','n_positions','n_records','record_equal_mass']:
                    compare(target+'_Full_preserves_v44_'+col,expected[col],trend[col])
            expected.insert(0,'method',method);expected.insert(0,'target',target);all_expected.append(expected)
    pd.concat(all_expected).to_csv(ROOT/'text_review/independent_four_method_trend.csv',index=False)

def within_type_audit():
    mapping=pd.read_csv(NEW/'within_type_codes.csv').sort_values('type_code')
    checks['exact_42108_unique_type_groups']=len(mapping)==42108 and np.array_equal(mapping.type_code,np.arange(42108))
    checks['all_3063_records_in_type_PR']=mapping.record_code.nunique()==3063
    groups_per_record=np.bincount(mapping.record_code,minlength=3063)
    plotted=pd.read_csv(NEW/'within_type_macro_pr.csv')
    endpoints=pd.read_csv(NEW/'within_type_all_selected_endpoint.csv')
    oldtype=pd.read_csv(OLD/'within_type_per_record_metrics.csv')
    expected_AP={(int(r.record_code),r.amino_acid,r.method):float(r.ap) for r in oldtype[oldtype.eligible].itertuples(index=False)}
    mapping_rows=list(mapping.itertuples(index=False))
    records=pd.read_csv(OLD/'record_codes.csv').sort_values('record_code')
    selected_records=[]
    for _,g in records.groupby(['dataset','fold']):
        vals=g.record_code.to_numpy();selected_records.extend([int(vals[0]),int(vals[len(vals)//2]),int(vals[-1])])
    selected_types=[]
    for rc in sorted(set(selected_records)):
        ids=mapping.loc[mapping.record_code==rc,'type_code'].to_numpy()
        selected_types.extend([int(ids[0]),int(ids[-1])])
    selected_types=set(selected_types)
    static={}
    for method,name in [('Static','static_sequence_complex_oof.npz'),('AA','aa_only_oof.npz')]:
        static[method]=np.load(CORE/'static_corrected'/name)['prediction']
    schema=json.loads((CORE/'static_corrected/training_schema.json').read_text(encoding='utf8'))['records']
    point=np.load(OLD/'prediction_reference_chi_entropy.npz')
    code=point['record_code'];amino=point['amino_acid'];ref=point['reference'];ri=point['residue_index']
    starts=np.r_[0,np.flatnonzero(np.diff(code))+1,len(code)]
    raw_total=0;spot_total=0
    for method in METHODS:
        path=NEW/f'raw_within_type_curves_{method}.npz'
        with np.load(path) as z:raw={k:z[k] for k in z.files}
        begin=np.r_[0,np.flatnonzero(np.diff(raw['type_code']))+1];end=np.r_[begin[1:],len(raw['type_code'])]
        checks[method+'_all_42108_raw_groups_retained']=len(begin)==42108 and np.array_equal(raw['type_code'][begin],np.arange(42108))
        sum_grid=np.zeros((3063,101));end_prev=np.zeros(3063);new_ap=[];old_ap=[];source_errors=[]
        for tc,(lo,hi) in enumerate(zip(begin,end)):
            row=mapping_rows[tc];rc=int(row.record_code)
            rr=raw['recall'][lo:hi];pp=raw['precision'][lo:hi]
            ix=np.searchsorted(rr,GRID,side='left')
            sampled=pp[ix].copy();sampled[0]=1
            sum_grid[rc]+=sampled
            end_prev[rc]+=pp[-1]
            a=float(np.sum(np.diff(rr)*pp[1:]))
            new_ap.append(a);old_ap.append(expected_AP[(rc,row.amino_acid,method)])
            if tc in selected_types:
                left,right=starts[rc],starts[rc+1]
                keep=amino[left:right]==row.amino_acid
                yy=ref[left:right][keep].astype(float)
                resid=ri[left:right][keep]
                if method in ['Full','Intra']:
                    variant='full' if method=='Full' else 'intrachain'
                    p=CORE/'neural_exports/evaluation'/variant/(str(row.sample_id)+'.npz');bind(p)
                    with np.load(p) as z:ss=z['mean'][resid,1].astype(float)
                else:ss=static[method][schema[rc]['start']+resid,1].astype(float)
                labels=yy>=np.quantile(yy,.9,method='linear')
                th=np.r_[np.inf,np.unique(ss)[::-1]]
                take=ss[:,None]>=th[None,:]
                tp=(take&labels[:,None]).sum(0)
                fp=(take&~labels[:,None]).sum(0)
                expected_rr=tp/labels.sum()
                expected_pp=np.divide(tp,tp+fp,out=np.ones(len(th),dtype=float),where=(tp+fp)>0)
                ecdf=(~labels[:,None]&take).sum(0)/(~labels).sum()
                source_errors.append(np.allclose(raw['threshold'][lo:hi],th,rtol=0,atol=0,equal_nan=True) and np.allclose(rr,expected_rr,atol=1e-14,rtol=0) and np.allclose(pp,expected_pp,atol=1e-14,rtol=0) and np.allclose(raw['fpr'][lo:hi],ecdf,atol=1e-14,rtol=0))
                spot_total+=1
        mean_grid=(sum_grid/groups_per_record[:,None]).mean(0)
        prev=float((end_prev/groups_per_record).mean())
        obs=plotted[plotted.method==method].sort_values('recall')
        compare(method+'_all_101_macro_PR_precision',mean_grid,obs.precision.to_numpy())
        compare(method+'_all_42108_raw_APs',new_ap,old_ap)
        compare(method+'_real_full_selection_endpoint',[prev],endpoints.loc[endpoints.method==method,'precision'])
        checks[method+'_original_prediction_type_curve_spotchecks_all_passed']=all(source_errors)
        checks[method+'_own_positive_prevalence_baseline']=np.allclose(obs.no_skill_precision,prev,atol=1e-14,rtol=0)
        if method=='AA':
            checks['AA_constant_score_AP_and_full_selection_precision_agree']=np.allclose(mean_grid[1:],prev,atol=1e-14,rtol=0)
        raw_total+=len(raw['recall'])
        print('within_type_method_completed',method,'groups',len(begin),flush=True)
    details['within_type_reconstruction']={'raw_threshold_points_all_methods':raw_total,'source_type_method_spotchecks':spot_total,'selected_records':len(set(selected_records)),'selection_rule':'first/middle/last record in each source/fold; first/last eligible amino-acid type; all four methods','aggregation':'each of 42,108 type grids equally within record, then 3,063 records equally'}

def main():
    source_audit();print('source_AP_audit_finished',flush=True)
    trend_audit();print('all_trend_audit_finished',flush=True)
    within_type_audit();print('all_within_type_PR_audit_finished',flush=True)
    report={'status':'passed' if all(checks.values()) else 'failed','checks':checks,'details':details,'inputs_sha256':inputs,'outputs_sha256':{str(p):sha(p) for p in NEW.iterdir() if p.is_file()},'audit_code_sha256':sha(__file__)}
    out=ROOT/'text_review/independent_metrics_audit.json';out.write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf8');print(json.dumps({'status':report['status'],'failed':[k for k,v in checks.items() if not v]},ensure_ascii=False))
if __name__=='__main__':main()
