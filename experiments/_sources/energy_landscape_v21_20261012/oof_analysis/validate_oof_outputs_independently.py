"""Read-only independent arithmetic/sample validation of completed OOF analysis."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import sys
import numpy as np
import pandas as pd
from scipy.stats import rankdata

BASE = Path('/root/flexppi_completion_20261010/energy_landscape_v21_20261012/oof_analysis')
SOURCE = Path('/root/flexppi_revision_20261002/corrected_training_v1')
def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def corr(x,y,rank=False):
    if rank:
        x,y = rankdata(x,method='average'),rankdata(y,method='average')
    return float(np.corrcoef(x,y)[0,1])

def category(values,aa):
    labels = np.full(len(values),'unclassified',dtype='<U24')
    for kind in np.unique(aa):
        selected = aa == kind
        if selected.sum() < 4:
            continue
        v = values[selected]
        bounds = np.percentile(v,[25,75],axis=0,method='linear')
        if np.any(bounds[0]==bounds[1]):
            continue
        for i in range(4):
            r,h = ('low','low') if i==0 else ('low','high') if i==1 else ('high','low') if i==2 else ('high','high')
            mask_r = v[:,0] <= bounds[0,0] if r=='low' else v[:,0]>=bounds[1,0]
            mask_h = v[:,1] <= bounds[0,1] if h=='low' else v[:,1]>=bounds[1,1]
            ix = np.flatnonzero(selected)[mask_r&mask_h]
            labels[ix] = r+'_RMSF_'+h+'_Hchi'
    return labels

def main():
    audit = json.loads((BASE/'run_audit.json').read_text())
    original_audit_hash = sha(BASE/'run_audit.json')
    assert audit['status']=='complete' and audit['records']==3063 and audit['paired_common_residues']==824988
    checks = []
    def check(name,value):
        checks.append(dict(name=name,passed=bool(value)))
        assert value,name
    check('All certified output hashes remain exact',all(sha(BASE/name)==digest for name,digest in audit['files'].items()))
    frozen = pd.read_csv(BASE/'frozen_input_records.csv',float_precision='round_trip')
    correlations = pd.read_csv(BASE/'per_record_correlations.csv',float_precision='round_trip')
    classes = pd.read_csv(BASE/'per_record_reference_classes.csv',float_precision='round_trip')
    confusion = pd.read_csv(BASE/'per_record_class_confusion.csv',float_precision='round_trip')
    source_counts = pd.read_csv(BASE/'raw_source_counts.csv').set_index('dataset')
    check('3063 distinct record inputs',frozen.sample_id.nunique()==len(frozen)==3063)
    check('Original source record counts',frozen.groupby('dataset').size().to_dict()=={'atlas':1929,'dynarepo':890,'itrace':244})
    check('Each correlation mode exactly same record IDs',all(set(group.sample_id)==set(frozen.sample_id) and group.sample_id.is_unique and group['count'].sum()==824988 for _,group in correlations.groupby('mode')))
    check('All record class counts partition exact common mask',classes.groupby('sample_id')['count'].sum().to_dict()==frozen.set_index('sample_id').common_valid.to_dict())
    check('All confusion row sums match corresponding reference classes',confusion.groupby(['sample_id','reference_class'])['count'].sum().to_dict()==classes.set_index(['sample_id','reference_class'])['count'].to_dict())
    check('Source common counts exact',source_counts.common_valid.to_dict()=={'atlas':367362,'dynarepo':282949,'itrace':174677})
    summaries = pd.read_csv(BASE/'correlation_record_macro.csv',float_precision='round_trip')
    for r in summaries.itertuples(index=False):
        selected = correlations[correlations['mode'].eq(r.mode)]
        if r.dataset != 'all':
            selected = selected[selected.dataset.eq(r.dataset)]
        for column in ['reference_spearman','prediction_spearman','reference_pearson','prediction_pearson','spearman_absolute_discrepancy']:
            value = getattr(r,column+'_record_macro')
            check('Record-macro coupling '+r.dataset+'/'+r.mode+'/'+column,abs(selected[column].mean()-value)<1e-12)
    group_summary = pd.read_csv(BASE/'reference_class_record_macro.csv',float_precision='round_trip')
    for r in group_summary.itertuples(index=False):
        selected = classes[classes.reference_class.eq(r.reference_class)]
        if r.dataset != 'all':
            selected = selected[selected.dataset.eq(r.dataset)]
        for column in ['reference_fraction','rmsf_mae','hchi_mae','rank_class_recall']:
            a,b = selected[column].mean(),getattr(r,column+'_record_macro')
            check('Record-macro class '+r.dataset+'/'+r.reference_class+'/'+column,(pd.isna(a) and pd.isna(b)) or abs(a-b)<1e-12)
    # Choice independent of result: first, middle and last alphabetic ID in
    # each source. Recalculate from saved source NPZ with different functions.
    chosen = []
    for dataset,group in frozen.groupby('dataset',sort=True):
        identifiers = sorted(group.sample_id)
        chosen += [identifiers[0],identifiers[len(identifiers)//2],identifiers[-1]]
    samples = []
    for sample in chosen:
        with np.load(SOURCE/'dataset/labels'/(sample+'.npz')) as label:
            names = label['metric_names'].astype(str).tolist()
            idx = [names.index('rmsf'),names.index('chi_entropy')]
        with np.load(SOURCE/'evaluation/full'/(sample+'.npz')) as cache:
            valid = cache['valid'][:,idx].all(axis=1)
            y,p = cache['target'][valid][:,idx].astype(float),cache['mean'][valid].astype(float)
            aa = cache['amino_acid'][valid].astype(str)
        y_center,p_center = y.copy(),p.copy()
        for kind in np.unique(aa):
            mask = aa==kind
            y_center[mask] -= np.mean(y[mask],axis=0)
            p_center[mask] -= np.mean(p[mask],axis=0)
        for mode,a,b in [('cross_type_raw',y,p),('within_record_AA_centered',y_center,p_center)]:
            row = correlations[correlations.sample_id.eq(sample)&correlations['mode'].eq(mode)].iloc[0]
            for prefix,values in [('reference',a),('prediction',b)]:
                for name,rank in [('spearman',True),('pearson',False)]:
                    check('NPZ '+sample+'/'+mode+'/'+prefix+'/'+name,abs(corr(values[:,0],values[:,1],rank)-row[prefix+'_'+name])<1e-12)
        reference_class,predicted_class = category(y,aa),category(p,aa)
        for name in sorted(set(classes.reference_class)):
            selected = reference_class==name
            row = classes[classes.sample_id.eq(sample)&classes.reference_class.eq(name)].iloc[0]
            check('NPZ class count '+sample+'/'+name,int(selected.sum())==int(row['count']))
            check('NPZ predicted class count '+sample+'/'+name,int((predicted_class==name).sum())==int(row.predicted_class_count))
            check('NPZ class intersection '+sample+'/'+name,int((selected&(predicted_class==name)).sum())==int(row.intersection_count))
            if selected.any():
                for k,output in enumerate(['rmsf','hchi']):
                    check('NPZ group MAE '+sample+'/'+name+'/'+output,abs(float(np.mean(abs(y[selected,k]-p[selected,k])))-row[output+'_mae'])<1e-12)
        samples.append(dict(sample_id=sample,count=len(y),same_AA_types=len(np.unique(aa))))
    check('No existing run certificate modified',sha(BASE/'run_audit.json')==original_audit_hash)
    result = dict(status='passed_readonly_independent_arithmetic_and_sample_recomputation',utc=datetime.now(timezone.utc).isoformat(),run_audit_sha256=original_audit_hash,verification_code_sha256=sha(__file__),checks=checks,checked_samples=samples,sample_selection='First, middle, last alphabetical ID in each source; no result-dependent choice',no_GPU=True,no_inference=True,no_training=True,no_existing_audit_modified=True)
    out = BASE/'validation'
    out.mkdir(exist_ok=True)
    truth = pd.read_csv('/root/flexppi_completion_20261010/baselines/cohort/truth_residues.csv.gz',keep_default_na=False)
    mask_inventory = truth.groupby(['dataset','amino_acid','valid_rmsf','valid_chi_entropy']).size().rename('residue_count').reset_index()
    mask_inventory.to_csv(out/'raw_mask_inventory_by_AA_and_source.csv',index=False)
    # No zero-position AA group is silently described as an undefined
    # correlation: record absence is disclosed separately from small groups.
    observed = pd.read_csv(BASE/'per_record_same_AA.csv',usecols=['sample_id','amino_acid'])
    observed_set = set(map(tuple,observed.to_numpy()))
    absent = [dict(sample_id=r.sample_id,dataset=r.dataset,amino_acid=kind,count=0,reason='no_common_valid_positions_of_this_AA') for r in frozen.itertuples(index=False) for kind in list('RNDCEQHILKMFPSTWYV') if (r.sample_id,kind) not in observed_set]
    pd.DataFrame(absent).to_csv(out/'absent_record_AA_groups.csv',index=False)
    result['raw_mask_inventory'] = mask_inventory.to_dict(orient='records')
    result['absent_record_AA_groups'] = len(absent)
    result['observed_record_AA_groups'] = len(observed)
    result['validation_output_sha256'] = {p.name:sha(p) for p in out.glob('*.csv')}
    (out/'independent_validation.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status=result['status'],checks=len(checks),samples=samples,run_audit_sha256=original_audit_hash),indent=2))

if __name__=='__main__':
    main()
