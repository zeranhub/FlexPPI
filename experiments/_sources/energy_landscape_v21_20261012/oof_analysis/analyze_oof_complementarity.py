"""Frozen Full OOF RMSF/Hchi complementarity; CPU only, no fitted mapping.

All model predictions are previously held-out predictions. Quartiles merely
define reference-conditioned descriptive subsets, not a learned calibration.
"""
import os
os.environ.setdefault('OMP_NUM_THREADS', '2')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '2')
os.environ.setdefault('MKL_NUM_THREADS', '2')
from pathlib import Path
from datetime import datetime, timezone
import argparse
import hashlib
import json
import time
import sys
import numpy as np
import pandas as pd
import yaml
from scipy.stats import spearmanr, pearsonr

VERSION = 'full_corrected3063_paired_complementarity_v1'
CLASSES = ['low_RMSF_low_Hchi','low_RMSF_high_Hchi','high_RMSF_low_Hchi','high_RMSF_high_Hchi','unclassified']
AA18 = list('RNDCEQHILKMFPSTWYV')
EXPECTED_SOURCE_FREEZE = '5b994bb59684e3a33c97f4a1801ee9a606e5dd856db9cfe23bf18515697707d3'
EXPECTED_PRED_MANIFEST = 'aa3bbeaa99cbaa35b278ed0b1abf97649d19f79d480c20b71dceb2cc9922fb97'

def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1048576), b''):
            h.update(block)
    return h.hexdigest()

def jsonwrite(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False)+'\n', encoding='utf-8')

def lines(path):
    return [json.loads(x) for x in Path(path).read_text(encoding='utf-8').splitlines() if x.strip()]

def require(condition, message):
    if not condition:
        raise ValueError(message)

def correlation(x, y):
    result = dict(spearman=np.nan, pearson=np.nan, correlation_reason='')
    if len(x) < 4:
        result['correlation_reason'] = 'fewer_than_4_common_positions'
    elif not np.isfinite(x).all() or not np.isfinite(y).all():
        raise ValueError('Finite validation failed inside analysis')
    elif np.ptp(x) == 0 and np.ptp(y) == 0:
        result['correlation_reason'] = 'both_variables_constant'
    elif np.ptp(x) == 0:
        result['correlation_reason'] = 'first_variable_constant'
    elif np.ptp(y) == 0:
        result['correlation_reason'] = 'second_variable_constant'
    else:
        result['spearman'] = float(spearmanr(x, y).statistic)
        result['pearson'] = float(pearsonr(x, y).statistic)
        require(np.isfinite(result['spearman']) and np.isfinite(result['pearson']), 'Unexpected undefined nonconstant correlation')
    return result

def errors(target, pred):
    r = dict(reference_mean=float(target.mean()),prediction_mean=float(pred.mean()),mae=float(np.abs(target-pred).mean()),rmse=float(np.sqrt(np.mean((target-pred)**2))),bias=float((pred-target).mean()))
    r.update(correlation(target,pred))
    return r

def tails(values):
    """Disjoint inclusive tails, or missing for a collapsed quartile interval."""
    result = np.full(len(values), '', dtype='<U4')
    if len(values) < 4:
        return result, np.nan, np.nan, 'fewer_than_4_same_AA_positions'
    q25, q75 = np.quantile(values,[.25,.75],method='linear')
    if q25 == q75:
        return result, float(q25), float(q75), 'collapsed_q25_q75'
    result[values <= q25] = 'low'
    result[values >= q75] = 'high'
    return result, float(q25), float(q75), ''

def joint_class(r, h):
    result = np.full(len(r),'unclassified',dtype='<U24')
    for a in ['low','high']:
        for b in ['low','high']:
            result[(r == a) & (h == b)] = a+'_RMSF_'+b+'_Hchi'
    return result

def plan():
    return dict(protocol=VERSION,variant='full',source='corrected3063 original system OOF, held-out corresponding fold only, seed7 epoch40',records=3063,
        no_training=True,no_GPU=True,no_MD=True,no_thermodynamic_free_energy_fit=True,no_fitting_on_OOF=True,
        units=dict(rmsf='angstrom',chi_entropy='nats'),cohort='Both original effective RMSF and Hchi masks, no finite/error-based subset selection after preflight',
        cross_output_correlations=['Spearman','Pearson'],correlation_minimum=4,
        centering='Within each record AND amino-acid type, independently subtract the reference and prediction mean for each output on exactly the common valid positions; all common positions retained, including groups below 4',
        same_AA_correlations='Each observed record/type group retained; correlations computed only with >=4 common positions and two nonconstant variables, all undefined reasons saved',
        quartiles='Within each record/type on common valid positions, NumPy linear 25th/75th percentiles; require >=4 and q25<q75 per output. Values <=q25 are low, >=q75 high; ties may enlarge tails. Collapsed intervals, middle values and short groups remain unclassified.',
        four_categories='Reference-defined low/low, low/high, high/low, high/high. Compute per-record reference and predicted means, errors and reference-vs-prediction ordering on the same positions. Separate predicted-own-quartile classification provides a 5x5 confusion including unclassified.',
        aggregation='Primary record-macro: every computable record contributes one vote within each source/subset. Supply separate raw residue counts and micro rates. No pooled correlation, no fold-mean SD or significance claims.',
        recovery='Paired absolute/signed discrepancy between reference and predicted cross-output correlation; four-category rank recovery, never exp(H) or a fitted energy mapping')

def preflight(args):
    root, cohort, out = Path(args.training_root), Path(args.cohort), Path(args.output)
    manifest_path = root/'evaluation/prediction_manifest.jsonl'
    freeze_path = root/'training_freeze.json'
    all_path = root/'splits/all_ready.jsonl'
    cohort_freeze_path = cohort/'cohort_freeze.json'
    eval_audit_path = root/'evaluation/evaluation_audit.json'
    require(sha(freeze_path) == EXPECTED_SOURCE_FREEZE,'Wrong corrected training freeze')
    require(sha(manifest_path) == EXPECTED_PRED_MANIFEST,'Wrong original OOF prediction manifest')
    cf = json.loads(cohort_freeze_path.read_text())
    require(cf['source_freeze_sha256'] == EXPECTED_SOURCE_FREEZE and cf['records'] == 3063,'Cohort is not corrected3063')
    require(sha(all_path) == cf['source_manifest_sha256'],'all_ready source manifest changed')
    for name, digest in cf['files'].items():
        require(sha(cohort/name) == digest, 'Cohort file hash changed: '+name)
    source_freeze = json.loads(freeze_path.read_text())
    frozen_hashes = {name.replace('\\','/'): digest for name,digest in source_freeze['input_hashes'].items()}
    originals = lines(all_path)
    original_map = {r['sample_id']: r for r in originals}
    predictions = [r for r in lines(manifest_path) if r['variant'] == 'full']
    pred_map = {r['sample_id']: r for r in predictions}
    require(len(originals) == len(original_map) == len(predictions) == len(pred_map) == 3063,'Duplicate, missing or wrong record count')
    require(set(original_map) == set(pred_map),'Original/prediction IDs differ')
    ea = json.loads(eval_audit_path.read_text())
    require(ea['status'] == 'complete' and ea['epoch'] == 40 and ea['seed'] == 7 and ea['input_freeze_sha256'] == EXPECTED_SOURCE_FREEZE,'OOF export audit protocol changed')
    truth = pd.read_csv(cohort/'truth_residues.csv.gz',keep_default_na=False,float_precision='round_trip',dtype={'chain':str,'author_residue_id':str,'amino_acid':str})
    require(set(truth.sample_id) == set(original_map),'Truth IDs differ')
    require(len(truth) == 961321 and int(truth.valid_rmsf.sum()) == 961321 and int(truth.valid_chi_entropy.sum()) == 824988,'Unexpected frozen valid counts')
    require(not truth.duplicated(['sample_id','residue_index']).any(),'Duplicate truth positions')
    structures = pd.read_csv(cohort/'structure_manifest.tsv',sep='\t',keep_default_na=False).set_index('sample_id')
    require(set(structures.index) == set(original_map),'Structure manifest IDs differ')
    split_hashes, checkpoint_hashes, config_hashes = {}, {}, {}
    for fold in range(5):
        test_path, train_path = root/f'splits/fold_{fold}_test.jsonl', root/f'splits/fold_{fold}_train.jsonl'
        test, train = lines(test_path), lines(train_path)
        test_ids, train_ids = {r['sample_id'] for r in test}, {r['sample_id'] for r in train}
        require(len(test_ids) == len(test) and len(train_ids) == len(train),'Duplicate split records')
        require(test_ids.isdisjoint(train_ids) and test_ids|train_ids == set(original_map),'Invalid OOF train/test membership')
        require(test_ids == {s for s,m in pred_map.items() if m['fold'] == fold},'Predictions do not correspond to heldout test fold')
        split_hashes[test_path.as_posix()],split_hashes[train_path.as_posix()] = sha(test_path),sha(train_path)
        weight = root/f'runs/full/fold_{fold}/epoch_040.pt'
        weight_hash = sha(weight)
        require({m['checkpoint_sha256'] for m in predictions if m['fold']==fold} == {weight_hash},'OOF checkpoint hash mismatch')
        checkpoint_hashes[weight.as_posix()] = weight_hash
        config_path = root/f'configs/full_fold{fold}.yaml'
        config = yaml.safe_load(config_path.read_text())
        require(config['model']['metrics']==['rmsf','chi_entropy'] and config['seed']==7 and config['train']['epochs']==40,'Unexpected model output channel order')
        config_hashes[config_path.as_posix()] = sha(config_path)
    matched, bundle, reproduced = [], [], []
    for position,(sample,rows) in enumerate(truth.groupby('sample_id',sort=True)):
        rows = rows.sort_values('residue_index')
        m, original, st = pred_map[sample], original_map[sample], structures.loc[sample]
        require(np.array_equal(rows.residue_index,np.arange(len(rows))),sample+': residue order')
        require(m['dataset']==original['dataset']==rows.dataset.iloc[0]==st.cohort and m['fold']==original['metadata']['fold']==int(rows.fold.iloc[0])==int(st.fold),sample+': source/fold mismatch')
        prediction_path = Path(m['cache_path'])
        require(prediction_path == root/'evaluation/full'/(sample+'.npz'),sample+': wrong cache root')
        label_name = m['record']['label_path'].replace('\\','/').split('/')[-1]
        label_path, structure_path = root/'dataset/labels'/label_name, root/'dataset/structures'/(sample+'.pdb')
        require(sha(label_path)==m['label_sha256']==st.label_sha256==frozen_hashes['dataset/labels/'+label_path.name],sample+': original label hash')
        require(sha(structure_path)==m['input_sha256']==st.original_input_sha256==frozen_hashes['dataset/structures/'+structure_path.name],sample+': original structure hash')
        with np.load(label_path,allow_pickle=False) as label:
            metric_names = label['metric_names'].astype(str).tolist()
            require(len(metric_names)==len(set(metric_names)) and 'rmsf' in metric_names and 'chi_entropy' in metric_names,sample+': unrecognized label metrics')
            columns = [metric_names.index('rmsf'),metric_names.index('chi_entropy')]
            label_values,label_valid = label['values'].astype(float),label['mask'].astype(bool)
        with np.load(prediction_path,allow_pickle=False) as cache:
            mean, raw_target, raw_valid = cache['mean'].astype(np.float64),cache['target'].astype(np.float64),cache['valid'].astype(bool)
            require(raw_target.shape==raw_valid.shape==label_values.shape==label_valid.shape==(len(rows),len(metric_names)),sample+': raw cache/label channel shape')
            require(np.array_equal(raw_target,label_values,equal_nan=True),sample+': all raw corrected label channels mismatch')
            target,valid = raw_target[:,columns],raw_valid[:,columns]
            aa, chain, author = cache['amino_acid'].astype(str),cache['chain'].astype(str),cache['author_residue_id'].astype(str)
        require(mean.shape==target.shape==valid.shape==(len(rows),2),sample+': two-output shape')
        require(np.isfinite(mean).all() and np.isfinite(target[valid]).all(),sample+': nonfinite outputs/valid labels')
        truth_target = rows[['target_rmsf','target_chi_entropy']].to_numpy(float)
        truth_valid = rows[['valid_rmsf','valid_chi_entropy']].to_numpy(bool)
        require(np.array_equal(target,truth_target,equal_nan=True) and np.array_equal(valid,truth_valid),sample+': exact cohort target/valid mismatch')
        require(np.array_equal(aa,rows.amino_acid.to_numpy(str)) and np.array_equal(chain,rows.chain.to_numpy(str)) and np.array_equal(author,rows.author_residue_id.to_numpy(str)),sample+': residue identities changed')
        require(np.array_equal(target,label_values[:,columns],equal_nan=True),sample+': raw corrected label values mismatch')
        require(np.all(label_valid[:,columns][valid]),sample+': effective mask not contained in original mask')
        common = valid.all(axis=1)
        require(set(aa[common]).issubset(AA18),sample+': unexpected Ala/Gly/unknown common Hchi type')
        matched.append(dict(sample_id=sample,dataset=m['dataset'],fold=m['fold'],residues=len(rows),valid_rmsf=int(valid[:,0].sum()),valid_hchi=int(valid[:,1].sum()),common_valid=int(common.sum()),common_excluded=int((~common).sum()),raw_label_metrics='|'.join(metric_names),analysis_label_columns='|'.join(map(str,columns)),prediction_path=prediction_path.as_posix(),prediction_sha256=sha(prediction_path),label_path=label_path.as_posix(),label_sha256=m['label_sha256'],structure_path=structure_path.as_posix(),structure_sha256=m['input_sha256'],checkpoint_sha256=m['checkpoint_sha256']))
        for k,name in enumerate(['rmsf','chi_entropy']):
            ix = valid[:,k]
            reproduced.append(dict(sample_id=sample,dataset=m['dataset'],fold=m['fold'],metric=name,count=int(ix.sum()),**errors(target[ix,k],mean[ix,k])))
        bundle.append((sample,m['dataset'],m['fold'],target[common],mean[common],aa[common],np.flatnonzero(common),chain[common],author[common]))
        if (position+1)%500 == 0:
            print('Preflight',position+1,'/3063',flush=True)
    matched = pd.DataFrame(matched)
    require(matched.valid_rmsf.sum()==961321 and matched.valid_hchi.sum()==matched.common_valid.sum()==824988,'Common-valid total mismatch')
    matched.to_csv(out/'frozen_input_records.csv',index=False)
    reproduced = pd.DataFrame(reproduced)
    reproduced.to_csv(out/'reproduced_main_record_metrics.csv',index=False)
    macro = reproduced.groupby('metric')[['mae','spearman']].mean()
    require(abs(macro.loc['rmsf','mae']-0.5771038628004094)<1e-6 and abs(macro.loc['chi_entropy','mae']-0.16586574902130044)<1e-6,'Main v20 mean MAE not reproduced')
    require(abs(macro.loc['rmsf','spearman']-0.7961867220004638)<1e-6 and abs(macro.loc['chi_entropy','spearman']-0.9411305304037728)<1e-6,'Main v20 ranking not reproduced')
    matched.groupby('dataset')[['residues','valid_rmsf','valid_hchi','common_valid','common_excluded']].sum().to_csv(out/'raw_source_counts.csv')
    frozen = dict(status='passed_input_freeze',protocol=VERSION,utc=datetime.now(timezone.utc).isoformat(),records=3063,common_valid=824988,source_counts=matched.groupby('dataset').size().astype(int).to_dict(),source_manifest_sha256=sha(all_path),source_training_freeze_sha256=sha(freeze_path),prediction_manifest_sha256=sha(manifest_path),cohort_freeze_sha256=sha(cohort_freeze_path),truth_sha256=sha(cohort/'truth_residues.csv.gz'),evaluation_audit_sha256=sha(eval_audit_path),frozen_records_sha256=sha(out/'frozen_input_records.csv'),split_hashes=split_hashes,checkpoint_hashes=checkpoint_hashes,config_hashes=config_hashes,reproduced_v20_record_macro=macro.to_dict(),analysis_code_sha256=sha(__file__),analysis_plan_sha256=sha(out/'analysis_plan.json'),input_assertions=['exact 3063 unique IDs/source/fold','OOF test record absent from matching train fold','same seed7 epoch40 hashes','label columns selected explicitly by names rmsf and chi_entropy; raw SASA channel is not an output','current corrected label and original PDB hashes match training freeze','exact both channel targets and effective masks match frozen cohort','exact chain/author number/AA order','all finite saved predictions','all common Hchi positions belong to 18 eligible types','v20 main recordmacro MAE/rho independently reproduced'],exclusion_policy='No failed files allowed, no accuracy/finite-driven residue removal. Hchi-ineligible positions are explicitly outside paired analysis, but remain in raw counts.')
    jsonwrite(out/'input_freeze.json',frozen)
    return bundle,frozen,matched

def summarize_numeric(frame, group_keys, numeric, out_path):
    rows = []
    for keys, group in frame.groupby(group_keys,dropna=False,sort=True):
        keys = (keys,) if not isinstance(keys,tuple) else keys
        row = dict(zip(group_keys,keys))
        row['record_rows'] = len(group)
        row['unique_records'] = int(group.sample_id.nunique())
        row['common_residues'] = int(group['count'].sum()) if 'count' in group else None
        for column in numeric:
            values = group[column].dropna().to_numpy(float)
            row[column+'_n_records'] = len(values)
            row[column+'_record_macro'] = float(values.mean()) if len(values) else np.nan
            row[column+'_median'] = float(np.median(values)) if len(values) else np.nan
            row[column+'_q25'] = float(np.quantile(values,.25)) if len(values) else np.nan
            row[column+'_q75'] = float(np.quantile(values,.75)) if len(values) else np.nan
        rows.append(row)
    pd.DataFrame(rows).to_csv(out_path,index=False)

def analyze(bundle,out):
    coupling, same_aa, groups, thresholds, confusion, residues = [], [], [], [], [], []
    missing = []
    for position,(sample,dataset,fold,target,pred,aa,residue_index,chain,author) in enumerate(bundle):
        base = dict(sample_id=sample,dataset=dataset,fold=fold)
        tc,pc = target.copy(),pred.copy()
        reference_class = np.full(len(target),'unclassified',dtype='<U24')
        predicted_class = reference_class.copy()
        residue_reason = np.full(len(target),'',dtype='<U110')
        for a in sorted(set(aa)):
            ix = np.flatnonzero(aa==a)
            tc[ix] -= target[ix].mean(axis=0)
            pc[ix] -= pred[ix].mean(axis=0)
            tr,th,pr,ph = [],[],[],[]
            row = dict(**base,amino_acid=a,count=len(ix))
            true_corr,pred_corr = correlation(target[ix,0],target[ix,1]),correlation(pred[ix,0],pred[ix,1])
            for prefix,c in [('reference',true_corr),('prediction',pred_corr)]:
                row.update({prefix+'_'+key:value for key,value in c.items()})
            for channel,label in enumerate(['rmsf','hchi']):
                row.update({label+'_'+key:value for key,value in errors(target[ix,channel],pred[ix,channel]).items()})
            same_aa.append(row)
            labels = []
            reasons = []
            for values, prefix in [(target[ix,0],'reference_rmsf'),(target[ix,1],'reference_hchi'),(pred[ix,0],'prediction_rmsf'),(pred[ix,1],'prediction_hchi')]:
                state,q25,q75,reason = tails(values)
                labels.append(state)
                reasons.append(reason)
                thresholds.append(dict(**base,amino_acid=a,output=prefix,count=len(ix),q25=q25,q75=q75,undefined_reason=reason,low_count=int((state=='low').sum()),high_count=int((state=='high').sum()),middle_or_undefined_count=int((state=='').sum())))
            reference_class[ix] = joint_class(labels[0],labels[1])
            predicted_class[ix] = joint_class(labels[2],labels[3])
            for j,k in enumerate(ix):
                if reference_class[k] == 'unclassified':
                    reason = '|'.join(sorted(set(r for r in reasons[:2] if r)))
                    if not reason:
                        reason = 'reference_RMSF_middle' if labels[0][j]=='' else ''
                        if labels[1][j]=='':
                            reason += ('|' if reason else '')+'reference_Hchi_middle'
                    residue_reason[k] = reason
        for mode,y,p in [('cross_type_raw',target,pred),('within_record_AA_centered',tc,pc)]:
            r = dict(**base,mode=mode,count=len(target),observed_AA_types=len(set(aa)))
            a,b = correlation(y[:,0],y[:,1]),correlation(p[:,0],p[:,1])
            for prefix,c in [('reference',a),('prediction',b)]:
                r.update({prefix+'_'+key:value for key,value in c.items()})
            for stat in ['spearman','pearson']:
                r[stat+'_prediction_minus_reference'] = b[stat]-a[stat]
                r[stat+'_absolute_discrepancy'] = abs(b[stat]-a[stat])
            coupling.append(r)
        for name in CLASSES:
            selected = reference_class==name
            n = int(selected.sum())
            n_pred = int((predicted_class==name).sum())
            n_match = int((selected & (predicted_class==name)).sum())
            r = dict(**base,reference_class=name,count=n,predicted_class_count=n_pred,intersection_count=n_match,reference_fraction=n/len(target) if len(target) else np.nan,predicted_fraction=n_pred/len(target) if len(target) else np.nan,rank_class_recall=n_match/n if n and name!='unclassified' else np.nan,rank_class_precision=n_match/n_pred if n_pred and name!='unclassified' else np.nan)
            for channel,label in enumerate(['rmsf','hchi']):
                if n:
                    r.update({label+'_'+key:value for key,value in errors(target[selected,channel],pred[selected,channel]).items()})
                else:
                    r.update({label+'_'+key:np.nan for key in ['reference_mean','prediction_mean','mae','rmse','bias','spearman','pearson']})
                    r[label+'_correlation_reason'] = 'no_reference_class_positions'
            groups.append(r)
            for pn in CLASSES:
                confusion.append(dict(**base,reference_class=name,predicted_class=pn,count=int((selected & (predicted_class==pn)).sum())))
        for reason in sorted(set(residue_reason[reference_class=='unclassified'])):
            missing.append(dict(**base,reason=reason,count=int((residue_reason==reason).sum())))
        # A complete remote row source supports regeneration of plots without
        # downloading or changing the original caches. It contains paired
        # valid positions only, and every omitted raw position is counted.
        residues.append(pd.DataFrame(dict(sample_id=sample,dataset=dataset,fold=fold,residue_index=residue_index,chain=chain,author_residue_id=author,amino_acid=aa,target_rmsf=target[:,0],target_hchi=target[:,1],pred_rmsf=pred[:,0],pred_hchi=pred[:,1],reference_class=reference_class,predicted_class=predicted_class,reference_unclassified_reason=residue_reason)))
        if (position+1)%500==0:
            print('Analysis',position+1,'/3063',flush=True)
    coupling, same_aa, groups, thresholds, confusion = map(pd.DataFrame,[coupling,same_aa,groups,thresholds,confusion])
    missing = pd.DataFrame(missing)
    all_pairs = pd.concat(residues,ignore_index=True)
    for name,frame in [('per_record_correlations',coupling),('per_record_same_AA',same_aa),('per_record_reference_classes',groups),('per_record_type_quartiles',thresholds),('per_record_class_confusion',confusion),('reference_unclassified_reasons',missing)]:
        frame.to_csv(out/(name+'.csv'),index=False)
    require(len(all_pairs)==824988,'Lost paired residue rows')
    require(confusion['count'].sum()==824988 and groups['count'].sum()==824988,'Class coverage is not exhaustive')
    all_pairs.to_csv(out/'paired_common_residues.csv.gz',index=False,compression='gzip')
    # Duplicate rows only for a separate all-source macro summary; no data
    # are duplicated within a source summary or individual record file.
    def with_all(frame):
        return pd.concat([frame,frame.assign(dataset='all')],ignore_index=True)
    numeric = ['reference_spearman','prediction_spearman','reference_pearson','prediction_pearson','spearman_prediction_minus_reference','spearman_absolute_discrepancy','pearson_prediction_minus_reference','pearson_absolute_discrepancy']
    summarize_numeric(with_all(coupling),['dataset','mode'],numeric,out/'correlation_record_macro.csv')
    numeric_aa = ['reference_spearman','prediction_spearman','reference_pearson','prediction_pearson']+[label+'_'+metric for label in ['rmsf','hchi'] for metric in ['reference_mean','prediction_mean','mae','spearman','pearson']]
    summarize_numeric(with_all(same_aa),['dataset','amino_acid'],numeric_aa,out/'same_AA_record_macro.csv')
    numeric_groups = ['reference_fraction','predicted_fraction','rank_class_recall','rank_class_precision']+[label+'_'+metric for label in ['rmsf','hchi'] for metric in ['reference_mean','prediction_mean','mae','rmse','bias','spearman','pearson']]
    # Empty class records count in fraction denominators but not amplitude/
    # error means. The column-specific number of computable records is kept.
    summarize_numeric(with_all(groups),['dataset','reference_class'],numeric_groups,out/'reference_class_record_macro.csv')
    with_all(confusion).groupby(['dataset','reference_class','predicted_class'],sort=True)['count'].sum().reset_index().to_csv(out/'class_confusion_residue_counts.csv',index=False)
    with_all(missing).groupby(['dataset','reason'],sort=True).agg(residue_count=('count','sum'),record_count=('sample_id','nunique')).reset_index().to_csv(out/'unclassified_reasons_summary.csv',index=False)
    with_all(thresholds).groupby(['dataset','amino_acid','output','undefined_reason'],sort=True,dropna=False).agg(record_type_groups=('sample_id','size'),residue_count=('count','sum'),low_count=('low_count','sum'),high_count=('high_count','sum'),middle_or_undefined_count=('middle_or_undefined_count','sum')).reset_index().to_csv(out/'type_quartile_coverage_summary.csv',index=False)
    reason_frames = []
    for frame,unit in [(coupling,'record_coupling'),(same_aa,'record_same_AA'),(groups,'reference_class')]:
        for column in [c for c in frame.columns if c.endswith('_correlation_reason')]:
            tmp = frame.loc[frame[column].ne(''),['sample_id','dataset','count',column]].rename(columns={column:'reason'}).copy()
            tmp['analysis_unit'],tmp['correlation'] = unit,column
            reason_frames.append(tmp)
    reason_frame = pd.concat(reason_frames,ignore_index=True)
    reason_frame.to_csv(out/'undefined_correlation_reasons.csv',index=False)
    with_all(reason_frame).groupby(['dataset','analysis_unit','correlation','reason']).agg(rows=('sample_id','size'),records=('sample_id','nunique'),residues=('count','sum')).reset_index().to_csv(out/'undefined_correlation_summary.csv',index=False)
    aa_count = all_pairs.groupby(['dataset','amino_acid','reference_class'])[['target_rmsf']].size().rename('residue_count').reset_index()
    aa_count.to_csv(out/'raw_AA_class_counts.csv',index=False)
    # Reproduce manuscript Full performance directly before interpreting any
    # new paired-subset calculation. RMSF common subset omits Ala/Gly and is
    # explicitly distinct from the 961321-position main RMSF score.
    return dict(records=3063,paired_common_residues=len(all_pairs),raw_AA_common_counts=all_pairs.amino_acid.value_counts().astype(int).to_dict(),reference_class_counts=all_pairs.reference_class.value_counts().astype(int).to_dict(),prediction_class_counts=all_pairs.predicted_class.value_counts().astype(int).to_dict(),correlation_record_rows=len(coupling),same_AA_record_rows=len(same_aa),no_failed_inputs=True,classification_sum_matches_all_common=True)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--training-root',required=True)
    parser.add_argument('--cohort',required=True)
    parser.add_argument('--output',required=True)
    args = parser.parse_args()
    out = Path(args.output)
    out.mkdir(parents=True,exist_ok=True)
    require(not (out/'run_audit.json').exists(),'Completed output audit exists; use a new directory rather than overwrite')
    start = time.monotonic()
    jsonwrite(out/'analysis_plan.json',plan())
    bundle,frozen,matched = preflight(args)
    analysis_start = time.monotonic()
    result = analyze(bundle,out)
    require(sha(__file__)==frozen['analysis_code_sha256'] and sha(out/'analysis_plan.json')==frozen['analysis_plan_sha256'],'Code/plan changed during analysis')
    require(sha(out/'frozen_input_records.csv')==frozen['frozen_records_sha256'],'Record input freeze changed')
    # Recheck every prediction, label and input file after statistics. No
    # downstream result may certify a concurrently changed source cache.
    for row in matched.itertuples(index=False):
        for path,digest in [(row.prediction_path,row.prediction_sha256),(row.label_path,row.label_sha256),(row.structure_path,row.structure_sha256)]:
            require(sha(path)==digest,'Input changed during analysis: '+path)
    result.update(status='complete',protocol=VERSION,utc=datetime.now(timezone.utc).isoformat(),CPU_only=True,existing_audits_modified=False,full_duration_seconds=time.monotonic()-start,analysis_duration_seconds=time.monotonic()-analysis_start,input_freeze_sha256=sha(out/'input_freeze.json'),analysis_code_sha256=sha(__file__),files={path.name:sha(path) for path in sorted(out.iterdir()) if path.is_file() and path.name!='run_audit.json'},interpretation='Descriptive dependence and paired recovery of two sampled observables; no physical free energy, kinetic rate, causal mechanism or independent-repeat significance is inferred.')
    jsonwrite(out/'run_audit.json',result)
    print(json.dumps(result,ensure_ascii=False,indent=2),flush=True)

if __name__=='__main__':
    try:
        main()
    except Exception as error:
        print(type(error).__name__+': '+str(error),file=sys.stderr,flush=True)
        raise
