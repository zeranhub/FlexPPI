"""Read-only diagnosis of source composition, OOF membership and exposure status.

Reads existing metrics only; no model, checkpoint tensors, network or inference.
Writes only this investigation folder. No source freeze or classification changes.
"""
from pathlib import Path
import hashlib,json,datetime,zipfile,xml.etree.ElementTree as ET
import numpy as np
import pandas as pd

OUT=Path(__file__).resolve().parent
ROOT=OUT.parents[2]
CORE=ROOT.parents[1]/'core_corrected_pipeline'

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def csv(p):return pd.read_csv(p,keep_default_na=False,na_values=[''])
def dump(p,x):Path(p).write_text(json.dumps(x,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
def boolean(x):return str(x).lower()=='true'

def main():
    sources={
      'membership':ROOT/'protein_complex_benchmark/classification/membership_frozen.csv',
      'cohort':ROOT/'protein_complex_benchmark/inputs/cohort_records.csv',
      'metrics':ROOT/'protein_complex_benchmark/evaluation/per_record_metrics.csv',
      'exposure':ROOT/'references/training_exposure_membership_frozen.csv',
      'exposure_protocol':ROOT/'references/TRAINING_EXPOSURE_SENSITIVITY_PROTOCOL.json',
      'exposure_audit':ROOT/'references/training_exposure_protocol_sources_sha256.json',
      'own_oof':ROOT/'references/training_exposure_own_oof_details.csv',
      'own_family':ROOT/'references/training_exposure_own_family_details.csv',
      'prediction_manifest':CORE/'neural_exports/evaluation/prediction_manifest.jsonl',
      'prediction_export_audit':CORE/'neural_exports/evaluation/evaluation_audit.json',
      'training_freeze':CORE/'training_payload/training_freeze.json',
      'training_config':CORE/'training_payload/configs/full_fold0.yaml',
      'training_sampler':CORE/'training_payload/src/flexppi/trainer.py',
      'original_aggregation':CORE/'type_analysis/analyze_corrected_main_comparison.py',
      'BF2_config':ROOT/'backflip/weights/backflip-2.1/config.yaml',
    }
    initial={k:sha(p) for k,p in sources.items()}
    member=csv(sources['membership']);cohort=csv(sources['cohort']);metrics=csv(sources['metrics']); exposure=csv(sources['exposure']); own=csv(sources['own_oof'])
    assert len(member)==len(cohort)==3063 and member.sample_id.is_unique and cohort.sample_id.is_unique
    records=cohort.merge(member[['sample_id','benchmark_group']],on='sample_id',validate='one_to_one')
    assert records.benchmark_group.value_counts().to_dict()=={'protein':2203,'protein_complex':860}
    # Cohort origin: do not use old single_chain/other_multichain labels as groups.
    composition=records.groupby(['benchmark_group','dataset','source_group'],dropna=False).agg(records=('sample_id','size'),valid_RMSF_positions=('n_valid_rmsf','sum')).reset_index()
    composition['within_group_fraction']=composition.records/composition.groupby('benchmark_group').records.transform('sum')
    composition.to_csv(OUT/'source_composition.csv',index=False)
    # Actual own-fold verification is repeated against saved prediction exports
    # and train/test manifests, rather than trusting the word OOF in a report.
    actual_train={};actual_test={}; fold_rows=[]; sampling_rows=[]
    for fold in range(5):
        for role,store in [('train',actual_train),('test',actual_test)]:
            p=CORE/'training_payload/splits'/('fold_%s_%s.jsonl'%(fold,role))
            sources['fold_%s_%s'%(fold,role)]=p;initial['fold_%s_%s'%(fold,role)]=sha(p)
            store[fold]={json.loads(s)['sample_id'] for s in p.read_text(encoding='utf-8').splitlines() if s.strip()}
        assert not actual_train[fold]&actual_test[fold]
        assert actual_train[fold]|actual_test[fold]==set(records.sample_id)
        assert actual_test[fold]==set(records.loc[records.original_fold==fold,'sample_id'])
        train=records[records.sample_id.isin(actual_train[fold])]
        for ds,g in train.groupby('dataset'):
            sampling_rows.append({'fold':fold,'dataset':ds,'training_records':len(g),'all_training_records':len(train),
                'training_record_fraction':len(g)/len(train),'configured_source_probability':1/train.dataset.nunique(),
                'expected_source_draws_per_epoch':len(train)/train.dataset.nunique(),
                'mean_expected_draws_per_record_per_epoch':len(train)/train.dataset.nunique()/len(g),
                'mean_expected_draws_per_record_over_40_epochs':40*len(train)/train.dataset.nunique()/len(g),
                'interpretation':'Expectation under source/system balanced sampling with replacement; actual sampled visits not reconstructed.'})
        for role,ids in [('train',actual_train[fold]),('test',actual_test[fold])]:
            part=records[records.sample_id.isin(ids)]
            for (group,dataset),g in part.groupby(['benchmark_group','dataset']):
                fold_rows.append({'fold':fold,'role':role,'group':group,'dataset':dataset,'records':len(g)})
    pd.DataFrame(fold_rows).to_csv(OUT/'actual_fold_source_counts.csv',index=False)
    pd.DataFrame(sampling_rows).to_csv(OUT/'configured_training_source_sampling.csv',index=False)
    manifest_count=0; own_full=[]
    for line in sources['prediction_manifest'].read_text(encoding='utf-8').splitlines():
        if not line.strip():continue
        row=json.loads(line);manifest_count+=1
        fold=int(row['fold']);sid=row['sample_id']
        assert sid in actual_test[fold] and sid not in actual_train[fold]
        assert int(row['record']['metadata']['fold'])==fold
        if row['variant']=='full':own_full.append({'sample_id':sid,'prediction_fold_verified':fold,'held_out_record':True})
    assert manifest_count==6126 and len(own_full)==3063
    own_full=pd.DataFrame(own_full).merge(own[own.variant=='full'],on='sample_id',validate='one_to_one').merge(records[['sample_id','benchmark_group']],on='sample_id',validate='one_to_one')
    own_full['same_family_in_train']=own_full.same_family_other_train_records.gt(0)
    family_table=own_full.groupby(['benchmark_group','dataset'],dropna=False).agg(records=('sample_id','size'),same_family_in_train=('same_family_in_train','sum'),same_biological_cluster_in_train=('same_biological_cluster_other_train_records',lambda x:int((x>0).sum())),matches_family_fold=('matches_family_fold',lambda x:int(sum(boolean(y) for y in x)))).reset_index()
    family_table.to_csv(OUT/'own_OOF_and_family_overlap.csv',index=False)
    assert not own_full.sample_in_corresponding_train_manifest.map(boolean).any()
    assert own_full.sample_in_corresponding_test_manifest.map(boolean).all()
    # Compare the exact same records, separately by method context and source.
    full=metrics[metrics.method=='Full'].set_index('sample_id')
    comparisons=[];paired=[]
    for method in ['backflip-2.1 (assembly)','backflip-2.1 (chains)']:
        comparator=metrics[metrics.method==method].set_index('sample_id')
        assert set(full.index)==set(comparator.index)==set(records.sample_id)
        p=records[['sample_id','benchmark_group','dataset','source_group','original_fold','family_fold','family_component_id']].merge(own_full[['sample_id','same_family_in_train']],on='sample_id',validate='one_to_one')
        for metric in ['spearman','top10_ap','z_mae','mae','pearson','rmse','bias']:
            p['Full_'+metric]=p.sample_id.map(full[metric]);p['BF2_'+metric]=p.sample_id.map(comparator[metric]);p['BF2_minus_Full_'+metric]=p['BF2_'+metric]-p['Full_'+metric]
        p['BF2_context']=method;paired.append(p)
        scopes=[('all',records.sample_id),('protein',records.loc[records.benchmark_group=='protein','sample_id']),('protein_complex',records.loc[records.benchmark_group=='protein_complex','sample_id'])]
        scopes.extend(('protein:'+ds,g.sample_id) for ds,g in records[records.benchmark_group=='protein'].groupby('dataset'))
        scopes.extend(('protein_source:'+sg,g.sample_id) for sg,g in records[records.benchmark_group=='protein'].groupby('source_group'))
        for scope,ids in scopes:
            part=p[p.sample_id.isin(ids)]
            for metric in ['spearman','top10_ap','z_mae','mae','pearson','rmse','bias']:
                valid=part[['Full_'+metric,'BF2_'+metric]].notna().all(axis=1);q=part[valid]
                comparisons.append({'comparator':method,'scope':scope,'metric':metric,'records':len(part),'paired_defined':len(q),
                    'Full_mean':q['Full_'+metric].mean(),'BF2_mean':q['BF2_'+metric].mean(),
                    'BF2_minus_Full':q['BF2_minus_Full_'+metric].mean(),
                    'BF2_better_records':int((q['BF2_minus_Full_'+metric]*(1 if metric in ['spearman','top10_ap','pearson'] else -1)>0).sum())})
    comparison=pd.DataFrame(comparisons);comparison.to_csv(OUT/'source_paired_comparison.csv',index=False)
    pd.concat(paired).to_csv(OUT/'paired_record_differences.csv',index=False)
    # Quantify each source's contribution to the protein macro gap exactly.
    contribution=[]
    for method in ['backflip-2.1 (assembly)','backflip-2.1 (chains)']:
        for metric in ['spearman','top10_ap','z_mae','mae']:
            allrow=comparison[(comparison.comparator==method)&(comparison.scope=='protein')&(comparison.metric==metric)].iloc[0]
            for ds in ['atlas','dynarepo']:
                row=comparison[(comparison.comparator==method)&(comparison.scope=='protein:'+ds)&(comparison.metric==metric)].iloc[0]
                contribution.append({'comparator':method,'metric':metric,'dataset':ds,'records':int(row.records),'protein_fraction':row.records/2203,'source_gap_BF2_minus_Full':row.BF2_minus_Full,'weighted_contribution':row.BF2_minus_Full*row.records/2203,'share_of_total_signed_gap':row.BF2_minus_Full*row.records/2203/allrow.BF2_minus_Full})
    pd.DataFrame(contribution).to_csv(OUT/'protein_gap_source_contribution.csv',index=False)
    # BF2's own training IDs are unknown. Other methods' published ID splits
    # cannot be relabelled as BF2 exposure or used to assert BF2 leakage.
    exp_summary=exposure.merge(member[['sample_id','benchmark_group']],on='sample_id',validate='many_to_one').groupby(['method','benchmark_group','dataset','documented_split','binding_strength'],dropna=False).size().reset_index(name='records')
    exp_summary.to_csv(OUT/'exposure_status_counts.csv',index=False)
    bf2_exp=exposure[exposure.method=='backflip-2.1']
    assert set(bf2_exp.documented_split)=={'unknown'} and len(bf2_exp)==3063
    own_strata=[]
    for p in paired:
        for (group,ds,fam),g in p.groupby(['benchmark_group','dataset','same_family_in_train']):
            if group!='protein':continue
            for metric in ['spearman','top10_ap','z_mae','mae']:
                own_strata.append({'comparator':g.BF2_context.iloc[0],'dataset':ds,'Full_same_family_other_record_in_train':bool(fam),'records':len(g),'metric':metric,'Full_mean':g['Full_'+metric].mean(),'BF2_mean':g['BF2_'+metric].mean(),'BF2_minus_Full':g['BF2_minus_Full_'+metric].mean()})
    pd.DataFrame(own_strata).to_csv(OUT/'protein_gap_by_Full_family_exposure.csv',index=False)
    # Existing external-method exposure analysis, paired only within each method.
    external_rows=[]
    for em,method in [('DynaProt-M','DynaProt-M'),('backflip-1.0','backflip-1.0 (assembly)'),('backflip-1.0-seq','backflip-1.0-seq (assembly)')]:
        comp=metrics[metrics.method==method].set_index('sample_id')
        for split,g in exposure[(exposure.method==em)&(exposure.dataset=='atlas')].groupby('documented_split'):
            ids=g.sample_id
            for metric in ['spearman','top10_ap','z_mae','mae']:
                f=full.loc[ids,metric];b=comp.loc[ids,metric]
                external_rows.append({'method':method,'exposure_manifest_method':em,'documented_split':split,'records':len(ids),'metric':metric,'Full_mean':f.mean(),'released_model_mean':b.mean(),'released_model_minus_Full':(b-f).mean(),'use_limit':'This is this named comparator exposure category; it is NOT a BF2 train/test category.'})
    pd.DataFrame(external_rows).to_csv(OUT/'other_model_published_exposure_strata.csv',index=False)
    manuscript=next(Path('E:/Project/FlexPPI/paper/revision_v14_20261010').glob('paper_v63*.docx'))
    sources['manuscript_v63']=manuscript;initial['manuscript_v63']=sha(manuscript)
    ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
    xml=ET.fromstring(zipfile.ZipFile(manuscript).read('word/document.xml'))
    paras=[''.join(t.text or '' for t in p.findall('.//w:t',ns)) for p in xml.findall('.//w:p',ns)]
    indices=[36,38,39,41,44,46,47,69,104,105,106,114,115,116,145,147,151,360]
    pd.DataFrame([{'paragraph_zero_based':i,'text':paras[i]} for i in indices]).to_csv(OUT/'manuscript_scope_evidence.csv',index=False)
    assert {k:sha(p) for k,p in sources.items()}==initial,'read-only source changed during review'
    audit={'generated_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source_hashes':{k:{'path':str(p),'sha256':initial[k]} for k,p in sources.items()},
        'status':'complete_read_only_review','source_files_modified':False,'network_used':False,'model_inference_run':False,'model_checkpoint_loaded':False,
        'records':3063,'protein_records':2203,'protein_atlas_records':1929,'protein_dynarepo_records':274,
        'actual_prediction_manifest_rows_verified':6126,'OOF_record_train_intersection':0,
        'BF2_own_training_ID_manifest_status':'unknown_all_3063','no_BF2_leakage_claim':True,
        'code_sha256':sha(Path(__file__))}
    dump(OUT/'review_audit.json',audit)
    print(comparison[(comparison.comparator=='backflip-2.1 (assembly)')&(comparison.scope.isin(['all','protein','protein_complex','protein:atlas','protein:dynarepo']))&(comparison.metric.isin(['spearman','mae','top10_ap','z_mae']))].to_string(index=False))
    print('\nFAMILY\n'+family_table.to_string(index=False))
    print('\nCONTRIBUTION\n'+pd.DataFrame(contribution).query("comparator == 'backflip-2.1 (assembly)'").to_string(index=False))

if __name__=='__main__':main()
