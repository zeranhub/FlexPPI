"""Independent verification, deliberately not importing evaluate_screening.py."""
from pathlib import Path
from collections import defaultdict,Counter
from datetime import datetime,timezone
import csv,json,hashlib,math,warnings
import numpy as np
from scipy.stats import spearmanr

HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
RESULTS=HERE/'results_fixed_protocol'
read=lambda p:json.loads(Path(p).read_text(encoding='utf8'))
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def csvrows(p):
    with Path(p).open(encoding='utf-8-sig',newline='') as handle:return list(csv.DictReader(handle))
def truthflag(x):return str(x).lower() in {'true','1','yes','passed'}

FEATURES={
 'zero_change':[],
 'sequence_only':['sequence_n_mutations','BLOSUM62_deficit'],
 'static_only':['delta_TCR_pMHC_BSA_A2','delta_peptide_TCR_contact_residue_pairs_4p5A'],
 'sequence_plus_static':['sequence_n_mutations','BLOSUM62_deficit','delta_TCR_pMHC_BSA_A2','delta_peptide_TCR_contact_residue_pairs_4p5A'],
 'sequence_plus_static_plus_R':['sequence_n_mutations','BLOSUM62_deficit','delta_TCR_pMHC_BSA_A2','delta_peptide_TCR_contact_residue_pairs_4p5A','R'],
 'sequence_plus_static_plus_H':['sequence_n_mutations','BLOSUM62_deficit','delta_TCR_pMHC_BSA_A2','delta_peptide_TCR_contact_residue_pairs_4p5A','H'],
 'sequence_plus_static_plus_RH':['sequence_n_mutations','BLOSUM62_deficit','delta_TCR_pMHC_BSA_A2','delta_peptide_TCR_contact_residue_pairs_4p5A','R','H'],
}
PRIMARY={'nyeso':'2BNQ','tax':'1AO7'}
fields=['analysis','label_variant','system','panel_id','template_pdb','source_system','source_template_pdb','model']

def main():
    assert sha(HERE/'evaluate_screening.py')=='f303f3fec546559ce7369e5570a81efa1d02cd51885a6fb9713de1d7be05b779'
    assert sha(HERE/'fixed_analysis_protocol.json')=='79e2bb8d848c5843619fdb58452d008e8637bb30ae1b26284969c11a2e2406dc'
    assert sha(ROOT/'feature_extraction/fixed_regions.json')=='8666f5fb4a112c14210bfd0231d6797a6a2356e5e652aee254c30afdff5cf8c4'
    assert sha(ROOT/'feature_extraction/candidate_features.csv')=='f28d14f1b549cea6bd2c21f883b8854336cd84d95cf36527636a2f259c36fe3e'
    audit=read(RESULTS/'evaluation_execution_audit.json')
    assert audit['status']=='completed' and audit['mode']=='full' and audit['prediction_rows']==2037
    assert audit['checkpoint_sha256']=='b3eed060919df15f075329bc892f959d9355232be5feb561417edd49ea82b39f'
    rules=read(ROOT/'experimental_labels/label_rule_freeze.json')
    split=read(HERE/'candidate_splits.json')
    definitions={(r['system'],r['sequence']):r for r in split['candidates']}
    labels={};metadata={}
    for variant,info in rules['variants'].items():
        assert sha(info['candidate_file'])==info['candidate_file_sha256']
        labels[variant]={}
        for row in csvrows(info['candidate_file']):
            system='nyeso' if row['tcr']=='1G4' else 'tax'
            definition=definitions[system,row['peptide_sequence']]
            cid=row['candidate_id']
            metadata[cid]=dict(definition,panel_id=row['panel_id'])
            labels[variant][cid]=float(row['delta_pKD_vs_parent'])
            assert row['censoring']=='none_reported'
    features={(r['candidate_id'],r['template_pdb']):r for r in csvrows(ROOT/'feature_extraction/candidate_features.csv')}
    assert len(features)==76
    for (cid,t),f in features.items():
        assert all(truthflag(f[k]) for k in ['preparation_qc_pass','dynamic_qc_pass','static_eligible'])
        assert int(f['sequence_n_mutations'])==metadata[cid]['mutation_count']
    observed=csvrows(RESULTS/'candidate_predictions.csv')
    coefs=csvrows(RESULTS/'ridge_coefficients_and_training_scales.csv')
    coef_index=defaultdict(dict)
    for c in coefs:
        k=tuple(c[f] for f in ['analysis','label_variant','system','template_pdb','source_system','source_template_pdb','fold_id','model'])
        assert c['feature'] not in coef_index[k]
        coef_index[k][c['feature']]=c
    grouped=defaultdict(list)
    for r in observed:
        key=tuple(r[f] for f in ['analysis','label_variant','system','template_pdb','source_system','source_template_pdb','fold_id','model'])
        grouped[key].append(r)
    max_prediction_difference=max_scale_difference=max_coefficient_difference=0.
    prediction_checks=[];train_fits=0;all_rebuilt=[]
    for key,rows in grouped.items():
        analysis,variant,system,template,source,source_template,fold,model=key
        testids=[r['candidate_id'] for r in rows]
        assert len(testids)==len(set(testids))
        if analysis=='within_system_leave_position_out':
            position=int(rows[0]['heldout_mutation_position'])
            assert all(metadata[c]['mutation_positions']==[position] for c in testids)
            expected=sorted(c for c,m in metadata.items() if m['system']==system and m['analysis_set']=='primary_single_substitution' and position not in m['mutation_positions'])
            assert source==system and source_template==template
        elif analysis=='secondary_cross_system_transport':
            assert source!=system and source_template==PRIMARY[source] and template==PRIMARY[system]
            assert all(metadata[c]['system']==system and metadata[c]['analysis_set']=='primary_single_substitution' for c in testids)
            expected=sorted(c for c,m in metadata.items() if m['system']==source and m['analysis_set']=='primary_single_substitution')
        else:
            assert analysis=='secondary_multiple_substitution'
            assert all(metadata[c]['analysis_set']=='secondary_multiple_substitution' for c in testids)
            expected=sorted(c for c,m in metadata.items() if m['system']==system and m['analysis_set']=='primary_single_substitution')
        assert not set(expected).intersection(testids)
        for r in rows:
            assert sorted(json.loads(r['training_candidate_ids']))==expected
            assert int(r['n_train_candidates'])==len(expected)
            assert float(r['true_delta_pKD'])==labels[variant][r['candidate_id']]
        cols=FEATURES[model]
        if cols:
            # Independent Python summation/RMS, then augmented least squares
            # rather than the evaluator's normal-equation solve.
            x=np.array([[float(features[c,source_template][f]) for f in cols] for c in expected])
            y=np.array([labels[variant][c] for c in expected])
            scales=np.array([math.sqrt(sum(float(v)**2 for v in x[:,j])/len(x)) for j in range(len(cols))])
            scales[scales==0]=1.
            xs=x/scales
            augmented=np.vstack([xs/math.sqrt(len(x)),np.eye(len(cols))])
            response=np.concatenate([y/math.sqrt(len(x)),np.zeros(len(cols))])
            beta=np.linalg.lstsq(augmented,response,rcond=None)[0]
            train_fits+=1
            saved=coef_index[key]
            assert set(saved)==set(cols)
            for i,f in enumerate(cols):
                c=saved[f]
                scale_delta=abs(float(c['training_RMS_scale'])-scales[i])
                coefficient_delta=abs(float(c['coefficient_scaled'])-beta[i])
                max_scale_difference=max(max_scale_difference,scale_delta)
                max_coefficient_difference=max(max_coefficient_difference,coefficient_delta)
                assert scale_delta<1e-11 and coefficient_delta<1e-10
                assert float(c['lambda_mean_MSE'])==1.0 and float(c['fitted_intercept'])==0.
                assert sorted(json.loads(c['training_candidate_ids']))==expected
            for row in rows:
                cid=row['candidate_id']
                xt=np.array([float(features[cid,template][f]) for f in cols])
                prediction=float(np.sum((xt/scales)*beta))
                delta=abs(prediction-float(row['predicted_delta_pKD']))
                max_prediction_difference=max(max_prediction_difference,delta)
                assert delta<1e-10,(key,cid,delta)
                all_rebuilt.append(dict(row,independent_prediction=prediction))
        else:
            assert key not in coef_index
            for row in rows:
                assert float(row['predicted_delta_pKD'])==0.
                all_rebuilt.append(dict(row,independent_prediction=0.))
        prediction_checks.append(dict(analysis=analysis,label_variant=variant,system=system,template_pdb=template,
            source_system=source,fold_id=fold,model=model,n_train=len(expected),n_test=len(testids),
            candidate_and_position_or_system_separation_passed=True,training_only_scale_and_target_recalibration_passed=True))
    assert len(all_rebuilt)==2037
    mgroups=defaultdict(list)
    for r in all_rebuilt:mgroups[tuple(r[f] for f in fields)].append(r)
    metric_index={tuple(r[f] for f in fields):r for r in csvrows(RESULTS/'candidate_level_metrics.csv')}
    metric_checks=[];main_table=[]
    for key,rows in mgroups.items():
        saved=metric_index[key]
        y=[float(r['true_delta_pKD']) for r in rows];p=[r['independent_prediction'] for r in rows]
        mae=sum(abs(a-b) for a,b in zip(y,p))/len(y)
        assert abs(mae-float(saved['MAE_delta_pKD']))<1e-10
        if len(rows)<3 or len(set(y))==1 or len(set(float(r['predicted_delta_pKD']) for r in rows))==1:
            rho=None;assert saved['Spearman_rho']==''
        else:
            # Recompute rho from archived decimal predictions to make exact
            # numerical ties invariant to a different linear solver's rounding.
            rho=float(spearmanr(y,[float(r['predicted_delta_pKD']) for r in rows]).statistic)
            assert abs(rho-float(saved['Spearman_rho']))<1e-12
        k=min(3,len(rows))
        best_true=sorted(rows,key=lambda r:(-float(r['true_delta_pKD']),r['candidate_id']))[:k]
        best_pred=sorted(rows,key=lambda r:(-float(r['predicted_delta_pKD']),r['candidate_id']))[:k]
        overlap=len({r['candidate_id'] for r in best_true}&{r['candidate_id'] for r in best_pred})
        assert overlap==int(saved['top_k_overlap'])
        item=dict(zip(fields,key));item.update(n_candidates=len(rows),independent_MAE=mae,independent_Spearman=rho,independent_top_k_overlap=overlap,
            MAE_rho_topk_match=True)
        metric_checks.append(item)
        if item['analysis']=='within_system_leave_position_out' and item['label_variant']=='reported_hybrid' and item['template_pdb']==PRIMARY[item['system']]:
            main_table.append(item)
    assert len(main_table)==14
    # Each candidate must be genuinely held out once for each model/template/
    # label version, never counted once per residue or retained as a parent0.
    count=Counter((r['label_variant'],r['system'],r['template_pdb'],r['model'],r['candidate_id']) for r in all_rebuilt if r['analysis']=='within_system_leave_position_out')
    assert all(v==1 for v in count.values()) and len(count)==1365
    result=dict(status='passed_independent_numerical_and_leakage_review',reviewed_at_utc=datetime.now(timezone.utc).isoformat(),
        method='Independent program; Python RMS summation, augmented penalized least-squares (not the evaluator normal equations), SciPy Spearman, independently sorted top3 and separately derived position/system training membership.',
        evaluator_sha256=sha(HERE/'evaluate_screening.py'),protocol_sha256=sha(HERE/'fixed_analysis_protocol.json'),
        candidate_features_sha256=sha(ROOT/'feature_extraction/candidate_features.csv'),
        source_results_audit_sha256=sha(RESULTS/'evaluation_execution_audit.json'),
        n_predictions_independently_rebuilt=2037,n_nonzero_ridge_training_fits_independently_rebuilt=train_fits,
        n_metric_groups_independently_checked=len(metric_checks),n_primary_system_model_groups=14,
        max_prediction_absolute_difference=max_prediction_difference,max_training_scale_absolute_difference=max_scale_difference,
        max_scaled_coefficient_absolute_difference=max_coefficient_difference,
        all_actual_training_candidate_lists_derived_from_frozen_position_or_system_rules=True,
        all_scaling_training_only_and_lambda1_mean_loss_and_zero_intercept_verified=True,
        no_parent_or_double_mutant_used_for_training=True,no_target_system_calibration_in_transport=True,
        primary_candidates_each_held_out_once_per_model_template_variant=True,
        main_results=sorted(main_table,key=lambda r:(r['system'],r['model'])),metric_checks=metric_checks,fold_checks=prediction_checks,
        frozen_protocol_or_evaluator_modified=False,any_post_result_feature_or_rule_selection=False,neural_model_training_or_inference_run_by_reviewer=False)
    (HERE/'independent_evaluation_review.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf8')
    print(json.dumps({k:v for k,v in result.items() if k not in ['metric_checks','fold_checks','main_results']},ensure_ascii=False))

if __name__=='__main__':main()
