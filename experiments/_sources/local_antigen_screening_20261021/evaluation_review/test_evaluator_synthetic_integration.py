"""Exercise all orchestration using synthetic labels/features only."""
from argparse import Namespace
from collections import Counter
from datetime import datetime,timezone
from pathlib import Path
import csv,hashlib,json,tempfile
import evaluate_screening as ev

HERE=Path(__file__).resolve().parent

def main():
    ev.self_test()
    protocol=ev.read_json(HERE/'fixed_analysis_protocol.json')
    split=ev.read_json(HERE/'candidate_splits.json')
    cases=split['candidates']
    modelsha=protocol['base_model']['sha256']
    with tempfile.TemporaryDirectory(prefix='_synthetic_evaluator_',dir=HERE) as temp:
        work=Path(temp)
        assert work.resolve().is_relative_to(HERE.resolve())
        feature_rows=[];byid={}
        for row in cases:
            cid=row['candidate_id'];byid[cid]=row
            pos=row['mutation_positions'];n=len(pos)
            changes=sum(ord(row['sequence'][p-1])-ord(row['parent_sequence'][p-1]) for p in pos)
            for tindex,template in enumerate(ev.TEMPLATES[row['system']]):
                factor=1+.1*tindex
                feature_rows.append(dict(candidate_id=cid,panel_id=row['system'],template_pdb=template,
                    sequence_n_mutations=n,BLOSUM62_deficit=n*2+sum(pos)/10,
                    delta_TCR_pMHC_BSA_A2=changes*.05,
                    delta_peptide_TCR_contact_residue_pairs_4p5A=(sum(pos)%3-1)*n,
                    R=.005*(sum(pos)+changes*.01)*factor,H=.02*changes*factor,
                    checkpoint_sha256=modelsha,template_is_primary=template==ev.PRIMARY[row['system']],
                    preparation_qc_pass=True,dynamic_qc_pass=True,static_eligible=True,
                    n_U_RMSF=10,n_U_Hchi=8,feature_region_id=template+'_SYNTHETIC_ONLY'))
        feature_path=work/'SYNTHETIC_NOT_REAL_candidate_features.csv'
        ev.write_csv(feature_path,feature_rows)
        variant_info={}
        for index,variant in enumerate(['reported_hybrid','all_constrained','all_fitted']):
            rows=[]
            for row in cases:
                pos=row['mutation_positions'];n=len(pos)
                truth=.7*n+.01*sum(pos)+.02*index*n
                study={'SYNTHETIC_study_A':truth}
                if row['system']=='nyeso':study={'SYNTHETIC_study_A':truth-.1*n,'SYNTHETIC_study_B':truth+.1*n}
                rows.append(dict(candidate_id=row['candidate_id'],panel_id=row['system'],
                    tcr='1G4' if row['system']=='nyeso' else 'A6',peptide_sequence=row['sequence'],
                    panel_reference_peptide=row['parent_sequence'],delta_pKD_vs_parent=truth,
                    censoring='none_reported',study_ids='|'.join(study),study_delta_pKD_json=json.dumps(study)))
            labels=work/f'SYNTHETIC_{variant}.csv';observations=work/f'SYNTHETIC_{variant}_observations.csv'
            ev.write_csv(labels,rows);ev.write_csv(observations,[{'KD_relation':'='}])
            variant_info[variant]=dict(candidate_file=str(labels),candidate_file_sha256=ev.digest(labels),
                observation_file=str(observations),primary_file_sha256=ev.digest(observations))
        rulepath=work/'SYNTHETIC_label_rule.json'
        ev.write_json(rulepath,dict(primary='reported_hybrid',frozen_before_model_results=True,variants=variant_info))
        checks=[]
        for mode,expected_rows in [('full',2037),('restricted',1455)]:
            output=work/('SYNTHETIC_outputs_'+mode)
            args=Namespace(features=str(feature_path),feature_manifest=None,label_rule=str(rulepath),
                model_checkpoint=protocol['base_model']['local_run_path'],mode=mode,output=str(output))
            ev.run_analysis(args)
            audit=ev.read_json(output/'evaluation_execution_audit.json')
            assert audit['prediction_rows']==expected_rows
            predictions=ev.read_csv(output/'candidate_predictions.csv')
            modes=Counter(r['analysis'] for r in predictions)
            assert modes['within_system_leave_position_out']==(1365 if mode=='full' else 975)
            assert modes['secondary_cross_system_transport']==(546 if mode=='full' else 390)
            assert modes['secondary_multiple_substitution']==(126 if mode=='full' else 90)
            for r in predictions:
                train=json.loads(r['training_candidate_ids'])
                assert r['candidate_id'] not in train
                if r['analysis']=='within_system_leave_position_out':
                    held=int(r['heldout_mutation_position'])
                    assert all(held not in byid[cid]['mutation_positions'] for cid in train)
                elif r['analysis']=='secondary_cross_system_transport':
                    assert all(byid[cid]['system']!=byid[r['candidate_id']]['system'] for cid in train)
                else:
                    assert byid[r['candidate_id']]['analysis_set']=='secondary_multiple_substitution'
                    assert all(byid[cid]['analysis_set']=='primary_single_substitution' for cid in train)
            assert ev.read_csv(output/'study_stratified_metrics.csv')
            assert audit['static_increment_claim_eligible']==(mode=='full')
            checks.append(dict(mode=mode,status='passed',prediction_rows=expected_rows,analysis_row_counts=dict(modes),
                all_candidate_and_position_and_system_exclusions_verified=True,
                parent_and_double_mutants_never_enter_fit=True,study_stratification_outputs_checked=True))
    result=dict(status='passed_synthetic_integration',tested_at_utc=datetime.now(timezone.utc).isoformat(),
        evaluator_path=str(HERE/'evaluate_screening.py'),evaluator_sha256=ev.digest(HERE/'evaluate_screening.py'),
        real_SPR_numeric_values_read=False,real_candidate_features_or_predictions_read=False,
        actual_neural_checkpoint_loaded=False,neural_inference_performed=False,neural_training_performed=False,
        complete_fixed_pipeline_tested_with_synthetic_labels_and_features=True,
        synthetic_fixtures_removed_after_test=True,checks=checks)
    ev.write_json(HERE/'evaluator_selftest_audit.json',result)
    print(json.dumps(result,ensure_ascii=False))

if __name__=='__main__':main()
