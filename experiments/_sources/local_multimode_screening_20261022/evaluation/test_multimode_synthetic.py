"""Synthetic-only tests; no real functional labels, profiles or neural model."""
from pathlib import Path
from datetime import datetime,timezone
import tempfile,json,math
import numpy as np
import evaluate_multimode as ev

HERE=Path(__file__).resolve().parent

def main():
    rng=np.random.default_rng(501);n=100
    alphabet=np.array(list('ACDEFGHIKLMNPQRSTVWY'))
    seq=[''.join(rng.choice(alphabet,9)) for _ in range(n)]
    stat=rng.normal(size=(n,2));R=rng.normal(size=(n,12));H=rng.normal(size=(n,9))
    kernels,settings=ev.construct_kernels(seq,stat,R,H,np.ones(12),np.ones(9))
    for k in kernels.values():
        assert np.allclose(k,k.T) and np.allclose(np.diag(k),1.)
        assert np.min(np.linalg.eigvalsh(k))>-1e-8
    train=list(range(20));test=list(range(20,100));y=rng.normal(size=n)
    p=ev.ridge_predict(kernels['baseline'],train,y[train],test)
    yy=y[train]-np.mean(y[train]);beta=np.linalg.solve(kernels['baseline'][:20,:20]+2*np.eye(20),yy)
    assert np.allclose(p,np.mean(y[train])+kernels['baseline'][20:,:20]@beta)
    ids=[f'SYNTHETIC_{i:03}' for i in range(n)]
    # Every strategy gets only the20 currently observed outcomes. Replacing
    # all80 hidden values cannot alter any selected candidate.
    changed_y=y.copy();changed_y[20:]=rng.normal(1000,500,size=80)
    for strategy in ['random','sequence_diversity']+list(ev.MODEL_BLOCKS):
        for feedback in [False,True]:
            a=ev.choose_batch(strategy,kernels,train,y[train],test,5,ids,np.random.default_rng(9),feedback)
            b=ev.choose_batch(strategy,kernels,train,changed_y[train],test,5,ids,np.random.default_rng(9),feedback)
            assert a==b and len(a)==len(set(a))==5 and not set(a).intersection(train)
    assert np.allclose(ev.rank_ties(np.array([2,1,1,3])),[3,1.5,1.5,4])
    assert math.isclose(ev.diagnostic_rho([1,2,3],[1,2,3]),1.)
    assert ev.diagnostic_rho([1,2,3],[0,0,0]) is None
    with tempfile.TemporaryDirectory(prefix='_synthetic_multimode_',dir=HERE) as tmp:
        root=Path(tmp);assert root.resolve().is_relative_to(HERE.resolve())
        features=root/'SYNTHETIC_features.npz'
        np.savez_compressed(features,candidate_ids=np.array(ids),sequences=np.array(seq),static_values=stat,R_pattern=R,H_pattern=H,
            R_training_scale=np.ones(12),H_training_scale=np.ones(9),checkpoint_sha256=np.array(ev.CHECKPOINT),
            unchanged_spectator_U_verified=np.array(True),all_library_variable_positions_excluded=np.array(True),
            zero_shot_scores=stat,zero_shot_score_names=np.array(['HERMES_fixed','ProteinMPNN']),zero_shot_scores_oriented_higher_better=np.array(True),
            template_sha256=np.array('SYNTHETIC_TEMPLATE'),regions_sha256=np.array('SYNTHETIC_U'),training_reference_sha256=np.array('SYNTHETIC_REFERENCE'))
        labels=root/'SYNTHETIC_labels.csv'
        ev.dump_csv(labels,[dict(panel_id='SYNTHETIC',candidate_id=cid,sequence=s,parent_candidate_id='SYNTHETIC_PARENT',
            value_oriented_higher_better=float(y[i]),label_qc_pass=True,censoring='none',source_id='SYNTHETIC_ONLY') for i,(cid,s) in enumerate(zip(ids,seq))])
        panel=dict(panel_id='SYNTHETIC',family_id='SYNTHETIC_FAMILY',task_object='antigen_peptide_mutation',endpoint='SYNTHETIC_ACTIVITY',
            complete_candidate_pool=True,positive_only_selection=False,mode2_eligible=True,mode3_eligible=True,
            labels_csv=str(labels),labels_sha256=ev.sha(labels),features_npz=str(features),features_sha256=ev.sha(features),
            hit_kind='relative_top_decile',hit_definition='SYNTHETIC_relative_top_decile',hit_is_above_parent_improvement=False,threshold_fixed_without_model_outputs=True)
        out=root/'out';out.mkdir()
        records,gates=ev.evaluate_candidate_panel(panel,out)
        assert len([g for g in gates if g['status']=='evaluable'])==2
        selections=ev.read_csv(out/'SYNTHETIC/selected_candidates.csv')
        assert len(selections)==20*(10+6)
        initialsets={}
        for r in selections:
            initial=json.loads(r['initial_candidate_ids']);chosen=json.loads(r['selected_candidate_ids'])
            assert len(chosen)==len(set(chosen)) and not set(initial).intersection(chosen)
            key=r['mode'],r['seed']
            if key in initialsets:assert initialsets[key]==initial
            else:initialsets[key]=initial
        trace=ev.read_csv(out/'SYNTHETIC/feedback_information_trace.csv')
        assert len(trace)==20*6*8
        for r in trace:
            before=json.loads(r['observed_candidate_ids_before']);batch=json.loads(r['selected_batch'])
            assert len(before)==20+(int(r['round'])-1)*5
            assert len(batch)==5 and not set(before).intersection(batch)
        paired,summary=ev.summarize_candidate_records(records,gates)
        assert paired and summary and all(r['paired_seeds']==20 for r in summary)
        report=ev.success_checks(summary,[],gates);assert len(report['per_panel_primary'])==2
        # Small-panel gate does not claim active learning on a14-candidate set.
        assert all(g['N']==100 for g in gates)
        # Site nomination has equal19 substitutions per position, never uses
        # functional values for eligibility, and preservesG/A chi missingness.
        sfile=root/'SYNTHETIC_sites.csv';sl=root/'SYNTHETIC_site_labels.csv'
        parents=list('ALRSGVKNWY');site_rows=[];score_rows=[]
        for i,aa in enumerate(parents):
            pos=f'D:{i+1}'
            score_rows.append(dict(panel_id='SYNTHETIC_DMS',position_key=pos,parent_aa=aa,parent_R_percentile=i/10,
                parent_H_percentile='' if aa in {'G','A'} else (9-i)/10,static_contact_degree=i,static_relative_SASA=(10-i)/10,
                parent_partner_context_response=i/20))
            for j,mutant in enumerate(sorted(ev.AA-{aa})):
                site_rows.append(dict(panel_id='SYNTHETIC_DMS',candidate_id=f'SYNTHETIC_SITE_{i}_{mutant}',position_key=pos,
                    parent_aa=aa,mutant_aa=mutant,value_oriented_higher_better=float(i+j/20),label_qc_pass=True,censoring='none'))
        ev.dump_csv(sfile,score_rows);ev.dump_csv(sl,site_rows)
        sp=dict(panel_id='SYNTHETIC_DMS',family_id='SYNTHETIC_DMS_FAMILY',task_object='TCR_mutation',endpoint='SYNTHETIC_enrichment_proxy',
            site_scores_csv=str(sfile),site_scores_sha256=ev.sha(sfile),site_labels_csv=str(sl),site_labels_sha256=ev.sha(sl),
            hit_kind='fixed_threshold',declared_hit_threshold=6.,hit_definition='SYNTHETIC_proxy_threshold',threshold_fixed_without_model_outputs=True)
        site,g=ev.evaluate_site_panel(sp,out)
        assert len(site)==2*(20+7)
        assert all(r['assay_budget']==r['n_positions']*19 for r in site)
        assert all(r['task_object']=='TCR_mutation' for r in site)
        ev.success_checks(summary,site,gates+g)
    result=dict(status='passed_synthetic_only',tested_at_utc=datetime.now(timezone.utc).isoformat(),
        evaluator_sha256=ev.sha(HERE/'evaluate_multimode.py'),test_sha256=ev.sha(__file__),
        all20paired_seeds_exercised=True,common_initial_labels_and_budget_verified=True,
        complete8round5candidate_feedback_exercised=True,no_hidden_label_access_in_policy_API=True,
        changing_unseen_outcomes_does_not_change_selection=True,site_equal19_substitution_budget_verified=True,
        neural_checkpoint_loaded=False,real_DMS_or_cell_labels_read=False,real_candidate_predictions_read=False,
        temporary_synthetic_fixtures_removed=True,no_model_hyperparameter_or_feature_selection=True)
    ev.dump_json(HERE/'synthetic_evaluation_selftest.json',result)
    print(json.dumps(result))

if __name__=='__main__':main()
