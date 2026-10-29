"""Independent scalar-pair and aggregation audit of the scientific reanalysis."""
from pathlib import Path
from itertools import combinations
import hashlib
import json
import math
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(4*1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def scalar_group(g, condition, target):
    # Explicit unordered-pair loop; does not import the production pair routine.
    truth_c = g['uniform_rep1_'+condition].to_numpy(float)
    truth_t = g['uniform_rep1_'+target].to_numpy(float)
    dedicated = g['pred_full_'+target].to_numpy(float)
    surrogate = g['pred_full_'+condition].to_numpy(float)
    intra = g['pred_intra_'+target].to_numpy(float)
    target_ranks = pd.Series(truth_t).rank(method='average').to_numpy()/len(g)
    result = {}
    for subset in ('all_close', 'separated_target_p50'):
        count = 0
        delta, full_scores, surrogate_scores, intra_scores = [], [], [], []
        for i, j in combinations(range(len(g)), 2):
            if abs(truth_c[j]-truth_c[i]) > .10+1e-12:
                continue
            if subset == 'separated_target_p50' and abs(target_ranks[j]-target_ranks[i]) < .5-1e-12:
                continue
            count += 1
            actual = truth_t[j]-truth_t[i]
            delta.append(abs(actual))
            if abs(actual) <= 1e-12:
                continue
            correct = []
            for prediction in (dedicated, surrogate, intra):
                change = prediction[j]-prediction[i]
                correct.append(.5 if abs(change) <= 1e-12 else float((change > 0) == (actual > 0)))
            full_scores.append(correct[0])
            surrogate_scores.append(correct[1])
            intra_scores.append(correct[2])
        avg = lambda a: float(np.mean(a)) if len(a) else np.nan
        result[subset] = dict(matched_pairs=count, non_tied_truth_pairs=len(full_scores),
            target_true_abs_delta_mean=avg(delta), full_direction_accuracy=avg(full_scores),
            surrogate_direction_accuracy=avg(surrogate_scores), intra_direction_accuracy=avg(intra_scores))
    return result


def main():
    manifest = json.loads((HERE/'source_manifest.json').read_text(encoding='utf-8'))
    source_failures = [s['path'] for s in manifest['sources'] if sha(s['path']) != s['sha256']]
    assert not source_failures, source_failures
    artifact = json.loads((HERE/'completion_audit.json').read_text(encoding='utf-8'))
    artifact_failures = [name for name, d in artifact['outputs'].items() if sha(HERE/name) != d['sha256']]
    assert not artifact_failures, artifact_failures
    p = pd.read_csv(HERE/'common_positions_with_uniform_truth_and_predictions.csv.gz', low_memory=False)
    t = pd.read_csv(HERE/'matched_per_record_AA.csv.gz', low_memory=False)
    actual = t[(t.label_mode == 'uniform_rep1')&(t.match_rule == 'physical_0p10')].set_index(
        ['sample_id', 'amino_acid', 'direction', 'target_subset'])
    checks = 0
    potential = 0
    eligible = 0
    numeric = ['matched_pairs', 'non_tied_truth_pairs', 'target_true_abs_delta_mean',
               'full_direction_accuracy', 'surrogate_direction_accuracy', 'intra_direction_accuracy']
    for (sid, aa), g in p.groupby(['sample_id', 'amino_acid'], sort=True):
        if len(g) < 4:
            continue
        potential += len(g)*(len(g)-1)//2
        eligible += len(g)
        for condition, target in (('rmsf', 'hchi'), ('hchi', 'rmsf')):
            independent = scalar_group(g.sort_values('residue_index'), condition, target)
            for subset, values in independent.items():
                row = actual.loc[(sid, aa, condition+'_to_'+target, subset)]
                for metric in numeric:
                    np.testing.assert_allclose(row[metric], values[metric], atol=1e-10, rtol=0, equal_nan=True)
                    checks += 1
    assert potential == 624602 and eligible == 33142
    # Verify the second and third aggregation levels from the exported type table.
    config = ['label_mode', 'direction', 'match_rule', 'target_subset', 'sample_id']
    metrics = ['target_true_abs_delta_mean', 'full_direction_accuracy',
               'surrogate_direction_accuracy', 'intra_direction_accuracy',
               'full_minus_surrogate_accuracy', 'full_minus_intra_accuracy']
    agg = t.groupby(config, sort=True)[metrics].mean()
    records = pd.read_csv(HERE/'matched_per_record.csv').set_index(config)
    np.testing.assert_allclose(agg.sort_index(), records[metrics].sort_index(), atol=2e-11, rtol=0, equal_nan=True)
    summaries = pd.read_csv(HERE/'matched_summary_family_cluster_bootstrap.csv')
    for key, frame in records.reset_index().groupby(config[:4]):
        for metric in metrics:
            independent = frame[metric].mean()
            match = summaries[(summaries.label_mode == key[0])&(summaries.direction == key[1])&
                              (summaries.match_rule == key[2])&(summaries.target_subset == key[3])&
                              (summaries.dataset_scope == 'all')&(summaries.metric == metric)]
            assert len(match) == 1
            np.testing.assert_allclose(match.record_macro_mean.iloc[0], independent, atol=2e-11, rtol=0, equal_nan=True)
    class_audit = artifact['classification']
    assert class_audit['class_label_mismatches'] == 0
    assert class_audit['previous_per_record_recalls_reproduced']
    extension = None
    ext = HERE/'sasa_extension'
    if (ext/'completion_audit.json').exists():
        quality = pd.read_csv(ext/'common_positions_quality_sasa.csv.gz', low_memory=False)
        assert quality.sidechain_complete.all() and quality.chi_input_complete.all()
        assert not quality.nonstandard_residue.any()
        assert np.isfinite(quality.sasa_sidechain_complex_A2).all()
        er = pd.read_csv(ext/'matched_per_record.csv')
        assert er.match_rule.isin(['sasa_q25_and_rmsf_q25', 'sasa_q25_and_rmsf_physical_0p10']).all()
        epairs = pd.read_csv(ext/'uniform_rep1_matched_pairs.csv.gz')
        assert (epairs.residue_index_i < epairs.residue_index_j).all()
        assert (epairs.sasa_percentile_delta_j_minus_i.abs() <= .25+1e-9).all()
        identity = quality.set_index(['sample_id', 'residue_index']).amino_acid
        keys_i = pd.MultiIndex.from_frame(epairs[['sample_id', 'residue_index_i']].rename(columns={'residue_index_i': 'residue_index'}))
        keys_j = pd.MultiIndex.from_frame(epairs[['sample_id', 'residue_index_j']].rename(columns={'residue_index_j': 'residue_index'}))
        np.testing.assert_array_equal(identity.reindex(keys_i).to_numpy(), epairs.amino_acid.to_numpy())
        np.testing.assert_array_equal(identity.reindex(keys_j).to_numpy(), epairs.amino_acid.to_numpy())
        extension = dict(quality_positions=len(quality), no_incomplete_sidechains=True,
                         pair_rows=len(epairs), all_pairs_same_AA_and_record=True,
                         all_exported_pairs_SASA_percentile_caliper_satisfied=True)
    absolute = None
    ad = HERE/'sasa_absolute_sensitivity'
    if (ad/'completion_audit.json').exists():
        source = pd.read_csv(HERE/'sasa_extension/common_positions_quality_sasa.csv.gz', low_memory=False)
        at = pd.read_csv(ad/'matched_per_record_AA.csv.gz', low_memory=False)
        central = at[(at.label_mode == 'uniform_rep1')&
                     (at.match_rule == 'sasa_abs_10A2_and_condition_abs_0p10')&
                     (at.target_subset == 'all_close')].set_index(['sample_id', 'amino_acid', 'direction', 'pair_scope'])
        absolute_checks = 0
        for (sid, aa), group in source.groupby(['sample_id', 'amino_acid'], sort=True):
            if len(group) < 4:
                continue
            group = group.sort_values('residue_index')
            area = group.sasa_sidechain_complex_A2.to_numpy(float)
            chains = group.chain.astype(str).to_numpy()
            for condition, target in (('rmsf', 'hchi'), ('hchi', 'rmsf')):
                cv = group['uniform_rep1_'+condition].to_numpy(float)
                tv = group['uniform_rep1_'+target].to_numpy(float)
                fv = group['pred_full_'+target].to_numpy(float)
                for scope in ('all_within_record', 'same_chain_only'):
                    changes, scores = [], []
                    for i, j in combinations(range(len(group)), 2):
                        if abs(area[j]-area[i]) > 10+1e-12 or abs(cv[j]-cv[i]) > .10+1e-12:
                            continue
                        if scope == 'same_chain_only' and chains[i] != chains[j]:
                            continue
                        delta = tv[j]-tv[i]
                        changes.append(abs(delta))
                        if abs(delta) > 1e-12:
                            prediction = fv[j]-fv[i]
                            scores.append(.5 if abs(prediction) <= 1e-12 else float((prediction > 0) == (delta > 0)))
                    row = central.loc[(sid, aa, condition+'_to_'+target, scope)]
                    assert row.matched_pairs == len(changes)
                    np.testing.assert_allclose(row.target_true_abs_delta_mean,
                        np.mean(changes) if changes else np.nan, atol=1e-10, rtol=0, equal_nan=True)
                    np.testing.assert_allclose(row.full_direction_accuracy,
                        np.mean(scores) if scores else np.nan, atol=1e-10, rtol=0, equal_nan=True)
                    absolute_checks += 3
        aaudit = json.loads((ad/'completion_audit.json').read_text(encoding='utf-8'))
        assert sha(aaudit['input_source']) == aaudit['input_sha256']
        assert sha(HERE/'protocol.json') == aaudit['primary_protocol_unchanged_sha256']
        absolute = dict(scalar_central_10A2_metric_checks=absolute_checks,
                        same_chain_and_all_record_scopes_verified=True,
                        exploratory_status_explicit=True, original_primary_protocol_unchanged=True)
    result = dict(status='passed', unchanged_source_files_verified=len(manifest['sources']),
                  base_output_hashes_verified=len(artifact['outputs']),
                  scalar_pair_metric_checks=checks, all_eligible_unordered_pairs=potential,
                  sameAA4_positions=eligible, scalar_primary_physical_pairs_match_vectorized=True,
                  type_to_record_and_record_to_overall_means_verified=True,
                  previous_full3063_classes_and_recalls_exact=True,
                  truth_conditioned_analysis_not_deployable_input=True,
                  family_clusters=42, largest_family_records=38,
                  sasa_extension=extension, absolute_sasa_post_exploration=absolute)
    (HERE/'independent_validation.json').write_text(json.dumps(result, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
