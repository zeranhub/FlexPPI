"""Explicit post-exploration robustness checks: absolute SASA and same-chain pairs.

Added after observing the rank-caliper SASA balance, not a replacement primary
analysis and not a retrospective preregistration. Thresholds fixed before this run.
"""
from pathlib import Path
import importlib.util
import gzip
import json
from datetime import datetime, timezone
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location('base_matching', HERE/'analyze_matched_complementarity.py')
BASE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(BASE)
OUT = HERE/'sasa_absolute_sensitivity'
SOURCE = HERE/'sasa_extension/common_positions_quality_sasa.csv.gz'
METRICS = BASE.METRICS+('sasa_true_abs_delta_mean', 'cross_chain_pair_fraction')
CONFIGS = ['label_mode', 'direction', 'match_rule', 'target_subset', 'pair_scope',
           'sample_id', 'dataset', 'fold', 'family_id']
RULE = {
    'version': 'absolute_sasa_robustness_post_exploration_v1',
    'timing': 'Exploratory sensitivity added after the main rank-SASA analysis showed a 14.00 A2 mean physical SASA difference. Main rank-SASA results are retained unchanged. User/root explicitly authorized symmetric absolute-area checks.',
    'eligibility': 'Same original record and AA; at least 4 quality-qualified common positions per record/AA, every unordered i<j pair. All positions in the source already satisfy finite SASA, complete input sidechain/chi atoms, and standard-residue checks.',
    'calipers': 'Absolute complex sidechain SASA difference <=10 A2 central sensitivity, <=5 and <=20 A2 additional sensitivities, each combined with absolute truth RMSF difference <=0.10 A for Hchi target, or absolute truth Hchi difference <=0.10 nats for RMSF target.',
    'pair_scopes': 'all_within_record and same_chain_only; same-chain subset retains the record/AA eligibility and target-rank definition, with pair chain identifiers required equal. A record is an assembly, not necessarily one protein chain.',
    'target_strata': 'all_close plus separated_target_p50 (target average-rank/n percentile difference >=0.50 within the original quality-qualified record/AA group). No prediction controls inclusion.',
    'aggregation_uncertainty': 'Same pair -> AA -> record equal-weight aggregation and 10000 whole-family-cluster percentile bootstrap as the original protocol; undefined records are retained. Family-equal and leave-largest-family-out summaries are descriptive sensitivities.',
    'truth_modes': list(BASE.MODES),
    'new_claim_limit': 'Robustness of complementary outputs, not evidence of causal effects, partner-induced changes, full thermodynamic entropy, or family-held-out generalization.',
}


def statistics(g, mode):
    n = len(g)
    i, j = np.triu_indices(n, 1)
    values = {k: g[mode+'_'+k].to_numpy(float) for k in ('rmsf', 'hchi')}
    td = {k: a[j]-a[i] for k, a in values.items()}
    rank = {k: pd.Series(a).rank(method='average').to_numpy()/n for k, a in values.items()}
    rd = {k: a[j]-a[i] for k, a in rank.items()}
    predictions = {k: g['pred_'+k].to_numpy(float) for k in ('full_rmsf', 'full_hchi', 'intra_rmsf', 'intra_hchi')}
    pdiff = {k: a[j]-a[i] for k, a in predictions.items()}
    sas = g.sasa_sidechain_complex_A2.to_numpy(float)
    ds = np.abs(sas[j]-sas[i])
    chains = g.chain.astype(str).to_numpy()
    same_chain = chains[i] == chains[j]
    meta = dict(sample_id=g.sample_id.iloc[0], dataset=g.dataset.iloc[0],
                fold=int(g.fold.iloc[0]), family_id=g.family_id.iloc[0],
                amino_acid=g.amino_acid.iloc[0], label_mode=mode)
    rows = []
    for condition, target in (('rmsf', 'hchi'), ('hchi', 'rmsf')):
        close = np.abs(td[condition]) <= .10+BASE.EPS
        for area in (5, 10, 20):
            area_close = ds <= area+BASE.EPS
            for scope in ('all_within_record', 'same_chain_only'):
                chain_ok = same_chain if scope == 'same_chain_only' else np.ones(len(i), bool)
                for subset in ('all_close', 'separated_target_p50'):
                    keep = close & area_close & chain_ok
                    if subset == 'separated_target_p50':
                        keep &= np.abs(rd[target]) >= .50-BASE.EPS
                    count = int(keep.sum())
                    row = meta | dict(direction=condition+'_to_'+target,
                        match_rule=f'sasa_abs_{area}A2_and_condition_abs_0p10',
                        target_subset=subset, pair_scope=scope,
                        eligible_positions=n, potential_pairs=len(i), matched_pairs=count,
                        matched_positions=len(np.unique(np.r_[i[keep], j[keep]])),
                        non_tied_truth_pairs=0, full_pred_tied_pairs=0)
                    for metric in METRICS:
                        row[metric] = np.nan
                    if count:
                        actual = td[target][keep]
                        row['target_true_abs_delta_mean'] = float(np.abs(actual).mean())
                        row['target_true_abs_delta_median'] = float(np.median(np.abs(actual)))
                        row['target_percentile_abs_delta_mean'] = float(np.abs(rd[target][keep]).mean())
                        row['condition_true_abs_delta_mean'] = float(np.abs(td[condition][keep]).mean())
                        row['condition_percentile_abs_delta_mean'] = float(np.abs(rd[condition][keep]).mean())
                        row['sasa_true_abs_delta_mean'] = float(ds[keep].mean())
                        row['cross_chain_pair_fraction'] = float((~same_chain[keep]).mean())
                        for model in ('full', 'intra'):
                            row[model+'_difference_mae'] = float(np.abs(pdiff[model+'_'+target][keep]-actual).mean())
                        row['full_minus_intra_difference_mae'] = row['full_difference_mae']-row['intra_difference_mae']
                        row['full_pred_abs_delta_mean'] = float(np.abs(pdiff['full_'+target][keep]).mean())
                        nt = np.abs(actual) > BASE.EPS
                        row['non_tied_truth_pairs'] = int(nt.sum())
                        if nt.any():
                            truth_sign = np.sign(actual[nt])
                            for model, key in (('full', 'full_'+target), ('surrogate', 'full_'+condition), ('intra', 'intra_'+target)):
                                change = pdiff[key][keep][nt]
                                score = np.where(np.abs(change) <= BASE.EPS, .5, (np.sign(change) == truth_sign).astype(float))
                                row[model+'_direction_accuracy'] = float(score.mean())
                            row['full_pred_tied_pairs'] = int((np.abs(pdiff['full_'+target][keep][nt]) <= BASE.EPS).sum())
                            row['full_minus_surrogate_accuracy'] = row['full_direction_accuracy']-row['surrogate_direction_accuracy']
                            row['full_minus_intra_accuracy'] = row['full_direction_accuracy']-row['intra_direction_accuracy']
                    rows.append(row)
    pairs = None
    if mode == 'uniform_rep1':
        # Includes every pair eligible for any absolute-area rule/direction.
        union = (ds <= 20+BASE.EPS)&((np.abs(td['rmsf']) <= .10+BASE.EPS)|(np.abs(td['hchi']) <= .10+BASE.EPS))
        pairs = pd.DataFrame(dict(sample_id=meta['sample_id'], dataset=meta['dataset'], amino_acid=meta['amino_acid'],
            residue_index_i=g.residue_index.to_numpy(int)[i[union]], residue_index_j=g.residue_index.to_numpy(int)[j[union]],
            same_chain=same_chain[union], sasa_abs_delta_A2=ds[union],
            **{k+'_truth_delta_j_minus_i': v[union] for k, v in td.items()},
            **{k+'_truth_percentile_delta_j_minus_i': v[union] for k, v in rd.items()},
            **{k+'_pred_delta_j_minus_i': v[union] for k, v in pdiff.items()}))
    return rows, pairs


def main():
    OUT.mkdir(exist_ok=True)
    path = OUT/'protocol_post_exploration.json'
    if path.exists():
        assert json.loads(path.read_text(encoding='utf-8'))['fixed_rules'] == RULE
    else:
        BASE.write_json(path, dict(created_utc=datetime.now(timezone.utc).isoformat(), fixed_rules=RULE,
                                  code_sha256_at_freeze=BASE.sha(__file__),
                                  unchanged_primary_protocol_sha256=BASE.sha(HERE/'protocol.json')))
    p = pd.read_csv(SOURCE, low_memory=False, dtype={'chain': str, 'author_residue_id': str})
    assert p.sidechain_complete.all() and p.chi_input_complete.all() and not p.nonstandard_residue.any()
    rows = []
    with gzip.open(OUT/'uniform_rep1_absolute_sasa_matched_pairs.csv.gz', 'wt', encoding='utf-8', newline='') as handle:
        first = True
        for mode in BASE.MODES:
            for _, g in p.groupby(['sample_id', 'amino_acid'], sort=True):
                if len(g) < 4:
                    continue
                r, pair = statistics(g.sort_values('residue_index'), mode)
                rows.extend(r)
                if pair is not None and len(pair):
                    pair.to_csv(handle, index=False, header=first, float_format='%.10g')
                    first = False
    per_type = pd.DataFrame(rows)
    BASE.csv(per_type, OUT/'matched_per_record_AA.csv.gz')
    records = []
    for key, g in per_type.groupby(CONFIGS, sort=True):
        row = dict(zip(CONFIGS, key))
        row['eligible_types'] = len(g)
        row['matched_types'] = int((g.matched_pairs > 0).sum())
        row['ordering_types'] = int(g.full_direction_accuracy.notna().sum())
        for col in BASE.COUNT_COLS:
            row[col] = int(g[col].sum())
        for col in METRICS:
            row[col] = BASE.mean(g[col])
        records.append(row)
    per_record = pd.DataFrame(records)
    BASE.csv(per_record, OUT/'matched_per_record.csv')
    BASE.csv(per_record[per_record.matched_types == 0], OUT/'no_matched_pairs_records.csv')
    BASE.csv(per_record[(per_record.full_direction_accuracy <= .5)|
                        (per_record.full_minus_surrogate_accuracy <= 0)|
                        per_record.full_direction_accuracy.isna()], OUT/'weak_or_negative_records.csv')
    summary = BASE.summarize(per_record, CONFIGS[:5], METRICS)
    BASE.csv(summary, OUT/'summary_family_cluster_bootstrap.csv')
    central = summary[(summary.label_mode == 'uniform_rep1')&(summary.dataset_scope == 'all')&
                      (summary.target_subset == 'all_close')&
                      (summary.match_rule == 'sasa_abs_10A2_and_condition_abs_0p10')]
    BASE.csv(central, OUT/'central_10A2_uniform_rep1_summary.csv')
    BASE.write_json(OUT/'completion_audit.json', dict(status='complete', exploratory_post_primary=True,
        input_source=str(SOURCE), input_sha256=BASE.sha(SOURCE), executed_code_sha256=BASE.sha(__file__),
        protocol_sha256=BASE.sha(path), primary_protocol_unchanged_sha256=BASE.sha(HERE/'protocol.json'),
        per_record_rows=len(per_record), per_type_rows=len(per_type), summary_rows=len(summary),
        original_sources_modified=False, word_modified=False, new_training=False, new_MD=False))
    print(central[central.metric.isin(['target_true_abs_delta_mean', 'full_direction_accuracy',
        'surrogate_direction_accuracy', 'full_minus_surrogate_accuracy', 'full_minus_intra_accuracy'])][
        ['direction', 'pair_scope', 'metric', 'record_macro_mean', 'family_cluster_ci95_low',
         'family_cluster_ci95_high', 'defined_records', 'matched_pairs']].to_string(index=False), flush=True)


if __name__ == '__main__':
    main()
