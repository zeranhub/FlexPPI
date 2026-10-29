"""Fixed candidate-level retrospective evaluation, without model inference.

This program only reads already prepared candidate features and curated SPR
labels. It never loads a neural network or updates base-model weights. Use the
pvb_cpu Python (NumPy required). Run --self-test without real data first.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from datetime import datetime, timezone
import csv
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
SEQUENCE = ['sequence_n_mutations', 'BLOSUM62_deficit']
STATIC = ['delta_TCR_pMHC_BSA_A2', 'delta_peptide_TCR_contact_residue_pairs_4p5A']
DYNAMIC = ['R', 'H']
FULL_MODELS = {
    'zero_change': [],
    'sequence_only': SEQUENCE,
    'static_only': STATIC,
    'sequence_plus_static': SEQUENCE + STATIC,
    'sequence_plus_static_plus_R': SEQUENCE + STATIC + ['R'],
    'sequence_plus_static_plus_H': SEQUENCE + STATIC + ['H'],
    'sequence_plus_static_plus_RH': SEQUENCE + STATIC + DYNAMIC,
}
RESTRICTED_MODELS = {
    'zero_change': [],
    'sequence_only': SEQUENCE,
    'sequence_plus_R': SEQUENCE + ['R'],
    'sequence_plus_H': SEQUENCE + ['H'],
    'sequence_plus_RH': SEQUENCE + DYNAMIC,
}
TEMPLATES = {'nyeso': ['2BNQ', '2BNR'], 'tax': ['1AO7', '1QSE', '1QSF']}
PRIMARY = {'nyeso': '2BNQ', 'tax': '1AO7'}
LAMBDA = 1.0
UNCENSORED = {'', 'none', 'none_reported', 'uncensored', 'not_censored', 'false', '0'}


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def write_json(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False), encoding='utf-8')


def read_csv(path):
    with Path(path).open(encoding='utf-8-sig', newline='') as handle:
        return list(csv.DictReader(handle))


def write_csv(path, rows):
    rows = list(rows)
    fields = list(dict.fromkeys(k for r in rows for k in r))
    with Path(path).open('w', encoding='utf-8-sig', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: json.dumps(v, ensure_ascii=False) if isinstance(v, (list, dict)) else v for k, v in row.items()})


def flag(value, name):
    text = str(value).strip().lower()
    if text in {'true', '1', 'yes', 'passed', 'pass'}:
        return True
    if text in {'false', '0', 'no', 'failed', 'fail'}:
        return False
    raise ValueError(f'Unrecognized Boolean {name}: {value!r}')


def finite_number(value):
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def rank_average_ties(values):
    values = np.asarray(values, dtype=np.float64)
    order = np.argsort(values, kind='stable')
    ranks = np.empty(len(values), dtype=np.float64)
    start = 0
    while start < len(values):
        stop = start + 1
        while stop < len(values) and values[order[stop]] == values[order[start]]:
            stop += 1
        ranks[order[start:stop]] = (start + stop - 1) / 2.0 + 1.0
        start = stop
    return ranks


def candidate_metrics(rows):
    ids = [r['candidate_id'] for r in rows]
    if len(ids) != len(set(ids)):
        raise ValueError('A candidate contributes more than once to a metric group')
    observed = np.asarray([r['true_delta_pKD'] for r in rows], dtype=np.float64)
    predicted = np.asarray([r['predicted_delta_pKD'] for r in rows], dtype=np.float64)
    n = len(rows)
    mae = float(np.mean(np.abs(predicted - observed))) if n else None
    zero = float(np.mean(np.abs(observed))) if n else None
    rho = None
    rho_reason = 'fewer_than_three_candidates'
    if n >= 3:
        ry, rp = rank_average_ties(observed), rank_average_ties(predicted)
        if np.std(ry) > 0 and np.std(rp) > 0:
            rho = float(np.corrcoef(ry, rp)[0, 1])
            rho_reason = ''
        else:
            rho_reason = 'constant_observed_or_predicted_ranking'
    comparable, credit, predicted_ties = 0, 0.0, 0
    for i in range(n):
        for j in range(i + 1, n):
            dy = observed[i] - observed[j]
            if abs(dy) <= 1e-12:
                continue
            comparable += 1
            dp = predicted[i] - predicted[j]
            if abs(dp) <= 1e-12:
                credit += 0.5
                predicted_ties += 1
            elif dy * dp > 0:
                credit += 1.0
    k = min(3, n)
    # Ties are always resolved lexically by biological candidate ID, never
    # chosen to maximize overlap. Cutoff ties are explicitly reported.
    truth_order = sorted(range(n), key=lambda i: (-observed[i], ids[i]))
    prediction_order = sorted(range(n), key=lambda i: (-predicted[i], ids[i]))
    top_truth = [ids[i] for i in truth_order[:k]]
    top_predicted = [ids[i] for i in prediction_order[:k]]
    truth_cutoff_tie = n > k and observed[truth_order[k - 1]] == observed[truth_order[k]]
    pred_cutoff_tie = n > k and predicted[prediction_order[k - 1]] == predicted[prediction_order[k]]
    return dict(n_candidates=n, MAE_delta_pKD=mae, zero_change_MAE_delta_pKD=zero,
                gain_fraction_vs_zero=(1 - mae / zero) if zero and mae is not None else None,
                Spearman_rho=rho, Spearman_unavailable_reason=rho_reason,
                candidate_pair_concordance=credit / comparable if comparable else None,
                comparable_candidate_pairs=comparable, predicted_tied_pairs=predicted_ties,
                top_k=k, top_k_overlap=len(set(top_truth) & set(top_predicted)),
                top_k_overlap_fraction=len(set(top_truth) & set(top_predicted)) / k if k else None,
                top_true_candidate_ids=top_truth, top_predicted_candidate_ids=top_predicted,
                top_true_cutoff_tie=bool(truth_cutoff_tie), top_predicted_cutoff_tie=bool(pred_cutoff_tie),
                independent_TCR_systems_in_this_group=1,
                formal_significance_test_performed=False)


def fit_fixed_ridge(x, y):
    """Source-only RMS scaling; solve mean-MSE ridge with lambda exactly1."""
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    if x.ndim != 2 or y.ndim != 1 or len(x) != len(y) or not len(y):
        raise ValueError('Invalid training arrays')
    if not np.isfinite(x).all() or not np.isfinite(y).all():
        raise ValueError('Non-finite training values')
    scale = np.sqrt(np.mean(x ** 2, axis=0))
    scale = np.where(scale == 0, 1.0, scale)
    xs = x / scale
    coefficient = np.linalg.solve(xs.T @ xs + len(y) * LAMBDA * np.eye(x.shape[1]), xs.T @ y)
    return coefficient, scale


def apply_fixed_ridge(x, coefficient, scale):
    return (np.asarray(x, dtype=np.float64) / scale) @ coefficient


def map_candidates(split, label_rows):
    definitions = {(r['system'], r['sequence']): r for r in split['candidates']}
    mapped = {}
    protocol_to_label = {}
    for row in label_rows:
        system = 'nyeso' if row['tcr'] == '1G4' else 'tax' if row['tcr'] == 'A6' else None
        if system is None:
            raise ValueError('Unexpected TCR: ' + row['tcr'])
        definition = definitions[system, row['peptide_sequence']]
        if row['panel_reference_peptide'] != definition['parent_sequence']:
            raise ValueError('SPR reference peptide differs from frozen sequence reference')
        cid = row['candidate_id']
        if cid in mapped:
            raise ValueError('Duplicate biological label candidate: ' + cid)
        mapped[cid] = dict(definition, candidate_id=cid, protocol_candidate_id=definition['candidate_id'], panel_id=row['panel_id'])
        protocol_to_label[definition['candidate_id']] = cid
    if len(mapped) != len(definitions):
        raise ValueError('Label candidates do not equal frozen candidate scope')
    return mapped, protocol_to_label


def load_labels(label_rule_path, split):
    rule = read_json(label_rule_path)
    if rule['primary'] != 'reported_hybrid' or not rule['frozen_before_model_results']:
        raise ValueError('Primary SPR fit rule does not match the frozen protocol')
    variants = {}
    source_files = []
    common_candidates = None
    for version, info in rule['variants'].items():
        path = Path(info['candidate_file'])
        if digest(path) != info['candidate_file_sha256']:
            raise ValueError('SPR candidate file changed since label freeze: ' + version)
        observation_path = Path(info['observation_file'])
        if digest(observation_path) != info['primary_file_sha256']:
            raise ValueError('SPR observation file changed since label freeze: ' + version)
        rows = read_csv(path)
        mapped, crosswalk = map_candidates(split, rows)
        scope = {(k, v['system'], v['sequence']) for k, v in mapped.items()}
        if common_candidates is None:
            common_candidates = scope
        elif scope != common_candidates:
            raise ValueError('SPR sensitivity variants have different candidate identities')
        byid = {}
        for row in rows:
            cid = row['candidate_id']
            censoring = str(row.get('censoring', '')).strip().lower()
            reasons = []
            value = finite_number(row.get('delta_pKD_vs_parent'))
            if censoring not in UNCENSORED:
                reasons.append('interval_censored_or_undefined_affinity')
            if value is None:
                reasons.append('nonfinite_delta_pKD')
            byid[cid] = dict(candidate_id=cid, delta_pKD=value, reasons=reasons,
                            censoring=row.get('censoring'), study_ids=row.get('study_ids', ''),
                            study_delta_pKD_json=row.get('study_delta_pKD_json', ''),
                            row_source=str(path))
        # The exact-equality relation is an independent source-level censoring
        # check; none_reported never means measurement uncertainty is absent.
        for observation in read_csv(observation_path):
            relation = str(observation.get('KD_relation', '=')).strip()
            if relation not in {'=', '=='}:
                raise ValueError('New censored observations require a versioned eligibility amendment, not midpoint substitution')
        variants[version] = byid
        source_files.extend([dict(role='SPR_candidates', variant=version, path=str(path), sha256=digest(path)),
                             dict(role='SPR_observations', variant=version, path=str(observation_path), sha256=digest(observation_path))])
    return rule, variants, mapped, crosswalk, source_files


def load_features(path, candidate_info, expected_checkpoint, mode, feature_manifest=None):
    rows = read_csv(path)
    if not rows:
        raise ValueError('Empty candidate feature table')
    required = {'candidate_id', 'panel_id', 'template_pdb'} | set(SEQUENCE + DYNAMIC)
    if mode == 'full':
        required |= set(STATIC)
    missing = required - set(rows[0])
    if missing:
        raise ValueError('Missing feature columns: ' + ', '.join(sorted(missing)))
    manifest_hash = None
    if feature_manifest:
        manifest = read_json(feature_manifest)
        manifest_hash = manifest.get('checkpoint_sha256', manifest.get('model_checkpoint_sha256'))
    parsed = {}
    checkpoint_rows = set()
    exclusions = []
    counts = defaultdict(set)
    for row in rows:
        cid = row['candidate_id']
        if cid not in candidate_info:
            raise ValueError('Feature candidate ID is not a biological SPR candidate: ' + cid)
        info = candidate_info[cid]
        template = row['template_pdb'].upper()
        if template not in TEMPLATES[info['system']]:
            raise ValueError('Unexpected template for candidate: ' + cid + '/' + template)
        if row['panel_id'] != info['panel_id']:
            raise ValueError('Feature and label panel identities disagree')
        key = cid, template
        if key in parsed:
            raise ValueError('Duplicate candidate/template row: ' + str(key))
        this_hash = row.get('checkpoint_sha256') or row.get('model_checkpoint_sha256') or manifest_hash
        if this_hash != expected_checkpoint:
            raise ValueError('Feature checkpoint SHA missing or differs from the fixed active model')
        checkpoint_rows.add(this_hash)
        primary_flag_name = 'template_is_primary' if 'template_is_primary' in row else 'is_primary_template'
        if primary_flag_name not in row or flag(row[primary_flag_name], primary_flag_name) != (template == PRIMARY[info['system']]):
            raise ValueError('Primary-template flag disagrees with fixed2BNQ/1AO7 choice')
        preparation_key = 'preparation_qc_pass' if 'preparation_qc_pass' in row else 'qc_pass'
        if preparation_key not in row:
            raise ValueError('Preparation QC flag is required')
        preparation_ok = flag(row[preparation_key], preparation_key)
        dynamic_ok = flag(row.get('dynamic_qc_pass', preparation_ok), 'dynamic_qc_pass')
        static_ok = flag(row.get('static_eligible', preparation_ok), 'static_eligible')
        needed = SEQUENCE + DYNAMIC + (STATIC if mode == 'full' else [])
        values = {f: finite_number(row.get(f)) for f in SEQUENCE + STATIC + DYNAMIC}
        reasons = []
        if not preparation_ok:
            reasons.append('preparation_qc_failed')
        if not dynamic_ok:
            reasons.append('dynamic_qc_failed')
        if mode == 'full' and not static_ok:
            reasons.append('static_baseline_ineligible')
        reasons.extend('nonfinite_' + f for f in needed if values[f] is None)
        if values['sequence_n_mutations'] != info['mutation_count']:
            raise ValueError('Mutation count disagrees with the frozen biological parent')
        region_counts = []
        for name in ['n_U_RMSF', 'n_U_Hchi']:
            if name in row:
                count = finite_number(row[name])
                if count is None or int(count) != count or count <= 0:
                    reasons.append('invalid_' + name)
                else:
                    counts[info['system'], template, name].add(int(count))
                region_counts.append(count)
        if 'feature_region_id' in row:
            counts[info['system'], template, 'feature_region_id'].add(row['feature_region_id'])
        if info['analysis_set'] == 'reference_only':
            bad_zero = [f for f in needed if values[f] is not None and abs(values[f]) > 1e-8]
            if bad_zero:
                raise ValueError('Parent features are not candidate-minus-the-same-template-parent zeros: ' + str(key) + str(bad_zero))
        parsed[key] = dict(candidate_id=cid, template_pdb=template, system=info['system'], panel_id=info['panel_id'],
                           values=values, feature_reasons=reasons, input_id=row.get('input_id'), region_counts=region_counts)
        exclusions.append(dict(candidate_id=cid, panel_id=info['panel_id'], system=info['system'], template_pdb=template,
                               analysis_set=info['analysis_set'], mutation_positions=info['mutation_positions'],
                               feature_eligible=not reasons, feature_exclusion_reasons=reasons,
                               preparation_qc_pass=preparation_ok, dynamic_qc_pass=dynamic_ok, static_eligible=static_ok))
    expected = {(cid, t) for cid, info in candidate_info.items() for t in TEMPLATES[info['system']]}
    if set(parsed) != expected:
        raise ValueError('The complete76-row candidate/template scope is required; missing=' + str(sorted(expected - set(parsed))))
    if len(parsed) != 76:
        raise ValueError('Unexpected template expansion count')
    for key, values in counts.items():
        if len(values) != 1:
            raise ValueError('Fixed environmental feature region changed across candidates: ' + str(key))
    return parsed, exclusions, sorted(checkpoint_rows)


def run_analysis(args):
    protocol_path, split_path, binding_path = HERE / 'fixed_analysis_protocol.json', HERE / 'candidate_splits.json', HERE / 'protocol_freeze_binding.json'
    protocol, split, binding = read_json(protocol_path), read_json(split_path), read_json(binding_path)
    if digest(protocol_path) != binding['files']['fixed_analysis_protocol.json'] or digest(split_path) != binding['files']['candidate_splits.json']:
        raise ValueError('Protocol or split file differs from its frozen hash')
    if float(protocol['calibrator']['lambda']) != LAMBDA or protocol['calibrator']['intercept']:
        raise ValueError('The fixed ridge contract changed')
    active_hash = protocol['base_model']['sha256']
    checkpoint = Path(args.model_checkpoint or protocol['base_model'].get('local_run_path', ''))
    if not checkpoint.is_file() or digest(checkpoint) != active_hash:
        raise ValueError('Actual local execution checkpoint is absent or has the wrong SHA')
    rule, variants, candidates, crosswalk, label_sources = load_labels(args.label_rule, split)
    features, feature_accounting, feature_hashes = load_features(args.features, candidates, active_hash, args.mode, args.feature_manifest)
    models = FULL_MODELS if args.mode == 'full' else RESTRICTED_MODELS
    baseline = 'sequence_plus_static' if args.mode == 'full' else 'sequence_only'
    joint = baseline + '_plus_RH' if args.mode == 'full' else 'sequence_plus_RH'
    single_r = baseline + '_plus_R' if args.mode == 'full' else 'sequence_plus_R'
    single_h = baseline + '_plus_H' if args.mode == 'full' else 'sequence_plus_H'
    predictions, coefficients, eligibility, skipped = [], [], [], []
    out = Path(args.output)
    if out.exists() and any(out.iterdir()):
        raise ValueError('Output folder is nonempty; use a new folder to preserve previous results')
    out.mkdir(parents=True, exist_ok=True)

    def eligible(version, cid, template):
        f, label = features[cid, template], variants[version][cid]
        return not f['feature_reasons'] and not label['reasons']

    def singles(version, system, template):
        return sorted(cid for cid, info in candidates.items() if info['system'] == system and info['analysis_set'] == 'primary_single_substitution' and eligible(version, cid, template))

    def sufficient(ids):
        positions = {candidates[cid]['mutation_positions'][0] for cid in ids}
        return len(ids) >= 8 and len(positions) >= 3

    for version in variants:
        for (cid, template), row in sorted(features.items()):
            label = variants[version][cid]
            info = candidates[cid]
            eligibility.append(dict(label_variant=version, candidate_id=cid, system=info['system'], panel_id=info['panel_id'],
                template_pdb=template, analysis_set=info['analysis_set'], mutation_positions=info['mutation_positions'],
                numeric_eligible=eligible(version, cid, template), exclusion_reasons=row['feature_reasons'] + label['reasons'],
                censoring=label['censoring'], study_ids=label['study_ids']))

    def fit_and_predict(version, analysis, system, template, train_ids, test_ids, fold_id,
                        source_system=None, source_template=None, held_position=None):
        if set(train_ids) & set(test_ids):
            raise ValueError('Biological candidate leaked from train to test')
        source_system = source_system or system
        source_template = source_template or template
        if any(candidates[cid]['analysis_set'] != 'primary_single_substitution' for cid in train_ids):
            raise ValueError('Parent or multiple mutant entered fitting')
        for model, columns in models.items():
            train_y = np.asarray([variants[version][cid]['delta_pKD'] for cid in train_ids], dtype=np.float64)
            if columns:
                train_x = np.asarray([[features[cid, source_template]['values'][c] for c in columns] for cid in train_ids], dtype=np.float64)
                test_x = np.asarray([[features[cid, template]['values'][c] for c in columns] for cid in test_ids], dtype=np.float64)
                beta, scale = fit_fixed_ridge(train_x, train_y)
                estimate = apply_fixed_ridge(test_x, beta, scale)
            else:
                beta, scale = np.asarray([]), np.asarray([])
                estimate = np.zeros(len(test_ids))
            for i, c in enumerate(columns):
                coefficients.append(dict(analysis=analysis, label_variant=version, system=system, template_pdb=template,
                    source_system=source_system, source_template_pdb=source_template, fold_id=fold_id, model=model,
                    feature=c, coefficient_scaled=float(beta[i]), training_RMS_scale=float(scale[i]),
                    coefficient_raw_units=float(beta[i] / scale[i]), n_train_candidates=len(train_ids),
                    lambda_mean_MSE=LAMBDA, fitted_intercept=0.0, training_candidate_ids=train_ids, test_candidate_ids=test_ids))
            for cid, value in zip(test_ids, estimate):
                if not math.isfinite(float(value)):
                    raise ValueError('Nonfinite heldout prediction')
                label = variants[version][cid]
                predictions.append(dict(analysis=analysis, label_variant=version, is_primary_label=version == rule['primary'],
                    system=system, panel_id=candidates[cid]['panel_id'], template_pdb=template,
                    is_primary_template=template == PRIMARY[system], source_system=source_system, source_template_pdb=source_template,
                    fold_id=fold_id, heldout_mutation_position=held_position, model=model, candidate_id=cid,
                    sequence=candidates[cid]['sequence'], mutation_positions=candidates[cid]['mutation_positions'],
                    analysis_set=candidates[cid]['analysis_set'], n_train_candidates=len(train_ids),
                    true_delta_pKD=float(label['delta_pKD']), predicted_delta_pKD=float(value),
                    absolute_error_delta_pKD=abs(float(value) - label['delta_pKD']), study_ids=label['study_ids'],
                    study_delta_pKD_json=label['study_delta_pKD_json'],
                    training_candidate_ids=train_ids,
                    numerical_unit='delta_pKD', independent_observation='biological_candidate'))

    # Every listed template is evaluated, with one test prediction per eligible
    # single mutant. The two primary templates were selected before outputs.
    for version in variants:
        for system, templates in TEMPLATES.items():
            for template in templates:
                ids = singles(version, system, template)
                if not sufficient(ids):
                    skipped.append(dict(analysis='within_system_leave_position_out', label_variant=version, system=system, template_pdb=template,
                        reason='Fewer than8 eligible single mutants or fewer than3 positions; fixed fallback is descriptive only.', n_eligible=len(ids)))
                    continue
                for fold in split['folds']:
                    if fold['system'] != system:
                        continue
                    train_ids = sorted(crosswalk[cid] for cid in fold['train_candidates'] if crosswalk[cid] in ids)
                    test_ids = sorted(crosswalk[cid] for cid in fold['test_candidates'] if crosswalk[cid] in ids)
                    if not test_ids:
                        continue
                    position = fold['heldout_position']
                    if any(position in candidates[cid]['mutation_positions'] for cid in train_ids):
                        raise ValueError('Heldout mutation position leaked into training')
                    fit_and_predict(version, 'within_system_leave_position_out', system, template, train_ids, test_ids,
                                    fold['fold_id'], held_position=position)
                # Multiple substitutions are an explicitly separate check,
                # never folded into the single-mutation primary sample size.
                doubles = sorted(cid for cid, info in candidates.items() if info['system'] == system and
                                 info['analysis_set'] == 'secondary_multiple_substitution' and eligible(version, cid, template))
                if doubles:
                    fit_and_predict(version, 'secondary_multiple_substitution', system, template, ids, doubles,
                                    system + '_train_all_single_test_multiple')

        # Exactly two transport directions, using fixed primary templates only.
        for source, target in [('nyeso', 'tax'), ('tax', 'nyeso')]:
            train_ids = singles(version, source, PRIMARY[source])
            test_ids = singles(version, target, PRIMARY[target])
            if not sufficient(train_ids) or not sufficient(test_ids):
                skipped.append(dict(analysis='secondary_cross_system_transport', label_variant=version, system=target, template_pdb=PRIMARY[target],
                    source_system=source, reason='Source or target primary cohort fails the predeclared minimum sample/position gate.'))
                continue
            fit_and_predict(version, 'secondary_cross_system_transport', target, PRIMARY[target], train_ids, test_ids,
                            source + '_to_' + target, source_system=source, source_template=PRIMARY[source])

    group_fields = ['analysis', 'label_variant', 'system', 'panel_id', 'template_pdb', 'source_system', 'source_template_pdb', 'model']
    grouped = defaultdict(list)
    for row in predictions:
        grouped[tuple(row[k] for k in group_fields)].append(row)
    metrics, rankings, study_metrics, study_errors = [], [], [], []
    for key, rows in sorted(grouped.items()):
        identity = dict(zip(group_fields, key))
        metrics.append(dict(identity, is_primary_label=identity['label_variant'] == rule['primary'],
                            is_primary_template=identity['template_pdb'] == PRIMARY[identity['system']], **candidate_metrics(rows)))
        ordered_pred = sorted(rows, key=lambda r: (-r['predicted_delta_pKD'], r['candidate_id']))
        ordered_true = sorted(rows, key=lambda r: (-r['true_delta_pKD'], r['candidate_id']))
        truth_ranks = {r['candidate_id']: i + 1 for i, r in enumerate(ordered_true)}
        for rank, row in enumerate(ordered_pred, 1):
            rankings.append(dict(identity, candidate_id=row['candidate_id'], sequence=row['sequence'],
                                 predicted_rank=rank, experimental_rank=truth_ranks[row['candidate_id']],
                                 true_delta_pKD=row['true_delta_pKD'], predicted_delta_pKD=row['predicted_delta_pKD'],
                                 cutoff_tie_rule='lexical_biological_candidate_ID', analysis_set=row['analysis_set']))
        # The already held-out candidate prediction is scored separately
        # against each contributing study. No study-specific model is chosen
        # and multiple studies never increase independent candidate n.
        per_study = defaultdict(list)
        for row in rows:
            values = json.loads(row['study_delta_pKD_json'] or '{}')
            if not isinstance(values, dict):
                raise ValueError('Study delta metadata must be a study_id->numeric_delta dictionary')
            for study, truth in values.items():
                truth = finite_number(truth)
                if truth is None:
                    raise ValueError('Invalid study-specific delta in a curated eligible label')
                item = dict(row, true_delta_pKD=truth)
                per_study[study].append(item)
                study_errors.append(dict(identity, study_id=study, candidate_id=row['candidate_id'],
                    predicted_delta_pKD=row['predicted_delta_pKD'], study_true_delta_pKD=truth,
                    absolute_error_delta_pKD=abs(row['predicted_delta_pKD']-truth),
                    prediction_source='Same held-out candidate prediction as the candidate-median analysis; not refitted by study.',
                    independent_unit='biological_candidate_within_study'))
        for study, study_rows in sorted(per_study.items()):
            study_metrics.append(dict(identity, study_id=study, **candidate_metrics(study_rows),
                prediction_source='Same held-out candidate predictions; no study-specific calibration or selection.',
                cannot_replace_primary_candidate_median_result=True))

    by_comparison = defaultdict(dict)
    for row in metrics:
        key = tuple(row[k] for k in group_fields if k != 'model')
        by_comparison[key][row['model']] = row
    increments = []
    comparison_fields = [k for k in group_fields if k != 'model']
    for key, rows in sorted(by_comparison.items()):
        base, together, only_r, only_h = [rows[m]['MAE_delta_pKD'] for m in [baseline, joint, single_r, single_h]]
        increments.append(dict(dict(zip(comparison_fields, key)), baseline_model=baseline, joint_model=joint,
            n_candidates=rows[baseline]['n_candidates'], baseline_MAE=base, joint_MAE=together,
            R_only_augmented_MAE=only_r, H_only_augmented_MAE=only_h,
            joint_increment_vs_baseline=base - together, joint_increment_vs_R_only=only_r - together,
            joint_increment_vs_H_only=only_h - together,
            static_increment_claim_eligible=args.mode == 'full'))

    position_groups = defaultdict(list)
    for row in predictions:
        if row['analysis'] == 'within_system_leave_position_out':
            key = tuple(row[k] for k in group_fields) + (row['heldout_mutation_position'],)
            position_groups[key].append(row)
    position_metrics = []
    for key, rows in sorted(position_groups.items()):
        position_metrics.append(dict(dict(zip(group_fields + ['heldout_mutation_position'], key)),
                                     n_train_candidates=rows[0]['n_train_candidates'], **candidate_metrics(rows)))

    claims = {}
    for version in variants:
        selected = [r for r in increments if r['analysis'] == 'within_system_leave_position_out' and
                    r['label_variant'] == version and r['template_pdb'] == PRIMARY[r['system']]]
        complete = {r['system'] for r in selected} == set(PRIMARY)
        claims[version] = dict(complete_two_system_primary_evaluation=complete,
            consistent_joint_increment_both_systems=complete and all(r['joint_increment_vs_baseline'] > 0 for r in selected),
            complementary_increment_over_each_readout_both_systems=complete and all(
                r['joint_increment_vs_baseline'] > 0 and r['joint_increment_vs_R_only'] > 0 and r['joint_increment_vs_H_only'] > 0 for r in selected),
            per_system=selected,
            interpretation='Descriptive finite-panel consistency only; no significance, broad-TCR generalization, new-antigen discovery or causal mechanism claim.')

    # Two-system macro averages retain the independent-system count explicitly.
    macro = []
    for version in variants:
        for model in models:
            selected = [r for r in metrics if r['analysis'] == 'within_system_leave_position_out' and
                        r['label_variant'] == version and r['model'] == model and r['is_primary_template']]
            if {r['system'] for r in selected} == set(PRIMARY):
                macro.append(dict(label_variant=version, model=model, analysis='within_system_leave_position_out_primary_templates',
                    systems=2, system_candidate_counts={r['system']: r['n_candidates'] for r in selected},
                    macro_MAE_delta_pKD=sum(r['MAE_delta_pKD'] for r in selected) / 2,
                    population_generalization_CI=None, significance_test=None))

    write_csv(out / 'candidate_predictions.csv', predictions)
    write_csv(out / 'candidate_rankings.csv', rankings)
    write_csv(out / 'candidate_eligibility.csv', eligibility)
    write_csv(out / 'feature_accounting.csv', feature_accounting)
    write_csv(out / 'ridge_coefficients_and_training_scales.csv', coefficients)
    write_csv(out / 'candidate_level_metrics.csv', metrics)
    write_csv(out / 'mutation_position_metrics.csv', position_metrics)
    write_csv(out / 'dynamic_increment_comparisons.csv', increments)
    write_csv(out / 'two_system_macro_metrics.csv', macro)
    write_csv(out / 'study_stratified_metrics.csv', study_metrics)
    write_csv(out / 'study_stratified_candidate_errors.csv', study_errors)
    write_json(out / 'predeclared_joint_claim_checks.json', claims)
    write_json(out / 'skipped_analyses.json', skipped)
    source_paths = [protocol_path, split_path, binding_path, Path(args.features), Path(args.label_rule), Path(__file__)]
    if args.feature_manifest:
        source_paths.append(Path(args.feature_manifest))
    audit = dict(status='completed' if predictions else 'descriptive_only_no_eligible_trained_comparison',
        generated_at_utc=datetime.now(timezone.utc).isoformat(), mode=args.mode, protocol_id=protocol['protocol_id'],
        checkpoint_sha256=active_hash, actual_checkpoint_path=str(checkpoint), actual_checkpoint_sha256=digest(checkpoint),
        training_records=protocol['base_model']['records'], same_weights_as_current_manuscript_2430=protocol['base_model'].get('same_weights_as_current_manuscript_2430', True),
        completely_transitive_family_disjoint=protocol['base_model'].get('fully_transitive_family_disjoint', True),
        source_files=[dict(path=str(p), sha256=digest(p)) for p in source_paths] + label_sources,
        n_biological_candidates=len(candidates), n_input_template_rows=len(features), independent_TCR_systems=2,
        prediction_rows=len(predictions), candidate_prediction_rows_are_not_independent_n=True,
        n_primary_single_substitution_candidates=sum(r['analysis_set'] == 'primary_single_substitution' for r in candidates.values()),
        all_prelisted_templates_evaluated=True, label_versions=list(variants), primary_label_version=rule['primary'],
        ridge_lambda_mean_MSE=LAMBDA, intercept=False, scaler_training_fold_only=True, candidate_and_position_leak_checks_passed=True,
        model_inference_performed=False, neural_weight_training_performed=False, remote_access=False,
        base_model_weight_updates=False, feature_selection=False, hyperparameter_search=False,
        primary_metrics='Candidate-level MAE for the13 prelisted single mutants per system after technical and label eligibility, with every candidate held out by mutation position.',
        static_increment_claim_eligible=args.mode == 'full',
        study_sensitivity='Each contributing study scores the same already-held-out candidate predictions. No study model is selected and these additional rows do not increase independent candidate n.',
        limitations=['Historical2431 pilot differs from the manuscript2430 weights and contains one disclosed indirect family-component member.',
                    'Only2 related HLA-A2 TCR systems; transport is exploratory.',
                    'Fixed-backbone sequence changes do not sample new conformations.',
                    'No candidate, position or structural feature is selected by the observed score.',
                    'Published study/assay-matched parent normalization is not evidence of a known same-batch paired experiment.'],
        undefined_rho_and_tie_rules='Spearman requires at least3 candidates and nonconstant ranks; ties use average ranks. Pairwise prediction ties receive0.5 credit. Top3 cutoff ties use lexical biological candidate ID and are flagged.',
        skipped_analyses=skipped)
    write_json(out / 'evaluation_execution_audit.json', audit)
    print(json.dumps(dict(status=audit['status'], output=str(out), mode=args.mode, prediction_rows=len(predictions), independent_TCR_systems=2), ensure_ascii=False))


def self_test():
    # Closed form: mean((2-beta)^2)+beta^2 has optimum beta=1.
    beta, scale = fit_fixed_ridge(np.ones((4, 1)), np.full(4, 2.0))
    assert np.allclose(beta, [1.0]) and np.allclose(scale, [1.0])
    # Heldout100 must not alter the scale learned from training[1,2].
    x = np.array([[1.0, 0.0], [2.0, 0.0]])
    y = np.array([2.0, 4.0])
    beta, scale = fit_fixed_ridge(x, y)
    assert np.allclose(scale, [math.sqrt(2.5), 1.0])
    assert np.allclose(apply_fixed_ridge([[100.0, 1000.0]], beta, scale), [100.0])
    # Direct objective derivative verifies the n_train*lambda convention.
    xs = x / scale
    derivative = 2 / len(y) * xs.T @ (xs @ beta - y) + 2 * LAMBDA * beta
    assert np.allclose(derivative, 0.0)
    assert np.allclose(rank_average_ties([3, 1, 1, 2]), [4, 1.5, 1.5, 3])
    metrics = candidate_metrics([dict(candidate_id='a', true_delta_pKD=1., predicted_delta_pKD=1.),
                                 dict(candidate_id='b', true_delta_pKD=2., predicted_delta_pKD=2.),
                                 dict(candidate_id='c', true_delta_pKD=3., predicted_delta_pKD=3.)])
    assert metrics['MAE_delta_pKD'] == 0 and math.isclose(metrics['Spearman_rho'], 1.)
    assert metrics['candidate_pair_concordance'] == 1. and metrics['top_k_overlap'] == 3
    ties = candidate_metrics([dict(candidate_id='a', true_delta_pKD=1., predicted_delta_pKD=0.),
                              dict(candidate_id='b', true_delta_pKD=2., predicted_delta_pKD=0.),
                              dict(candidate_id='c', true_delta_pKD=3., predicted_delta_pKD=0.)])
    assert ties['Spearman_rho'] is None and ties['candidate_pair_concordance'] == .5
    # Frozen positions: every primary candidate is tested exactly once, and no
    # held-out position is used to calibrate its own prediction.
    split = read_json(HERE / 'candidate_splits.json')
    byid = {r['candidate_id']: r for r in split['candidates']}
    tested = []
    for fold in split['folds']:
        assert set(fold['train_candidates']).isdisjoint(fold['test_candidates'])
        assert all(fold['heldout_position'] not in byid[cid]['mutation_positions'] for cid in fold['train_candidates'])
        assert all(byid[cid]['mutation_positions'] == [fold['heldout_position']] for cid in fold['test_candidates'])
        tested.extend(fold['test_candidates'])
    assert Counter(tested) == Counter(r['candidate_id'] for r in split['candidates'] if r['analysis_set'] == 'primary_single_substitution')
    print(json.dumps(dict(status='passed_synthetic_and_split_tests', real_affinity_values_read=False,
                         real_candidate_features_read=False, ridge_mean_loss_scaling_verified=True,
                         source_only_scaling_verified=True, frozen_candidate_position_splits_verified=True)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--features', help='Complete76-row biological candidate/template feature CSV')
    parser.add_argument('--feature-manifest', help='Optional feature provenance JSON; checkpoint_sha256 required here or in every CSV row')
    parser.add_argument('--label-rule', default=str(HERE.parent / 'experimental_labels/label_rule_freeze.json'))
    parser.add_argument('--model-checkpoint', help='Actual local execution checkpoint to byte-verify; default active protocol path')
    parser.add_argument('--output', default=str(HERE / 'results_fixed_protocol'))
    parser.add_argument('--mode', choices=['full', 'restricted'], default='full', help='Full static baseline by default; restricted explicitly forbids a static-increment claim')
    parser.add_argument('--self-test', action='store_true', help='Run synthetic numerical and frozen-split checks only, without real labels/features')
    args = parser.parse_args()
    if args.self_test:
        self_test()
    else:
        if not args.features:
            parser.error('--features is required unless --self-test is used')
        run_analysis(args)


if __name__ == '__main__':
    main()
