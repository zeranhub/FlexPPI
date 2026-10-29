"""Audit exact-sequence isolated/complex pairs using existing family OOF exports.

No model is loaded or run. All residues are retained in the profile export;
scoring uses the intersection of the two states' label-valid masks.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import warnings

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--directory', type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    directory = args.directory.resolve()
    inventory = json.loads((directory / 'paired_inventory.json').read_text(encoding='utf-8'))
    target_pairs = [p for p in inventory['target_chain_pairs'] if p['sequence_match'] == 'exact']
    sample_ids = {p[k] for p in target_pairs for k in ['isolated_sample_id', 'complex_sample_id']}
    source = directory / 'archive_extract/family_folds/evaluation/oof_per_residue.csv.gz'
    chunks = []
    for chunk in pd.read_csv(source, chunksize=250_000, keep_default_na=True, dtype={'author_residue_id': str}):
        selected = chunk[chunk.sample_id.isin(sample_ids)]
        if not selected.empty:
            chunks.append(selected)
    frame = pd.concat(chunks, ignore_index=True)
    assert set(frame.variant.unique()) == {'full', 'intrachain'}
    assert not frame.duplicated(['sample_id', 'variant', 'residue_index']).any()
    for mask_column in ['valid_rmsf', 'valid_chi_entropy']:
        assert pd.api.types.is_bool_dtype(frame[mask_column]), (mask_column, frame[mask_column].dtype)

    records = inventory['records']
    original_pairs = {
        (p['isolated_sample_id'], p['complex_sample_id']): p
        for p in inventory['pairs']
    }
    cache: dict[tuple[str, str, str], pd.DataFrame] = {}

    def chain_table(sample_id: str, variant: str, chain: str) -> pd.DataFrame:
        key = sample_id, variant, chain
        if key in cache:
            return cache[key]
        table = frame[(frame.sample_id == sample_id) & (frame.variant == variant) & (frame.chain == chain)]
        table = table.sort_values('residue_index').reset_index(drop=True)
        expected = [c for c in records[sample_id]['actual_input_chains'] if c['chain'] == chain]
        assert len(expected) == 1, (sample_id, chain, 'ambiguous input chain')
        expected = expected[0]
        assert len(table) == expected['n_residues'], (sample_id, variant, chain, len(table), expected['n_residues'])
        observed_sequence = ''.join(table.amino_acid)
        assert observed_sequence == expected['sequence'], (sample_id, variant, chain, 'OOF/input sequence mismatch')
        assert table.fold.nunique() == 1 and int(table.fold.iloc[0]) == records[sample_id]['fold']
        cache[key] = table
        return table

    metric_rows = []
    profile_rows = []
    for pair in target_pairs:
        isolated_id, complex_id = pair['isolated_sample_id'], pair['complex_sample_id']
        isolated_chain, complex_chain = pair['isolated_chain'], pair['complex_chain']
        isolated, complex_record = records[isolated_id], records[complex_id]
        assert isolated['fold'] == complex_record['fold'] == pair['fold'] == pair['complex_fold']
        assert isolated['family_component_id'] == complex_record['family_component_id'] == pair['family_component_id']
        tables = {
            (state, variant): chain_table(sample_id, variant, chain)
            for state, sample_id, chain in [('isolated', isolated_id, isolated_chain), ('complex', complex_id, complex_chain)]
            for variant in ['full', 'intrachain']
        }
        length = len(tables['isolated', 'full'])
        expected_sequence = ''.join(tables['isolated', 'full'].amino_acid)
        for table in tables.values():
            assert len(table) == length
            assert ''.join(table.amino_acid) == expected_sequence
        for state in ['isolated', 'complex']:
            full, intra = tables[state, 'full'], tables[state, 'intrachain']
            assert np.array_equal(full.residue_index, intra.residue_index)
            assert np.array_equal(full.author_residue_id.astype(str), intra.author_residue_id.astype(str))
            for metric in ['rmsf', 'chi_entropy']:
                assert np.array_equal(full[f'valid_{metric}'], intra[f'valid_{metric}'])
                assert np.allclose(full[f'target_{metric}'], intra[f'target_{metric}'], equal_nan=True, rtol=0, atol=1e-7)
        original_pair = original_pairs[isolated_id, complex_id]
        matched_copies = sum(m['sequence_match'] == 'exact' for m in original_pair['target_chain_mappings'])
        same_protocol = (
            pair['same_temperature']
            and pair['source_group_pair'] == ['test53_single', 'test53_complex']
        )
        pair_id = f'{complex_id}:{complex_chain}__minus__{isolated_id}:{isolated_chain}'
        meta = {
            'pair_id': pair_id,
            'pdb': pair['pdb'],
            'complex_sample_id': complex_id,
            'isolated_sample_id': isolated_id,
            'complex_chain': complex_chain,
            'isolated_chain': isolated_chain,
            'description': '; '.join(pair['descriptions']),
            'fold': pair['fold'],
            'family_component_id': pair['family_component_id'],
            'sequence_match': 'exact',
            'n_sequence_positions': length,
            'same_temperature': bool(pair['same_temperature']),
            'isolated_temperature_K': isolated['temperature_K'],
            'complex_temperature_K': complex_record['temperature_K'],
            'isolated_source_group': pair['source_group_pair'][0],
            'complex_source_group': pair['source_group_pair'][1],
            'same_protocol_test53': bool(same_protocol),
            'virus_related_exclude': bool(pair['virus_related_exclude']),
            'primary_analysis': bool(same_protocol and not pair['virus_related_exclude']),
            'analysis_subset': 'primary_test53_same_temperature' if same_protocol and not pair['virus_related_exclude'] else ('excluded_virus' if pair['virus_related_exclude'] else 'exploratory_source_or_temperature_confounding'),
            'same_pdb_record_pair_count': original_pair['same_pdb_pair_count'],
            'same_pdb_complex_record_count': len(original_pair['same_pdb_distinct_complex_records']),
            'same_record_exact_target_copy_count': matched_copies,
            'independent_system_id': pair['pdb'],
            'independence_warning': 'Same PDB and matching target copies are repeated views of one source system.',
        }
        profile = pd.DataFrame({
            **{k: [v] * length for k, v in meta.items()},
            'target_position_1based': np.arange(1, length + 1),
            'amino_acid': list(expected_sequence),
            'isolated_residue_index': tables['isolated', 'full'].residue_index,
            'complex_residue_index': tables['complex', 'full'].residue_index,
            'isolated_author_residue_id': tables['isolated', 'full'].author_residue_id,
            'complex_author_residue_id': tables['complex', 'full'].author_residue_id,
        })
        for metric, threshold in [('rmsf', 0.1), ('chi_entropy', 0.1)]:
            singles, complexes = tables['isolated', 'full'], tables['complex', 'full']
            valid_single = singles[f'valid_{metric}'].to_numpy(dtype=bool)
            valid_complex = complexes[f'valid_{metric}'].to_numpy(dtype=bool)
            valid = valid_single & valid_complex
            truth_single = singles[f'target_{metric}'].to_numpy(dtype=float)
            truth_complex = complexes[f'target_{metric}'].to_numpy(dtype=float)
            assert np.isfinite(truth_single[valid_single]).all()
            assert np.isfinite(truth_complex[valid_complex]).all()
            truth_delta = truth_complex - truth_single
            truth_delta_masked = np.where(valid, truth_delta, np.nan)
            responsive = valid & (np.abs(truth_delta) >= threshold)
            zero_mae = float(np.abs(truth_delta[valid]).mean()) if valid.any() else np.nan
            profile[f'valid_isolated_{metric}'] = valid_single
            profile[f'valid_complex_{metric}'] = valid_complex
            profile[f'valid_intersection_{metric}'] = valid
            profile[f'md_isolated_{metric}'] = truth_single
            profile[f'md_complex_{metric}'] = truth_complex
            profile[f'md_delta_{metric}'] = truth_delta_masked
            profile[f'descriptive_response_threshold_{metric}'] = threshold
            profile[f'md_response_ge_threshold_{metric}'] = responsive
            for variant in ['full', 'intrachain']:
                pred_single = tables['isolated', variant][f'pred_{metric}'].to_numpy(dtype=float)
                pred_complex = tables['complex', variant][f'pred_{metric}'].to_numpy(dtype=float)
                assert np.isfinite(pred_single[valid]).all() and np.isfinite(pred_complex[valid]).all()
                pred_delta = pred_complex - pred_single
                errors = np.abs(pred_delta[valid] - truth_delta[valid])
                mae = float(errors.mean()) if valid.any() else np.nan
                if valid.sum() >= 2 and np.ptp(pred_delta[valid]) > 0 and np.ptp(truth_delta[valid]) > 0:
                    with warnings.catch_warnings():
                        warnings.simplefilter('error')
                        rho = float(spearmanr(truth_delta[valid], pred_delta[valid]).statistic)
                else:
                    rho = np.nan
                direction_agreement = float(np.mean(np.sign(pred_delta[responsive]) == np.sign(truth_delta[responsive]))) if responsive.any() else np.nan
                row = {
                    **meta,
                    'metric': metric,
                    'variant': variant,
                    'delta_definition': 'complex_minus_isolated',
                    'unit': 'angstrom' if metric == 'rmsf' else 'nats',
                    'n_valid_intersection': int(valid.sum()),
                    'delta_mae': mae,
                    'zero_change_delta_mae': zero_mae,
                    'mae_improvement_over_zero': zero_mae - mae,
                    'relative_mae_improvement_over_zero': (zero_mae - mae) / zero_mae if zero_mae > 0 else np.nan,
                    'beats_zero_change': bool(mae < zero_mae) if valid.any() else False,
                    'delta_spearman': rho,
                    'md_delta_mean': float(truth_delta[valid].mean()) if valid.any() else np.nan,
                    'md_delta_median': float(np.median(truth_delta[valid])) if valid.any() else np.nan,
                    'pred_delta_mean': float(pred_delta[valid].mean()) if valid.any() else np.nan,
                    'pred_delta_median': float(np.median(pred_delta[valid])) if valid.any() else np.nan,
                    'descriptive_response_abs_threshold': threshold,
                    'n_md_response_ge_threshold': int(responsive.sum()),
                    'direction_agreement_on_thresholded_md': direction_agreement,
                    'direction_is_descriptive_only': True,
                }
                metric_rows.append(row)
                profile[f'{variant}_isolated_{metric}'] = pred_single
                profile[f'{variant}_complex_{metric}'] = pred_complex
                profile[f'{variant}_delta_{metric}'] = pred_delta
                profile[f'{variant}_delta_abs_error_{metric}'] = np.where(valid, np.abs(pred_delta - truth_delta), np.nan)
                profile[f'{variant}_direction_matches_{metric}'] = np.where(responsive, np.sign(pred_delta) == np.sign(truth_delta), np.nan)
            profile[f'zero_change_delta_abs_error_{metric}'] = np.where(valid, np.abs(truth_delta), np.nan)
        profile_rows.append(profile)
    metrics = pd.DataFrame(metric_rows)
    profiles = pd.concat(profile_rows, ignore_index=True)
    assert len(metrics) == len(target_pairs) * 4
    assert len(profiles) == sum(p['shared_length'] for p in target_pairs)
    metrics.to_csv(directory / 'pair_delta_metrics.csv', index=False, float_format='%.10g')
    profiles.to_csv(directory / 'pair_delta_profiles.csv.gz', index=False, compression='gzip', float_format='%.10g')

    full = metrics[(metrics.variant == 'full') & metrics.primary_analysis]
    summary = []
    for pair_id, group in full.groupby('pair_id', sort=True):
        assert set(group.metric) == {'rmsf', 'chi_entropy'}
        result = group.set_index('metric')
        first = group.iloc[0]
        summary.append({
            'pair_id': pair_id, 'pdb': first.pdb, 'description': first.description,
            'full_beats_zero_both_metrics': bool(group.beats_zero_change.all()),
            'rmsf': {k: float(result.loc['rmsf', k]) for k in ['delta_mae', 'zero_change_delta_mae', 'relative_mae_improvement_over_zero', 'delta_spearman']},
            'chi_entropy': {k: float(result.loc['chi_entropy', k]) for k in ['delta_mae', 'zero_change_delta_mae', 'relative_mae_improvement_over_zero', 'delta_spearman']},
        })
    audit = {
        'source': str(source),
        'exact_target_chain_pairs_including_virus': len(target_pairs),
        'metric_rows': len(metrics), 'profile_rows_all_positions': len(profiles),
        'primary_nonvirus_pair_count': int(metrics[metrics.primary_analysis].pair_id.nunique()),
        'primary_nonvirus_pdb_system_count': int(metrics[metrics.primary_analysis].pdb.nunique()),
        'assertions_passed': ['input sequence matches OOF residue order', 'all four paired tables share the identical target sequence', 'family folds/components match', 'Full/Intra labels and validity masks agree in each state', 'no duplicated OOF residue rows', 'all exact target positions retained'],
        'notes': ['All score improvements are descriptive comparisons on existing OOF outputs, not evidence of statistical significance.', 'Response thresholds are used only for descriptive direction agreement, never to choose positions for MAE/Spearman.', 'Virus-related pairs are exported and flagged but excluded from primary_analysis.', 'Temperature/source metadata do not prove identical force fields, equilibration, alignment, or independent replicate availability.'],
        'primary_full_summary': summary,
    }
    (directory / 'pair_delta_audit.json').write_text(json.dumps(audit, ensure_ascii=False, indent=2, allow_nan=True), encoding='utf-8')
    print(json.dumps({k: audit[k] for k in ['exact_target_chain_pairs_including_virus', 'metric_rows', 'profile_rows_all_positions', 'primary_nonvirus_pair_count', 'primary_nonvirus_pdb_system_count']}, indent=2))
    for entry in summary:
        print(entry['pdb'], entry['pair_id'], entry['description'], 'Full both:', entry['full_beats_zero_both_metrics'], 'RMSF rel:', round(entry['rmsf']['relative_mae_improvement_over_zero'], 4), 'Hchi rel:', round(entry['chi_entropy']['relative_mae_improvement_over_zero'], 4))


if __name__ == '__main__':
    main()
