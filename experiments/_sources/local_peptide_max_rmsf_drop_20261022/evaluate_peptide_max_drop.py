"""Frozen, zero-label test of the user-specified any-peptide-site RMSF-drop rule."""
from pathlib import Path
from collections import Counter
from datetime import datetime, timezone
import csv, hashlib, json
import numpy as np
import pandas as pd
from scipy.stats import hypergeom, spearmanr
from sklearn.metrics import roc_auc_score, average_precision_score
from Bio.Align import substitution_matrices
from Bio.PDB import PDBParser

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT.parent / 'local_multimode_screening_20261022'
OUT = ROOT / 'functional_screening'
OUT.mkdir(exist_ok=True)
PANELS = [f'MSKCC_TCR{i}' for i in range(1, 8)]
CHECKPOINT = 'b3eed060919df15f075329bc892f959d9355232be5feb561417edd49ea82b39f'

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def read(path):
    with Path(path).open(encoding='utf-8-sig', newline='') as handle:
        return list(csv.DictReader(handle))

def write(rows, name):
    pd.DataFrame(rows).to_csv(OUT / name, index=False, encoding='utf-8-sig')

def top_weights(score, budget):
    """Fractional membership is exact expectation under random ordering of ties."""
    score = np.asarray(score, dtype=float)
    assert np.isfinite(score).all() and 0 < budget <= len(score)
    cutoff = np.sort(score)[-budget]
    higher = score > cutoff
    equal = score == cutoff
    weights = higher.astype(float)
    weights[equal] = (budget - higher.sum()) / equal.sum()
    assert abs(weights.sum() - budget) < 1e-10
    return weights

def matched_expectation(weights, positions, hits):
    answer = 0.0
    for position in sorted(set(positions)):
        mask = positions == position
        answer += float(weights[mask].sum() * hits[mask].mean())
    return answer

def parent_contacts(path):
    model = PDBParser(QUIET=True).get_structure('parent', path)[0]
    output = {}
    partners = {'hla': ['A'], 'tcr': ['C', 'D']}
    for residue in model['B']:
        if residue.id[0] != ' ':
            continue
        residue_id = str(residue.id[1]) + residue.id[2].strip()
        peptide_atoms = np.array([a.coord for a in residue if a.element not in ('H', 'D')])
        record = {}
        for role, chains in partners.items():
            close = []
            for chain in chains:
                for partner in model[chain]:
                    if partner.id[0] != ' ':
                        continue
                    atoms = np.array([a.coord for a in partner if a.element not in ('H', 'D')])
                    if not len(atoms):
                        continue
                    distance = np.sqrt(((peptide_atoms[:, None] - atoms[None])**2).sum(-1)).min()
                    if distance <= 4.5:
                        close.append(f'{chain}:{partner.id[1]}{partner.id[2].strip()}')
            record[f'parent_{role}_contacts_4p5A'] = len(close)
            record[f'parent_{role}_contact_residues'] = ';'.join(close)
        output[residue_id] = record
    assert len(output) == 9
    return output

def main():
    protocol = ROOT / 'protocol.json'
    protocol_hash = sha(protocol)
    assert protocol_hash == 'b41817fb7b7536db8644124fc5b5b29450c33bac1fb992f509852a15e5f748f4'
    label_path = SOURCE / 'antigen_labels/curated/curated_candidates.csv'
    assert sha(label_path) == json.loads(protocol.read_text(encoding='utf-8'))['source_labels_sha256']
    labels = read(label_path)
    static_path = SOURCE / 'structures/static_features.csv'
    static = {r['candidate_id']: r for r in read(static_path)}
    manifest_path = SOURCE / 'structures/input_manifest.json'
    templates = {r['parent_candidate_id'].removesuffix('__parent'): r for r in json.loads(manifest_path.read_text())['templates'] if r['template_id'].endswith('__AF3')}
    assert set(templates) == set(PANELS)
    blosum = substitution_matrices.load('BLOSUM62')
    candidates, profiles, source_hashes, contacts_table = [], [], [], []
    max_delta_residual = 0.0
    for panel in PANELS:
        parent_path = SOURCE / 'predictions/full_pool' / f'{panel}__parent/prediction.npz'
        parent = np.load(parent_path)
        assert parent['metrics'].tolist() == ['rmsf', 'chi_entropy']
        assert str(parent['checkpoint_sha256']) == CHECKPOINT
        mask = parent['chain'] == 'B'
        assert mask.sum() == 9 and parent['geometry_valid'][mask].all()
        parent_r = parent['mean'][mask, 0].astype(float)
        parent_aa = parent['amino_acid'][mask]
        residue_ids = parent['author_residue_id'][mask]
        assert residue_ids.tolist() == [str(i) for i in range(1, 10)]
        template = templates[panel]
        assert sha(template['input_pdb']) == template['input_sha256']
        contacts = parent_contacts(template['input_pdb'])
        for j, residue_id in enumerate(residue_ids):
            contacts_table.append(dict(panel_id=panel, peptide_position=j+1, author_residue_id=residue_id, parent_aa=parent_aa[j], parent_rmsf_A=parent_r[j], **contacts[residue_id]))
        source_hashes.append(dict(panel_id=panel, candidate_id=f'{panel}__parent', path=str(parent_path), sha256=sha(parent_path)))
        panel_labels = [r for r in labels if r['panel_id'] == panel]
        assert len(panel_labels) == 171
        for label in panel_labels:
            path = SOURCE / 'predictions/full_pool' / label['candidate_id'] / 'prediction.npz'
            d = np.load(path)
            assert d['metrics'].tolist() == ['rmsf', 'chi_entropy']
            assert str(d['checkpoint_sha256']) == CHECKPOINT
            assert np.array_equal(d['chain'], parent['chain'])
            assert np.array_equal(d['author_residue_id'], parent['author_residue_id'])
            assert d['geometry_valid'][mask].all()
            assert str(d['parent_candidate_id']) == f'{panel}__parent'
            mutant_aa = d['amino_acid'][mask]
            mutated = np.flatnonzero(mutant_aa != parent_aa)
            assert len(mutated) == 1
            mutation_index = int(mutated[0])
            assert int(label['position']) == mutation_index + 1
            assert label['parent_aa'] == parent_aa[mutation_index]
            assert label['mutant_aa'] == mutant_aa[mutation_index]
            assert label['sequence'] == ''.join(mutant_aa)
            mutant_r = d['mean'][mask, 0].astype(float)
            drops = parent_r - mutant_r
            residual = float(abs(drops + d['delta_mean_vs_parent'][mask, 0]).max())
            max_delta_residual = max(max_delta_residual, residual)
            assert residual < 1e-6
            assert np.array_equal(d['parent_mean'][mask, 0], parent['mean'][mask, 0])
            maximum = float(drops.max())
            sites = np.flatnonzero(drops == maximum)
            # All tied maximum sites are retained, never selected according to labels.
            site = int(sites[0])
            s = static[label['candidate_id']]
            record = dict(panel_id=panel, family_id=label['family_id'], candidate_id=label['candidate_id'], sequence=label['sequence'], mutation=label['mutation'], mutation_position=mutation_index+1, parent_aa=label['parent_aa'], mutant_aa=label['mutant_aa'], functional_response=float(label['value_oriented_higher_better']), is_hit=int(label['is_hit']), max_rmsf_drop_A=maximum, second_largest_drop_A=float(np.sort(drops)[-2]), mutated_site_drop_A=float(drops[mutation_index]), max_drop_position=site+1, all_tied_max_positions=';'.join(str(i+1) for i in sites), max_at_mutated_site=bool(mutation_index in sites), parent_rmsf_at_max_A=float(parent_r[site]), raw_BLOSUM62=float(blosum[label['parent_aa'], label['mutant_aa']]), HERMES=float(s['HERMES_fixed_raw']), ProteinMPNN020=float(s['ProteinMPNN_raw']), ProteinMPNN002=float(s['ProteinMPNN002_raw']), **contacts[str(site+1)])
            assert record['is_hit'] == int(record['functional_response'] >= .3)
            for j in range(9):
                record[f'drop_P{j+1}_A'] = float(drops[j])
                profiles.append(dict(panel_id=panel, candidate_id=label['candidate_id'], mutation=label['mutation'], peptide_position=j+1, parent_aa=parent_aa[j], mutant_aa=mutant_aa[j], parent_rmsf_A=parent_r[j], mutant_rmsf_A=mutant_r[j], rmsf_drop_A=drops[j], is_mutated_site=j == mutation_index, is_max_drop_site=j in sites, functional_response=record['functional_response'], is_hit=record['is_hit']))
            candidates.append(record)
            source_hashes.append(dict(panel_id=panel, candidate_id=label['candidate_id'], path=str(path), sha256=sha(path)))
    assert len(candidates) == 1197 and len(profiles) == 10773
    write(candidates, 'all_candidate_scores.csv')
    write(profiles, 'all_peptide_residue_profiles.csv')
    write(contacts_table, 'parent_site_contacts.csv')
    write(source_hashes, 'prediction_source_hashes.csv')
    frame = pd.DataFrame(candidates)
    score_columns = ['max_rmsf_drop_A', 'raw_BLOSUM62', 'HERMES', 'ProteinMPNN020', 'ProteinMPNN002', 'second_largest_drop_A']
    metrics, curves, thresholds, position_rows, selections = [], [], [], [], []
    for panel in PANELS:
        pool = frame[frame.panel_id == panel].sort_values('candidate_id').reset_index(drop=True)
        positions = pool.mutation_position.to_numpy()
        response = pool.functional_response.to_numpy()
        primary_hits = pool.is_hit.to_numpy()
        n = len(pool)
        for endpoint in ['functional_hit_0p3', 'relative_top_decile']:
            cutoff = .3 if endpoint == 'functional_hit_0p3' else float(np.quantile(response, .9))
            hits = (response >= cutoff).astype(int)
            positives = int(hits.sum())
            for score_column in score_columns:
                score = pool[score_column].to_numpy()
                association = dict(AUROC=float(roc_auc_score(hits, score)), average_precision=float(average_precision_score(hits, score)), spearman_response=float(spearmanr(score, response).statistic))
                for budget in range(1, n+1):
                    weights = top_weights(score, budget)
                    found = float(weights @ hits)
                    expected = budget * positives / n
                    matched = matched_expectation(weights, positions, hits)
                    row = dict(panel_id=panel, family_id=pool.family_id.iloc[0], endpoint=endpoint, endpoint_cutoff=cutoff, score=score_column, budget=budget, n_candidates=n, pool_hits=positives, expected_hits=found, precision=found/budget, pool_prevalence=positives/n, enrichment=found/expected if expected else np.nan, random_expected_hits=expected, random_hits_p10=float(hypergeom.ppf(.1, n, positives, budget)), random_hits_p90=float(hypergeom.ppf(.9, n, positives, budget)), mutation_position_matched_random_hits=matched, gain_vs_position_matched_random=found-matched, score_cutoff=float(np.sort(score)[-budget]))
                    curves.append(row)
                    if budget in (10, 20, 30, 60):
                        metrics.append(dict(**row, **association))
                    if endpoint == 'functional_hit_0p3' and score_column == 'max_rmsf_drop_A' and budget == 30:
                        order = np.lexsort((pool.candidate_id.to_numpy(), -score))
                        for rank, index in enumerate(order, 1):
                            selections.append(dict(panel_id=panel, rank_lexical_within_ties=rank, candidate_id=pool.iloc[index].candidate_id, max_rmsf_drop_A=score[index], top30_expected_membership=weights[index], functional_response=response[index], is_hit=int(hits[index])))
                        for position in range(1, 10):
                            mask = positions == position
                            position_rows.append(dict(panel_id=panel, mutation_position=position, pool_n=int(mask.sum()), pool_hits=int(hits[mask].sum()), selected_expected_n=float(weights[mask].sum()), selected_expected_hits=float(weights[mask] @ hits[mask]), matched_random_expected_hits=float(weights[mask].sum()*hits[mask].mean())))
        # Complete threshold curve, including a pre-specified strict D>0 descriptive gate.
        scores = pool.max_rmsf_drop_A.to_numpy()
        for threshold in sorted(set(scores.tolist())):
            selected = scores >= threshold
            k = int(selected.sum())
            thresholds.append(dict(panel_id=panel, threshold_A=threshold, condition='D>=threshold', selected_n=k, hits=int(primary_hits[selected].sum()), precision=float(primary_hits[selected].mean()), random_expected_hits=k*float(primary_hits.mean())))
        selected = scores > 0
        k = int(selected.sum())
        thresholds.append(dict(panel_id=panel, threshold_A=0., condition='D>0', selected_n=k, hits=int(primary_hits[selected].sum()), precision=float(primary_hits[selected].mean()) if k else np.nan, random_expected_hits=k*float(primary_hits.mean())))
    write(metrics, 'fixed_budget_metrics.csv')
    write(curves, 'all_budget_curves.csv')
    write(thresholds, 'all_thresholds.csv')
    write(position_rows, 'top30_mutation_position_diagnostic.csv')
    write(selections, 'max_drop_rankings_all_candidates.csv')
    # Optional pre-declared context with identical observed seeds; this is NOT zero-label validation.
    observed_rows = []
    for panel in PANELS:
        pool = frame[frame.panel_id == panel].sort_values('candidate_id').reset_index(drop=True)
        records = read(SOURCE / 'evaluation/results_stronger_baselines' / panel / 'selected_candidates.csv')
        seeds = sorted({int(r['seed']) for r in records if r['mode'] == 'mode2'})
        assert seeds == list(range(20261022, 20261042))
        by_id = pool.set_index('candidate_id')
        for seed in seeds:
            rows = [r for r in records if r['mode'] == 'mode2' and int(r['seed']) == seed]
            initial_sets = {tuple(json.loads(r['initial_candidate_ids'])) for r in rows}
            assert len(initial_sets) == 1
            initial = set(next(iter(initial_sets)))
            assert len(initial) == 18
            remaining = pool[~pool.candidate_id.isin(initial)].reset_index(drop=True)
            weights = top_weights(remaining.max_rmsf_drop_A, 30)
            hits = remaining.is_hit.to_numpy()
            random_expected = 30*float(hits.mean())
            observed_rows.append(dict(panel_id=panel, seed=seed, strategy='max_rmsf_drop_A_no_calibration', observed_initial_n=18, extra_budget=30, remaining_n=len(remaining), new_hits=float(weights @ hits), random_expected_hits=random_expected))
            for method in ('sequence_only_control', 'BLOSUM_only_control', 'baseline_R', 'baseline_RH'):
                r, = [r for r in rows if r['strategy'] == method]
                chosen = json.loads(r['selected_candidate_ids'])
                assert len(chosen) == len(set(chosen)) == 30 and not (set(chosen) & initial)
                observed_rows.append(dict(panel_id=panel, seed=seed, strategy=method, observed_initial_n=18, extra_budget=30, remaining_n=len(remaining), new_hits=float(by_id.loc[chosen, 'is_hit'].sum()), random_expected_hits=random_expected))
    write(observed_rows, 'observed18_context_seed_metrics.csv')
    observed_summary = pd.DataFrame(observed_rows).groupby(['panel_id','strategy'], as_index=False).agg(mean_new_hits=('new_hits','mean'), min_new_hits=('new_hits','min'), max_new_hits=('new_hits','max'), seeds=('seed','nunique'))
    observed_summary.to_csv(OUT/'observed18_context_summary.csv', index=False, encoding='utf-8-sig')
    summary = pd.DataFrame(metrics)
    summary = summary[(summary.budget == 30) & (summary.endpoint == 'functional_hit_0p3')]
    summary.to_csv(OUT/'primary_top30_metrics.csv', index=False, encoding='utf-8-sig')
    diagnostics = []
    for panel, group in frame.groupby('panel_id', sort=False):
        scores = group.max_rmsf_drop_A
        diagnostics.append(dict(panel_id=panel, total_n=len(group), max_drop_min_A=float(scores.min()), max_drop_median_A=float(scores.median()), max_drop_max_A=float(scores.max()), positive_max_drop_n=int((scores>0).sum()), max_at_mutated_site_n=int(group.max_at_mutated_site.sum()), max_position_counts=dict(Counter(group.max_drop_position.astype(int))), top30_max_position_counts=dict(Counter(group.sort_values(['max_rmsf_drop_A','candidate_id'],ascending=[False,True]).head(30).max_drop_position.astype(int)))))
    audit = dict(status='passed', finished_at_utc=datetime.now(timezone.utc).isoformat(), frozen_protocol_sha256=protocol_hash, checkpoint_sha256=CHECKPOINT, candidates=1197, peptide_residue_comparisons=10773, parent_predictions=7, prediction_source_files=len(source_hashes), max_abs_delta_sign_audit_residual_A=max_delta_residual, labels_sha256=sha(label_path), static_scores_sha256=sha(static_path), input_manifest_sha256=sha(manifest_path), primary_rule='signed maximum over all nine peptide positions', score_columns=score_columns, primary_budget=30, all_budgets_archived=True, ties='exact expected random tie order', no_new_model_training=True, no_new_neural_inference=True, no_new_MD=True, manuscript_modified=False, diagnostics=diagnostics)
    (OUT/'input_and_calculation_audit.json').write_text(json.dumps(audit, ensure_ascii=False, indent=2),encoding='utf-8')
    print(summary[['panel_id','score','expected_hits','random_expected_hits','mutation_position_matched_random_hits','AUROC','spearman_response']].to_string(index=False))
    print('INPUT_AND_CALCULATION_AUDIT: passed')

if __name__ == '__main__':
    main()
