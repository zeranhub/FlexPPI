"""Exploratory case discovery across exact target sequences and partner contexts.

No population performance claim is made from ranked cases. Raw labels and the
frozen, same-fold predictor are used without calibrating the observed effects.
"""
from pathlib import Path
import hashlib
import itertools
import json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parent
OUT = ROOT/'same_chain_context_candidates'
VARIANTS = ['full', 'intra', 'aa_sasa', 'aa_sasa_contacts', 'static47_sasa']


def rho(x, y):
    return float(spearmanr(x, y).statistic) if len(x) >= 4 and np.ptp(x) > 0 and np.ptp(y) > 0 else np.nan


def main():
    OUT.mkdir(exist_ok=True)
    cols = ['sample_id', 'residue_index', 'dataset', 'chain', 'author_residue_id', 'amino_acid',
            'valid_rmsf', 'valid_hchi', 'md_rmsf', 'md_hchi', 'pred_full_rmsf', 'pred_full_hchi',
            'sidechain_sasa_complex_A2', 'sidechain_sasa_isolated_A2', 'delta_sidechain_sasa_A2',
            'intra_ca_contacts_8A', 'partner_ca_contacts_8A', 'distance_to_chain_end', 'fold']
    d = pd.read_parquet(ROOT/'feature_baseline/features_all.parquet', columns=cols,
                        filters=[('dataset', '==', 'dynarepo')])
    predcols = ['sample_id', 'residue_index'] + [f'pred_{v}_{m}' for v in VARIANTS[1:] for m in ['rmsf', 'hchi']]
    p = pd.read_parquet(ROOT/'feature_baseline/oof_predictions_all.parquet', columns=predcols,
                        filters=[('dataset', '==', 'dynarepo')])
    d = d.merge(p, on=['sample_id', 'residue_index'], validate='one_to_one')
    inventory = pd.read_csv(ROOT/'feature_baseline/cohort_records.csv').set_index('sample_id')
    family = pd.read_csv('E:/Project/FlexPPI/experiments/_sources/family_folds/record_family_folds.csv').set_index('sample_id')
    clean = family.groupby('component_id').original_fold.nunique().eq(1)
    source = [json.loads(s) for s in Path('E:/Project/FlexPPI/experiments/core_corrected_pipeline/training_payload/splits/all_ready.jsonl').read_text().splitlines()]
    temperatures = {s['sample_id']: s.get('temperature') for s in source}
    chains = []
    seqgroups = {}
    for sid, record in d.groupby('sample_id', sort=False):
        seqs = {chain: ''.join(c.sort_values('residue_index').amino_acid) for chain, c in record.groupby('chain', sort=False)}
        for chain, data in record.groupby('chain', sort=False):
            seq = seqs[chain]
            if len(seq) < 50:
                continue
            others = sorted(seq2 for chain2, seq2 in seqs.items() if chain2 != chain)
            partner_signature = hashlib.sha256('|'.join(others).encode()).hexdigest()[:16]
            rec = dict(sample_id=sid, chain=chain, sequence=seq, length=len(seq),
                       n_chains=len(seqs), partner_signature=partner_signature,
                       fold=int(data.fold.iloc[0]), temperature=temperatures[sid],
                       pdb_id=inventory.loc[sid, 'pdb_id'], component_id=family.loc[sid, 'component_id'],
                       family_clean=bool(clean[family.loc[sid, 'component_id']]),
                       structure_path=inventory.loc[sid, 'structure_path'])
            chains.append(rec)
            seqgroups.setdefault(seq, []).append((rec, data.sort_values('residue_index').reset_index(drop=True)))
    pd.DataFrame(chains).to_csv(OUT/'all_chain_inventory.csv', index=False)
    results = []
    for sequence, states in seqgroups.items():
        if len(states) < 2:
            continue
        seqhash = hashlib.sha256(sequence.encode()).hexdigest()[:12]
        for (a, da), (b, db) in itertools.combinations(states, 2):
            if a['sample_id'] == b['sample_id'] or a['partner_signature'] == b['partner_signature']:
                continue
            if a['fold'] != b['fold'] or a['temperature'] != b['temperature']:
                continue
            meta = dict(target_sequence_sha256_prefix=seqhash, target_length=len(sequence),
                        state_a=a['sample_id'], chain_a=a['chain'], pdb_a=a['pdb_id'],
                        state_b=b['sample_id'], chain_b=b['chain'], pdb_b=b['pdb_id'],
                        fold=a['fold'], temperature=a['temperature'],
                        family_clean=a['family_clean'] and b['family_clean'],
                        component_a=a['component_id'], component_b=b['component_id'],
                        same_pdb=a['pdb_id'] == b['pdb_id'], n_chains_a=a['n_chains'], n_chains_b=b['n_chains'])
            for metric in ['rmsf', 'hchi']:
                good = da[f'valid_{metric}'].to_numpy() & db[f'valid_{metric}'].to_numpy()
                y = db[f'md_{metric}'].to_numpy()[good]-da[f'md_{metric}'].to_numpy()[good]
                if len(y) < 30:
                    continue
                strong = np.abs(y) >= .10
                row = dict(meta, metric=metric, n=len(y), true_abs_mean=float(np.abs(y).mean()),
                           true_mean=float(y.mean()), n_strong=int(strong.sum()),
                           n_true_increase=int((y >= .1).sum()), n_true_decrease=int((y <= -.1).sum()),
                           zero_mae=float(np.abs(y).mean()),
                           always_decrease_accuracy=float((y[strong] < 0).mean()) if strong.any() else np.nan)
                for variant in VARIANTS:
                    z = db[f'pred_{variant}_{metric}'].to_numpy()[good]-da[f'pred_{variant}_{metric}'].to_numpy()[good]
                    row[f'{variant}_mae'] = float(np.abs(z-y).mean())
                    row[f'{variant}_rho'] = rho(y, z)
                    row[f'{variant}_gain_over_zero'] = float(np.abs(y).mean()-np.abs(z-y).mean())
                    row[f'{variant}_direction_strong'] = float(((np.sign(z[strong]) == np.sign(y[strong])).astype(float) + .5*(z[strong] == 0)).mean()) if strong.any() else np.nan
                    row[f'{variant}_pred_mean'] = float(z.mean())
                row['full_gain_over_static47_sasa'] = row['static47_sasa_mae']-row['full_mae']
                row['full_gain_over_sasa_contacts'] = row['aa_sasa_contacts_mae']-row['full_mae']
                results.append(row)
    table = pd.DataFrame(results)
    table.to_csv(OUT/'all_samefold_exact_chain_context_pairs.csv', index=False)
    if len(table):
        eligible = table[(table.full_gain_over_zero > 0) & (table.full_gain_over_static47_sasa > 0) &
                         (table.full_gain_over_sasa_contacts > 0) & (table.full_rho >= .25)].copy()
        eligible = eligible.sort_values(['family_clean', 'full_gain_over_zero'], ascending=False)
        eligible.to_csv(OUT/'positive_delta_candidates.csv', index=False)
        print(eligible[['state_a', 'chain_a', 'state_b', 'chain_b', 'pdb_a', 'pdb_b', 'metric', 'target_length',
                        'family_clean', 'zero_mae', 'full_mae', 'static47_sasa_mae', 'aa_sasa_contacts_mae',
                        'full_rho', 'full_direction_strong', 'always_decrease_accuracy']].head(30).to_string(index=False), flush=True)
    (OUT/'audit.json').write_text(json.dumps(dict(status='complete', target_chain_count=len(chains),
        sequence_groups=len(seqgroups), state_pair_metrics=len(results),
        fixed_model_comparison='Exact target sequence, same OOF fold, same reported temperature, different partner sequence set',
        limitation='Exploratory screened cases; different PDBs can have different preparations and protocols. '
                   'OOF labels can share source or homologs. family_clean disclosed. No rank-based case is a population performance estimate.'), indent=2), encoding='utf-8')


if __name__ == '__main__':
    main()
