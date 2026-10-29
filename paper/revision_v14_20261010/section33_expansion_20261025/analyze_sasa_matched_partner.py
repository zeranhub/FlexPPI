"""Exploratory robustness of interface gains at similar initial exposure.

Absolute sidechain SASA bins 5/10/20 A^2 were added after rank-SASA inspection;
they are not a prospective confirmatory analysis. No prediction selects groups.
"""
from pathlib import Path
import argparse
import json
import numpy as np
import pandas as pd
from analyze_partner_cohort import cluster_summary

ROOT = Path(__file__).resolve().parent
FAMILY = Path('E:/Project/FlexPPI/experiments/_sources/family_folds/record_family_folds.csv')


def main(priority=False):
    tag = 'uniform83_initial_structure' if priority else 'full3063'
    path = ROOT/'feature_baseline'/('features_priority83.parquet' if priority else 'features_all.parquet')
    out = ROOT/'partner_cohort'/f'sasa_matched_{tag}'
    out.mkdir(parents=True, exist_ok=True)
    cols = ['sample_id', 'dataset', 'chain', 'amino_acid', 'valid_rmsf', 'valid_hchi',
            'md_rmsf', 'md_hchi', 'pred_full_rmsf', 'pred_full_hchi', 'pred_intra_rmsf',
            'pred_intra_hchi', 'sidechain_complete', 'chi_input_complete', 'nonstandard_residue',
            'sidechain_sasa_complex_A2', 'intra_ca_contacts_8A', 'interface_ca8']
    d = pd.read_parquet(path, columns=cols)
    d = d[d.sidechain_complete & ~d.nonstandard_residue & np.isfinite(d.sidechain_sasa_complex_A2)].copy()
    d['contact_bin'] = d.intra_ca_contacts_8A // 4
    family = pd.read_csv(FAMILY).set_index('sample_id')
    clean = family.groupby('component_id').original_fold.nunique().eq(1)
    rows, rejected = [], []
    for width in (5, 10, 20):
        d['sasa_bin'] = np.floor(d.sidechain_sasa_complex_A2/width).astype(int)
        for sid, record in d.groupby('sample_id', sort=False):
            if not record.interface_ca8.any():
                continue
            meta = dict(sample_id=sid, dataset=record.dataset.iloc[0], component_id=family.loc[sid, 'component_id'],
                        family_clean=bool(clean[family.loc[sid, 'component_id']]), sasa_bin_A2=width)
            for metric in ('rmsf', 'hchi'):
                g = record[record[f'valid_{metric}']].copy()
                if metric == 'hchi':
                    g = g[g.chi_input_complete]
                g['gain'] = np.abs(g[f'pred_intra_{metric}']-g[f'md_{metric}'])-np.abs(g[f'pred_full_{metric}']-g[f'md_{metric}'])
                matched = []
                for (chain, aa, contact_bin, sasa_bin), block in g.groupby(
                        ['chain', 'amino_acid', 'contact_bin', 'sasa_bin'], sort=False):
                    near, away = block[block.interface_ca8], block[~block.interface_ca8]
                    if len(near) and len(away):
                        a, b = float(near.gain.mean()), float(away.gain.mean())
                        matched.append(dict(amino_acid=aa, interface_gain=a, outside_gain=b,
                                            difference_in_gain=a-b, n_near=len(near), n_away=len(away)))
                if matched:
                    t = pd.DataFrame(matched)
                    values = t.groupby('amino_acid')[['interface_gain', 'outside_gain', 'difference_in_gain']].mean().mean()
                    rows.append(dict(meta, metric=metric, **values, n_context_groups=len(t),
                                     n_types=int(t.amino_acid.nunique()), n_positions=int(t.n_near.sum()+t.n_away.sum())))
                else:
                    rejected.append(dict(meta, metric=metric, reason='No same-chain/AA/intra-contact/SASA bin with both interface and outside positions'))
    r = pd.DataFrame(rows)
    r.to_csv(out/'matched_per_record.csv', index=False)
    pd.DataFrame(rejected).to_csv(out/'no_qualifying_record.csv', index=False)
    summaries = []
    for (metric, width), g in r.groupby(['metric', 'sasa_bin_A2']):
        for pop, population in [('all', g), ('family_clean', g[g.family_clean])]:
            for source, subset in [('all', population)] + list(population.groupby('dataset')):
                for value in ('interface_gain', 'outside_gain', 'difference_in_gain'):
                    summaries.append(dict(metric=metric, sasa_bin_A2=width, population=pop, dataset=source,
                                          statistic=value, **cluster_summary(subset, value)))
    pd.DataFrame(summaries).to_csv(out/'family_cluster_summary.csv', index=False)
    (out/'audit.json').write_text(json.dumps(dict(status='complete', source=str(path),
        label_scope='original frozen labels even for 83 subset',
        grouping='same record, chain, AA, intra-CA8 count bin width4, SASA bin width5/10/20 A2',
        definition='abs error(Intra)-abs error(Full), and interface minus outside gain',
        prospective=False, reason='Absolute-SASA bins are exploratory robustness after SASA rank analysis',
        interpretation='Conditional observational error comparison of separately trained models, not causal binding response'), indent=2), encoding='utf-8')
    print('SASA matched partner complete', tag, len(r), flush=True)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--priority83', action='store_true')
    main(ap.parse_args().priority83)
