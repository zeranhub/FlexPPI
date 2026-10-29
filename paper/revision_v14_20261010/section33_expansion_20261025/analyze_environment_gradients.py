"""All-site, same-amino-acid environment associations; frozen labels and OOF.

The SASA grouping uses initial structures only. Shared-trajectory quantities
are related summaries, not independent experimental measurements.
"""
from pathlib import Path
import argparse
import json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, rankdata
from analyze_partner_cohort import cluster_summary

ROOT = Path(__file__).resolve().parent
FAM = Path('E:/Project/FlexPPI/experiments/_sources/family_folds/record_family_folds.csv')
UNIFORM = Path('E:/Project/FlexPPI/experiments/_sources/label_sensitivity/labels/uniform_rep1')


def corr(x, y):
    return float(spearmanr(x, y).statistic) if len(x) >= 4 and np.ptp(x) > 0 and np.ptp(y) > 0 else np.nan


def analyze(path, name, uniform=False):
    out = ROOT/'environment_gradients'/name
    out.mkdir(parents=True, exist_ok=True)
    columns = ['sample_id', 'residue_index', 'dataset', 'chain', 'author_residue_id', 'amino_acid',
               'valid_hchi', 'md_hchi', 'pred_full_hchi', 'pred_intra_hchi', 'pred_full_rmsf',
               'sidechain_sasa_complex_A2', 'sidechain_sasa_isolated_A2', 'delta_sidechain_sasa_A2',
               'sidechain_complete', 'chi_input_complete', 'nonstandard_residue',
               'intra_ca_contacts_8A', 'partner_ca_contacts_8A', 'interface_ca8']
    d = pd.read_parquet(path, columns=columns)
    initial_rows = len(d)
    d = d[d.valid_hchi & d.sidechain_complete & d.chi_input_complete & ~d.nonstandard_residue].copy()
    d = d[np.isfinite(d.sidechain_sasa_complex_A2)]
    if uniform:
        values, valid = {}, {}
        for sid, g in d.groupby('sample_id', sort=False):
            with np.load(UNIFORM/f'{sid}.npz') as z:
                idx = g.residue_index.to_numpy(int)
                values[sid] = pd.Series(z['values'][idx, 1], index=g.index)
                valid[sid] = pd.Series(z['mask'][idx, 1], index=g.index)
        d['md_hchi'] = pd.concat(values.values()).reindex(d.index)
        d = d[pd.concat(valid.values()).reindex(d.index).astype(bool)]
    d = d[np.isfinite(d.md_hchi) & np.isfinite(d.pred_full_hchi)]
    fam = pd.read_csv(FAM).set_index('sample_id')
    d['component_id'] = d.sample_id.map(fam.component_id)
    assert d.component_id.notna().all()
    family_clean = fam.groupby('component_id').original_fold.nunique().eq(1)
    d['family_clean'] = d.component_id.map(family_clean)
    all_valid = len(d)
    n = d.groupby(['sample_id', 'amino_acid']).sample_id.transform('size')
    d = d[n >= 4].copy()
    rows, gradient, type_rows, exception_rows = [], [], [], []
    for i, (sid, group) in enumerate(d.groupby('sample_id', sort=False)):
        record = group.copy()
        for source, target in [('sidechain_sasa_complex_A2', 'sasa_c'), ('md_hchi', 'truth_c'),
                               ('pred_full_hchi', 'full_c'), ('pred_intra_hchi', 'intra_c'),
                               ('pred_full_rmsf', 'rmsf_c')]:
            record[target] = record[source] - record.groupby('amino_acid')[source].transform('mean')
        meta = dict(sample_id=sid, dataset=group.dataset.iloc[0], component_id=group.component_id.iloc[0],
                    family_clean=bool(group.family_clean.iloc[0]), n_positions=len(group))
        record_types = []
        for aa, t in record.groupby('amino_acid', sort=False):
            sasa, truth = t.sidechain_sasa_complex_A2.to_numpy(), t.md_hchi.to_numpy()
            full = t.pred_full_hchi.to_numpy()
            raw_rank = rankdata(sasa, method='average')
            pct = (raw_rank-1)/(len(t)-1)
            classes = np.minimum((pct*4).astype(int), 3)
            eligible_gradient = np.ptp(sasa) > 0
            low, high = pct <= .25, pct >= .75
            buried, exposed = sasa <= 1., sasa >= 50.
            r = dict(meta, amino_acid=aa, n_type=len(t), rho_sasa_truth=corr(sasa, truth),
                     rho_full_truth=corr(full, truth),
                     high_minus_low_truth=float(truth[high].mean()-truth[low].mean())
                     if eligible_gradient and low.any() and high.any() else np.nan,
                     high_minus_low_full=float(full[high].mean()-full[low].mean())
                     if eligible_gradient and low.any() and high.any() else np.nan,
                     exposed_minus_buried_truth=float(truth[exposed].mean()-truth[buried].mean())
                     if exposed.any() and buried.any() else np.nan,
                     n_buried=int(buried.sum()), n_exposed=int(exposed.sum()))
            type_rows.append(r)
            record_types.append(r)
            if eligible_gradient:
                for quartile in range(4):
                    take = classes == quartile
                    if take.any():
                        gradient.append(dict(meta, amino_acid=aa, exposure_quartile=quartile+1,
                                             n_type_bin=int(take.sum()),
                                             mean_sasa=float(sasa[take].mean()),
                                             centered_md_hchi=float(t.truth_c.to_numpy()[take].mean()),
                                             centered_full_hchi=float(t.full_c.to_numpy()[take].mean())))
            hlo, hhi = np.quantile(truth, [.25, .75])
            if hlo < hhi:
                for state, mask in [('buried_high_hchi', buried & (truth >= hhi)),
                                    ('exposed_low_hchi', exposed & (truth <= hlo))]:
                    for pos in np.flatnonzero(mask):
                        v = t.iloc[pos]
                        exception_rows.append(dict(meta, amino_acid=aa, state=state,
                                                   residue_index=int(v.residue_index), chain=v.chain,
                                                   author_residue_id=v.author_residue_id,
                                                   sasa=float(sasa[pos]), md_hchi=float(truth[pos]),
                                                   pred_full_hchi=float(full[pos]),
                                                   intra_contact8=int(v.intra_ca_contacts_8A),
                                                   partner_contact8=int(v.partner_ca_contacts_8A)))
        tt = pd.DataFrame(record_types)
        rows.append(dict(meta, rho_sasa_truth_centered=corr(record.sasa_c, record.truth_c),
                         rho_full_truth_centered=corr(record.full_c, record.truth_c),
                         rho_intra_truth_centered=corr(record.intra_c, record.truth_c),
                         rho_rmsf_truth_centered=corr(record.rmsf_c, record.truth_c),
                         n_types=len(tt), **{k: float(tt[k].mean()) for k in
                             ['high_minus_low_truth', 'high_minus_low_full', 'exposed_minus_buried_truth']}))
        if (i+1) % 300 == 0:
            print('Environment', name, i+1, flush=True)
    per = pd.DataFrame(rows)
    per.to_csv(out/'per_record_environment_metrics.csv', index=False)
    pd.DataFrame(type_rows).to_csv(out/'per_record_type_environment_metrics.csv.gz', index=False, compression='gzip')
    gr = pd.DataFrame(gradient)
    gr.to_csv(out/'exposure_quartile_type_source.csv.gz', index=False, compression='gzip')
    gr_record = gr.groupby(['sample_id', 'dataset', 'component_id', 'exposure_quartile'], as_index=False)[
        ['centered_md_hchi', 'centered_full_hchi', 'mean_sasa']].mean()
    gr_record.to_csv(out/'exposure_quartile_record_source.csv', index=False)
    pd.DataFrame(exception_rows).to_csv(out/'all_structural_exceptions.csv.gz', index=False, compression='gzip')
    summaries = []
    metrics = ['rho_sasa_truth_centered', 'rho_full_truth_centered', 'rho_intra_truth_centered',
               'rho_rmsf_truth_centered', 'high_minus_low_truth', 'high_minus_low_full',
               'exposed_minus_buried_truth']
    for pop, pp in [('all', per), ('family_clean', per[per.family_clean])]:
        for source, ss in [('all', pp)] + list(pp.groupby('dataset')):
            for metric in metrics:
                summaries.append(dict(population=pop, dataset=source, metric=metric,
                                      **cluster_summary(ss, metric)))
    for quartile, ss in gr_record.groupby('exposure_quartile'):
        for metric in ['centered_md_hchi', 'centered_full_hchi']:
            summaries.append(dict(population='exposure_quartile', dataset='all', metric=metric,
                                  exposure_quartile=int(quartile), **cluster_summary(ss, metric)))
    pd.DataFrame(summaries).to_csv(out/'family_cluster_summary.csv', index=False)
    audit = dict(status='complete', cohort=name, input=str(path), rows_before_filter=initial_rows,
                 quality_and_label_valid=all_valid, same_type_n_ge4=len(d), records=len(per),
                 families=int(per.component_id.nunique()),
                 truth='uniform first 0-100ns window' if uniform else 'frozen original labels',
                 structure='initial complete protein assembly, not trajectory averaged',
                 quartiles='same record and AA, average rank; constant SASA types excluded from gradient',
                 buried='sidechain SASA <= 1 A^2', exposed='sidechain SASA >= 50 A^2',
                 limitation='Correlations are descriptive; not an independently calibrated SASA baseline. '
                            'Quartile, buried/exposed and exception selections do not use predictions.')
    (out/'audit.json').write_text(json.dumps(audit, indent=2), encoding='utf-8')
    print(json.dumps(audit, indent=2), flush=True)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--all', action='store_true')
    args = ap.parse_args()
    if args.all:
        analyze(ROOT/'feature_baseline/features_all.parquet', 'full3063', uniform=False)
    else:
        analyze(ROOT/'feature_baseline/features_priority83.parquet', 'uniform83', uniform=True)
