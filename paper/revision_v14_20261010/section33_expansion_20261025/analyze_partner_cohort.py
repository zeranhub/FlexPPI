"""Frozen OOF Full/Intra comparison by partner contact and matched chain context.

This is a comparison of separately trained predictors, not apo/holo response.
Run with the existing local scientific Python runtime. Original inputs are read-only.
"""
from pathlib import Path
import json
import os
os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('MKL_NUM_THREADS', '1')
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'partner_cohort'
CORE = Path('E:/Project/FlexPPI/experiments/core_corrected_pipeline')
FAMILY = Path('E:/Project/FlexPPI/experiments/_sources/family_folds/record_family_folds.csv')
SEED = 20261025


def rho(y, p):
    if len(y) < 4 or np.ptp(y) == 0 or np.ptp(p) == 0:
        return np.nan
    return float(spearmanr(y, p).statistic)


def typed_rho(y, p, aa):
    keep = np.zeros(len(y), bool)
    y, p = y.copy(), p.copy()
    for a in np.unique(aa):
        ix = aa == a
        if ix.sum() >= 4:
            keep |= ix
            y[ix] -= y[ix].mean()
            p[ix] -= p[ix].mean()
    return rho(y[keep], p[keep]), int(keep.sum())


def cluster_summary(frame, value):
    d = frame.dropna(subset=[value])
    if len(d) == 0:
        return dict(n_records=0, n_families=0, mean=None, ci_low=None, ci_high=None)
    groups = d.groupby('component_id')[value].agg(['sum', 'count', 'mean'])
    sums, counts = groups['sum'].to_numpy(), groups['count'].to_numpy()
    rng = np.random.default_rng(SEED)
    boot = []
    if len(groups) >= 2:
        for start in range(0, 10000, 128):
            idx = rng.integers(0, len(groups), size=(min(128, 10000-start), len(groups)))
            boot.extend((sums[idx].sum(1) / counts[idx].sum(1)).tolist())
    biggest = groups['count'].idxmax()
    reduced = d[d.component_id != biggest]
    return dict(n_records=len(d), n_families=len(groups), mean=float(d[value].mean()),
                ci_low=float(np.quantile(boot, .025)) if boot else None,
                ci_high=float(np.quantile(boot, .975)) if boot else None,
                positive_fraction=float((d[value] > 0).mean()),
                family_equal_mean=float(groups['mean'].mean()),
                largest_family_records=int(groups['count'].max()),
                excluding_largest_family_mean=float(reduced[value].mean()) if len(reduced) else None)


def matched_error_gain(y, full, intra, aa, chain, inner, interface):
    # Exact chain and type; contact bins 0–3,4–7,... . No response-derived selection.
    t = pd.DataFrame(dict(aa=aa, chain=chain, inner_bin=inner//4, interface=interface,
                          full_abs=np.abs(full-y), intra_abs=np.abs(intra-y)))
    rows = []
    for (ch, a, b), g in t.groupby(['chain', 'aa', 'inner_bin'], sort=False):
        near, away = g[g.interface], g[~g.interface]
        if len(near) and len(away):
            ng = float((near.intra_abs-near.full_abs).mean())
            ag = float((away.intra_abs-away.full_abs).mean())
            rows.append(dict(aa=a, interface_gain=ng, outside_gain=ag,
                             difference_in_gain=ng-ag, n_near=len(near), n_away=len(away)))
    if not rows:
        return None
    g = pd.DataFrame(rows)
    means = g.groupby('aa')[['interface_gain', 'outside_gain', 'difference_in_gain']].mean().mean()
    return dict(means, n_types=int(g.aa.nunique()), n_context_groups=len(g),
                n_positions=int(g.n_near.sum()+g.n_away.sum()))


def summarize_saved():
    d = pd.read_csv(OUT/'per_record_contact_metrics.csv')
    match = pd.read_csv(OUT/'matched_interface_outside_gains.csv')
    summaries = []
    for (metric, scope), group in d.groupby(['metric', 'scope']):
        for pop, population in [('all', group), ('multichain', group[group.n_chains > 1]),
                                ('family_clean', group[group.family_clean])]:
            for source, subset in [('all', population)] + list(population.groupby('dataset')):
                if len(subset):
                    for value in ['mae_gain', 'typed_rho_gain']:
                        summaries.append(dict(metric=metric, scope=scope, population=pop,
                                              dataset=source, statistic=value,
                                              full_mae=float(subset.full_mae.mean()),
                                              intra_mae=float(subset.intra_mae.mean()),
                                              **cluster_summary(subset, value)))
    for metric, group in match.groupby('metric'):
        for source, subset in [('all', group)] + list(group.groupby('dataset')):
            for value in ['interface_gain', 'outside_gain', 'difference_in_gain']:
                summaries.append(dict(metric=metric, scope='matched_chain_type_intra_contact',
                                      population='multichain', dataset=source, statistic=value,
                                      **cluster_summary(subset, value)))
    pd.DataFrame(summaries).to_csv(OUT/'family_cluster_summary.csv', index=False)


def main():
    OUT.mkdir(exist_ok=True)
    schema = json.loads((CORE/'static_corrected/training_schema.json').read_text())
    names = schema['names']
    with np.load(CORE/'static_corrected/training_features.npz') as z:
        x = z['X']
        inner = {r: np.rint(np.expm1(x[:, names.index(f'intra_contacts_{r}A')])).astype(np.int16)
                 for r in (6, 8, 10, 12)}
        partner = {r: np.rint(np.expm1(x[:, names.index(f'partner_contacts_{r}A')])).astype(np.int16)
                   for r in (6, 8, 10, 12)}
        del x
    family = pd.read_csv(FAMILY)
    assert family.sample_id.is_unique
    clean = family.groupby('component_id').original_fold.nunique().eq(1)
    fam = family.set_index('sample_id').component_id.to_dict()
    rows, matched = [], []
    total = np.zeros(2, int)
    multi_count = 0
    for index, record in enumerate(schema['records']):
        sid, start, stop = record['sample_id'], record['start'], record['stop']
        with np.load(CORE/f'neural_exports/evaluation/full/{sid}.npz') as z:
            y, p, valid = z['target'][:, :2], z['mean'][:, :2], z['valid'][:, :2]
            aa, chain = z['amino_acid'], z['chain']
        with np.load(CORE/f'neural_exports/evaluation/intrachain/{sid}.npz') as z:
            q = z['mean'][:, :2]
            assert np.array_equal(valid, z['valid'][:, :2]), sid
            assert np.allclose(y[valid], z['target'][:, :2][valid]), sid
            assert np.array_equal(aa, z['amino_acid']), sid
        assert len(y) == stop-start
        total += valid.sum(0)
        nchains = len(np.unique(chain))
        multi_count += nchains > 1
        meta = dict(sample_id=sid, dataset=record['dataset'], fold=record['fold'],
                    component_id=fam[sid], family_clean=bool(clean[fam[sid]]), n_chains=nchains)
        for m, metric in enumerate(['rmsf', 'hchi']):
            good = valid[:, m]
            scopes = {'all': np.ones(len(y), bool)}
            if nchains > 1:
                c8 = partner[8][start:stop]
                scopes.update(interface8=c8 > 0, outside8=c8 == 0,
                              partner_count_1_2=(c8 >= 1)&(c8 <= 2),
                              partner_count_3_5=(c8 >= 3)&(c8 <= 5), partner_count_ge6=c8 >= 6)
                for radius in (6, 10, 12):
                    scopes[f'interface{radius}'] = partner[radius][start:stop] > 0
            for scope, selected in scopes.items():
                ix = good & selected
                if not ix.any():
                    continue
                a, b, c = y[ix, m], p[ix, m], q[ix, m]
                trp, tn = typed_rho(a, b, aa[ix])
                trq, _ = typed_rho(a, c, aa[ix])
                fm, im = float(np.abs(b-a).mean()), float(np.abs(c-a).mean())
                rows.append(dict(meta, metric=metric, scope=scope, n_positions=int(ix.sum()),
                                 full_mae=fm, intra_mae=im, mae_gain=im-fm,
                                 full_rho=rho(a, b), intra_rho=rho(a, c),
                                 full_typed_rho=trp, intra_typed_rho=trq, n_typed=tn,
                                 typed_rho_gain=trp-trq))
            if nchains > 1:
                c8 = partner[8][start:stop]
                v = matched_error_gain(y[good, m], p[good, m], q[good, m], aa[good],
                                       chain[good], inner[8][start:stop][good], c8[good] > 0)
                if v:
                    matched.append(dict(meta, metric=metric, **v))
        if (index+1) % 300 == 0:
            print('Partner cohort', index+1, len(schema['records']), flush=True)
    d, match = pd.DataFrame(rows), pd.DataFrame(matched)
    d.to_csv(OUT/'per_record_contact_metrics.csv', index=False)
    match.to_csv(OUT/'matched_interface_outside_gains.csv', index=False)
    summarize_saved()
    audit = dict(status='complete', records=len(schema['records']), multichain_records=int(multi_count),
                 effective_rmsf=int(total[0]), effective_hchi=int(total[1]),
                 family_components=int(family.component_id.nunique()), bootstrap_repetitions=10000,
                 seed=SEED, unit='record means, families resampled as intact clusters',
                 interface_definition='other-chain CA count within 8A > 0; also 6,10,12A',
                 match='same record, chain, amino-acid; intra-CA8 count binned by 4; then AA and record equal',
                 limitation='Full and Intra separately trained. Gains are not a true binding response. '
                            'Cluster CI is not a family-isolated training experiment.')
    (OUT/'audit.json').write_text(json.dumps(audit, indent=2), encoding='utf-8')
    print(json.dumps(audit, indent=2), flush=True)


if __name__ == '__main__':
    import sys
    if '--summary-only' in sys.argv:
        summarize_saved()
    else:
        main()
