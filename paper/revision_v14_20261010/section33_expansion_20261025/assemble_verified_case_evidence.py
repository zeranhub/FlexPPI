"""Review a frozen candidate using every eligible residue and matched pair.

Discovery is exploratory. All angles come from the same trajectory as labels.
No manuscript or source data is changed.
"""
from pathlib import Path
import itertools
import json
import sys
import os
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['OPENBLAS_NUM_THREADS'] = '1'
print('Loading numerical libraries', flush=True)
import numpy as np
import pandas as pd
import requests
from types import SimpleNamespace


def spearmanr(a, b):
    ra = pd.Series(np.asarray(a)).rank(method='average').to_numpy()
    rb = pd.Series(np.asarray(b)).rank(method='average').to_numpy()
    value = np.corrcoef(ra, rb)[0, 1] if len(ra) >= 3 and np.ptp(ra) and np.ptp(rb) else np.nan
    return SimpleNamespace(statistic=float(value))

ROOT = Path(__file__).resolve().parent
CORE = Path('E:/Project/FlexPPI/experiments/core_corrected_pipeline')
sys.path.insert(0, str(CORE/'training_payload/src'))
from flexppi.constants import THREE_TO_ONE
from flexppi.data.pairing import global_sequence_alignment

OUT = ROOT/'verified_case_evidence'
VARIANTS = ['full', 'intra', 'aa_sasa', 'aa_sasa_contacts', 'static47_sasa']


def read_residues(path):
    rows = {}
    for line in path.read_text(errors='replace').splitlines():
        if line.startswith('ENDMDL'):
            break
        if line[:6].strip() not in ['ATOM', 'HETATM'] or line[17:20].strip() not in THREE_TO_ONE:
            continue
        if line[16:17] not in [' ', 'A', '1']:
            continue
        key = (line[21:22].strip() or '_', line[22:26].strip()+line[26:27].strip())
        if key not in rows:
            rows[key] = dict(chain=key[0], author=key[1], aa=THREE_TO_ONE[line[17:20].strip()])
    return list(rows.values())


def number_mapping(sid, pid):
    path = OUT/(pid+'_public.pdb')
    if not path.exists():
        r = requests.get('https://files.rcsb.org/download/'+pid+'.pdb', timeout=50)
        r.raise_for_status()
        path.write_bytes(r.content)
    public = read_residues(path)
    local = read_residues(CORE/f'training_payload/dataset/structures/{sid}.pdb')
    rows, audit = [], []
    for chain in dict.fromkeys(r['chain'] for r in local):
        ll = [r for r in local if r['chain'] == chain]
        candidates = []
        for pc in dict.fromkeys(r['chain'] for r in public):
            pp = [r for r in public if r['chain'] == pc]
            a = global_sequence_alignment(''.join(r['aa'] for r in ll), ''.join(r['aa'] for r in pp))
            candidates.append((a.coverage, a.identity, a.score, pc == chain, pc, a, pp))
        best = max(candidates, key=lambda x: x[:4])
        coverage, identity, score, _, pc, alignment, pp = best
        assert coverage > .9 and identity > .99, (sid, chain, coverage, identity)
        audit.append(dict(sample_id=sid, local_chain=chain, public_chain=pc, coverage=coverage,
                          identity=identity, equal_best_public_chains=';'.join(c[4] for c in candidates if c[:3] == best[:3]),
                          rule='Sequence correspondence; same chain name breaks sequence-identical homomer ties'))
        for i, j in zip(alignment.indices_a, alignment.indices_b):
            rows.append(dict(sample_id=sid, chain=chain, author_residue_id=ll[i]['author'],
                             public_chain=pc, public_residue_id=pp[j]['author'], amino_acid=ll[i]['aa']))
    return pd.DataFrame(rows), audit


def main():
    OUT.mkdir(exist_ok=True)
    print('Reading frozen candidate tables', flush=True)
    sid = 'A003W'
    d = pd.read_parquet(ROOT/'feature_baseline/features_all.parquet', filters=[('sample_id', '==', sid)])
    p = pd.read_parquet(ROOT/'feature_baseline/oof_predictions_all.parquet', filters=[('sample_id', '==', sid)])
    additional = [c for c in p if c not in d and c.startswith('pred_')]
    d = d.merge(p[['sample_id', 'residue_index']+additional], on=['sample_id', 'residue_index'], validate='one_to_one')
    d['author_residue_id'] = d.author_residue_id.astype(str)
    occ = pd.read_csv(ROOT/f'verified_case_trajectories/{sid}/all_sidechains_occupancy.csv', dtype={'author_residue_id':str})
    windows = pd.read_csv(ROOT/f'verified_case_trajectories/{sid}/all_sidechains_per_window.csv')
    extra = [c for c in occ if c not in d]
    d = d.merge(occ[['residue_index']+extra], on='residue_index', validate='one_to_one', how='left')
    mapping, map_audit = number_mapping(sid, '2V5Q')
    d = d.merge(mapping, on=['sample_id', 'chain', 'author_residue_id', 'amino_acid'], how='left', validate='one_to_one')
    d.to_csv(OUT/'2V5Q_all_residues.csv', index=False)
    scope = d[d.chain.eq('A') & d.valid_hchi & d.amino_acid.ne('P')].copy()
    scopes = {'whole_chain_nonPro':scope, 'interface_nonPro':scope[scope.interface_ca8],
              'outside_nonPro':scope[~scope.interface_ca8]}
    metrics = []
    for name, g in scopes.items():
        aa_n = g.groupby('amino_acid').amino_acid.transform('size')
        same = g[aa_n >= 4]
        for variant in VARIANTS:
            y, z = g.md_hchi, g[f'pred_{variant}_hchi']
            center_y = same.md_hchi-same.groupby('amino_acid').md_hchi.transform('mean')
            center_z = same[f'pred_{variant}_hchi']-same.groupby('amino_acid')[f'pred_{variant}_hchi'].transform('mean')
            metrics.append(dict(scope=name, model=variant, n=len(g), mae=float(np.abs(y-z).mean()),
                                typed_n=len(same), typed_rho=float(spearmanr(center_y, center_z).statistic)))
    pd.DataFrame(metrics).to_csv(OUT/'2V5Q_whole_chain_metrics.csv', index=False)
    eligible = scope[scope.distance_to_chain_end >= 3].copy()
    eligible['sasa_rank'] = eligible.groupby('amino_acid').sidechain_sasa_complex_A2.rank(pct=True)
    wlookup = {i:g.sort_values('block') for i,g in windows.groupby('residue_index')}
    pairs = []
    for aa, g in eligible.groupby('amino_acid'):
        if (scope.amino_acid == aa).sum() < 4:
            continue
        for (_, a), (_, b) in itertools.combinations(g.iterrows(), 2):
            if abs(a.sidechain_sasa_complex_A2-b.sidechain_sasa_complex_A2)>10 or abs(a.sasa_rank-b.sasa_rank)>.25 or abs(a.md_rmsf-b.md_rmsf)>.1:
                continue
            if a.md_hchi > b.md_hchi:
                a, b = b, a
            wa, wb = wlookup[a.residue_index], wlookup[b.residue_index]
            row = dict(aa=aa, low_index=a.residue_index, high_index=b.residue_index,
                       low_label=aa+a.author_residue_id, high_label=aa+b.author_residue_id,
                       low_public_label=aa+str(a.public_residue_id), high_public_label=aa+str(b.public_residue_id),
                       true_gap=b.md_hchi-a.md_hchi, at_least_one_interface=bool(a.interface_ca8 or b.interface_ca8),
                       low_sasa=a.sidechain_sasa_complex_A2, high_sasa=b.sidechain_sasa_complex_A2,
                       low_rmsf=a.md_rmsf, high_rmsf=b.md_rmsf,
                       low_true_hchi=a.md_hchi, high_true_hchi=b.md_hchi)
            for variant in VARIANTS:
                row[variant+'_gap'] = b[f'pred_{variant}_hchi']-a[f'pred_{variant}_hchi']
                row[variant+'_site_mae'] = (abs(a[f'pred_{variant}_hchi']-a.md_hchi)+abs(b[f'pred_{variant}_hchi']-b.md_hchi))/2
                row[variant+'_low'] = a[f'pred_{variant}_hchi']
                row[variant+'_high'] = b[f'pred_{variant}_hchi']
            for col in ['hchi', 'fine_entropy_mean_per_chi', 'width_mean_per_chi', 'shifted60_joint_entropy']:
                gaps = wb[col].to_numpy()-wa[col].to_numpy()
                row[col+'_block_positive_count'] = int((gaps>0).sum())
                row[col+'_gap_mean'] = float(gaps.mean())
                row[col+'_gaps'] = ';'.join(f'{v:.6f}' for v in gaps)
            pairs.append(row)
    pairs = pd.DataFrame(pairs)
    pairs.to_csv(OUT/'2V5Q_all_similar_environment_pairs.csv', index=False)
    summaries = []
    for name, g in [('all', pairs), ('interface', pairs[pairs.at_least_one_interface]),
                    ('gap_ge_0p4', pairs[pairs.true_gap >= .4])]:
        for variant in VARIANTS:
            summaries.append(dict(scope=name, model=variant, n=len(g),
                                  direction=float(((g[variant+'_gap']>0)+.5*(g[variant+'_gap']==0)).mean()),
                                  gap_mae=float(np.abs(g[variant+'_gap']-g.true_gap).mean())))
    pd.DataFrame(summaries).to_csv(OUT/'2V5Q_all_pairs_summary.csv', index=False)
    special = pairs[((pairs.low_label=='E102') & (pairs.high_label=='E151')) |
                    ((pairs.low_label=='K105') & (pairs.high_label=='K162'))]
    special.to_csv(OUT/'2V5Q_two_disjoint_pairs_physical_check.csv', index=False)
    maps = [mapping]
    for sample in ['A00C7', 'A00DH']:
        mp, ma = number_mapping(sample, '3H6G')
        maps.append(mp)
        map_audit.extend(ma)
    pd.concat(maps).to_csv(OUT/'native_pdb_residue_mapping.csv', index=False)
    (OUT/'mapping_audit.json').write_text(json.dumps(map_audit, indent=2), encoding='utf-8')
    print(pd.DataFrame(metrics).to_string(index=False))
    print(pd.DataFrame(summaries).to_string(index=False))
    print(special.to_string(index=False))


if __name__ == '__main__':
    main()
