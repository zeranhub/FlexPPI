"""Alignment-free MD diagnostic for all prespecified local C-alpha pairs.

Distances are invariant to translation/rotation, unlike aligned RMSF. This is
one additional summary of the same trajectories, not independent evidence.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
OUT = ROOT/'glur6_independent_audit'


def pair_fluctuation(coords, ii, jj):
    result = np.empty((5, len(ii)), dtype=np.float32)
    for start in range(0, len(ii), 256):
        stop = min(start+256, len(ii))
        distance = np.linalg.norm(coords[:, ii[start:stop]]-coords[:, jj[start:stop]], axis=-1)
        for block in range(5):
            result[block, start:stop] = distance[block*1000:(block+1)*1000].std(axis=0)
    return result


def main():
    m = np.load(ROOT/'verified_case_trajectories/A00DH/ca_coordinates.npz')
    d = np.load(ROOT/'verified_case_trajectories/A00C7/ca_coordinates.npz')
    mono = m['coords_A']
    assert mono.shape == (5000, 383, 3)
    ii, jj = np.triu_indices(383, k=4)
    initial_distance = np.linalg.norm(mono[0, ii]-mono[0, jj], axis=1)
    within = initial_distance <= 12.
    ii, jj, initial_distance = ii[within], jj[within], initial_distance[within]
    msd = pair_fluctuation(mono, ii, jj)
    rows = []
    dimers = {}
    for chain in ['A', 'B']:
        mask = d['chain'] == chain
        assert np.array_equal(d['amino_acid'][mask], m['amino_acid'])
        dimers[chain] = pair_fluctuation(d['coords_A'][:, mask], ii, jj)
    dimers['AB_mean'] = (dimers['A']+dimers['B'])/2
    pairs = pd.DataFrame(dict(sequence_position_i=ii+1, sequence_position_j=jj+1,
                              initial_monomer_distance_A=initial_distance,
                              monomer_distance_sd_A=msd.mean(axis=0)))
    for chain, ds in dimers.items():
        delta = ds-msd
        pairs['dimer_'+chain+'_distance_sd_A'] = ds.mean(axis=0)
        pairs['delta_'+chain+'_distance_sd_A'] = delta.mean(axis=0)
        for block in range(5):
            rows.append(dict(target=chain, block=block+1, n_pairs=len(ii),
                             monomer_mean_distance_sd_A=float(msd[block].mean()),
                             dimer_mean_distance_sd_A=float(ds[block].mean()),
                             mean_delta_sd_A=float(delta[block].mean()),
                             fraction_pairs_lower=float((delta[block] < 0).mean())))
        rows.append(dict(target=chain, block='mean_of_5', n_pairs=len(ii),
                         monomer_mean_distance_sd_A=float(msd.mean()),
                         dimer_mean_distance_sd_A=float(ds.mean()),
                         mean_delta_sd_A=float(delta.mean()),
                         fraction_pairs_lower=float((delta.mean(axis=0) < 0).mean())))
    pairs.to_csv(OUT/'alignment_free_all_local_pairs.csv', index=False)
    pd.DataFrame(rows).to_csv(OUT/'alignment_free_local_pair_summary.csv', index=False)
    (OUT/'alignment_free_diagnostic_contract.json').write_text(json.dumps(dict(
        status='complete', n_pairs=len(ii), units='angstrom',
        pairs='All sequence-separated-by-at-least-four C-alpha pairs within 12 A in the first monomer frame; fixed without examining dynamic response or model predictions',
        quantity='Mean over five 100 ns blocks of within-block standard deviation of pairwise C-alpha distance',
        limitation='Same MD trajectories; no independent experiment, causal mechanism, or model prediction of this distance statistic. Not a replacement for RMSF evaluation.'
    ), indent=2), encoding='utf-8')
    print(pd.DataFrame(rows).to_string(index=False), flush=True)


if __name__ == '__main__':
    main()
