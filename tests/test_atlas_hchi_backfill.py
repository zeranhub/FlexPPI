"""Check scientifically material invariants of the ATLAS H-chi backfill."""
import importlib.util
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

spec = importlib.util.spec_from_file_location(
    'atlas_hchi_backfill', Path(__file__).resolve().parents[1] / 'scripts/backfill_atlas_hchi.py'
)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def inputs():
    structure = SimpleNamespace(chain_ids=['_', '_'], author_residue_ids=['1', '2'],
                                residue_names=['SER', 'ALA'])
    values = np.asarray([[1.25, np.nan, 8], [2.5, np.nan, np.nan]], dtype=np.float32)
    variances = np.asarray([[.25, 0, .5], [.75, 0, 0]], dtype=np.float32)
    mask = np.asarray([[True, False, True], [True, False, False]])
    replicas = [[{'chain': 'SYSTEM', 'residue_number': '1', 'insertion_code': '',
                  'amino_acid_3': 'SER', 'n_chi_angles': 1, 'chi_entropy': value,
                  'n_frames': 1001}] for value in [0., .5, 1.]]
    return structure, replicas, values, variances, mask


def test_replica_mean_variance_and_undefined_chi_mask():
    s, rows, old_v, old_var, old_m = inputs()
    v, var, mask, reps, nchi, frames = module.merge_replica_rows(s, rows, old_v, old_var, old_m)
    assert v[0, 1] == .5 and var[0, 1] == .25
    assert mask[0, 1] and not mask[1, 1] and np.isnan(v[1, 1])
    for before, after in [(old_v, v), (old_var, var), (old_m, mask)]:
        np.testing.assert_array_equal(before[:, [0, 2]], after[:, [0, 2]])
    np.testing.assert_array_equal(frames, [1001] * 3)
    assert nchi.tolist() == [1, 0]


def test_duplicate_residue_rejected():
    args = list(inputs()); args[1][0].append(dict(args[1][0][0]))
    with pytest.raises(ValueError, match='duplicate'):
        module.merge_replica_rows(*args)


def test_missing_replica_not_silently_averaged():
    args = list(inputs()); args[1] = args[1][:2]
    with pytest.raises(ValueError, match='three replicas'):
        module.merge_replica_rows(*args)


def test_wrong_residue_identity_rejected():
    args = list(inputs()); args[1][0][0]['amino_acid_3'] = 'THR'
    with pytest.raises(ValueError, match='Unmapped'):
        module.merge_replica_rows(*args)


def test_entropy_above_state_space_bound_rejected():
    args = list(inputs()); args[1][0][0]['chi_entropy'] = 1.5
    with pytest.raises(ValueError, match='Invalid H-chi'):
        module.merge_replica_rows(*args)
