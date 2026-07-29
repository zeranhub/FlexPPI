import math

import numpy as np

from flexppi.data.labels import _shrake_rupley_residue_sasa, _sphere_points


def test_shrake_rupley_isolated_atom_has_full_sphere_area():
    radius = 1.70 + 1.40
    area = _shrake_rupley_residue_sasa(
        positions=np.asarray([[0.0, 0.0, 0.0]], dtype=np.float32),
        expanded_radii=np.asarray([radius], dtype=np.float32),
        atom_to_residue=np.asarray([0], dtype=np.int64),
        n_residues=1,
        sphere=_sphere_points(48),
    )
    assert np.isclose(area[0], 4.0 * math.pi * radius * radius, rtol=1e-6)


def test_shrake_rupley_overlapping_atoms_are_partly_occluded():
    radius = 1.70 + 1.40
    area = _shrake_rupley_residue_sasa(
        positions=np.asarray([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]], dtype=np.float32),
        expanded_radii=np.asarray([radius, radius], dtype=np.float32),
        atom_to_residue=np.asarray([0, 1], dtype=np.int64),
        n_residues=2,
        sphere=_sphere_points(96),
    )
    full_area = 4.0 * math.pi * radius * radius
    assert 0.0 < area[0] < full_area
    assert 0.0 < area[1] < full_area
