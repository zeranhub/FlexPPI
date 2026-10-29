"""Tests for leakage-free priors and descriptive uncertainty diagnostics."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from analyze_atlas_hchi_v1_oof import attach_training_priors, regression, uncertainty_curves


def test_training_priors_exclude_held_out_fold():
    data = pd.DataFrame(dict(fold=[0, 1, 2, 3, 4], amino_acid=["K"] * 5,
                             target_rmsf=[100., 1., 1., 1., 1.],
                             target_chi_entropy=[100., 1., 1., 1., 1.],
                             valid_rmsf=True, valid_chi_entropy=True))
    attach_training_priors(data)
    assert data.loc[0, "aa_prior_rmsf"] == 1.
    assert data.loc[1, "aa_prior_chi_entropy"] == 25.75


def test_residual_correction_preserves_error():
    p, y, prior = np.array([1., 4., 7.]), np.array([2., 5., 8.]), np.array([0., 2., 6.])
    assert regression(p, y)["mae"] == regression(p - prior, y - prior)["mae"]
    assert regression(p, y)["spearman"] == 1.


def test_uncertainty_risk_uses_all_at_full_coverage():
    risk, bins, coverage = uncertainty_curves([0., 1., 2., 3.], [0., 0., 0., 0.],
                                             [1., 2., 3., 4.], [1., 1., 1., 1.])
    assert risk[-1]["retained_count"] == 4
    assert risk[-1]["relative_mae"] == 1.
    assert sum(b["count"] for b in bins) == 4
    model = [r["observed_coverage"] for r in coverage if r["variance_kind"] == "model_only"]
    total = [r["observed_coverage"] for r in coverage if r["variance_kind"] == "model_plus_label"]
    assert np.all(np.diff(model) >= 0)
    assert np.all(np.array(total) >= model)


def test_nonfinite_values_not_silently_dropped():
    with pytest.raises(ValueError):
        regression([1., np.nan], [1., 2.])
    with pytest.raises(ValueError):
        uncertainty_curves([1.], [2.], [-1.], [0.])
