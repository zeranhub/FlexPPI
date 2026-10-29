"""Regression checks for a layout change that must not change the statistics."""
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd

SCRIPT = Path(__file__).resolve().parents[1] / "scripts/plot_atlas_hchi_v1_amino_acid_integrated.py"
SPEC = importlib.util.spec_from_file_location("integrated_plot", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_all_unique_observations_retained_and_raw_hchi_not_repeated():
    data, stats, missing, _ = MODULE.prepare_data()
    assert len(data) == 330
    assert len(stats) == 66
    assert len(data[data.block == "system_macro"]) == 75
    assert len(data[data.block == "within_aa_pooled"]) == 255
    raw = data[(data.metric == "chi_entropy") & (data.condition == "FlexPPI")]
    assert len(raw) == 15
    assert len(missing) == 3
    assert set(data.loc[data.block == "within_aa_pooled", "condition"]) == set(MODULE.AA_ORDER)


def test_means_sd_match_source_and_missing_strata_are_not_zero_filled():
    data, stats, missing, _ = MODULE.prepare_data()
    original = pd.read_csv(MODULE.DATA / "source_fold_metrics.csv")
    for row in stats.loc[stats.block == "system_macro"].itertuples():
        group = original[(original.fold >= 0) & (original.dataset == row.dataset)
                         & (original.metric == row.metric) & (original.predictor == row.condition)
                         & (original.scope == "macro_by_system")]
        assert np.isclose(row.mean, group.spearman.mean(), rtol=0, atol=1e-14)
        assert np.isclose(row.sd, group.spearman.std(ddof=1), rtol=0, atol=1e-14)
    for entry in missing:
        assert data[(data.block == "within_aa_pooled") & (data.dataset == entry["dataset"])
                    & (data.condition == entry["amino_acid"])].empty
    # Proline is a real near-zero result and must not be mistaken for missing data.
    assert len(data[(data.block == "within_aa_pooled") & (data.condition == "P")]) == 15
