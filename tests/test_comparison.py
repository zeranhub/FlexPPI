import numpy as np
from types import SimpleNamespace

from flexppi.comparison import (
    _paired_comparison,
    _system_metrics,
    _top_fraction_metrics,
)
from flexppi.context_comparison import _aggregate_chains_to_complexes, _interface_mask
from flexppi.pymol_visualization import _chain_percentiles


def test_system_metrics_reward_correct_rank_and_scale():
    target = np.asarray([0.2, 0.4, 0.8, 1.6, 3.2])
    perfect = _system_metrics(target.copy(), target, absolute_scale=True)
    reversed_prediction = _system_metrics(target[::-1], target, absolute_scale=False)

    assert np.isclose(perfect["spearman"], 1.0)
    assert perfect["z_mae"] == 0.0
    assert perfect["absolute_mae_angstrom"] == 0.0
    assert np.isclose(reversed_prediction["spearman"], -1.0)
    assert "absolute_mae_angstrom" not in reversed_prediction


def test_top_fraction_precision_is_one_for_exact_order():
    values = np.arange(10, dtype=float)
    metrics = _top_fraction_metrics(values, values, fraction=0.2)
    assert metrics["auroc"] == 1.0
    assert metrics["average_precision"] == 1.0
    assert metrics["precision_at_k"] == 1.0


def test_paired_comparison_flips_lower_is_better_error_metrics():
    rows = [
        {
            "candidate": {"spearman": 0.8, "z_mae": 0.3},
            "reference": {"spearman": 0.6, "z_mae": 0.5},
        },
        {
            "candidate": {"spearman": 0.7, "z_mae": 0.4},
            "reference": {"spearman": 0.5, "z_mae": 0.6},
        },
    ]
    report = _paired_comparison(rows, "candidate", "reference", 100, 7)
    assert np.isclose(report["spearman"]["mean_advantage"], 0.2)
    assert np.isclose(report["z_mae"]["mean_advantage"], 0.2)
    assert report["z_mae"]["win_rate"] == 1.0


def test_interface_mask_requires_a_different_chain_within_cutoff():
    backbone = np.zeros((3, 4, 3), dtype=np.float32)
    backbone[:, 1] = np.asarray([[0.0, 0.0, 0.0], [20.0, 0.0, 0.0], [5.0, 0.0, 0.0]])
    backbone_mask = np.ones((3, 4), dtype=bool)
    structure = SimpleNamespace(
        num_residues=3,
        backbone=backbone,
        backbone_mask=backbone_mask,
        chain_ids=["A", "A", "B"],
    )
    assert _interface_mask(structure, 8.0).tolist() == [True, False, True]


def test_chain_metrics_are_averaged_before_complex_macro_evaluation():
    rows = [
        {"sample_id": "complex1", "residue_count": 10, "method": {"spearman": 0.9}},
        {"sample_id": "complex1", "residue_count": 100, "method": {"spearman": 0.3}},
    ]
    aggregated = _aggregate_chains_to_complexes(rows, ("method",))
    assert len(aggregated) == 1
    assert aggregated[0]["residue_count"] == 110
    assert np.isclose(aggregated[0]["method"]["spearman"], 0.6)


def test_visual_percentiles_are_normalized_independently_per_chain():
    values = np.asarray([1.0, 2.0, 3.0, 100.0, 200.0])
    percentiles = _chain_percentiles(values, ["A", "A", "A", "B", "B"])
    assert percentiles.tolist() == [0.0, 50.0, 100.0, 0.0, 100.0]
