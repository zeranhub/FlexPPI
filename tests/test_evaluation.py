import numpy as np

from flexppi.data.records import remap_path
from flexppi.evaluation import (
    _direction_report,
    _regression_report,
    _retrieval_report,
)
from flexppi.tcr3d import _cdr_regions


def test_linux_manifest_path_can_be_remapped_to_local_prefix(tmp_path):
    mapped = remap_path(
        "/remote/data/labels/example.npz",
        {"/remote/data": str(tmp_path)},
    )
    assert mapped == str(tmp_path / "labels" / "example.npz")


def test_regression_reports_pooled_and_equal_system_metrics():
    target = np.asarray([0.0, 1.0, 2.0, 4.0])
    prediction = target + np.asarray([0.1, -0.1, 0.2, -0.2])
    groups = np.asarray(["a", "a", "b", "b"])
    report = _regression_report(prediction, target, groups)
    assert report["count"] == 4
    assert report["system_count"] == 2
    assert report["mae"] == np.mean(np.abs(prediction - target))
    assert report["macro_by_system"]["correlation_system_count"] == 2


def test_retrieval_and_dead_zone_direction_are_not_trivial_sign_metrics():
    target = np.asarray([-1.0, -0.01, 0.02, 2.0, 3.0])
    prediction = np.asarray([-0.8, 0.0, 0.0, 1.5, 2.5])
    groups = np.asarray(["a"] * len(target))
    direction = _direction_report(prediction, target, threshold=0.1)
    retrieval = _retrieval_report(prediction, target, groups, fraction=0.2)
    assert direction["active_count"] == 3
    assert direction["active_sign_accuracy"] == 1.0
    assert retrieval["macro_top_fraction_recall"] == 1.0


def test_tcr_cdr_sequences_are_mapped_without_using_residue_numbers():
    regions = _cdr_regions(
        "AAAACDEFGHIKLMNPQRSTVWYAAAA",
        {
            "cdr1_sequences": "CDEF",
            "cdr2_sequences": "KLMN",
            "cdr3_sequences": "STVW",
        },
    )
    assert regions.tolist().count(1) == 4
    assert regions.tolist().count(2) == 4
    assert regions.tolist().count(3) == 4
