"""The reviewed gate changes only RMSF RMSE tolerance, never correlations."""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import resume_atlas_hchi_v1_platform_audit as audit


def result(statistic, difference):
    return dict(passed=False, comparisons=[dict(metric="rmsf", statistic=statistic,
                                                cpu_minus_server=difference, passed=False)])


def test_review_records_original_gate(monkeypatch, tmp_path):
    monkeypatch.setattr(audit, "STRICT_VERIFY", lambda *a: result("rmse", .0002761))
    reviewed = audit.reviewed_verify(2, {}, tmp_path, tmp_path)
    assert reviewed["passed"]
    assert not reviewed["strict_2e_4_passed"]
    assert reviewed["comparisons"][0]["reviewed_tolerance"] == .001
    assert (tmp_path / "fold_2_platform_review.json").exists()


def test_correlation_gate_not_relaxed(monkeypatch, tmp_path):
    monkeypatch.setattr(audit, "STRICT_VERIFY", lambda *a: result("spearman", .0002761))
    with pytest.raises(ValueError, match="exceeds reviewed"):
        audit.reviewed_verify(2, {}, tmp_path, tmp_path)


def test_larger_rmse_difference_still_stops(monkeypatch, tmp_path):
    monkeypatch.setattr(audit, "STRICT_VERIFY", lambda *a: result("rmse", .0011))
    with pytest.raises(ValueError, match="exceeds reviewed"):
        audit.reviewed_verify(2, {}, tmp_path, tmp_path)
