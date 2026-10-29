"""OOF cache and split safeguards, independent of real checkpoints and datasets."""

import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

SPEC = importlib.util.spec_from_file_location(
    "v2_oof_export", Path(__file__).resolve().parents[1] / "scripts/export_v2_oof_predictions.py"
)
exporter = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(exporter)


def fixture_inventory(tmp_path):
    records = [{"sample_id": f"sample_{i}", "dataset": "atlas",
                "metadata": {"fold": i % 5, "biological_system": f"pdb:{i}",
                             "cluster_id": f"cluster_{i}"}} for i in range(10)]
    def write(name, values):
        (tmp_path / name).write_text("\n".join(json.dumps(r) for r in values), encoding="utf-8")
    write("all_ready.jsonl", records)
    for fold in range(5):
        write(f"fold_{fold}_train.jsonl", [r for r in records if r["metadata"]["fold"] != fold])
        write(f"fold_{fold}_test.jsonl", [r for r in records if r["metadata"]["fold"] == fold])
    return records, write


def test_inventory_all_samples_once(tmp_path):
    fixture_inventory(tmp_path)
    audit = exporter.audit_inventory(tmp_path, expected_count=10)
    assert len(audit["assignments"]) == 10
    assert all(len(fold) == 2 for fold in audit["folds"].values())


def test_inventory_rejects_overlap(tmp_path):
    records, write = fixture_inventory(tmp_path)
    write("fold_0_test.jsonl", [records[0], records[1], records[5]])
    with pytest.raises(ValueError, match="overlap"):
        exporter.audit_inventory(tmp_path, expected_count=10)


def test_inventory_rejects_missing_and_changed_records(tmp_path):
    records, write = fixture_inventory(tmp_path)
    write("fold_0_test.jsonl", [records[0]])
    with pytest.raises(ValueError, match="Incomplete"):
        exporter.audit_inventory(tmp_path, expected_count=10)
    records[0]["dataset"] = "changed"
    write("fold_0_test.jsonl", [records[0], records[5]])
    with pytest.raises(ValueError, match="differs from all_ready"):
        exporter.audit_inventory(tmp_path, expected_count=10)


def fixture_arrays(identity):
    n = 3
    arrays = {key: np.zeros(n, dtype=np.int64) for key in ("residue_index", "aatype")}
    for key in ("chain", "author_residue_id", "amino_acid"):
        arrays[key] = np.asarray(["A", "B", "C"])
    for key in ("target_mask", "geometry_mask", "output_mask"):
        arrays[key] = np.ones(n, dtype=bool)
    for key in ("label_mask", "valid"):
        arrays[key] = np.ones((n, 2), dtype=bool)
    for key in ("target", "label_variance", "mean", "model_variance"):
        arrays[key] = np.ones((n, 2), dtype=np.float32)
    arrays["metric_names"] = np.asarray(exporter.METRICS)
    arrays["metadata_json"] = np.asarray(json.dumps(
        {"identity": identity, "status": "complete", "residue_count": n}))
    return arrays


def test_cache_roundtrip_and_identity_refusal(tmp_path):
    identity = {"checkpoint_sha256": "model_a", "label_sha256": "label_a"}
    arrays = fixture_arrays(identity)
    destination = tmp_path / "cache.npz"
    exporter.save_cache(destination, arrays, identity)
    loaded = exporter.load_cache(destination, identity)
    assert np.array_equal(loaded["target"], arrays["target"])
    with pytest.raises(ValueError, match="identity/protocol"):
        exporter.load_cache(destination, {**identity, "checkpoint_sha256": "model_b"})
    with pytest.raises(ValueError, match="identity/protocol"):
        exporter.load_cache(destination, {**identity, "label_sha256": "label_b"})
    assert not list(tmp_path.glob("*.tmp"))


@pytest.mark.parametrize("field", ["mean", "target", "model_variance", "label_variance"])
def test_cache_rejects_nonfinite_values(field):
    identity = {"id": "test"}
    arrays = fixture_arrays(identity)
    arrays[field][0, 0] = np.nan
    with pytest.raises(ValueError, match="Non-finite"):
        exporter.validate_cache_arrays(arrays, identity)


def test_cache_rejects_mask_and_shape_drift():
    identity = {"id": "test"}
    arrays = fixture_arrays(identity)
    arrays["valid"][0, 0] = False
    with pytest.raises(ValueError, match="masks disagree"):
        exporter.validate_cache_arrays(arrays, identity)
    arrays["mean"] = np.zeros((3, 3))
    with pytest.raises(ValueError, match="Malformed"):
        exporter.validate_cache_arrays(arrays, identity)


def test_path_mapping_prefers_existing_and_supports_reverse(tmp_path):
    actual = tmp_path / "source" / "a.npz"
    actual.parent.mkdir()
    actual.write_bytes(b"labels")
    mappings = {str(actual.parent): "/remote/labels"}
    assert exporter.resolve_existing("/remote/labels/a.npz", mappings) == actual.resolve()
    assert exporter.resolve_existing(str(actual), {}) == actual.resolve()
    with pytest.raises(FileNotFoundError):
        exporter.resolve_existing("/remote/missing.npz", mappings)


def test_cache_filename_cannot_escape_output():
    name = exporter.cache_name("../../outside/file")
    assert "/" not in name and "\\" not in name
    assert name.endswith(".npz")
