"""Prevent v3 OOF runs from accepting changed splits or the wrong model protocol."""

import copy
import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import pytest
import torch


SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))
SPEC = importlib.util.spec_from_file_location("v3_oof_export", SCRIPTS / "export_v3_oof_predictions.py")
exporter = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(exporter)


@pytest.fixture
def frozen_split(tmp_path):
    split = tmp_path / "frozen"
    split.mkdir()
    names = ["all_ready.jsonl", "summary.json", "leakage_audit.json"]
    names += [f"fold_{fold}_{role}.jsonl" for fold in range(5) for role in ("train", "test")]
    for name in names:
        payload = {"passed": True, "violations": {"pdb_id": []}} if name == "leakage_audit.json" else {}
        (split / name).write_text(json.dumps(payload), encoding="utf-8")

    def write_checksums(omit=None):
        lines = [f"{hashlib.sha256((split / name).read_bytes()).hexdigest()}  {name}"
                 for name in names if name != omit]
        (split / "SHA256SUMS").write_text("\n".join(lines) + "\n", encoding="utf-8")

    write_checksums()
    return split, write_checksums


def test_changed_frozen_fold_is_rejected_before_membership_audit(frozen_split):
    split, _ = frozen_split
    exporter.audit_frozen_checksums(split)
    with (split / "fold_2_test.jsonl").open("a", encoding="utf-8") as handle:
        handle.write('\n{"sample_id":"unexpected"}')
    with pytest.raises(ValueError, match="checksum mismatch"):
        exporter.audit_frozen_checksums(split)


def test_removing_a_fold_checksum_does_not_bypass_freeze(frozen_split):
    split, write_checksums = frozen_split
    write_checksums(omit="fold_4_train.jsonl")
    with pytest.raises(ValueError, match="Missing frozen checksum entries"):
        exporter.audit_frozen_checksums(split)


@pytest.mark.parametrize("audit", [
    {"passed": False, "violations": {}},
    {"passed": True, "violations": {"pdb_id": ["1ABC"]}},
])
def test_valid_checksums_do_not_override_failed_leakage_audit(frozen_split, audit):
    split, write_checksums = frozen_split
    (split / "leakage_audit.json").write_text(json.dumps(audit), encoding="utf-8")
    write_checksums()
    with pytest.raises(ValueError, match="leakage audit does not pass"):
        exporter.audit_frozen_checksums(split)


@pytest.fixture
def checkpoint_bundle(tmp_path):
    folder = tmp_path / exporter.CHECKPOINT_DIRECTORY / "fold_0"
    folder.mkdir(parents=True)
    config = {
        "data": {
            "state_manifest": f"{exporter.SPLIT_DIRECTORY}/fold_0_train.jsonl",
            "test_state_manifest": f"{exporter.SPLIT_DIRECTORY}/fold_0_test.jsonl",
            "validation_state_manifest": None,
            "path_mappings": {"/original": "E:/original"},
        },
        "train": {"epochs": 40, "val_fraction": 0.0},
        "model": {"metrics": ["rmsf", "chi_entropy"], "use_condition_features": False,
                  "enable_response_head": False, "sequence_backend": "learned"},
    }

    def write(config_value=None, epoch=40, resolved=None):
        config_value = copy.deepcopy(config if config_value is None else config_value)
        resolved = copy.deepcopy(config_value if resolved is None else resolved)
        torch.save({"epoch": epoch, "config": config_value, "model": {}}, folder / "epoch_040.pt")
        (folder / "resolved_config.json").write_text(json.dumps(resolved), encoding="utf-8")

    write()
    return tmp_path, folder, config, write


def load_validation(root, mappings=None):
    return exporter.load_fold(root, 0, torch.device("cpu"), False, mappings or {})


def test_final_checkpoint_metadata_is_checked_even_if_filename_says_epoch40(checkpoint_bundle):
    root, _, _, write = checkpoint_bundle
    write(epoch=39)
    with pytest.raises(ValueError, match="not epoch 40"):
        load_validation(root)


def test_resolved_config_cannot_disguise_checkpoint_training_settings(checkpoint_bundle):
    root, _, config, write = checkpoint_bundle
    resolved = copy.deepcopy(config)
    resolved["train"]["epochs"] = 41
    write(resolved=resolved)
    with pytest.raises(ValueError, match="configs disagree"):
        load_validation(root)


@pytest.mark.parametrize("field,bad_path", [
    ("state_manifest", f"{exporter.SPLIT_DIRECTORY}/fold_1_train.jsonl"),
    ("test_state_manifest", f"{exporter.SPLIT_DIRECTORY}/fold_1_test.jsonl"),
    ("state_manifest", "data/splits/state_rmsf_hchi_v2_performance_first/fold_0_train.jsonl"),
])
def test_self_consistent_checkpoint_for_wrong_fold_or_version_is_rejected(
        checkpoint_bundle, field, bad_path):
    root, _, config, write = checkpoint_bundle
    config["data"][field] = bad_path
    write(config)
    with pytest.raises(ValueError, match="unexpected " + field):
        load_validation(root)


@pytest.mark.parametrize("section,field,bad_value", [
    ("data", "validation_state_manifest", "validation.jsonl"),
    ("train", "val_fraction", 0.2),
    ("train", "epochs", 60),
])
def test_epoch40_from_a_different_training_budget_is_not_the_frozen_protocol(
        checkpoint_bundle, section, field, bad_value):
    root, _, config, write = checkpoint_bundle
    config[section][field] = bad_value
    write(config)
    with pytest.raises(ValueError, match="no-validation protocol"):
        load_validation(root)


@pytest.mark.parametrize("field,bad_value", [
    ("metrics", ["chi_entropy", "rmsf"]),
    ("use_condition_features", True),
    ("enable_response_head", True),
])
def test_same_fold_with_different_prediction_task_is_rejected(checkpoint_bundle, field, bad_value):
    root, _, config, write = checkpoint_bundle
    config["model"][field] = bad_value
    write(config)
    with pytest.raises(ValueError, match="structure-only two-output model"):
        load_validation(root)


def test_remote_prefix_and_local_mapping_preserve_the_frozen_checkpoint(checkpoint_bundle):
    root, folder, config, write = checkpoint_bundle
    for field in ("state_manifest", "test_state_manifest"):
        config["data"][field] = "/ranze/repos/FlexPPI/" + config["data"][field]
    write(config)
    before = hashlib.sha256((folder / "epoch_040.pt").read_bytes()).hexdigest()
    model, resolved, path, digest = load_validation(root, {"/ranze/data": "D:/data"})
    assert model is None
    assert path == folder / "epoch_040.pt"
    assert digest == before
    assert resolved["data"]["path_mappings"] == {
        "/original": "E:/original", "/ranze/data": "D:/data"}
    saved = json.loads((folder / "resolved_config.json").read_text(encoding="utf-8"))
    assert saved == config
    assert hashlib.sha256(path.read_bytes()).hexdigest() == before
