"""Guard new-run identity, safe loading and smoke coverage."""
import copy
import importlib.util
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pytest
import torch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
spec = importlib.util.spec_from_file_location("atlas_oof", ROOT / "scripts/export_atlas_hchi_v1_oof_predictions.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def config():
    return dict(stage="pretrain", init_checkpoint=None,
                data=dict(state_manifest=module.SPLIT_DIRECTORY + "/fold_0_train.jsonl",
                          test_state_manifest=module.SPLIT_DIRECTORY + "/fold_0_test.jsonl",
                          validation_state_manifest=None),
                train=dict(epochs=40, val_fraction=0),
                model=dict(metrics=["rmsf", "chi_entropy"], sequence_backend="learned",
                           use_condition_features=False, enable_response_head=False))


def test_new_split_configuration_and_wrong_fold():
    c = config()
    module.validate_configuration(dict(epoch=40, config=c), c, 0)
    with pytest.raises(ValueError, match="fold/split"):
        module.validate_configuration(dict(epoch=40, config=c), c, 1)


def test_old_checkpoint_or_epoch_cannot_pass():
    c = config()
    with pytest.raises(ValueError, match="epoch 40"):
        module.validate_configuration(dict(epoch=39, config=c), c, 0)
    old = copy.deepcopy(c)
    old["data"]["state_manifest"] = "data/splits/state_rmsf_hchi_v2_performance_first/fold_0_train.jsonl"
    with pytest.raises(ValueError, match="fold/split"):
        module.validate_configuration(dict(epoch=40, config=old), old, 0)


def test_standard_numpy_rng_safe_load(tmp_path):
    path = tmp_path / "rng.pt"
    torch.save(dict(model={"w": torch.ones(2)}, numpy_rng=np.random.get_state()), path)
    result = module.safe_load_checkpoint(path)
    assert result["model"]["w"].tolist() == [1, 1]


def test_unexpected_checkpoint_global_rejected(tmp_path):
    path = tmp_path / "unsupported.pt"
    torch.save(dict(unexpected=datetime(2026, 1, 1)), path)
    with pytest.raises(ValueError, match="Unexpected checkpoint globals"):
        module.safe_load_checkpoint(path)


def test_smoke_covers_sources_and_full_preserves_order():
    records = [(str(i), dict(dataset=s)) for i, s in enumerate(["atlas", "atlas", "dynarepo", "itrace"])]
    assert module.select_records(records, None) == records
    assert [r[1]["dataset"] for r in module.select_records(records, 3)] == ["atlas", "dynarepo", "itrace"]


def test_protocol_and_label_counts_are_new():
    assert module.EXPECTED_COUNTS == dict(rmsf=959347, chi_entropy=798385)
    assert module.PROTOCOL["name"] != module.shared.PROTOCOL["name"]
    assert "atlas_hchi_v1" in module.CHECKPOINT_DIRECTORY
