import pytest
import torch

from flexppi.config import DataConfig, ExperimentConfig
from flexppi.losses import state_loss
from flexppi.model.network import FlexPPI
from flexppi.trainer import _fixed_state_validation_sets


class RecordDataset:
    def __init__(self, groups):
        self.records = []
        for index, group in enumerate(groups):
            record = type("Record", (), {})()
            record.sample_id = "sample_{}".format(index)
            record.metadata = {"cluster_id": group}
            self.records.append(record)

    def __len__(self):
        return len(self.records)


def test_two_metric_state_loss_selects_rmsf_and_chi_columns():
    output = {
        "mean": torch.ones(1, 3, 2),
        "log_variance": torch.zeros(1, 3, 2),
        "mask": torch.ones(1, 3, dtype=torch.bool),
    }
    batch = {
        "label_values": torch.ones(1, 3, 3),
        "label_variances": torch.zeros(1, 3, 3),
        "label_mask": torch.ones(1, 3, 3, dtype=torch.bool),
        "target_mask": torch.ones(1, 3, dtype=torch.bool),
    }
    loss = state_loss(
        output,
        batch,
        [1.0, 1.0],
        metric_names=["rmsf", "chi_entropy"],
    )
    assert torch.isfinite(loss)


def test_fixed_validation_rejects_shared_cluster(monkeypatch):
    train = RecordDataset(["cluster_a", "cluster_b"])
    validation = RecordDataset(["cluster_b"])
    monkeypatch.setattr("flexppi.trainer.StateDataset", lambda *args, **kwargs: validation)
    config = ExperimentConfig(
        data=DataConfig(validation_state_manifest="validation.jsonl")
    )
    with pytest.raises(ValueError, match="share biological groups"):
        _fixed_state_validation_sets(config, train)


def test_state_only_model_has_no_response_head():
    config = ExperimentConfig()
    config.model.enable_response_head = False
    model = FlexPPI(config.model)
    assert model.response_head is None
