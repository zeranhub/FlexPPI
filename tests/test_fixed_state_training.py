import pytest
import torch

from flexppi.config import DataConfig, ExperimentConfig
from flexppi.losses import heteroscedastic_nll, state_loss
from flexppi.model.network import FlexPPI
from flexppi.trainer import _fixed_state_validation_sets, _source_system_balanced_weights


class RecordDataset:
    def __init__(self, groups):
        self.records = []
        for index, group in enumerate(groups):
            record = type("Record", (), {})()
            record.sample_id = "sample_{}".format(index)
            record.metadata = {"cluster_id": group, "biological_system": group}
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
    config = ExperimentConfig(data=DataConfig(validation_state_manifest="validation.jsonl"))
    with pytest.raises(ValueError, match="share biological groups"):
        _fixed_state_validation_sets(config, train)


def test_state_only_model_has_no_response_head():
    config = ExperimentConfig()
    config.model.enable_response_head = False
    model = FlexPPI(config.model)
    assert model.response_head is None


def test_state_loss_averages_each_system_before_batch():
    predicted = torch.zeros(2, 100, 1)
    target = torch.zeros_like(predicted)
    target[0, :, 0] = 1.0
    target[1, 0, 0] = 3.0
    valid = torch.zeros_like(predicted, dtype=torch.bool)
    valid[0, :, 0] = True
    valid[1, 0, 0] = True
    loss = heteroscedastic_nll(
        predicted,
        torch.zeros_like(predicted),
        target,
        torch.zeros_like(predicted),
        valid,
        [1.0],
    )
    assert loss.item() == pytest.approx(2.5)


def test_source_system_weights_do_not_reward_duplicate_records():
    dataset = RecordDataset(["system_a", "system_a", "system_b", "system_c"])
    dataset.records[0].dataset = "source_one"
    dataset.records[1].dataset = "source_one"
    dataset.records[2].dataset = "source_one"
    dataset.records[3].dataset = "source_two"
    weights = _source_system_balanced_weights(dataset)
    assert weights == pytest.approx([0.125, 0.125, 0.25, 0.5])


def test_structure_only_model_disables_condition_path():
    config = ExperimentConfig()
    config.model.use_condition_features = False
    model = FlexPPI(config.model)
    assert model.condition_mlp is None
