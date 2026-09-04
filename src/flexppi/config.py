"""Typed YAML configuration with conservative defaults."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Union

import yaml


@dataclass
class DataConfig:
    state_manifest: Optional[Union[str, List[str]]] = None
    validation_state_manifest: Optional[Union[str, List[str]]] = None
    test_state_manifest: Optional[Union[str, List[str]]] = None
    pair_manifest: Optional[str] = None
    max_residues: int = 768
    batch_size: int = 1
    num_workers: int = 0
    path_mappings: Dict[str, str] = field(default_factory=dict)


@dataclass
class ModelConfig:
    sequence_backend: str = "learned"
    esm_model_name: str = "facebook/esm2_t6_8M_UR50D"
    freeze_esm: bool = True
    single_dim: int = 128
    pair_dim: int = 48
    sequence_dim: int = 128
    ipa_blocks: int = 4
    ipa_heads: int = 8
    ipa_scalar_dim: int = 16
    ipa_qk_points: int = 4
    ipa_v_points: int = 8
    rbf_bins: int = 24
    dropout: float = 0.1
    use_condition_features: bool = True
    enable_response_head: bool = True
    metrics: List[str] = field(default_factory=lambda: ["rmsf", "chi_entropy"])


@dataclass
class TrainConfig:
    epochs: int = 20
    val_fraction: float = 0.1
    learning_rate: float = 3e-4
    weight_decay: float = 1e-4
    grad_clip: float = 1.0
    log_every: int = 25
    checkpoint_every: int = 1
    sampling_strategy: str = "uniform_records"
    samples_per_epoch: Optional[int] = None
    state_weight: float = 1.0
    response_weight: float = 1.0
    consistency_weight: float = 0.25
    metric_weights: Dict[str, float] = field(
        default_factory=lambda: {"rmsf": 1.0, "chi_entropy": 1.0}
    )


@dataclass
class ExperimentConfig:
    stage: str = "pretrain"
    seed: int = 7
    device: str = "auto"
    output_dir: str = "artifacts/run"
    init_checkpoint: Optional[str] = None
    data: DataConfig = field(default_factory=DataConfig)
    model: ModelConfig = field(default_factory=ModelConfig)
    train: TrainConfig = field(default_factory=TrainConfig)

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


def _construct_dataclass(cls: Any, value: Optional[Mapping[str, Any]]) -> Any:
    value = {} if value is None else dict(value)
    allowed = set(cls.__dataclass_fields__)
    unknown = sorted(set(value) - allowed)
    if unknown:
        raise ValueError("Unknown {} keys: {}".format(cls.__name__, ", ".join(unknown)))
    return cls(**value)


def load_config(path: str) -> ExperimentConfig:
    """Load a YAML file and reject misspelled keys instead of ignoring them."""

    config_path = Path(path)
    with config_path.open("r", encoding="utf-8") as handle:
        raw = yaml.safe_load(handle) or {}

    nested = dict(raw)
    data = _construct_dataclass(DataConfig, nested.pop("data", None))
    model = _construct_dataclass(ModelConfig, nested.pop("model", None))
    train = _construct_dataclass(TrainConfig, nested.pop("train", None))
    allowed = {"stage", "seed", "device", "output_dir", "init_checkpoint"}
    unknown = sorted(set(nested) - allowed)
    if unknown:
        raise ValueError("Unknown experiment keys: {}".format(", ".join(unknown)))
    return ExperimentConfig(data=data, model=model, train=train, **nested)
