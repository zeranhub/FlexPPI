"""PyTorch datasets for absolute states and paired responses."""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Union

import numpy as np
import torch
from torch.utils.data import Dataset

from flexppi.constants import METRIC_NAMES, UNKNOWN_INDEX
from flexppi.data.labels import load_label_file
from flexppi.data.records import (
    PairRecord,
    StateRecord,
    load_pair_records,
    load_state_records,
    remap_pair_record,
    remap_state_record,
)
from flexppi.data.structures import parse_pdb, target_mask


def _record_has_labels(record: StateRecord) -> bool:
    return bool(record.label_path and Path(record.label_path).exists())


def _fits_limit(record: StateRecord, max_residues: Optional[int]) -> bool:
    if max_residues is None:
        return True
    count = int(record.metadata.get("residue_count", 0))
    return count == 0 or count <= max_residues


def _condition_vector(record: StateRecord) -> np.ndarray:
    dataset_code = 0.0 if record.dataset.lower() == "atlas" else 1.0
    return np.asarray(
        [
            (float(record.temperature) - 300.0) / 20.0,
            math.log1p(max(float(record.window_ns), 0.0)) / math.log1p(500.0),
            dataset_code,
        ],
        dtype=np.float32,
    )


def load_state_example(record: StateRecord) -> Dict[str, Any]:
    structure = parse_pdb(record.structure_path)
    if not _record_has_labels(record):
        raise FileNotFoundError("Missing labels for {}".format(record.sample_id))
    values, variances, label_mask = load_label_file(str(record.label_path))
    if values.shape[0] != structure.num_residues:
        raise ValueError(
            "Structure/label length mismatch for {}: {} versus {}".format(
                record.sample_id, structure.num_residues, values.shape[0]
            )
        )
    geometry_mask = structure.frame_mask
    target = target_mask(structure, record.target_chains)
    label_mask = label_mask & geometry_mask[:, None]
    return {
        "sample_id": record.sample_id,
        "dataset": record.dataset,
        "source_group": record.source_group,
        "metadata": record.metadata,
        "aatype": torch.from_numpy(structure.aatype),
        "residue_index": torch.from_numpy(structure.residue_index),
        "chain_index": torch.from_numpy(structure.chain_index),
        "backbone": torch.from_numpy(structure.backbone),
        "backbone_mask": torch.from_numpy(structure.backbone_mask),
        "geometry_mask": torch.from_numpy(geometry_mask),
        "target_mask": torch.from_numpy(target),
        "condition": torch.from_numpy(_condition_vector(record)),
        "label_values": torch.from_numpy(values),
        "label_variances": torch.from_numpy(variances),
        "label_mask": torch.from_numpy(label_mask),
        "metric_names": METRIC_NAMES,
    }


class StateDataset(Dataset):
    def __init__(
        self,
        manifest_path: Union[str, Sequence[str]],
        max_residues: Optional[int] = None,
        path_mappings: Optional[Mapping[str, str]] = None,
    ):
        manifest_paths = [manifest_path] if isinstance(manifest_path, str) else list(manifest_path)
        records: List[StateRecord] = []
        for path in manifest_paths:
            records.extend(
                remap_state_record(record, path_mappings)
                for record in load_state_records(path)
            )
        self.records = [
            record for record in records
            if _record_has_labels(record) and _fits_limit(record, max_residues)
        ]

    def __len__(self) -> int:
        return len(self.records)

    def __getitem__(self, index: int) -> Dict[str, Any]:
        return load_state_example(self.records[index])


class PairDataset(Dataset):
    def __init__(
        self,
        manifest_path: str,
        max_residues: Optional[int] = None,
        path_mappings: Optional[Mapping[str, str]] = None,
    ):
        records = [
            remap_pair_record(record, path_mappings)
            for record in load_pair_records(manifest_path)
        ]
        self.records = [
            record for record in records
            if _record_has_labels(record.apo)
            and _record_has_labels(record.holo)
            and _fits_limit(record.apo, max_residues)
            and _fits_limit(record.holo, max_residues)
        ]

    def __len__(self) -> int:
        return len(self.records)

    def __getitem__(self, index: int) -> Dict[str, Any]:
        record = self.records[index]
        return {
            "pair_id": record.pair_id,
            "apo": load_state_example(record.apo),
            "holo": load_state_example(record.holo),
            "apo_indices": torch.tensor(record.apo_indices, dtype=torch.long),
            "holo_indices": torch.tensor(record.holo_indices, dtype=torch.long),
            "sequence_identity": record.sequence_identity,
            "coverage": record.coverage,
            "metadata": record.metadata,
        }


def collate_states(examples: Sequence[Dict[str, Any]]) -> Dict[str, Any]:
    if not examples:
        raise ValueError("Cannot collate an empty batch")
    batch_size = len(examples)
    max_length = max(int(example["aatype"].shape[0]) for example in examples)
    metric_count = len(METRIC_NAMES)
    output: Dict[str, Any] = {
        "sample_id": [example["sample_id"] for example in examples],
        "dataset": [example["dataset"] for example in examples],
        "source_group": [example.get("source_group", "unknown") for example in examples],
        "metadata": [example.get("metadata", {}) for example in examples],
        "aatype": torch.full((batch_size, max_length), UNKNOWN_INDEX, dtype=torch.long),
        "residue_index": torch.zeros((batch_size, max_length), dtype=torch.long),
        "chain_index": torch.zeros((batch_size, max_length), dtype=torch.long),
        "backbone": torch.zeros((batch_size, max_length, 4, 3), dtype=torch.float32),
        "backbone_mask": torch.zeros((batch_size, max_length, 4), dtype=torch.bool),
        "geometry_mask": torch.zeros((batch_size, max_length), dtype=torch.bool),
        "target_mask": torch.zeros((batch_size, max_length), dtype=torch.bool),
        "padding_mask": torch.zeros((batch_size, max_length), dtype=torch.bool),
        "condition": torch.stack([example["condition"] for example in examples]),
        "label_values": torch.full(
            (batch_size, max_length, metric_count), float("nan"), dtype=torch.float32
        ),
        "label_variances": torch.zeros(
            (batch_size, max_length, metric_count), dtype=torch.float32
        ),
        "label_mask": torch.zeros(
            (batch_size, max_length, metric_count), dtype=torch.bool
        ),
    }
    tensor_keys = [
        "aatype", "residue_index", "chain_index", "backbone", "backbone_mask",
        "geometry_mask", "target_mask", "label_values", "label_variances", "label_mask",
    ]
    for batch_index, example in enumerate(examples):
        length = int(example["aatype"].shape[0])
        output["padding_mask"][batch_index, :length] = True
        for key in tensor_keys:
            output[key][batch_index, :length] = example[key]
    return output


def collate_pairs(examples: Sequence[Dict[str, Any]]) -> Dict[str, Any]:
    if not examples:
        raise ValueError("Cannot collate an empty pair batch")
    apo = collate_states([example["apo"] for example in examples])
    holo = collate_states([example["holo"] for example in examples])
    max_pairs = max(int(example["apo_indices"].numel()) for example in examples)
    apo_indices = torch.zeros((len(examples), max_pairs), dtype=torch.long)
    holo_indices = torch.zeros((len(examples), max_pairs), dtype=torch.long)
    pair_mask = torch.zeros((len(examples), max_pairs), dtype=torch.bool)
    for batch_index, example in enumerate(examples):
        length = int(example["apo_indices"].numel())
        apo_indices[batch_index, :length] = example["apo_indices"]
        holo_indices[batch_index, :length] = example["holo_indices"]
        pair_mask[batch_index, :length] = True
    return {
        "pair_id": [example["pair_id"] for example in examples],
        "apo": apo,
        "holo": holo,
        "apo_indices": apo_indices,
        "holo_indices": holo_indices,
        "pair_mask": pair_mask,
        "sequence_identity": torch.tensor(
            [example["sequence_identity"] for example in examples], dtype=torch.float32
        ),
        "coverage": torch.tensor(
            [example["coverage"] for example in examples], dtype=torch.float32
        ),
        "metadata": [example["metadata"] for example in examples],
    }
