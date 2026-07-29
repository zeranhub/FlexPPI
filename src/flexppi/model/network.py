"""Shared state encoder and paired response model."""

from __future__ import annotations

from typing import Any, Dict

import torch
from torch import nn

from flexppi.config import ModelConfig
from flexppi.constants import UNKNOWN_INDEX
from flexppi.model.geometry import GeometryPairEncoder, frames_from_backbone
from flexppi.model.heads import ResponseHead, StateHead
from flexppi.model.ipa import IPABlock
from flexppi.model.sequence import build_sequence_encoder


def gather_residues(tensor: torch.Tensor, indices: torch.Tensor) -> torch.Tensor:
    batch_indices = torch.arange(tensor.shape[0], device=tensor.device)[:, None]
    return tensor[batch_indices, indices]


class FlexPPI(nn.Module):
    def __init__(self, config: ModelConfig):
        super().__init__()
        self.config = config
        self.metric_count = len(config.metrics)
        self.sequence_encoder = build_sequence_encoder(
            config.sequence_backend,
            config.sequence_dim,
            config.esm_model_name,
            config.freeze_esm,
        )
        self.aa_embedding = nn.Embedding(UNKNOWN_INDEX + 1, 32)
        self.single_input = nn.Linear(self.sequence_encoder.output_dim + 32, config.single_dim)
        self.condition_mlp = nn.Sequential(
            nn.Linear(3, config.single_dim),
            nn.SiLU(),
            nn.Linear(config.single_dim, config.single_dim),
        )
        self.pair_encoder = GeometryPairEncoder(config.pair_dim, config.rbf_bins)
        self.blocks = nn.ModuleList(
            [
                IPABlock(
                    config.single_dim,
                    config.pair_dim,
                    config.ipa_heads,
                    config.ipa_scalar_dim,
                    config.ipa_qk_points,
                    config.ipa_v_points,
                    config.dropout,
                )
                for _ in range(config.ipa_blocks)
            ]
        )
        self.final_norm = nn.LayerNorm(config.single_dim)
        self.state_head = StateHead(config.single_dim, self.metric_count)
        self.response_head = ResponseHead(config.single_dim, self.metric_count)

    def forward_state(self, batch: Dict[str, Any]) -> Dict[str, torch.Tensor]:
        padding_mask = batch["padding_mask"].bool()
        geometry_mask = batch["geometry_mask"].bool() & padding_mask
        sequence = self.sequence_encoder(
            batch["aatype"], batch["chain_index"], padding_mask
        )
        single = self.single_input(
            torch.cat((sequence, self.aa_embedding(batch["aatype"])), dim=-1)
        )
        single = single + self.condition_mlp(batch["condition"])[:, None, :]
        rotation, translation, frame_mask = frames_from_backbone(
            batch["backbone"], batch["backbone_mask"]
        )
        geometry_mask = geometry_mask & frame_mask
        pair = self.pair_encoder(
            batch["backbone"], rotation, batch["chain_index"],
            batch["residue_index"], geometry_mask,
        )
        for block in self.blocks:
            single = block(single, pair, rotation, translation, geometry_mask)
        representation = self.final_norm(single)
        mean, log_variance = self.state_head(representation)
        return {
            "representation": representation,
            "mean": mean,
            "log_variance": log_variance,
            "mask": geometry_mask,
        }

    def forward_pair(self, batch: Dict[str, Any]) -> Dict[str, Any]:
        apo = self.forward_state(batch["apo"])
        holo = self.forward_state(batch["holo"])
        apo_indices = batch["apo_indices"]
        holo_indices = batch["holo_indices"]
        apo_representation = gather_residues(apo["representation"], apo_indices)
        holo_representation = gather_residues(holo["representation"], holo_indices)
        apo_mean = gather_residues(apo["mean"], apo_indices)
        holo_mean = gather_residues(holo["mean"], holo_indices)
        state_difference = holo_mean - apo_mean
        delta_mean, delta_log_variance = self.response_head(
            apo_representation, holo_representation, state_difference
        )
        pair_mask = batch["pair_mask"].bool()
        pair_mask = pair_mask & gather_residues(apo["mask"], apo_indices)
        pair_mask = pair_mask & gather_residues(holo["mask"], holo_indices)
        return {
            "apo": apo,
            "holo": holo,
            "delta_mean": delta_mean,
            "delta_log_variance": delta_log_variance,
            "state_difference": state_difference,
            "pair_mask": pair_mask,
        }

    def forward(self, batch: Dict[str, Any]) -> Dict[str, Any]:
        if "apo" in batch and "holo" in batch:
            return self.forward_pair(batch)
        return self.forward_state(batch)

