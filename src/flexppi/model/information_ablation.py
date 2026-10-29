"""Retrainable models with explicit information restrictions.

These models are for *retrained* ablations, not for replacing an input in an
already trained checkpoint.  ``sequence`` never reads coordinates, atom masks,
geometry masks, or PDB residue numbering.  Its attention is restricted to each
chain and uses ordinal positions in that chain.  ``geometry`` removes explicit
amino-acid identity from both embedding paths, but retains backbone geometry,
chain boundaries and sequence separation; these can still correlate with amino
acid identity.  ``intrachain`` evaluates the original encoder one chain at a
time, sharing all weights, so no IPA pathway can reach another chain.

For information-source experiments, use ``use_condition_features=False`` in
every arm, match the label masks externally, and report parameter counts.  The
sequence architecture has fewer parameters than IPA and is not a capacity-
matched claim by itself.  The existing StateDataset already applies the same
frame-validity restriction to label_mask for every model.
"""

from __future__ import annotations

import math
from typing import Any, Dict

import torch
from torch import nn

from flexppi.config import ModelConfig
from flexppi.constants import UNKNOWN_INDEX
from flexppi.model.heads import ResponseHead, StateHead
from flexppi.model.network import FlexPPI
from flexppi.model.sequence import build_sequence_encoder


class ChainScalarAttention(nn.Module):
    """Scalar attention with learned, signed chain-local position biases."""

    def __init__(self, single_dim: int, heads: int, scalar_dim: int):
        super().__init__()
        self.heads = heads
        self.scalar_dim = scalar_dim
        self.linear_q = nn.Linear(single_dim, heads * scalar_dim)
        self.linear_kv = nn.Linear(single_dim, heads * scalar_dim * 2)
        self.relative_bias = nn.Embedding(65, heads)
        self.linear_out = nn.Linear(heads * scalar_dim, single_dim)
        # Match the residual attention initialization used by IPABlock.
        nn.init.zeros_(self.linear_out.weight)
        nn.init.zeros_(self.linear_out.bias)

    def forward(
        self, single: torch.Tensor, chain_index: torch.Tensor, mask: torch.Tensor
    ) -> torch.Tensor:
        batch_size, residues, _ = single.shape
        query = self.linear_q(single).view(
            batch_size, residues, self.heads, self.scalar_dim
        )
        key_value = self.linear_kv(single).view(
            batch_size, residues, self.heads, self.scalar_dim * 2
        )
        key, value = torch.split(key_value, self.scalar_dim, dim=-1)
        same_chain = chain_index[:, :, None] == chain_index[:, None, :]
        # Count preceding valid tokens belonging to the same chain.  Absolute
        # PDB numbering and the number/order of other chains cannot affect it.
        preceding = torch.ones(
            residues, residues, dtype=torch.bool, device=single.device
        ).tril()
        position = (
            same_chain & preceding[None] & mask[:, None, :]
        ).sum(dim=-1) - 1
        separation = position[:, None, :] - position[:, :, None]
        position_bias = self.relative_bias(
            separation.clamp(min=-32, max=32).long() + 32
        ).permute(0, 3, 1, 2)
        logits = torch.einsum("bihc,bjhc->bhij", query, key)
        logits = logits / math.sqrt(self.scalar_dim) + position_bias
        allowed = same_chain & mask[:, :, None] & mask[:, None, :]
        # Finite masking plus explicit zeroing also handles an all-padding row.
        attention = torch.softmax(logits.masked_fill(~allowed[:, None], -1e4), dim=-1)
        attention = attention * allowed[:, None].to(attention.dtype)
        result = torch.einsum("bhij,bjhc->bihc", attention, value)
        result = self.linear_out(result.reshape(batch_size, residues, -1))
        return result * mask[..., None].to(result.dtype)


class ChainScalarBlock(nn.Module):
    def __init__(self, config: ModelConfig):
        super().__init__()
        self.input_norm = nn.LayerNorm(config.single_dim)
        self.attention = ChainScalarAttention(
            config.single_dim, config.ipa_heads, config.ipa_scalar_dim
        )
        self.attention_dropout = nn.Dropout(config.dropout)
        self.transition_norm = nn.LayerNorm(config.single_dim)
        self.transition = nn.Sequential(
            nn.Linear(config.single_dim, config.single_dim * 4),
            nn.SiLU(),
            nn.Linear(config.single_dim * 4, config.single_dim),
        )
        self.transition_dropout = nn.Dropout(config.dropout)

    def forward(
        self, single: torch.Tensor, chain_index: torch.Tensor, mask: torch.Tensor
    ) -> torch.Tensor:
        single = single + self.attention_dropout(
            self.attention(self.input_norm(single), chain_index, mask)
        )
        single = single + self.transition_dropout(
            self.transition(self.transition_norm(single))
        )
        return single * mask[..., None].to(single.dtype)


class SequenceOnlyFlexPPI(FlexPPI):
    """Sequence-context control sharing the state/loss/paired-forward interface."""

    def __init__(self, config: ModelConfig):
        # Do not instantiate IPA or a geometry encoder with unused parameters.
        nn.Module.__init__(self)
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
        self.condition_mlp = (
            nn.Sequential(
                nn.Linear(3, config.single_dim),
                nn.SiLU(),
                nn.Linear(config.single_dim, config.single_dim),
            )
            if config.use_condition_features
            else None
        )
        self.blocks = nn.ModuleList(
            [ChainScalarBlock(config) for _ in range(config.ipa_blocks)]
        )
        self.final_norm = nn.LayerNorm(config.single_dim)
        self.state_head = StateHead(config.single_dim, self.metric_count)
        self.response_head = (
            ResponseHead(config.single_dim, self.metric_count)
            if config.enable_response_head
            else None
        )

    def forward_state(self, batch: Dict[str, Any]) -> Dict[str, torch.Tensor]:
        mask = batch["padding_mask"].bool()
        sequence = self.sequence_encoder(batch["aatype"], batch["chain_index"], mask)
        single = self.single_input(
            torch.cat((sequence, self.aa_embedding(batch["aatype"])), dim=-1)
        )
        if self.condition_mlp is not None:
            single = single + self.condition_mlp(batch["condition"])[:, None, :]
        for block in self.blocks:
            single = block(single, batch["chain_index"], mask)
        representation = self.final_norm(single)
        mean, log_variance = self.state_head(representation)
        return {
            "representation": representation,
            "mean": mean,
            "log_variance": log_variance,
            "mask": mask,
        }


class GeometryOnlyFlexPPI(FlexPPI):
    """Keep structural/topological inputs but hide explicit residue identity."""

    def forward_state(self, batch: Dict[str, Any]) -> Dict[str, torch.Tensor]:
        restricted_batch = dict(batch)
        # Both the sequence encoder and aa_embedding consume this tensor.
        restricted_batch["aatype"] = torch.full_like(batch["aatype"], UNKNOWN_INDEX)
        return super().forward_state(restricted_batch)


class IntrachainFlexPPI(FlexPPI):
    """Run every chain independently through one shared, unchanged IPA encoder."""

    def forward_state(self, batch: Dict[str, Any]) -> Dict[str, torch.Tensor]:
        padding = batch["padding_mask"].bool()
        batch_size, residues = padding.shape
        residue_keys = (
            "aatype", "chain_index", "residue_index", "backbone", "backbone_mask",
            "geometry_mask", "padding_mask",
        )
        output_rows: Dict[str, list] = {
            "representation": [], "mean": [], "log_variance": [], "mask": []
        }
        reference = self.single_input.weight
        for batch_index in range(batch_size):
            row = {
                "representation": reference.new_zeros(residues, self.config.single_dim),
                "mean": reference.new_zeros(residues, self.metric_count),
                "log_variance": reference.new_zeros(residues, self.metric_count),
                "mask": torch.zeros(residues, dtype=torch.bool, device=padding.device),
            }
            chains = torch.unique(batch["chain_index"][batch_index][padding[batch_index]])
            for chain in chains:
                indices = torch.nonzero(
                    padding[batch_index] & (batch["chain_index"][batch_index] == chain),
                    as_tuple=False,
                ).squeeze(-1)
                chain_batch = {
                    key: batch[key][batch_index:batch_index + 1].index_select(1, indices)
                    for key in residue_keys
                }
                if self.condition_mlp is not None:
                    chain_batch["condition"] = batch["condition"][batch_index:batch_index + 1]
                prediction = super().forward_state(chain_batch)
                for key in row:
                    # Out-of-place index_copy preserves gradients through all chains.
                    row[key] = row[key].index_copy(0, indices, prediction[key][0])
            for key in row:
                output_rows[key].append(row[key])
        return {key: torch.stack(rows, dim=0) for key, rows in output_rows.items()}


def build_information_ablation(config: ModelConfig, variant: str) -> FlexPPI:
    """Build ``sequence``, ``geometry``, ``intrachain`` or unchanged ``full``.

    Existing full checkpoints load directly into ``full``, ``geometry`` and
    ``intrachain`` because their parameter names are unchanged.  Loading such a
    checkpoint into a restricted arm is only a sensitivity experiment; claims
    about predictive information must use independently retrained models.
    """
    models = {
        "sequence": SequenceOnlyFlexPPI,
        "geometry": GeometryOnlyFlexPPI,
        "intrachain": IntrachainFlexPPI,
        "full": FlexPPI,
    }
    if variant not in models:
        raise ValueError("Unknown information ablation: {}. Choose {}".format(
            variant, ", ".join(models)
        ))
    return models[variant](config)
