"""Fixed-frame invariant point attention for static structures.

AlphaFold updates frames to generate coordinates. FlexPPI keeps the observed
backbone frames fixed and only updates residue representations, because its
output is a dynamic state descriptor rather than a new structure.
"""

from __future__ import annotations

import math

import torch
from torch import nn

from flexppi.model.geometry import global_to_local, local_to_global


class InvariantPointAttention(nn.Module):
    def __init__(
        self,
        single_dim: int,
        pair_dim: int,
        heads: int,
        scalar_dim: int,
        qk_points: int,
        v_points: int,
        eps: float = 1e-8,
    ):
        super().__init__()
        self.heads = heads
        self.scalar_dim = scalar_dim
        self.qk_points = qk_points
        self.v_points = v_points
        self.eps = eps
        self.linear_q = nn.Linear(single_dim, heads * scalar_dim)
        self.linear_kv = nn.Linear(single_dim, heads * scalar_dim * 2)
        self.linear_q_points = nn.Linear(single_dim, heads * qk_points * 3)
        self.linear_kv_points = nn.Linear(
            single_dim, heads * (qk_points + v_points) * 3
        )
        self.linear_pair_bias = nn.Linear(pair_dim, heads)
        self.point_weights = nn.Parameter(torch.zeros(heads))
        output_dim = heads * (scalar_dim + v_points * 4 + pair_dim)
        self.linear_out = nn.Linear(output_dim, single_dim)
        nn.init.zeros_(self.linear_out.weight)
        nn.init.zeros_(self.linear_out.bias)

    def forward(
        self,
        single: torch.Tensor,
        pair: torch.Tensor,
        rotation: torch.Tensor,
        translation: torch.Tensor,
        mask: torch.Tensor,
    ) -> torch.Tensor:
        batch, residues, _ = single.shape
        q = self.linear_q(single).view(batch, residues, self.heads, self.scalar_dim)
        kv = self.linear_kv(single).view(batch, residues, self.heads, self.scalar_dim * 2)
        k, v = torch.split(kv, self.scalar_dim, dim=-1)
        q_points_local = self.linear_q_points(single).view(
            batch, residues, self.heads, self.qk_points, 3
        )
        kv_points_local = self.linear_kv_points(single).view(
            batch, residues, self.heads, self.qk_points + self.v_points, 3
        )
        k_points_local, v_points_local = torch.split(
            kv_points_local, [self.qk_points, self.v_points], dim=-2
        )
        q_points = local_to_global(q_points_local, rotation, translation)
        k_points = local_to_global(k_points_local, rotation, translation)
        v_points = local_to_global(v_points_local, rotation, translation)

        scalar_logits = torch.einsum("bihc,bjhc->bhij", q, k)
        scalar_logits = scalar_logits * math.sqrt(1.0 / (3.0 * self.scalar_dim))
        pair_logits = self.linear_pair_bias(pair).permute(0, 3, 1, 2)
        pair_logits = pair_logits * math.sqrt(1.0 / 3.0)
        point_difference = q_points[:, :, None] - k_points[:, None, :]
        point_distance2 = torch.sum(point_difference * point_difference, dim=-1)
        point_scale = torch.nn.functional.softplus(self.point_weights)
        point_scale = point_scale * math.sqrt(
            1.0 / (3.0 * (self.qk_points * 9.0 / 2.0))
        )
        point_logits = -0.5 * torch.sum(
            point_distance2 * point_scale[None, None, None, :, None], dim=-1
        )
        point_logits = point_logits.permute(0, 3, 1, 2)
        logits = scalar_logits + pair_logits + point_logits
        pair_mask = mask[:, None, :, None] & mask[:, None, None, :]
        logits = logits.masked_fill(~pair_mask, -1e5)
        attention = torch.softmax(logits, dim=-1)

        scalar_output = torch.einsum("bhij,bjhc->bihc", attention, v)
        point_output_global = torch.einsum("bhij,bjhpc->bihpc", attention, v_points)
        point_output_local = global_to_local(point_output_global, rotation, translation)
        point_norm = torch.sqrt(
            torch.sum(point_output_local * point_output_local, dim=-1) + self.eps
        )
        pair_output = torch.einsum("bhij,bijc->bihc", attention, pair)
        combined = torch.cat(
            (
                scalar_output.reshape(batch, residues, -1),
                point_output_local.reshape(batch, residues, -1),
                point_norm.reshape(batch, residues, -1),
                pair_output.reshape(batch, residues, -1),
            ),
            dim=-1,
        )
        return self.linear_out(combined) * mask[..., None].to(single.dtype)


class IPABlock(nn.Module):
    def __init__(
        self,
        single_dim: int,
        pair_dim: int,
        heads: int,
        scalar_dim: int,
        qk_points: int,
        v_points: int,
        dropout: float,
    ):
        super().__init__()
        self.input_norm = nn.LayerNorm(single_dim)
        self.ipa = InvariantPointAttention(
            single_dim, pair_dim, heads, scalar_dim, qk_points, v_points
        )
        self.ipa_dropout = nn.Dropout(dropout)
        self.transition_norm = nn.LayerNorm(single_dim)
        self.transition = nn.Sequential(
            nn.Linear(single_dim, single_dim * 4),
            nn.SiLU(),
            nn.Linear(single_dim * 4, single_dim),
        )
        self.transition_dropout = nn.Dropout(dropout)

    def forward(
        self,
        single: torch.Tensor,
        pair: torch.Tensor,
        rotation: torch.Tensor,
        translation: torch.Tensor,
        mask: torch.Tensor,
    ) -> torch.Tensor:
        single = single + self.ipa_dropout(
            self.ipa(self.input_norm(single), pair, rotation, translation, mask)
        )
        single = single + self.transition_dropout(self.transition(self.transition_norm(single)))
        return single * mask[..., None].to(single.dtype)

