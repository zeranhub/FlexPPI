"""Absolute-state and signed-response prediction heads."""

from __future__ import annotations

import torch
from torch import nn


class StateHead(nn.Module):
    def __init__(self, single_dim: int, metric_count: int):
        super().__init__()
        self.metric_count = metric_count
        self.network = nn.Sequential(
            nn.LayerNorm(single_dim),
            nn.Linear(single_dim, single_dim),
            nn.SiLU(),
            nn.Linear(single_dim, metric_count * 2),
        )

    def forward(self, representation: torch.Tensor) -> tuple:
        raw = self.network(representation).view(
            representation.shape[0], representation.shape[1], self.metric_count, 2
        )
        mean = torch.nn.functional.softplus(raw[..., 0])
        log_variance = torch.clamp(raw[..., 1], min=-8.0, max=6.0)
        return mean, log_variance


class ResponseHead(nn.Module):
    def __init__(self, single_dim: int, metric_count: int):
        super().__init__()
        self.metric_count = metric_count
        input_dim = single_dim * 4 + metric_count
        self.network = nn.Sequential(
            nn.LayerNorm(input_dim),
            nn.Linear(input_dim, single_dim * 2),
            nn.SiLU(),
            nn.Linear(single_dim * 2, metric_count * 2),
        )

    def forward(
        self,
        apo_representation: torch.Tensor,
        holo_representation: torch.Tensor,
        state_difference: torch.Tensor,
    ) -> tuple:
        features = torch.cat(
            (
                apo_representation,
                holo_representation,
                holo_representation - apo_representation,
                holo_representation * apo_representation,
                state_difference,
            ),
            dim=-1,
        )
        raw = self.network(features).view(
            features.shape[0], features.shape[1], self.metric_count, 2
        )
        mean = raw[..., 0]
        log_variance = torch.clamp(raw[..., 1], min=-8.0, max=6.0)
        return mean, log_variance

