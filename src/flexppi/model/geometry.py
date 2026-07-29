"""Rigid residue frames and invariant pair features."""

from __future__ import annotations

import torch
from torch import nn


def safe_normalize(vector: torch.Tensor, eps: float = 1e-8) -> torch.Tensor:
    return vector / torch.sqrt(torch.sum(vector * vector, dim=-1, keepdim=True) + eps)


def frames_from_backbone(
    backbone: torch.Tensor,
    backbone_mask: torch.Tensor,
) -> tuple:
    """Build local-to-global frames from N, CA, C atoms.

    The returned rotation stores local basis vectors as columns. Therefore a
    local point p is transformed as R @ p + t, where t is the CA coordinate.
    """

    n_atom = backbone[..., 0, :]
    ca_atom = backbone[..., 1, :]
    c_atom = backbone[..., 2, :]
    x_axis = safe_normalize(c_atom - ca_atom)
    n_direction = n_atom - ca_atom
    y_axis = n_direction - torch.sum(n_direction * x_axis, dim=-1, keepdim=True) * x_axis
    y_axis = safe_normalize(y_axis)
    z_axis = safe_normalize(torch.cross(x_axis, y_axis, dim=-1))
    rotation = torch.stack((x_axis, y_axis, z_axis), dim=-1)
    frame_mask = torch.all(backbone_mask[..., :3], dim=-1)
    identity = torch.eye(3, dtype=rotation.dtype, device=rotation.device)
    rotation = torch.where(frame_mask[..., None, None], rotation, identity)
    translation = torch.where(frame_mask[..., None], ca_atom, torch.zeros_like(ca_atom))
    return rotation, translation, frame_mask


def local_to_global(
    points: torch.Tensor,
    rotation: torch.Tensor,
    translation: torch.Tensor,
) -> torch.Tensor:
    transformed = torch.einsum("bnij,bnhpj->bnhpi", rotation, points)
    return transformed + translation[:, :, None, None, :]


def global_to_local(
    points: torch.Tensor,
    rotation: torch.Tensor,
    translation: torch.Tensor,
) -> torch.Tensor:
    shifted = points - translation[:, :, None, None, :]
    inverse = rotation.transpose(-1, -2)[:, :, None, None, :, :]
    return torch.matmul(inverse, shifted.unsqueeze(-1)).squeeze(-1)


class GeometryPairEncoder(nn.Module):
    """Encode distances, chain relation, sequence separation, and direction."""

    def __init__(self, pair_dim: int, rbf_bins: int = 24, max_distance: float = 24.0):
        super().__init__()
        self.rbf_bins = rbf_bins
        self.max_distance = max_distance
        centers = torch.linspace(0.0, max_distance, rbf_bins)
        self.register_buffer("centers", centers)
        width = float(centers[1] - centers[0]) if rbf_bins > 1 else max_distance
        self.register_buffer("width", torch.tensor(width))
        raw_dim = rbf_bins + 2 + 3
        self.projection = nn.Sequential(
            nn.Linear(raw_dim, pair_dim),
            nn.SiLU(),
            nn.Linear(pair_dim, pair_dim),
        )

    def forward(
        self,
        backbone: torch.Tensor,
        rotation: torch.Tensor,
        chain_index: torch.Tensor,
        residue_index: torch.Tensor,
        mask: torch.Tensor,
    ) -> torch.Tensor:
        ca = backbone[..., 1, :]
        difference = ca[:, None, :, :] - ca[:, :, None, :]
        distance = torch.sqrt(torch.sum(difference * difference, dim=-1) + 1e-8)
        rbf = torch.exp(-((distance[..., None] - self.centers) / self.width) ** 2)
        same_chain = (chain_index[:, :, None] == chain_index[:, None, :]).to(ca.dtype)
        separation = torch.abs(residue_index[:, :, None] - residue_index[:, None, :]).to(ca.dtype)
        separation = torch.clamp(separation, max=32.0) / 32.0
        inverse_i = rotation.transpose(-1, -2)[:, :, None, :, :]
        local_difference = torch.matmul(inverse_i, difference.unsqueeze(-1)).squeeze(-1)
        direction = local_difference / torch.sqrt(
            torch.sum(local_difference * local_difference, dim=-1, keepdim=True) + 1e-8
        )
        raw = torch.cat((rbf, same_chain[..., None], separation[..., None], direction), dim=-1)
        pair = self.projection(raw)
        pair_mask = mask[:, :, None] & mask[:, None, :]
        return pair * pair_mask[..., None].to(pair.dtype)

