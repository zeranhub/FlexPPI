"""Noise-aware absolute and paired objectives."""

from __future__ import annotations

from typing import Dict, Optional, Sequence

import torch
import torch.nn.functional as functional

from flexppi.constants import METRIC_NAMES
from flexppi.model.network import gather_residues


def _metric_weight_tensor(weights: Sequence[float], reference: torch.Tensor) -> torch.Tensor:
    return torch.as_tensor(weights, dtype=reference.dtype, device=reference.device)


def heteroscedastic_nll(
    predicted_mean: torch.Tensor,
    predicted_log_variance: torch.Tensor,
    target: torch.Tensor,
    label_variance: torch.Tensor,
    mask: torch.Tensor,
    metric_weights: Sequence[float],
) -> torch.Tensor:
    valid = mask.bool() & torch.isfinite(target)
    if not torch.any(valid):
        return predicted_mean.sum() * 0.0
    safe_target = torch.where(valid, target, predicted_mean.detach())
    total_variance = torch.exp(predicted_log_variance) + torch.clamp(label_variance, min=0.0)
    total_variance = torch.clamp(total_variance, min=1e-6, max=1e6)
    nll = 0.5 * (
        (predicted_mean - safe_target) ** 2 / total_variance + torch.log(total_variance)
    )
    weights = _metric_weight_tensor(metric_weights, nll)
    weighted_valid = valid.to(nll.dtype) * weights
    return torch.sum(nll * weighted_valid) / torch.clamp(weighted_valid.sum(), min=1.0)


def state_loss(
    output: Dict[str, torch.Tensor],
    batch: Dict[str, torch.Tensor],
    metric_weights: Sequence[float],
    metric_names: Optional[Sequence[str]] = None,
) -> torch.Tensor:
    label_values = batch["label_values"]
    label_variances = batch["label_variances"]
    label_mask = batch["label_mask"]
    if metric_names is not None:
        unknown = sorted(set(metric_names) - set(METRIC_NAMES))
        if unknown:
            raise ValueError("Unknown state metrics: {}".format(", ".join(unknown)))
        indices = [METRIC_NAMES.index(name) for name in metric_names]
        label_values = label_values[..., indices]
        label_variances = label_variances[..., indices]
        label_mask = label_mask[..., indices]
    if output["mean"].shape[-1] != label_values.shape[-1]:
        raise ValueError(
            "Model predicts {} metrics but the selected labels contain {}".format(
                output["mean"].shape[-1], label_values.shape[-1]
            )
        )
    mask = label_mask & batch["target_mask"][..., None] & output["mask"][..., None]
    return heteroscedastic_nll(
        output["mean"], output["log_variance"], label_values,
        label_variances, mask, metric_weights,
    )


def paired_loss(
    output: Dict[str, torch.Tensor],
    batch: Dict[str, torch.Tensor],
    metric_weights: Sequence[float],
    state_weight: float = 0.5,
    response_weight: float = 1.0,
    consistency_weight: float = 0.25,
    metric_names: Optional[Sequence[str]] = None,
) -> Dict[str, torch.Tensor]:
    apo_state = state_loss(output["apo"], batch["apo"], metric_weights, metric_names)
    holo_state = state_loss(output["holo"], batch["holo"], metric_weights, metric_names)
    apo_indices = batch["apo_indices"]
    holo_indices = batch["holo_indices"]
    apo_values = gather_residues(batch["apo"]["label_values"], apo_indices)
    holo_values = gather_residues(batch["holo"]["label_values"], holo_indices)
    apo_variances = gather_residues(batch["apo"]["label_variances"], apo_indices)
    holo_variances = gather_residues(batch["holo"]["label_variances"], holo_indices)
    apo_mask = gather_residues(batch["apo"]["label_mask"], apo_indices)
    holo_mask = gather_residues(batch["holo"]["label_mask"], holo_indices)
    if metric_names is not None:
        indices = [METRIC_NAMES.index(name) for name in metric_names]
        apo_values = apo_values[..., indices]
        holo_values = holo_values[..., indices]
        apo_variances = apo_variances[..., indices]
        holo_variances = holo_variances[..., indices]
        apo_mask = apo_mask[..., indices]
        holo_mask = holo_mask[..., indices]
    delta_target = holo_values - apo_values
    delta_variance = holo_variances + apo_variances
    delta_mask = apo_mask & holo_mask & output["pair_mask"][..., None]
    response = heteroscedastic_nll(
        output["delta_mean"], output["delta_log_variance"], delta_target,
        delta_variance, delta_mask, metric_weights,
    )
    consistency_elements = functional.smooth_l1_loss(
        output["delta_mean"], output["state_difference"].detach(), reduction="none"
    )
    consistency_mask = output["pair_mask"][..., None].to(consistency_elements.dtype)
    consistency = torch.sum(consistency_elements * consistency_mask) / torch.clamp(
        consistency_mask.sum() * consistency_elements.shape[-1], min=1.0
    )
    total = (
        state_weight * 0.5 * (apo_state + holo_state)
        + response_weight * response
        + consistency_weight * consistency
    )
    return {
        "loss": total,
        "apo_state_loss": apo_state,
        "holo_state_loss": holo_state,
        "response_loss": response,
        "consistency_loss": consistency,
    }

