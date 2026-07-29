import math

import torch

from flexppi.config import ModelConfig
from flexppi.losses import paired_loss
from flexppi.model.network import FlexPPI


def make_state(backbone):
    batch, residues = backbone.shape[:2]
    return {
        "aatype": torch.arange(residues).unsqueeze(0).repeat(batch, 1) % 20,
        "residue_index": torch.arange(residues).unsqueeze(0).repeat(batch, 1),
        "chain_index": torch.zeros(batch, residues, dtype=torch.long),
        "backbone": backbone,
        "backbone_mask": torch.ones(batch, residues, 4, dtype=torch.bool),
        "geometry_mask": torch.ones(batch, residues, dtype=torch.bool),
        "target_mask": torch.ones(batch, residues, dtype=torch.bool),
        "padding_mask": torch.ones(batch, residues, dtype=torch.bool),
        "condition": torch.zeros(batch, 3),
        "label_values": torch.rand(batch, residues, 3),
        "label_variances": torch.full((batch, residues, 3), 0.01),
        "label_mask": torch.ones(batch, residues, 3, dtype=torch.bool),
    }


def backbone_coordinates(residues=6):
    ca = torch.stack(
        (torch.arange(residues, dtype=torch.float32) * 3.8, torch.zeros(residues), torch.zeros(residues)),
        dim=-1,
    )
    n = ca + torch.tensor([-1.2, 0.7, 0.1])
    c = ca + torch.tensor([1.3, 0.5, -0.1])
    o = ca + torch.tensor([2.0, 0.8, -0.1])
    return torch.stack((n, ca, c, o), dim=1).unsqueeze(0)


def small_config():
    return ModelConfig(
        sequence_backend="learned",
        single_dim=32,
        pair_dim=16,
        sequence_dim=24,
        ipa_blocks=2,
        ipa_heads=4,
        ipa_scalar_dim=8,
        ipa_qk_points=2,
        ipa_v_points=3,
        rbf_bins=8,
        dropout=0.0,
    )


def test_state_prediction_is_rigid_transform_invariant():
    torch.manual_seed(3)
    model = FlexPPI(small_config()).eval()
    for block in model.blocks:
        torch.nn.init.normal_(block.ipa.linear_out.weight, std=0.03)
    backbone = backbone_coordinates()
    angle = 0.7
    rotation = torch.tensor(
        [[math.cos(angle), -math.sin(angle), 0.0],
         [math.sin(angle), math.cos(angle), 0.0],
         [0.0, 0.0, 1.0]],
        dtype=torch.float32,
    )
    transformed = torch.einsum("ij,bnkj->bnki", rotation, backbone)
    transformed = transformed + torch.tensor([10.0, -4.0, 7.0])
    first = model.forward_state(make_state(backbone))["mean"]
    second = model.forward_state(make_state(transformed))["mean"]
    assert torch.allclose(first, second, atol=2e-5, rtol=2e-5)


def test_pair_forward_loss_and_backward():
    torch.manual_seed(5)
    model = FlexPPI(small_config())
    apo = make_state(backbone_coordinates())
    holo_backbone = backbone_coordinates() + torch.tensor([0.0, 0.2, 0.1])
    holo = make_state(holo_backbone)
    residues = apo["aatype"].shape[1]
    batch = {
        "apo": apo,
        "holo": holo,
        "apo_indices": torch.arange(residues).unsqueeze(0),
        "holo_indices": torch.arange(residues).unsqueeze(0),
        "pair_mask": torch.ones(1, residues, dtype=torch.bool),
    }
    output = model.forward_pair(batch)
    assert output["delta_mean"].shape == (1, residues, 3)
    values = paired_loss(output, batch, [1.0, 0.5, 0.25])
    assert torch.isfinite(values["loss"])
    values["loss"].backward()
    assert any(parameter.grad is not None for parameter in model.parameters() if parameter.requires_grad)

