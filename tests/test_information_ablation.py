"""Information-path tests use nonzero attention outputs to avoid vacuous passes."""

import copy
import math

import pytest
import torch

from flexppi.config import ModelConfig
from flexppi.losses import state_loss
from flexppi.model.information_ablation import build_information_ablation
from flexppi.model.network import FlexPPI


def config():
    return ModelConfig(
        single_dim=32, pair_dim=16, sequence_dim=24, ipa_blocks=2,
        ipa_heads=4, ipa_scalar_dim=8, ipa_qk_points=2, ipa_v_points=3,
        rbf_bins=8, dropout=0.0, use_condition_features=False,
        enable_response_head=False,
    )


def state():
    generator = torch.Generator().manual_seed(72)
    ca = torch.randn(2, 7, 3, generator=generator) * 3.0
    backbone = torch.stack(
        (ca + torch.tensor([-1.2, 0.7, 0.1]), ca,
         ca + torch.tensor([1.3, 0.5, -0.1]), ca + torch.tensor([2.0, 0.8, -0.1])),
        dim=2,
    )
    return {
        "aatype": torch.arange(14).view(2, 7) % 20,
        "chain_index": torch.tensor([[2, 2, 2, 2, 9, 9, 9], [2, 2, 2, 9, 9, 9, 9]]),
        "residue_index": torch.tensor([[1, 2, 3, 4, 101, 102, 103], [8, 9, 10, 21, 22, 23, 24]]),
        "backbone": backbone,
        "backbone_mask": torch.ones(2, 7, 4, dtype=torch.bool),
        "geometry_mask": torch.ones(2, 7, dtype=torch.bool),
        "padding_mask": torch.ones(2, 7, dtype=torch.bool),
        "target_mask": torch.ones(2, 7, dtype=torch.bool),
        "label_values": torch.rand(2, 7, 3, generator=generator),
        "label_variances": torch.full((2, 7, 3), 0.01),
        "label_mask": torch.ones(2, 7, 3, dtype=torch.bool),
    }


def model(variant):
    torch.manual_seed(83)
    result = build_information_ablation(config(), variant).eval()
    for block in result.blocks:
        attention = getattr(block, "ipa", getattr(block, "attention", None))
        torch.nn.init.normal_(attention.linear_out.weight, std=0.07)
    return result


def assert_predictions_equal(first, second, indices=None):
    for key in ("representation", "mean", "log_variance", "mask"):
        left, right = first[key], second[key]
        if indices is not None:
            left, right = left[indices], right[indices]
        torch.testing.assert_close(left, right, atol=3e-5, rtol=3e-5)


def test_sequence_ignores_all_geometry_inputs_and_pdb_numbering():
    network = model("sequence")
    original = state()
    changed = copy.deepcopy(original)
    changed["backbone"] = torch.full_like(changed["backbone"], float("nan"))
    changed["backbone_mask"][:] = False
    changed["geometry_mask"][:] = False
    changed["residue_index"] = torch.flip(changed["residue_index"] * 37, dims=[1])
    assert_predictions_equal(network(original), network(changed))
    for key in ("backbone", "backbone_mask", "geometry_mask", "residue_index"):
        del changed[key]
    assert_predictions_equal(network(original), network(changed))


def test_geometry_ignores_amino_acid_identity_but_responds_to_geometry():
    network = model("geometry")
    original = state()
    changed = copy.deepcopy(original)
    changed["aatype"] = (changed["aatype"] + 8) % 20
    first = network(original)
    assert_predictions_equal(first, network(changed))
    changed["backbone"][:, 0] += torch.tensor([13.0, -8.0, 7.0])
    assert not torch.allclose(first["mean"], network(changed)["mean"], atol=1e-5)


@pytest.mark.parametrize("variant", ["sequence", "intrachain"])
def test_chain_restricted_models_ignore_partner_changes_and_removal(variant):
    network = model(variant)
    original = state()
    target = original["chain_index"] == 2
    changed = copy.deepcopy(original)
    changed["aatype"][~target] = (changed["aatype"][~target] + 7) % 20
    changed["backbone"][~target] += torch.tensor([23.0, -14.0, 17.0])
    first = network(original)
    assert_predictions_equal(first, network(changed), target)
    changed["padding_mask"][~target] = False
    assert_predictions_equal(first, network(changed), target)
    # A positive control establishes that the original full encoder can use partners.
    full = model("full")
    assert not torch.allclose(
        full(original)["mean"][target], full(changed)["mean"][target], atol=1e-5
    )


def test_sequence_uses_context_and_local_order_without_absolute_chain_ids():
    network = model("sequence")
    original = state()
    changed = copy.deepcopy(original)
    changed["aatype"][:, 1] = (changed["aatype"][:, 1] + 10) % 20
    first = network(original)
    assert not torch.allclose(first["mean"][:, 0], network(changed)["mean"][:, 0], atol=1e-5)
    changed = copy.deepcopy(original)
    changed["chain_index"] = changed["chain_index"] * 100 + 47
    assert_predictions_equal(first, network(changed))
    # The samples have different chain boundaries (4+3 and 3+4), so each
    # needs its own interleaving to preserve the internal order of both chains.
    permutations = (
        torch.tensor([0, 4, 1, 5, 2, 6, 3]),
        torch.tensor([0, 3, 1, 4, 2, 5, 6]),
    )
    for sample_index, permutation in enumerate(permutations):
        changed = {
            key: value[sample_index:sample_index + 1, permutation]
            if value.ndim >= 2 and value.shape[1] == 7
            else value[sample_index:sample_index + 1]
            for key, value in original.items()
        }
        reordered = network(changed)
        inverse = permutation.argsort()
        expected = {key: value[sample_index:sample_index + 1] for key, value in first.items()}
        assert_predictions_equal(
            expected, {key: value[:, inverse] for key, value in reordered.items()}
        )


@pytest.mark.parametrize("variant", ["full", "geometry", "intrachain"])
def test_geometry_models_are_rigid_transform_invariant(variant):
    network = model(variant)
    original = state()
    changed = copy.deepcopy(original)
    angle = 0.7
    rotation = torch.tensor([
        [math.cos(angle), -math.sin(angle), 0.0],
        [math.sin(angle), math.cos(angle), 0.0], [0.0, 0.0, 1.0],
    ])
    changed["backbone"] = torch.einsum("ij,bnkj->bnki", rotation, original["backbone"])
    changed["backbone"] += torch.tensor([10.0, -4.0, 7.0])
    assert_predictions_equal(network(original), network(changed))


@pytest.mark.parametrize("variant", ["sequence", "geometry", "intrachain", "full"])
def test_variants_support_masked_state_loss_and_gradients(variant):
    network = model(variant)
    batch = state()
    batch["padding_mask"][1, 5:] = False
    batch["label_mask"][0, 0, 1:] = False
    output = network(batch)
    loss = state_loss(output, batch, [1.0, 1.0], config().metrics)
    assert torch.isfinite(loss)
    loss.backward()
    for name, parameter in network.named_parameters():
        assert parameter.grad is not None, "No gradient for {}".format(name)
        assert torch.isfinite(parameter.grad).all(), "Invalid gradient for {}".format(name)
    assert not output["mask"][1, 5:].any()


def test_full_factory_preserves_checkpoint_contract():
    torch.manual_seed(17)
    original = FlexPPI(config()).eval()
    for variant in ("full", "geometry", "intrachain"):
        network = build_information_ablation(config(), variant).eval()
        network.load_state_dict(original.state_dict(), strict=True)
    full = build_information_ablation(config(), "full").eval()
    full.load_state_dict(original.state_dict(), strict=True)
    assert_predictions_equal(original(state()), full(state()))


def test_invalid_variant_is_rejected():
    with pytest.raises(ValueError, match="Unknown information ablation"):
        build_information_ablation(config(), "zero_distances")
