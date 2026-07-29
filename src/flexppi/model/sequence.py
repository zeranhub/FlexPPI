"""Learned and optional ESM sequence encoders."""

from __future__ import annotations

from contextlib import nullcontext
from typing import List

import torch
from torch import nn

from flexppi.constants import RESTYPES, UNKNOWN_INDEX


class LearnedSequenceEncoder(nn.Module):
    def __init__(self, output_dim: int):
        super().__init__()
        self.output_dim = output_dim
        self.embedding = nn.Embedding(UNKNOWN_INDEX + 1, output_dim)

    def forward(
        self,
        aatype: torch.Tensor,
        chain_index: torch.Tensor,
        padding_mask: torch.Tensor,
    ) -> torch.Tensor:
        del chain_index
        return self.embedding(aatype) * padding_mask[..., None].to(torch.float32)


class HuggingFaceESMEncoder(nn.Module):
    """Encode each chain separately, matching the chain treatment used by Ibex."""

    def __init__(self, model_name: str, freeze: bool = True):
        super().__init__()
        try:
            from transformers import AutoModel, AutoTokenizer
        except ImportError as exc:
            raise RuntimeError(
                "ESM backend requires: pip install -e .[esm]"
            ) from exc
        self.tokenizer = AutoTokenizer.from_pretrained(model_name)
        self.model = AutoModel.from_pretrained(model_name)
        self.output_dim = int(self.model.config.hidden_size)
        self.freeze = freeze
        if freeze:
            for parameter in self.model.parameters():
                parameter.requires_grad_(False)

    def forward(
        self,
        aatype: torch.Tensor,
        chain_index: torch.Tensor,
        padding_mask: torch.Tensor,
    ) -> torch.Tensor:
        device = aatype.device
        chains: List[str] = []
        locations = []
        alphabet = RESTYPES + "X"
        for batch_index in range(aatype.shape[0]):
            valid_chains = torch.unique(chain_index[batch_index][padding_mask[batch_index]])
            for chain in valid_chains.tolist():
                indices = torch.nonzero(
                    padding_mask[batch_index] & (chain_index[batch_index] == int(chain)),
                    as_tuple=False,
                ).squeeze(-1)
                sequence = "".join(
                    alphabet[min(int(value), UNKNOWN_INDEX)]
                    for value in aatype[batch_index, indices].tolist()
                )
                chains.append(sequence)
                locations.append((batch_index, indices))
        if not chains:
            return torch.zeros(
                aatype.shape[0], aatype.shape[1], self.output_dim, device=device
            )
        tokenized = self.tokenizer(chains, padding=True, return_tensors="pt")
        tokenized = {key: value.to(device) for key, value in tokenized.items()}
        context = torch.no_grad() if self.freeze else nullcontext()
        with context:
            hidden = self.model(**tokenized).last_hidden_state
        output = torch.zeros(
            aatype.shape[0], aatype.shape[1], self.output_dim,
            dtype=hidden.dtype, device=device,
        )
        for chain_number, (batch_index, indices) in enumerate(locations):
            length = int(indices.numel())
            output[batch_index, indices] = hidden[chain_number, 1 : length + 1]
        return output


def build_sequence_encoder(backend: str, output_dim: int, model_name: str, freeze: bool) -> nn.Module:
    backend = backend.lower()
    if backend == "learned":
        return LearnedSequenceEncoder(output_dim)
    if backend in {"esm", "huggingface_esm"}:
        return HuggingFaceESMEncoder(model_name, freeze=freeze)
    raise ValueError("Unknown sequence backend: {}".format(backend))

