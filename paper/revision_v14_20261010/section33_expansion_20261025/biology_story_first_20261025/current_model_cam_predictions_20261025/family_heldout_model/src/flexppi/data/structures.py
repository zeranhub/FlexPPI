"""Minimal, dependency-free PDB parsing for residue-level learning.

Only the first model is read. Hydrogen atoms are ignored naturally because the
model only requires N, CA, C, and O. Residues with an incomplete N-CA-C frame
remain in the sequence but receive a false geometry mask.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import numpy as np

from flexppi.constants import (
    BACKBONE_ATOM_INDEX,
    BACKBONE_ATOMS,
    RESTYPE_TO_INDEX,
    THREE_TO_ONE,
    UNKNOWN_INDEX,
)


@dataclass
class ProteinStructure:
    sequence: str
    residue_names: List[str]
    chain_ids: List[str]
    author_residue_ids: List[str]
    aatype: np.ndarray
    residue_index: np.ndarray
    chain_index: np.ndarray
    backbone: np.ndarray
    backbone_mask: np.ndarray

    @property
    def num_residues(self) -> int:
        return len(self.sequence)

    @property
    def frame_mask(self) -> np.ndarray:
        return self.backbone_mask[:, :3].all(axis=-1)

    def unique_chains(self) -> List[str]:
        return list(dict.fromkeys(self.chain_ids))

    def indices_for_chain(self, chain_id: str) -> np.ndarray:
        return np.asarray([i for i, value in enumerate(self.chain_ids) if value == chain_id])

    def chain_sequence(self, chain_id: str) -> str:
        indices = self.indices_for_chain(chain_id)
        return "".join(self.sequence[int(i)] for i in indices)


def _normalise_chain_id(value: str) -> str:
    value = value.strip()
    return value if value else "_"


def parse_pdb(path: str) -> ProteinStructure:
    """Parse standard amino-acid residues and backbone coordinates from a PDB."""

    residue_order: List[Tuple[str, str, str]] = []
    residues: Dict[Tuple[str, str, str], Dict[str, object]] = {}
    saw_model = False

    with Path(path).open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            record = line[:6].strip()
            if record == "MODEL":
                if saw_model:
                    break
                saw_model = True
                continue
            if record == "ENDMDL":
                break
            if record not in {"ATOM", "HETATM"} or len(line) < 54:
                continue
            altloc = line[16:17]
            if altloc not in {" ", "A", "1"}:
                continue
            residue_name = line[17:20].strip().upper()
            if residue_name not in THREE_TO_ONE:
                continue
            atom_name = line[12:16].strip().upper()
            chain_id = _normalise_chain_id(line[21:22])
            residue_number = line[22:26].strip()
            insertion_code = line[26:27].strip()
            key = (chain_id, residue_number, insertion_code)
            try:
                xyz = np.asarray(
                    [float(line[30:38]), float(line[38:46]), float(line[46:54])],
                    dtype=np.float32,
                )
            except ValueError:
                continue
            if key not in residues:
                residues[key] = {"name": residue_name, "atoms": {}}
                residue_order.append(key)
            atom_dict = residues[key]["atoms"]
            if atom_name not in atom_dict:
                atom_dict[atom_name] = xyz

    if not residue_order:
        raise ValueError("No protein residues found in {}".format(path))

    n_res = len(residue_order)
    backbone = np.zeros((n_res, len(BACKBONE_ATOMS), 3), dtype=np.float32)
    backbone_mask = np.zeros((n_res, len(BACKBONE_ATOMS)), dtype=bool)
    residue_names: List[str] = []
    chain_ids: List[str] = []
    author_ids: List[str] = []
    sequence_chars: List[str] = []
    residue_index = np.zeros(n_res, dtype=np.int64)
    chain_index = np.zeros(n_res, dtype=np.int64)
    chain_to_number: Dict[str, int] = {}
    chain_positions: Dict[str, int] = {}

    for i, key in enumerate(residue_order):
        chain_id, residue_number, insertion_code = key
        residue = residues[key]
        residue_name = str(residue["name"])
        atoms = residue["atoms"]
        one_letter = THREE_TO_ONE.get(residue_name, "X")
        residue_names.append(residue_name)
        chain_ids.append(chain_id)
        author_ids.append(residue_number + insertion_code)
        sequence_chars.append(one_letter)
        if chain_id not in chain_to_number:
            chain_to_number[chain_id] = len(chain_to_number)
            chain_positions[chain_id] = 0
        chain_index[i] = chain_to_number[chain_id]
        residue_index[i] = chain_positions[chain_id]
        chain_positions[chain_id] += 1
        for atom_name, atom_idx in BACKBONE_ATOM_INDEX.items():
            if atom_name in atoms:
                backbone[i, atom_idx] = atoms[atom_name]
                backbone_mask[i, atom_idx] = True

    sequence = "".join(sequence_chars)
    aatype = np.asarray(
        [RESTYPE_TO_INDEX.get(letter, UNKNOWN_INDEX) for letter in sequence], dtype=np.int64
    )
    return ProteinStructure(
        sequence=sequence,
        residue_names=residue_names,
        chain_ids=chain_ids,
        author_residue_ids=author_ids,
        aatype=aatype,
        residue_index=residue_index,
        chain_index=chain_index,
        backbone=backbone,
        backbone_mask=backbone_mask,
    )


def target_mask(structure: ProteinStructure, target_chains: Sequence[str]) -> np.ndarray:
    if not target_chains:
        return np.ones(structure.num_residues, dtype=bool)
    wanted = {_normalise_chain_id(value) for value in target_chains}
    return np.asarray([chain_id in wanted for chain_id in structure.chain_ids], dtype=bool)

