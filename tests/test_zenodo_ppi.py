import pytest

from flexppi.data.zenodo_ppi import (
    chain_ids_from_atom_counts,
    normalise_residue_name,
    protein_residue_names,
)


def test_force_field_residue_aliases_preserve_standard_sequence():
    names = ["ACE", "ASH", "CYX", "GLH", "HIE", "NME"]
    assert protein_residue_names(names) == ["ASP", "CYS", "GLU", "HIS"]
    assert normalise_residue_name("ALA") == "ALA"


def test_chain_ids_are_expanded_at_exact_atom_boundary():
    assert chain_ids_from_atom_counts(5, [2, 3], ["A", "B"]) == [
        "A", "A", "B", "B", "B"
    ]


def test_chain_atom_count_mismatch_is_rejected():
    with pytest.raises(ValueError, match="do not sum"):
        chain_ids_from_atom_counts(6, [2, 3], ["A", "B"])
