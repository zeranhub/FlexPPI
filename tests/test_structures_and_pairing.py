from pathlib import Path

from flexppi.data.pairing import global_sequence_alignment, match_apo_to_holo
from flexppi.data.structures import parse_pdb


def atom_line(serial, name, residue, chain, residue_id, x, y, z):
    return (
        "ATOM  {serial:5d} {name:^4s} {residue:>3s} {chain:1s}{residue_id:4d}    "
        "{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           {element:>2s}\n"
    ).format(
        serial=serial,
        name=name,
        residue=residue,
        chain=chain,
        residue_id=residue_id,
        x=x,
        y=y,
        z=z,
        element=name[0],
    )


def write_chain(path: Path, chains):
    lines = []
    serial = 1
    residue_names = {"A": "ALA", "G": "GLY", "S": "SER", "V": "VAL"}
    for chain, sequence, y_offset in chains:
        for index, amino_acid in enumerate(sequence, start=1):
            x = float(index * 3)
            for atom, dx, dy in [("N", -1.0, 0.2), ("CA", 0.0, 0.0), ("C", 1.0, 0.1), ("O", 1.5, 0.5)]:
                lines.append(
                    atom_line(
                        serial, atom, residue_names[amino_acid], chain, index,
                        x + dx, y_offset + dy, 0.0,
                    )
                )
                serial += 1
    lines.append("END\n")
    path.write_text("".join(lines), encoding="utf-8")


def test_parse_and_chain_matching(tmp_path):
    apo_path = tmp_path / "apo.pdb"
    holo_path = tmp_path / "holo.pdb"
    write_chain(apo_path, [("A", "AGSV", 0.0)])
    write_chain(holo_path, [("X", "VVVV", 20.0), ("R", "AGSV", 0.0)])
    apo = parse_pdb(str(apo_path))
    holo = parse_pdb(str(holo_path))
    assert apo.sequence == "AGSV"
    assert holo.unique_chains() == ["X", "R"]
    apo_indices, holo_indices, chains, identity, coverage = match_apo_to_holo(apo, holo)
    assert apo_indices == [0, 1, 2, 3]
    assert holo_indices == [4, 5, 6, 7]
    assert chains == ["R"]
    assert identity == 1.0
    assert coverage == 1.0


def test_alignment_keeps_exact_correspondence_across_gap():
    alignment = global_sequence_alignment("AGSV", "AGGSV")
    assert alignment.indices_a == [0, 1, 2, 3]
    # The two consecutive G residues make the gap position mathematically
    # ambiguous; either exact-G correspondence is a valid global alignment.
    assert alignment.indices_b in ([0, 1, 3, 4], [0, 2, 3, 4])
    assert alignment.coverage == 1.0
