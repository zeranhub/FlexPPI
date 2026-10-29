"""Structural pre-screen; these contact counts are not model predictions."""
from pathlib import Path
import csv
import hashlib
import json
import sys
import types

import numpy as np

SOURCE = Path(r"E:\Project\FlexPPI\paper\revision_v14_20261010\section33_expansion_20261025\biology_story_first_20261025\current_model_cam_predictions_20261025\family_heldout_model\src")
# Load the exact frozen parser without importing the neural network or Torch.
for package_name, package_path in (
    ("flexppi", SOURCE / "flexppi"),
    ("flexppi.data", SOURCE / "flexppi" / "data"),
):
    package = types.ModuleType(package_name)
    package.__path__ = [str(package_path)]
    sys.modules[package_name] = package
from flexppi.constants import THREE_TO_ONE
from flexppi.data.structures import parse_pdb

ROOT = Path(__file__).resolve().parent
JOBS = [
    ("1MXE", 1, "E", "A", None),
    ("5T35", 1, "A", "D", "759"),
    ("5T35", 2, "E", "H", "759"),
    ("6W8I", 1, "A", "D", "TKY"),
    ("6W8I", 2, "B", "E", "TKY"),
    ("6W8I", 3, "C", "F", "TKY"),
    ("6W7O", 1, "A", "C", "TL7"),
    ("6W7O", 2, "B", "D", "TL7"),
]


def read_atoms(path):
    atoms = []
    seen = set()
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("ENDMDL"):
            break
        if line[:6].strip() not in {"ATOM", "HETATM"}:
            continue
        if line[16] not in {" ", "A", "1"}:
            continue
        element = line[76:78].strip()
        if element in {"H", "D"}:
            continue
        key = (line[21], line[22:27].strip(), line[12:16].strip())
        if key in seen:
            continue
        seen.add(key)
        atoms.append({
            "chain": line[21], "residue": line[22:27].strip(),
            "name": line[17:20].strip(), "atom": line[12:16].strip(),
            "protein": line[17:20].strip() in THREE_TO_ONE,
            "xyz": [float(line[30:38]), float(line[38:46]), float(line[46:54])],
        })
    return atoms


def nearest_residue_distances(source, other):
    if not other:
        return {atom["residue"]: float("inf") for atom in source}
    other_xyz = np.asarray([atom["xyz"] for atom in other])
    distances = {}
    for atom in source:
        distance = float(np.sqrt(((other_xyz - atom["xyz"]) ** 2).sum(axis=1)).min())
        residue = atom["residue"]
        distances[residue] = min(distance, distances.get(residue, float("inf")))
    return distances


def input_signature(path):
    structure = parse_pdb(str(path))
    digest = hashlib.sha256()
    for field in ("aatype", "residue_index", "chain_index", "backbone", "backbone_mask"):
        digest.update(getattr(structure, field).tobytes())
    return digest.hexdigest(), structure


rows = []
details = []
for pdb, assembly, target_chain, partner_chain, ligand_name in JOBS:
    path = ROOT / f"{pdb}_assembly{assembly}.pdb"
    atoms = read_atoms(path)
    target = [a for a in atoms if a["protein"] and a["chain"] == target_chain]
    partner = [a for a in atoms if a["protein"] and a["chain"] == partner_chain]
    ligand = [a for a in atoms if not a["protein"] and a["name"] == ligand_name]
    protein_only_path = ROOT / f"{pdb}_assembly{assembly}_protein_only.pdb"
    protein_only_path.write_text("\n".join(
        line for line in path.read_text(encoding="utf-8").splitlines()
        if line[:6].strip() in {"ATOM", "HETATM"}
        and line[17:20].strip() in THREE_TO_ONE
    ) + "\nEND\n", encoding="utf-8")
    original_signature, parsed = input_signature(path)
    stripped_signature, _ = input_signature(protein_only_path)
    target_distance = nearest_residue_distances(target, partner)
    target_ligand = nearest_residue_distances(target, ligand)
    partner_distance = nearest_residue_distances(partner, target)
    partner_ligand = nearest_residue_distances(partner, ligand)
    target_indices = parsed.indices_for_chain(target_chain)
    partner_indices = parsed.indices_for_chain(partner_chain)
    row = {
        "pdb": pdb, "assembly": assembly,
        "target_chain": target_chain, "partner_chain": partner_chain,
        "target_modeled_residues": len(target_indices),
        "partner_modeled_residues": len(partner_indices),
        "target_incomplete_N_CA_C": int((~parsed.frame_mask[target_indices]).sum()),
        "partner_incomplete_N_CA_C": int((~parsed.frame_mask[partner_indices]).sum()),
        "ligand_heavy_atoms": len(ligand),
        "protein_input_identical_after_hetero_removal": original_signature == stripped_signature,
    }
    for threshold in (4.0, 4.5, 5.0):
        suffix = str(threshold).replace(".", "p")
        row[f"target_protein_contact_{suffix}A"] = sum(d <= threshold for d in target_distance.values())
        row[f"partner_protein_contact_{suffix}A"] = sum(d <= threshold for d in partner_distance.values())
        row[f"target_ligand_contact_{suffix}A"] = sum(d <= threshold for d in target_ligand.values())
        row[f"partner_ligand_contact_{suffix}A"] = sum(d <= threshold for d in partner_ligand.values())
        row[f"target_ligand_without_protein_contact_{suffix}A"] = sum(
            d <= threshold and target_distance[r] > threshold for r, d in target_ligand.items())
        row[f"partner_ligand_without_protein_contact_{suffix}A"] = sum(
            d <= threshold and partner_distance[r] > threshold for r, d in partner_ligand.items())
    rows.append(row)
    for side, protein_atoms, pd, ld in (
        ("target", target, target_distance, target_ligand),
        ("partner", partner, partner_distance, partner_ligand),
    ):
        names = {a["residue"]: a["name"] for a in protein_atoms}
        for residue, distance in pd.items():
            details.append({
                "pdb": pdb, "assembly": assembly, "side": side,
                "residue": residue, "amino_acid": names[residue],
                "nearest_partner_heavy_atom_A": distance,
                "nearest_PROTAC_heavy_atom_A": ld[residue] if ligand else None,
            })

with (ROOT / "structural_input_audit.csv").open("w", newline="", encoding="utf-8-sig") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
with (ROOT / "residue_contact_details.csv").open("w", newline="", encoding="utf-8-sig") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(details[0]))
    writer.writeheader()
    writer.writerows(details)
(ROOT / "structural_input_audit.json").write_text(json.dumps(rows, indent=2), encoding="utf-8")
for row in rows:
    print(json.dumps({k: v for k, v in row.items() if "4p0" not in k and "5p0" not in k}))
