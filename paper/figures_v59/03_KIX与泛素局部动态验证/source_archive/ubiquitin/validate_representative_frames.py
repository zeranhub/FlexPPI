"""Verify bundled PDB side-chain angles against the archived raw MD angles."""
import csv
from pathlib import Path

import MDAnalysis as mda
from MDAnalysis.lib.distances import calc_dihedrals
import numpy as np

BUNDLE = Path(__file__).resolve().parent
with (BUNDLE / "representative_frames.csv").open(encoding="utf-8", newline="") as handle:
    sources = list(csv.DictReader(handle))
definitions = [
    ["N", "CA", "CB", "CG"],
    ["CA", "CB", "CG", "CD"],
    ["CB", "CG", "CD", "OE1"],
]
checks = []
for row in sources:
    protein = mda.Universe(str(BUNDLE / row["pdb_path"]))
    residue = protein.select_atoms(f'segid A and resid {row["input_residue_id"]}')
    assert len(residue.residues) == 1 and residue.residues[0].resname == "GLN"
    measured = []
    differences = []
    for i, definition in enumerate(definitions, 1):
        points = [residue.select_atoms(f"name {name}").positions for name in definition]
        assert all(point.shape == (1, 3) for point in points)
        angle = float(np.rad2deg(calc_dihedrals(*points))[0])
        expected = float(row[f"chi{i}_degrees"])
        difference = float(abs((angle - expected + 180) % 360 - 180))
        measured.append(angle)
        differences.append(difference)
    assert max(differences) < 0.15, (row["pdb_path"], differences)
    checks.append(dict(residue=row["residue"], joint_code=row["joint_code"],
                       frame_zero_based=row["trajectory_frame_zero_based"],
                       pdb_path=row["pdb_path"],
                       measured_chi1_degrees=measured[0], measured_chi2_degrees=measured[1],
                       measured_chi3_degrees=measured[2],
                       max_circular_angle_difference_degrees=max(differences),
                       tolerance_degrees=0.15, status="PASS",
                       note="small differences are due to PDB 0.001 A coordinate rounding"))
with (BUNDLE / "frame_geometry_audit.csv").open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(checks[0]))
    writer.writeheader()
    writer.writerows(checks)
print(f'PASS: all {len(checks)} real-frame side chains preserve their archived MD angles; '
      f'max error {max(row["max_circular_angle_difference_degrees"] for row in checks):.4f} degrees.')
