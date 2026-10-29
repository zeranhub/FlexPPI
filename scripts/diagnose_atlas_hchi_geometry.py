"""Read-only diagnosis of atom pairs rejected by ATLAS Hchi extraction."""
import argparse
import json
from pathlib import Path
import shutil
import sys
import tempfile

import MDAnalysis as mda
import numpy as np
from MDAnalysis.lib.distances import minimize_vectors

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from flexppi.constants import CHI_ATOMS


def diagnose(sid, root):
    reports = []
    pdb = root / sid / (sid + ".pdb")
    topology = mda.Universe(str(pdb))
    reference = topology.atoms.positions.copy()
    pairs = set()
    for residue in topology.residues:
        lookup = {str(atom.name): int(atom.index) for atom in residue.atoms}
        if residue.resname == "ILE" and "CD1" not in lookup and "CD" in lookup:
            lookup["CD1"] = lookup["CD"]
        for definition in CHI_ATOMS.get(residue.resname, []):
            if all(name in lookup for name in definition):
                indices = [lookup[name] for name in definition]
                pairs.update(zip(indices, indices[1:]))
    indices = np.asarray(sorted(pairs))
    for replica in range(1, 4):
        original = root / sid / f"{sid}_R{replica}.xtc"
        with tempfile.TemporaryDirectory(prefix="flexppi_hchi_diagnostic_") as directory:
            temporary = Path(directory) / original.name
            shutil.copyfile(original, temporary)
            universe = mda.Universe(str(pdb), str(temporary))
            bad = {}
            maximum = 0.0
            minimum = float("inf")
            for ts in universe.trajectory:
                vectors = ts.positions[indices[:, 1]] - ts.positions[indices[:, 0]]
                distances = np.linalg.norm(vectors, axis=1)
                minimum = min(minimum, float(distances.min()))
                maximum = max(maximum, float(distances.max()))
                for index in np.flatnonzero((distances < 0.5) | (distances > 3.0)):
                    pair = tuple(indices[index])
                    if pair not in bad:
                        atoms = universe.atoms[list(pair)]
                        bad[pair] = {"atoms": [{"index": int(a.index), "resid": int(a.resid),
                                                "resname": str(a.resname), "name": str(a.name)} for a in atoms],
                                     "pdb_distance": float(np.linalg.norm(reference[pair[1]] - reference[pair[0]])),
                                     "frames_bad": 0, "first_frame": int(ts.frame), "min_distance": float("inf"),
                                     "max_distance": 0.0, "min_image_min": float("inf"), "min_image_max": 0.0,
                                     "box": ts.dimensions.tolist() if ts.dimensions is not None else None}
                    row = bad[pair]
                    distance = float(distances[index])
                    mic = float(np.linalg.norm(minimize_vectors(vectors[index:index + 1], ts.dimensions))) if ts.dimensions is not None else distance
                    row["frames_bad"] += 1
                    row["min_distance"] = min(row["min_distance"], distance)
                    row["max_distance"] = max(row["max_distance"], distance)
                    row["min_image_min"] = min(row["min_image_min"], mic)
                    row["min_image_max"] = max(row["min_image_max"], mic)
            frames = len(universe.trajectory)
            universe.trajectory.close()
            report = {"sample_id": sid, "replica": replica, "frames": frames,
                      "min_distance": minimum, "max_distance": maximum, "bad_pairs": list(bad.values())}
            reports.append(report)
            print(json.dumps(report), flush=True)
    return reports


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sample_ids", nargs="+")
    parser.add_argument("--atlas-root", type=Path, default=Path("E:/data/protein_md_datasets/ATLAS_extracted"))
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    reports = []
    for sample_id in args.sample_ids:
        reports.extend(diagnose(sample_id, args.atlas_root))
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(reports, indent=2) + "\n", encoding="utf-8")
