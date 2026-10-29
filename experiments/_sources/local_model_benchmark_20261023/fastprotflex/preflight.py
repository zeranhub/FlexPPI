"""Size-only, label-free preflight for the frozen complete structures."""
from pathlib import Path
import hashlib
import json
import time

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
AA = set("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split())
records = pd.read_csv(HERE.parent / "cohort/records.csv", keep_default_na=False)
rows = []
started = time.perf_counter()
for row in records.itertuples(index=False):
    data = Path(row.structure_path).read_bytes()
    assert hashlib.sha256(data).hexdigest() == row.structure_sha256
    n = 0
    ca = 0
    alternate = 0
    for line in data.decode("ascii").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        if line[17:20].strip() not in AA or line[76:78].strip().upper() in {"H", "D"}:
            continue
        n += 1
        ca += line[12:16].strip() == "CA"
        alternate += bool(line[16:17].strip())
    assert ca == row.n_residues, (row.sample_id, ca, row.n_residues)
    rows.append({
        "sample_id": row.sample_id, "group": row.main_group_candidate,
        "n_residues": row.n_residues, "n_chains": row.n_chains,
        "protein_heavy_atoms": n, "alternate_location_atoms": alternate,
        "official_branch": "GDVallatonce" if n <= 10000 else "GDVperpartes_new2",
    })
sizes = pd.DataFrame(rows)
sizes.to_csv(HERE / "preflight_sizes.csv", index=False)
smoke = pd.read_csv(HERE / "smoke_results.csv")
coef = np.linalg.lstsq(
    np.column_stack([np.ones(len(smoke)), smoke.input_protein_heavy_atoms, smoke.input_protein_heavy_atoms**2]),
    smoke.wall_seconds, rcond=None,
)[0]
estimate = np.maximum(1.0, np.column_stack([
    np.ones(len(sizes)), sizes.protein_heavy_atoms, sizes.protein_heavy_atoms**2,
]) @ coef)
audit = {
    "records": len(rows), "total_ca_residues": int(sizes.n_residues.sum()),
    "structures_sha256_verified": len(rows),
    "max_protein_heavy_atoms": int(sizes.protein_heavy_atoms.max()),
    "official_branch_counts": sizes.official_branch.value_counts().to_dict(),
    "alternate_location_atoms": int(sizes.alternate_location_atoms.sum()),
    "smoke_based_estimate_wall_hours": float(estimate.sum() / 3600),
    "estimate_limit": "Three smoke structures only; excludes initial import and includes R startup per structure. Extrapolated, not measured full-run time.",
    "largest_structure_ids": sizes.nlargest(8, "protein_heavy_atoms").to_dict(orient="records"),
    "preflight_wall_seconds": time.perf_counter() - started,
}
(HERE / "preflight_audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
print(json.dumps(audit, indent=2), flush=True)
