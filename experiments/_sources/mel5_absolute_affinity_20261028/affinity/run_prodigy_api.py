"""Extract unrounded predictions from the installed, unmodified PRODIGY API."""
from pathlib import Path
import csv
from datetime import datetime, timezone
from importlib.metadata import version
import json
import math
import sys

from prodigy_prot.modules.parsers import parse_structure
from prodigy_prot.modules.prodigy import Prodigy

ROOT = Path(__file__).resolve().parent
EXPERIMENT = {
    "3hg1": ("MEL5 / ELA", 18e-6),
    "4jff": ("alpha24beta17 / ELA", 0.6e-9),
    "4jfd": ("alpha24beta17 / G4A", 36e-6),
    "4jfe": ("alpha24beta17 / I7A", 31e-6),
}

rows = []
raw = []
audit = json.loads((ROOT / "structure_audit.json").read_text(encoding="utf-8"))
for pdb, (state, experimental_kd) in EXPERIMENT.items():
    path = ROOT / "structures" / f"{pdb}_crystal_protein.pdb"
    models, _, _ = parse_structure(str(path))
    assert len(models) == 1
    expected = next(v for v in audit["records"] if v["pdb"] == pdb)
    observed = {c.id: len(list(c.get_residues())) for c in models[0].get_chains()}
    assert observed == {c: v["residues"] for c, v in expected["chains"].items()}
    model = Prodigy(model=models[0], name=pdb, selection=["D,E", "A,B,C"], temp=25.0)
    model.predict(distance_cutoff=5.5, acc_threshold=0.05)
    prediction = model.as_dict()
    raw.append({"pdb": pdb, "parsed_chain_residue_counts": observed, **prediction})
    rows.append({
        "pdb": pdb, "state": state, "model": "PRODIGY 2.4.0", "input_kind": "crystal",
        "temperature_C": 25.0, "TCR_chains": "DE", "pMHC_chains": "ABC",
        "experimental_KD_M": experimental_kd,
        "experimental_pKD": -math.log10(experimental_kd),
        "predicted_DG_kcal_mol": prediction["ba_val"],
        "predicted_KD_M": prediction["kd_val"],
        "predicted_pKD": -math.log10(prediction["kd_val"]),
        "pKD_error_pred_minus_exp": -math.log10(prediction["kd_val"]) + math.log10(experimental_kd),
        "intermolecular_contacts": prediction["ICs"],
    })
output = ROOT / "absolute_affinity_results.csv"
with output.open("w", encoding="utf-8-sig", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
(ROOT / "prodigy_full_precision.json").write_text(json.dumps({
    "calculated_at_UTC": datetime.now(timezone.utc).isoformat(),
    "python": sys.version,
    "versions": {p: version(p) for p in ["prodigy-prot", "freesasa", "biopython", "numpy"]},
    "raw_predictions": raw, "results": rows,
    "experimental_source": "Madura et al. JBC 2013 Table 2, https://doi.org/10.1074/jbc.M113.464560; SPR at 25 C",
}, indent=2), encoding="utf-8")
print(json.dumps(rows, indent=2))
