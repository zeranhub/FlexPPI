"""Join published SPR values with both independently computed absolute outputs."""
from pathlib import Path
import csv
import json
import math
from hashlib import sha256

ROOT = Path(__file__).resolve().parent
path = ROOT / "absolute_affinity_results.csv"
with path.open(encoding="utf-8-sig", newline="") as handle:
    reader = csv.DictReader(handle)
    columns = reader.fieldnames
    prodigy = [row for row in reader if row["model"].startswith("PRODIGY")]
assert len(prodigy) == 4
records = list(prodigy)
audits = []
for base in prodigy:
    pdb = base["pdb"]
    audit_path = ROOT / "raw_outputs/proaffinity" / pdb / "audit.json"
    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    assert audit["all_chain_sequences_match_input"]
    assert len(audit["graph_checks"]) == 3
    assert all(x["valid"] for x in audit["graph_checks"])
    assert audit["state_dict_loading"] == [{"missing_keys": [], "unexpected_keys": []}]
    prediction = float(audit["prediction_as_reported_pKa"])
    row = dict(base)
    row.update(model="ProAffinity-GNN (simplified, corrected)", temperature_C="",
               predicted_DG_kcal_mol="", predicted_KD_M=10 ** (-prediction),
               predicted_pKD=prediction,
               pKD_error_pred_minus_exp=prediction - float(base["experimental_pKD"]),
               intermolecular_contacts="")
    records.append(row)
    audits.append({"pdb": pdb, "audit_sha256": sha256(audit_path.read_bytes()).hexdigest(),
                   "upstream_bug_fix": audit["upstream_bug_fix"],
                   "GNN_checkpoint_sha256": audit["checkpoint_sha256"]})
with path.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=columns)
    writer.writeheader()
    writer.writerows(records)
summary = []
for model in ["PRODIGY 2.4.0", "ProAffinity-GNN (simplified, corrected)"]:
    lookup = {row["pdb"]: row for row in records if row["model"] == model}
    for before, after, label in [("3hg1", "4jff", "TCR maturation"),
                                  ("4jff", "4jfd", "G4A peptide substitution"),
                                  ("4jff", "4jfe", "I7A peptide substitution")]:
        exp = float(lookup[after]["experimental_pKD"]) - float(lookup[before]["experimental_pKD"])
        predicted = float(lookup[after]["predicted_pKD"]) - float(lookup[before]["predicted_pKD"])
        summary.append({"model": model, "contrast": after + "-" + before,
                        "description": label, "experimental_pKD_after_minus_before": exp,
                        "predicted_pKD_after_minus_before": predicted,
                        "experimental_KD_after_over_before": 10 ** (-exp),
                        "predicted_KD_after_over_before": 10 ** (-predicted),
                        "same_direction": exp * predicted > 0,
                        "note": "Differences derived only after separate absolute predictions; no dedicated ddG model."})
(ROOT / "absolute_prediction_comparison.json").write_text(json.dumps({
    "rows": records, "contrasts": summary, "GNN_audits": audits,
    "PRODIGY": "Official prodigy-prot 2.4.0 / FreeSASA 2.2.1, no approximation substitutions.",
    "ProAffinity_scope": "Authors simplified inference checkpoint, single upstream variable correction and explicit PDBQT preparation differences. Report implementation conditions; do not assert general model or family-held-out performance.",
}, ensure_ascii=False, indent=2), encoding="utf-8")
print(json.dumps(summary, ensure_ascii=False, indent=2))
