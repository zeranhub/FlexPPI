"""Download public source records and audit candidate PDBs; no model inference.

This intentionally does not construct missing residues or pretend a stripped
complex is an experimental apo structure. Generated files are provenance data.
"""
from __future__ import annotations

import concurrent.futures
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "artifacts/apo_holo_mechanistic_literature_v2_20260913"
PDBS = ["1BFE", "1BE9", "5HEB", "6QJD", "3QH3", "1QRN", "4FTV", "2AGH", "6MHF", "1AO7"]
DOIS = [
    "10.1073/pnas.0904492106", "10.1074/jbc.RA118.003832",
    "10.1371/journal.pone.0090030", "10.1101/2025.05.28.656549",
    "10.1038/ncomms15260", "10.1371/journal.pcbi.1009977",
    "10.1073/pnas.1906658116", "10.1038/ncomms15163",
    "10.1038/s41598-021-88432-z", "10.1021/acsomega.4c03375",
    "10.1038/s41467-026-73600-4",
    "10.1371/journal.pcbi.1002429", "10.1073/pnas.1705311114",
]


def fetch(item):
    name, url = item
    target = OUT / "sources" / name
    target.parent.mkdir(parents=True, exist_ok=True)
    record = {"name": name, "url": url, "retrieved_utc": datetime.now(timezone.utc).isoformat()}
    try:
        if not target.exists():
            response = requests.get(url, timeout=(10, 35))
            response.raise_for_status()
            content = response.content
            if name.endswith(".json"):
                json.loads(content)
            elif not any(line.startswith(b"ATOM  ") for line in content.splitlines()):
                raise ValueError("Response has no PDB ATOM records")
            target.write_bytes(content)
        record.update(status="ok", path=str(target), sha256=hashlib.sha256(target.read_bytes()).hexdigest())
    except Exception as exc:
        record.update(status="failed", error=str(exc))
    return record


def audit_pdb(pdb_id):
    path = OUT / "sources" / f"{pdb_id}.pdb"
    if not path.exists():
        return {"pdb_id": pdb_id, "status": "not_downloaded"}
    lines = path.read_text(encoding="utf-8").splitlines()
    chains, heterogens = {}, set()
    for line in lines:
        if line.startswith("ENDMDL"):
            break
        if line.startswith("HETATM"):
            heterogens.add(line[17:20].strip())
        if not line.startswith("ATOM  ") or line[16] not in " A1":
            continue
        chain = line[21].strip() or "_"
        residues = chains.setdefault(chain, {})
        key = line[22:27].strip()
        residue = residues.setdefault(key, {"name": line[17:20].strip(), "atoms": set()})
        residue["atoms"].add(line[12:16].strip())
    result = {"pdb_id": pdb_id,
              "title": " ".join(line[10:].strip() for line in lines if line.startswith("TITLE ")),
              "model_count": sum(line.startswith("MODEL ") for line in lines) or 1,
              "first_model_heterogens": sorted(heterogens),
              "missing_residue_records": [line for line in lines if line.startswith("REMARK 465")],
              "sequence_discrepancies": [line for line in lines if line.startswith("SEQADV")],
              "chains": []}
    for chain, residues in chains.items():
        result["chains"].append({"author_chain": chain, "observed_residues": len(residues),
            "first_residue": next(iter(residues)), "last_residue": list(residues)[-1],
            "residue_names": [r["name"] for r in residues.values()],
            "incomplete_backbone": [k for k, r in residues.items() if not {"N", "CA", "C", "O"} <= r["atoms"]]})
    return result


def main():
    jobs = [(f"{p}.pdb", f"https://files.rcsb.org/download/{p}.pdb") for p in PDBS]
    jobs += [(f"{p}.entry.json", f"https://data.rcsb.org/rest/v1/core/entry/{p}") for p in PDBS]
    jobs += [(f"doi_{i:02d}.json", f"https://api.crossref.org/works/{d}") for i, d in enumerate(DOIS, 1)]
    with concurrent.futures.ThreadPoolExecutor(max_workers=6) as pool:
        records = list(pool.map(fetch, jobs))
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "download_manifest.json").write_text(json.dumps(records, indent=2, ensure_ascii=False), encoding="utf-8")
    audits = [audit_pdb(p) for p in PDBS]
    (OUT / "structure_audit.json").write_text(json.dumps(audits, indent=2, ensure_ascii=False), encoding="utf-8")
    for audit in audits:
        print(audit["pdb_id"], audit.get("title", audit.get("status")))
        for chain in audit.get("chains", []):
            print(" ", chain["author_chain"], chain["observed_residues"], chain["first_residue"], chain["last_residue"],
                  "missing backbone:", len(chain["incomplete_backbone"]))
    for i, doi in enumerate(DOIS, 1):
        path = OUT / "sources" / f"doi_{i:02d}.json"
        if path.exists():
            metadata = json.loads(path.read_text(encoding="utf-8"))["message"]
            print(doi, metadata.get("type"), metadata.get("title"), metadata.get("published"), metadata.get("relation"))
    print("Downloads:", sum(r["status"] == "ok" for r in records), "/", len(records))


if __name__ == "__main__":
    main()
