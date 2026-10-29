"""Attach cached official protein names to every evaluated protein record."""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import pandas as pd

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[1]
ANNOTATIONS = ROOT.parent / "functional_annotations_complete"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


records = pd.read_csv(OUT / "all_2203_proteins_ranked.csv", encoding="utf-8-sig", float_precision="round_trip")
cached = pd.read_csv(ANNOTATIONS / "record_annotations.csv", encoding="utf-8-sig").set_index("sample_id")
entries, sources = {}, {}
for path in sorted(ANNOTATIONS.glob("batch_*.json")):
    raw = json.loads(path.read_text(encoding="utf-8"))
    for entry in raw.get("response", {}).get("data", {}).get("entries", []) or []:
        if entry:
            entries[entry["rcsb_id"]] = entry
            sources[entry["rcsb_id"]] = path
rows = []
used_paths = {ANNOTATIONS / "record_annotations.csv"}
for record in records.itertuples(index=False):
    annotation = cached.loc[record.sample_id]
    actual_name = ""
    uniprot, entity_ids = [], []
    mapping = "entry_metadata_only_not_certified_chain_identity"
    source = ANNOTATIONS / "record_annotations.csv"
    if record.dataset == "atlas":
        pdb_id, chain = record.sample_id.rsplit("_", 1)
        entry = entries.get(pdb_id.upper(), {})
        matches = [entity for entity in entry.get("polymer_entities", []) or []
            if chain in (entity.get("rcsb_polymer_entity_container_identifiers") or {}).get("auth_asym_ids", [])]
        actual_name = ";".join(sorted({(entity.get("rcsb_polymer_entity") or {}).get("pdbx_description", "") for entity in matches}))
        uniprot = sorted({accession for entity in matches for accession in
            (entity.get("rcsb_polymer_entity_container_identifiers") or {}).get("uniprot_ids", []) or []})
        entity_ids = sorted(entity["rcsb_id"] for entity in matches)
        if len(matches) == 1:
            mapping = "ATLAS_original_PDB_author_chain_to_official_entity"
        if entry:
            source = sources[pdb_id.upper()]
            used_paths.add(source)
    rows.append(dict(sample_id=record.sample_id, dataset=record.dataset, n_residues=record.n_residues,
        original_fold=record.original_fold, actual_input_chain_count=record.n_chains,
        protein_name=actual_name, name_mapping=mapping, official_entity_ids=";".join(entity_ids),
        uniprot_ids=";".join(uniprot), source_entry_descriptions=annotation.protein_descriptions,
        source_pdb_classification=annotation.pdb_classification, annotation_source=str(source),
        full_rho=getattr(record, "_0", None)))
names = pd.DataFrame(rows).drop(columns="full_rho")
selected = ["sample_id", "spearman::Full", "spearman::backflip-2.1 (chains)", "mae::Full",
    "mae::backflip-2.1 (chains)", "delta_rho_full_minus_bf2"]
names = names.merge(records[selected], on="sample_id", validate="one_to_one")
names.to_csv(OUT / "all_2203_test_proteins_with_names.csv", index=False)
names.head(30).to_csv(OUT / "worst30_protein_identities.csv", index=False)
audit = {"created_utc": datetime.now(timezone.utc).isoformat(), "script_sha256": sha(Path(__file__)),
    "input_sha256": {str(OUT / "all_2203_proteins_ranked.csv"): sha(OUT / "all_2203_proteins_ranked.csv"),
                     **{str(path): sha(path) for path in sorted(used_paths)}},
    "records": len(names), "name_mapping_counts": names.name_mapping.value_counts().to_dict(),
    "scope": "Names describe actual ATLAS input author chains; source entry summaries are retained separately. Non-ATLAS chain identities rely on the frozen biological grouping evidence, not this entry-level name lookup.",
    "output_sha256": {name: sha(OUT / name) for name in ("all_2203_test_proteins_with_names.csv", "worst30_protein_identities.csv")}}
(OUT / "test_protein_annotation_audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
print(names.head(8)[["sample_id", "n_residues", "protein_name", "spearman::Full", "spearman::backflip-2.1 (chains)"]].to_string(index=False))
