"""Validate external chains have no search/exact hits in the purged training set."""
from pathlib import Path
import json
from collections import defaultdict
HERE=Path(__file__).resolve().parent
OUT=HERE/"sequence_audit"
chains=json.loads((OUT/"chain_inventory.json").read_text(encoding="utf-8"))
allowed={json.loads(x)["sample_id"] for x in (OUT/"purged_training_local.jsonl").read_text(encoding="utf-8").splitlines()}
members=defaultdict(set)
neighbors=defaultdict(set)
for r in chains:
    if r["sample_id"] in allowed:members[r["sequence_id"]].add(r["sample_id"])
    neighbors[r["sequence_id"]].add(r["sequence_id"])
with (OUT/"all_chain_similarity.tsv").open() as f:
    for line in f:
        a=line.rstrip().split("\t")
        if float(a[2])>=.3 and min((int(a[5])-int(a[4])+1)/int(a[6]),(int(a[8])-int(a[7])+1)/int(a[9]))>=.7:
            neighbors[a[0]].add(a[1]);neighbors[a[1]].add(a[0])
rows=[]
for r in chains:
    if r["role"]!="external":continue
    candidates=neighbors[r["sequence_id"]] if r["length"]>=20 else {r["sequence_id"]}
    hits=set().union(*(members[key] for key in candidates))
    assert not hits,(r["sample_id"],r["chain"],hits)
    rows.append(dict(sample_id=r["sample_id"],chain=r["chain"],length=r["length"],purged_training_hits=len(hits)))
report={"status":"passed","strict_training_records":len(allowed),"external_chain_count":len(rows),
        "maximum_external_chain_hits_after_purge":0,"chains":rows,
        "scope":"Absence under the declared all-chain MMseqs search and exact-short rule; not proof of absence of every evolutionary relationship."}
(OUT/"strict_overlap_validation.json").write_text(json.dumps(report,indent=2),encoding="utf-8")
print("Strict external chain overlap verified:",len(rows),"chains, no retained training hits")
