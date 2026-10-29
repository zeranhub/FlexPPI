"""Deduplicated all-chain search inputs, without reading model predictions."""
from pathlib import Path
import sys
import json
import hashlib
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parents[1]/"reproducibility/src"))
from flexppi.data.structures import parse_pdb
from run_corrected_static import corrected_records

def main():
    out=HERE/"sequence_audit"
    out.mkdir(exist_ok=True)
    primary=corrected_records()
    external=[json.loads(x) for x in (HERE/"external_thermopcd/ready_test.jsonl").read_text(encoding="utf-8").splitlines()]
    rows=[]
    sequences={}
    for role,records in [("training",primary),("external",external)]:
        for r in records:
            s=parse_pdb(r["structure_path"])
            for chain in s.unique_chains():
                sequence=s.chain_sequence(chain)
                key="seq_"+hashlib.sha256(sequence.encode()).hexdigest()
                sequences[key]=sequence
                rows.append(dict(role=role,sample_id=r["sample_id"],chain=chain,sequence_id=key,
                    length=len(sequence),sequence=sequence,dataset=r["dataset"],fold=r["metadata"].get("fold",-1),
                    biological_system=r["metadata"].get("biological_system",r["source_group"])))
    long={key:seq for key,seq in sequences.items() if len(seq)>=20}
    (out/"long_unique.fasta").write_text("".join(f">{key}\n{seq}\n" for key,seq in sorted(long.items())),encoding="ascii")
    (out/"chain_inventory.json").write_text(json.dumps(rows,indent=2),encoding="utf-8")
    metadata={"training_records":len(primary),"external_records":len(external),"chain_observations":len(rows),
        "unique_long_sequences":len(long),"unique_short_sequences":sum(len(x)<20 for x in sequences.values()),
        "fasta_sha256":hashlib.sha256((out/"long_unique.fasta").read_bytes()).hexdigest(),
        "inputs":"Frozen corrected training and prespecified ThermoPCD; no labels or predictions used for membership decisions."}
    (out/"search_input_audit.json").write_text(json.dumps(metadata,indent=2),encoding="utf-8")
    print(json.dumps(metadata))

if __name__=="__main__":main()
