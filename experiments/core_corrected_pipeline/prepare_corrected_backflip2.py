"""Same corrected static structures for official per-chain BackFlip2 inference."""
from pathlib import Path
import json
import csv
import hashlib
import shutil
import tarfile
import numpy as np
from run_corrected_static import corrected_records
from prepare_corrected_labels import parse_pdb,sha
HERE=Path(__file__).resolve().parent

def main():
    out=HERE/"backflip2_corrected_inputs"
    if (out/"input_freeze.json").exists():raise RuntimeError("Corrected BackFlip2 inputs already frozen")
    (out/"inputs").mkdir(parents=True,exist_ok=True)
    manifest=[]
    records=[]
    external=[json.loads(x) for x in (HERE/"external_thermopcd/ready_test.jsonl").read_text(encoding="utf-8").splitlines()]
    mapping_file=(out/"residue_mapping.csv").open("w",encoding="utf-8",newline="")
    writer=csv.DictWriter(mapping_file,fieldnames=["role","sample_id","chain","input_id","input_index","record_residue_index","author_residue_id","amino_acid"])
    writer.writeheader()
    for role,cohort in [("training",corrected_records()),("external",external)]:
        for record in cohort:
            s=parse_pdb(record["structure_path"])
            raw=Path(record["structure_path"]).read_text(encoding="utf-8").splitlines(True)
            stop=next((i for i,line in enumerate(raw) if line.startswith("ENDMDL")),len(raw))
            raw=raw[:stop]
            for chain in s.unique_chains():
                ids=s.indices_for_chain(chain)
                number={(s.author_residue_ids[i],s.residue_names[i]):j+1 for j,i in enumerate(ids)}
                key="bf2_"+hashlib.sha256((record["sample_id"]+"|"+chain).encode()).hexdigest()[:24]
                path=out/"inputs"/(key+".pdb")
                lines=[]
                seen=set()
                for line in raw:
                    if not line.startswith(("ATOM  ","HETATM")) or len(line)<54 or line[16] not in [" ","A"]:continue
                    observed_chain=line[21].strip() or "_"
                    author=line[22:26].strip()+line[26].strip()
                    residue_key=(author,line[17:20])
                    if observed_chain!=chain or residue_key not in number:continue
                    atom_key=(residue_key,line[12:16].strip())
                    if atom_key in seen:continue
                    seen.add(atom_key)
                    lines.append("ATOM  "+f"{len(lines)+1:5d}"+line[11:16]+" "+line[17:21]+"A"+f"{number[residue_key]:4d}"+" "+line[27:])
                path.write_text("".join(lines)+"END\n",encoding="utf-8")
                isolated=parse_pdb(str(path))
                assert isolated.sequence==s.chain_sequence(chain)
                assert np.array_equal(isolated.backbone_mask,s.backbone_mask[ids])
                assert np.allclose(isolated.backbone,s.backbone[ids],atol=0,rtol=0)
                manifest.append(dict(input_id=key,sample_id=record["sample_id"],chain=chain,role=role,
                    n_residues=len(ids),pdb_file="inputs/"+path.name,input_sha256=sha(path)))
                for j,i in enumerate(ids):writer.writerow(dict(role=role,sample_id=record["sample_id"],chain=chain,
                    input_id=key,input_index=j,record_residue_index=int(i),author_residue_id=s.author_residue_ids[i],amino_acid=s.sequence[i]))
            records.append(dict(role=role,record=record,input_sha256=sha(record["structure_path"]),label_sha256=sha(record["label_path"])))
    mapping_file.close()
    with (out/"input_manifest.csv").open("w",encoding="utf-8",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(manifest[0]));writer.writeheader();writer.writerows(manifest)
    (out/"record_manifest.json").write_text(json.dumps(records,indent=2),encoding="utf-8")
    shutil.copy2(HERE.parent/"backflip2/run_inference.py",out/"run_inference.py")
    audit=dict(records=len(records),chains=len(manifest),training_records=3063,external_records=32,
        input_manifest_sha256=sha(out/"input_manifest.csv"),mapping_sha256=sha(out/"residue_mapping.csv"),
        model="Official BackFlip2 tag2.1, unmodified weights, independent chain inference on CPU, covariance nm2 to RMSF Angstrom",
        input_check="Each isolated chain preserves corrected sequence, backbone masks and all backbone coordinates exactly; identity mapping retained.",
        failure_rule="Keep and disclose every inference failure; common-mask comparisons require matching actual prediction coverage. No selection on errors.")
    (out/"input_freeze.json").write_text(json.dumps(audit,indent=2),encoding="utf-8")
    with tarfile.open(HERE/"backflip2_corrected_inputs.tar.gz","w:gz") as tar:
        for p in out.iterdir():tar.add(p,arcname=p.name)
    print(json.dumps(audit,indent=2))

if __name__=="__main__":main()
