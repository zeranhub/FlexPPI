# -*- coding: utf-8 -*-
"""Freeze an untrained-source cohort BEFORE extraction or model scoring."""
from pathlib import Path
import sys
import json
import csv
import time
import ctypes
import argparse
import re
from concurrent.futures import ProcessPoolExecutor
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from prepare_corrected_labels import canonical_pdb, parse_pdb, extract, replica_map, sha, _summarise_samples, save_label_file, METRIC_NAMES, staged_trajectory, ALIASES, THREE_TO_ONE

OUT = HERE / "external_thermopcd"
SOURCE = HERE.parents[1] / "artifacts/data/thermopcd_performance_first_candidates.jsonl"
INVENTORY = HERE.parents[1] / "artifacts/thermopcd_audit/local_entries.csv"

def inventory_key(path):
    # The historical CSV contains U+FFFD in place of the degree symbol. Match
    # the otherwise identical ASCII pathname, rejecting any key collision.
    return re.sub(r"[^\x00-\x7f]+", "_", str(Path(path))).lower()

def freeze():
    OUT.mkdir(exist_ok=True)
    if (OUT / "frozen_input_manifest.jsonl").exists():
        raise RuntimeError("External cohort already frozen; refusing overwrite")
    candidates = [json.loads(x) for x in SOURCE.read_text(encoding="utf-8").splitlines()]
    inventory_rows = list(csv.DictReader(INVENTORY.open(encoding="utf-8-sig",newline="")))
    entries = {inventory_key(r["path"]):r for r in inventory_rows}
    if len(entries) != len(inventory_rows):
        raise ValueError("Non-ASCII inventory path keys collide")
    accepted, exclusions = [], []
    for r in candidates:
        try:
            source_dirs = r["metadata"]["source_directories"]
            resolved_dirs = {}
            for d in source_dirs:
                row = entries[inventory_key(d)]
                if "|" in row["temperature_k"] or "|" in row["seed"]:
                    raise ValueError("merged temperature/seed provenance")
                resolved_dirs[inventory_key(d)] = Path(row["path"])
            topologies = []
            resolved_trajectories = []
            for t in r["trajectory_paths"]:
                folder = resolved_dirs[inventory_key(Path(t).parent)]
                trajectory = folder / Path(t).name
                if not trajectory.is_file():
                    raise ValueError("missing trajectory")
                choices = list(folder.glob("*.pdb"))
                preferred = [p for p in choices if p.name.lower().endswith("-m.pdb")]
                preferred += [p for p in choices if "0ns" in p.name.lower() and p not in preferred]
                if not preferred:
                    raise ValueError("no traceable initial topology")
                topologies.append(str(preferred[0]))
                resolved_trajectories.append(str(trajectory))
            original_pdb = resolved_dirs[inventory_key(Path(r["structure_path"]).parent)] / Path(r["structure_path"]).name
            if not original_pdb.exists():
                raise ValueError("original initial structure not resolved")
            target = OUT / "structures" / (r["sample_id"]+".pdb")
            canonical_pdb(original_pdb,target)
            structure = parse_pdb(str(target))
            if structure.num_residues > 1024:
                raise ValueError("corrected length exceeds 1024")
            r["original_structure_path"] = r["structure_path"]
            r["structure_path"] = str(target.resolve())
            r["trajectory_paths"] = resolved_trajectories
            r["replica_topologies"] = topologies
            r["metadata"] = {**r["metadata"],"residue_count":structure.num_residues,
                             "input_sha256":sha(target),"test_role":"external_source_never_used_for_training",
                             "evaluation_version":"thermopcd_external_frozen_20261002"}
            accepted.append(r)
        except Exception as e:
            exclusions.append({"sample_id":r["sample_id"],"reason":str(e)})
    manifest = OUT / "frozen_input_manifest.jsonl"
    manifest.write_text("\n".join(json.dumps(r,ensure_ascii=False) for r in accepted)+"\n",encoding="utf-8")
    protocol = {"timestamp_utc":time.strftime("%Y-%m-%dT%H:%M:%SZ",time.gmtime()),
                "candidate_count":len(candidates),"accepted_count":len(accepted),"excluded_before_labels_or_scores":exclusions,
                "manifest_sha256":sha(manifest),"source_manifest_sha256":sha(SOURCE),"inventory_sha256":sha(INVENTORY),
                "selection":"Existing 310K preferred, else311K production-complex candidate inventory; no merged temperatures/seeds; <=1024 corrected residues; all qualifying records.",
                "label_protocol":"All readable production frames per replica; per-chain CA fit to first frame; joint 3x120deg chi states; replica mean/sample variance; canonical chemical names and sequence identity.",
                "models":"New corrected fivefold full and intrachain, fixed epoch40, seed7; external predictions equally averaged over five models; static S+g+p fit on all corrected training records.",
                "primary_metrics":"System-equal MAE and Spearman for RMSF/Hchi, source-cohort summary; AA-only and strong static controls; all failed samples disclosed.",
                "homology":"All-chain identities/long-chain homology and short exact peptides to training audited; external source does not imply unseen family.",
                "statistics":"Family/system dependence units where estimable; no residue-wise p values; confidence intervals conditional on fitted models.",
                "not_blind_new_data":"Some ThermoPCD structures/paired RMSF were inspected in older exploratory work. Newly extracted common labels/model source-out scores were not used to select this cohort.",
                "no_tuning":"No test-driven epoch, thresholds, temperature or architecture selection; input/label QC rules fixed before scores."}
    (OUT / "prespecified_protocol.json").write_text(json.dumps(protocol,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps(protocol,ensure_ascii=False,indent=2))

def shortpath(p):
    buffer = ctypes.create_unicode_buffer(32768)
    n = ctypes.windll.kernel32.GetShortPathNameW(str(p),buffer,len(buffer))
    return buffer.value if n else str(p)

def derive_matching_gro(gro,target,canonical,expected_atoms):
    """Recover atom-complete topology only when identity/order/coordinates agree.

    Two source folders have hydrogen-complete complex.gro but only a heavy-atom
    PDB. Their model input remains the frozen PDB. This conversion supplies the
    raw XTC atom order, with chain/author identities verified against that PDB.
    """
    import MDAnalysis as mda
    u=mda.Universe(str(gro))
    if len(u.atoms)!=expected_atoms:
        raise ValueError("GRO atom count differs from XTC")
    residues=[r for r in u.residues if ALIASES.get(r.resname,r.resname) in THREE_TO_ONE and "CA" in r.atoms.names]
    sequence="".join(THREE_TO_ONE[ALIASES.get(r.resname,r.resname)] for r in residues)
    if sequence!=canonical.sequence:
        raise ValueError("GRO ordered protein sequence differs from frozen PDB")
    assignments={}
    max_delta=0.
    for i,r in enumerate(residues):
        atoms={a.name:a for a in r.atoms}
        for k,name in enumerate(["N","CA","C"]):
            if canonical.backbone_mask[i,k]:
                if name not in atoms: raise ValueError("GRO lacks frozen backbone atom")
                delta=float(np.max(np.abs(atoms[name].position-canonical.backbone[i,k])))
                max_delta=max(max_delta,delta)
        assignments[r.ix]=i
    # GRO coordinates have 0.01-A rounding; no fitting/chain permutation is
    # accepted here. Matching positions resolve repeated identical chains.
    if max_delta>.03:
        raise ValueError("GRO initial backbone coordinates differ from frozen PDB")
    lines=[]
    for serial,atom in enumerate(u.atoms,1):
        name=atom.name
        residue_name=ALIASES.get(atom.resname,atom.resname)
        if residue_name=="ILE" and name=="CD":name="CD1"
        if atom.residue.ix in assignments:
            i=assignments[atom.residue.ix]
            chain=canonical.chain_ids[i]
            author=canonical.author_residue_ids[i]
            match=re.fullmatch(r"(-?\d+)([A-Za-z]?)",author)
            if not match:raise ValueError("Unsupported author residue identity in frozen PDB")
            residue_id=int(match.group(1));insertion=match.group(2) or " "
            prefix="ATOM  "
        else:
            chain="Z";residue_id=int(atom.resid)%10000;insertion=" ";prefix="HETATM"
        x,y,z=atom.position
        lines.append(f"{prefix}{serial:5d} {name:>4s} {residue_name:>3s} {chain:1s}{residue_id:4d}{insertion:1s}   {x:8.3f}{y:8.3f}{z:8.3f}{1.:6.2f}{0.:6.2f}\n")
    target.parent.mkdir(exist_ok=True)
    target.write_text("".join(lines)+"END\n",encoding="utf-8")
    return {"raw_gro":str(gro),"raw_gro_sha256":sha(gro),"atom_count":expected_atoms,
            "ordered_sequence_identity":1.,"max_initial_backbone_coordinate_delta_A":max_delta,
            "rule":"Exact protein sequence/order and initial N/CA/C coordinate agreement <=0.03 A; GRO rounding tolerance. Frozen model input unchanged."}

def process(r):
    auditfile = OUT / "audit" / (r["sample_id"]+".json")
    auditfile.parent.mkdir(exist_ok=True)
    if auditfile.exists():
        old = json.loads(auditfile.read_text(encoding="utf-8"))
        if old.get("status") == "complete": return old
        attempts=OUT/"attempts"
        attempts.mkdir(exist_ok=True)
        (attempts/(r["sample_id"]+"_"+str(time.time_ns())+".json")).write_text(
            json.dumps(old,ensure_ascii=False,indent=2),encoding="utf-8")
    started = time.monotonic()
    audit = {"sample_id":r["sample_id"]}
    try:
        canonical = parse_pdb(r["structure_path"])
        samples, traces, sources = [],[],[]
        for k,(top,xtc) in enumerate(zip(r["replica_topologies"],r["trajectory_paths"])):
            # A -m.pdb can omit hydrogens while the XTC retains them. The raw
            # topology must match atom count; this is QC, never score selection.
            from MDAnalysis.coordinates.XTC import XTCReader
            with staged_trajectory(xtc) as (staged,trajectory_hash):
                reader=XTCReader(str(staged))
                expected_atoms=reader.n_atoms
                reader.close()
                candidates=[Path(top)]+[p for p in Path(xtc).parent.glob("*.pdb") if str(p)!=top]
                matched=[]
                for candidate in candidates:
                    count=0
                    for line in candidate.read_text(encoding="utf-8",errors="replace").splitlines():
                        if line.startswith("ENDMDL"):break
                        if line.startswith(("ATOM  ","HETATM")):count+=1
                    if count==expected_atoms:matched.append(candidate)
                derivation=None
                if not matched:
                    gro=Path(xtc).parent/"complex.gro"
                    if not gro.exists():raise ValueError(f"No source topology matches XTC {expected_atoms} atoms")
                    derived=OUT/"topologies"/f"{r['sample_id']}_rep{k}_from_gro.pdb"
                    derivation=derive_matching_gro(gro,derived,canonical,expected_atoms)
                    top=str(derived)
                else:top=str(matched[0])
                tmp = OUT / "topologies" / f"{r['sample_id']}_rep{k}.pdb"
                canonical_pdb(top,tmp)
                replica, xs, stats = extract(tmp, staged, 0)
            mapping, mapped_chains = replica_map(canonical,replica)
            for x in xs:
                y = np.full((canonical.num_residues,len(METRIC_NAMES)),np.nan)
                for i,j in mapping: y[i] = x[j]
                samples.append(y)
            traces.append({**stats,"chain_mapping":mapped_chains})
            sources.append({"topology":top,"topology_sha256":sha(top),"xtc":xtc,"xtc_sha256":trajectory_hash,"topology_derivation":derivation})
        values, variances, mask = _summarise_samples(samples,canonical.num_residues)
        label = OUT / "labels" / (r["sample_id"]+".npz")
        save_label_file(str(label),values,variances,mask)
        record = dict(r)
        record["label_path"] = str(label.resolve())
        record["trajectory_paths"] = []
        record["metadata"] = {**r["metadata"],"label_status":"ready","label_sha256":sha(label)}
        for key in ["replica_topologies","original_structure_path"]: record.pop(key,None)
        audit.update(status="complete",record=record,valid_counts=mask[:,:2].sum(0).tolist(),
                     replicas=traces,source_hashes=sources,seconds=time.monotonic()-started)
    except Exception as e:
        import traceback
        audit.update(status="failed",error=repr(e),traceback=traceback.format_exc(),seconds=time.monotonic()-started)
    auditfile.write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding="utf-8")
    return audit

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode",choices=["freeze","extract"])
    ap.add_argument("--workers",type=int,default=2)
    args = ap.parse_args()
    if args.mode == "freeze": return freeze()
    protocol = json.loads((OUT/"prespecified_protocol.json").read_text(encoding="utf-8"))
    if sha(OUT/"frozen_input_manifest.jsonl") != protocol["manifest_sha256"]:
        raise RuntimeError("Frozen test membership changed")
    records = [json.loads(x) for x in (OUT/"frozen_input_manifest.jsonl").read_text(encoding="utf-8").splitlines()]
    audits = []
    with ProcessPoolExecutor(args.workers) as pool:
        for a in pool.map(process,records):
            audits.append(a)
            print(a["sample_id"],a["status"],round(a["seconds"],1),flush=True)
    (OUT/"ready_test.jsonl").write_text("\n".join(json.dumps(a["record"],ensure_ascii=False) for a in audits if a["status"]=="complete")+"\n",encoding="utf-8")
    (OUT/"extraction_summary.json").write_text(json.dumps({"prespecified":len(records),
        "completed":sum(a["status"]=="complete" for a in audits),"failed":[a for a in audits if a["status"]!="complete"],
        "code_sha256":sha(__file__)},ensure_ascii=False,indent=2),encoding="utf-8")

if __name__ == "__main__": main()
