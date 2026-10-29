# -*- coding: utf-8 -*-
"""Versioned chemical-name and replica-identity corrections from real trajectories.

Original topology and XTC files are read only. ATLAS's already corrected labels
are retained. DynaRepo keeps its 1000-frame windows; ITrace keeps whole-replica
statistics, with exact sequence alignment rather than chain-letter merging.
"""
from __future__ import annotations
from pathlib import Path
import sys
import os
import json
import hashlib
import shutil
import time
import traceback
import argparse
import contextlib
import tempfile
from concurrent.futures import ProcessPoolExecutor
from collections import Counter
import numpy as np

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
HERE = Path(__file__).resolve().parent
REV = HERE
ORIG = HERE.parents[1]
sys.path.insert(0, str(REV / "reproducibility/src"))
from flexppi.constants import THREE_TO_ONE, CHI_ATOMS, METRIC_NAMES
from flexppi.data.structures import parse_pdb
from flexppi.data.labels import _dihedral_degrees, _entropy, _summarise_samples, save_label_file

ALIASES = {name: "HIS" for name in ["HSD", "HSE", "HSP", "HID", "HIE", "HIP"]}
PATHS = {
    "/ranze/repos/FlexPPI/artifacts/data": str(ORIG / "artifacts/data"),
    "/ranze/data/ATLAS_extracted": "E:/data/protein_md_datasets/ATLAS_extracted",
    "/ranze/data/dynarepo/dynarepo_non_nucleic_full": "D:/data/dynarepo_non_nucleic_full",
    "/ranze/data/ITrace": "E:/data/ITrace",
}

def resolve(path):
    if not path:
        return None
    p = str(path).replace("\\", "/")
    for old, new in PATHS.items():
        if p.startswith(old):
            return Path(new + p[len(old):])
    return Path(p)

def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda: f.read(1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()

@contextlib.contextmanager
def staged_trajectory(source):
    """Bounded SSD staging; hash original bytes during one sequential read."""
    import shutil
    root=HERE/"_trajectory_staging"
    root.mkdir(exist_ok=True)
    folder=Path(tempfile.mkdtemp(prefix="xtc_",dir=root))
    target=folder/"trajectory.xtc"
    h=hashlib.sha256()
    try:
        with Path(source).open("rb") as src,target.open("wb") as dst:
            for block in iter(lambda:src.read(4*1024*1024),b""):
                dst.write(block)
                h.update(block)
        if target.stat().st_size!=Path(source).stat().st_size:
            raise ValueError("Staged trajectory byte count differs")
        yield target,h.hexdigest()
    finally:
        # This unique temporary directory is generated inside the verified
        # workspace staging root, never from a source path or sample name.
        if not folder.resolve().is_relative_to(root.resolve()):
            raise RuntimeError("Temporary cleanup path escaped staging root")
        shutil.rmtree(folder)

def canonical_pdb(source, target):
    lines = Path(source).read_text(encoding="utf-8", errors="replace").splitlines(True)
    changes = Counter()
    ile_atoms = {}
    for line in lines:
        if line.startswith(("ATOM  ", "HETATM")) and len(line) >= 54 and line[17:20] == "ILE":
            ile_atoms.setdefault(line[21:27], set()).add(line[12:16].strip())
    if any({"CD", "CD1"}.issubset(atoms) for atoms in ile_atoms.values()):
        raise ValueError("Ambiguous ILE CD and CD1 coexist in " + str(source))
    out = []
    for line in lines:
        if line.startswith(("ATOM  ", "HETATM")) and len(line) >= 54:
            old = line[17:20]
            if old in ALIASES:
                line = line[:17] + ALIASES[old] + line[20:]
                changes[old + "->HIS_atoms"] += 1
            if old == "ILE" and line[12:16].strip() == "CD":
                line = line[:12] + " CD1" + line[16:]
                changes["ILE_CD->CD1_atoms"] += 1
        out.append(line)
    Path(target).parent.mkdir(parents=True, exist_ok=True)
    Path(target).write_text("".join(out), encoding="utf-8")
    return dict(changes)

def chain_fit(coordinates, chain_ids):
    """All CA atoms of each chain, fitted independently to its first frame."""
    values = np.full(coordinates.shape[1], np.nan, dtype=np.float64)
    chain_ids = np.asarray(chain_ids)
    for chain in dict.fromkeys(chain_ids.tolist()):
        ids = np.flatnonzero(chain_ids == chain)
        if len(ids) < 3:
            continue
        mobile = coordinates[:, ids].astype(np.float64)
        centered = mobile - mobile.mean(axis=1, keepdims=True)
        ref = mobile[0] - mobile[0].mean(axis=0)
        covariance = np.einsum("tni,nj->tij", centered, ref, optimize=True)
        u, _, vt = np.linalg.svd(covariance)
        rotation = u @ vt
        neg = np.linalg.det(rotation) < 0
        u[neg, :, -1] *= -1
        rotation[neg] = u[neg] @ vt[neg]
        aligned = np.einsum("tni,tij->tnj", centered, rotation, optimize=True)
        values[ids] = np.sqrt(np.mean(np.sum((aligned-aligned.mean(0))**2, axis=-1), axis=0))
    return values

def extract(topology, trajectory, window_frames):
    import MDAnalysis as mda
    structure = parse_pdb(str(topology))
    u = mda.Universe(str(topology), str(trajectory))
    try:
        residues = [r for r in u.residues if r.resname in THREE_TO_ONE and "CA" in r.atoms.names]
        keys = [(str(r.atoms.chainIDs[0]).strip() or "_", str(r.resid)+str(r.icode).strip(), str(r.resname)) for r in residues]
        expected = list(zip(structure.chain_ids, structure.author_residue_ids, structure.residue_names))
        if keys != expected:
            raise ValueError("MD topology/PDB residue identity mismatch")
        ca, quads, chis = [], [], {}
        for i, r in enumerate(residues):
            lookup = {a.name: a.index for a in r.atoms}
            ca.append(lookup["CA"])
            definitions = CHI_ATOMS.get(r.resname, [])
            if definitions and all(all(a in lookup for a in quad) for quad in definitions):
                chis[i] = list(range(len(quads), len(quads)+len(definitions)))
                quads += [[lookup[a] for a in quad] for quad in definitions]
        coords, states = [], []
        quad_array = np.asarray(quads, dtype=int)
        for ts in u.trajectory:
            coords.append(ts.positions[ca].copy())
            if len(quads):
                angles = _dihedral_degrees(ts.positions[quad_array])
                states.append(np.clip(np.floor((angles+180.0)/120.0).astype(np.int8), 0, 2))
        coords = np.asarray(coords)
        states = np.asarray(states)
        nframes = len(coords)
        ranges = [(0, nframes)] if window_frames <= 0 else [
            (start, start+window_frames) for start in range(0, nframes-window_frames+1, window_frames)]
        if not ranges:
            ranges = [(0, nframes)]
        samples = []
        for start, end in ranges:
            sample = np.full((structure.num_residues, len(METRIC_NAMES)), np.nan)
            sample[:, 0] = chain_fit(coords[start:end], structure.chain_ids)
            for i, columns in chis.items():
                codes = sum(states[start:end, col].astype(np.int32) * 3**p for p, col in enumerate(columns))
                sample[i, 1] = _entropy(codes)
            samples.append(sample)
        return structure, samples, {"frames": nframes, "windows": len(samples), "chi_residues": len(chis)}
    finally:
        # MDAnalysis retains cyclic references; close the XTC handle explicitly
        # before the temporary SSD copy is deleted on Windows.
        u.trajectory.close()


def replica_map(canonical, replica):
    """Reject sequence substitutions; allow only verified terminal/internal gaps."""
    from Bio.Align import PairwiseAligner
    from scipy.optimize import linear_sum_assignment
    a, b = canonical.unique_chains(), replica.unique_chains()
    if len(a) != len(b):
        raise ValueError("Replica chain counts differ")
    aligner = PairwiseAligner(mode="global", match_score=2, mismatch_score=-8,
                              open_gap_score=-3, extend_gap_score=-.3)
    options = {}
    scores = np.full((len(a), len(b)), -1e6)
    for i, ac in enumerate(a):
        ai = canonical.indices_for_chain(ac)
        for j, bc in enumerate(b):
            bi = replica.indices_for_chain(bc)
            seq_a, seq_b = canonical.chain_sequence(ac), replica.chain_sequence(bc)
            alignment = aligner.align(seq_a, seq_b)[0]
            pairs = []
            for aa, bb in zip(*alignment.aligned):
                pairs += list(zip(range(*aa), range(*bb)))
            identity = sum(seq_a[x] == seq_b[y] for x, y in pairs) / max(len(pairs), 1)
            coverage = len(pairs) / max(len(ai), len(bi))
            if identity == 1 and coverage >= .95:
                options[(i, j)] = ([(int(ai[x]), int(bi[y])) for x,y in pairs], coverage)
                scores[i,j] = coverage + (.00001 if ac == bc else 0)
    ii, jj = linear_sum_assignment(-scores)
    mapping, stats = [], []
    for i,j in zip(ii,jj):
        if (i,j) not in options:
            raise ValueError("No exact-identity replica chain mapping")
        pairs, coverage = options[(i,j)]
        mapping += pairs
        stats.append({"canonical_chain": a[i], "replica_chain": b[j], "identity": 1, "coverage": coverage})
    if len({x for x,_ in mapping}) != len(mapping) or len({y for _,y in mapping}) != len(mapping):
        raise ValueError("Nonunique residue mapping")
    return mapping, stats

def process(record):
    sample = record["sample_id"]
    outdir = HERE / "dataset"
    auditpath = outdir / "audit" / (sample+".json")
    if auditpath.exists():
        existing = json.loads(auditpath.read_text(encoding="utf-8"))
        if existing.get("status") == "complete":
            return existing
    started = time.monotonic()
    report = {"sample_id": sample, "dataset": record["dataset"], "fold": record["metadata"]["fold"]}
    try:
        original_pdb = resolve(record["structure_path"])
        original_label = resolve(record["label_path"])
        newpdb = outdir / "structures" / (sample+".pdb")
        changes = canonical_pdb(original_pdb, newpdb)
        canonical = parse_pdb(str(newpdb))
        old_structure = parse_pdb(str(original_pdb))
        if canonical.num_residues > 1024:
            raise ValueError("Corrected structure exceeds predeclared 1024-residue limit")
        label = outdir / "labels" / (sample+".npz")
        protocols, mappings, trajs = [], [], []
        if record["dataset"] == "atlas":
            # ATLAS extraction already handles ILE CD/CD1. Renaming a side-chain
            # atom changes neither its sequence nor backbone-only model input.
            if canonical.sequence != old_structure.sequence:
                raise ValueError("Unexpected ATLAS sequence correction; cannot reuse labels")
            label.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(original_label, label)
        else:
            samples = []
            if record["dataset"] == "itrace":
                pdbid = record["metadata"]["pdb_id"].lower()
                pairs = [(Path(f"E:/data/ITrace/trajectories/{pdbid}_run{k}_topology.pdb"),
                          Path(f"E:/data/ITrace/trajectories/{pdbid}_run{k}_md_processed.xtc")) for k in [1,2,3]]
                windows = 0
            elif record["dataset"] == "dynarepo":
                pairs = [(original_pdb, resolve(t)) for t in record["trajectory_paths"]]
                windows = 1000
            else:
                raise ValueError("Unexpected training source")
            for k, (top, xtc) in enumerate(pairs):
                tmp = outdir / "topologies" / f"{sample}_rep{k}.pdb"
                canonical_pdb(top, tmp)
                with staged_trajectory(xtc) as (staged,trajectory_hash):
                    replica, extracted, protocol = extract(tmp, staged, windows)
                mapping, ms = replica_map(canonical, replica)
                for x in extracted:
                    mapped = np.full((canonical.num_residues, len(METRIC_NAMES)), np.nan)
                    for target, source in mapping:
                        mapped[target] = x[source]
                    samples.append(mapped)
                protocols.append(protocol)
                mappings.append(ms)
                trajs.append({"topology": str(top), "topology_sha256": sha(top),
                              "trajectory": str(xtc), "trajectory_sha256": trajectory_hash, "bytes": xtc.stat().st_size,
                              "staging":"One sequential source read with SHA256, then SSD frame reads; unique temporary copy removed"})
            values, variances, mask = _summarise_samples(samples, canonical.num_residues)
            save_label_file(str(label), values, variances, mask)
        with np.load(label) as z:
            mask = z["mask"] & canonical.frame_mask[:,None]
            newvalues = z["values"]
        with np.load(original_label) as z:
            oldmask, oldvalues = z["mask"], z["values"]
        oldkeys = list(zip(old_structure.chain_ids, old_structure.author_residue_ids, old_structure.sequence))
        newkeys = list(zip(canonical.chain_ids, canonical.author_residue_ids, canonical.sequence))
        lookup = {key: i for i,key in enumerate(newkeys)}
        matched = [(i,lookup[key]) for i,key in enumerate(oldkeys) if key in lookup]
        differences = {}
        for m, name in enumerate(METRIC_NAMES[:2]):
            diffs = [abs(float(oldvalues[i,m])-float(newvalues[j,m])) for i,j in matched if oldmask[i,m] and mask[j,m]]
            differences[name] = {"common": len(diffs), "max_abs_delta": max(diffs,default=0),
                                 "changed_gt_1e_5": sum(d>1e-5 for d in diffs),
                                 "old_valid": int(oldmask[:,m].sum()), "new_valid": int(mask[:,m].sum())}
        newrecord = dict(record)
        newrecord.update(structure_path=str(newpdb.resolve()), label_path=str(label.resolve()), trajectory_paths=[])
        newrecord["metadata"] = {**record["metadata"], "residue_count": canonical.num_residues,
                                 "split_version": "chemical_identity_corrected_v1", "label_status": "ready",
                                 "label_protocol": "canonical_His_ILE_CD1_sequence_mapped_replicas_v1",
                                 "original_record_sha256": hashlib.sha256(json.dumps(record,sort_keys=True).encode()).hexdigest()}
        report.update(status="complete", record=newrecord, residue_count=canonical.num_residues,
                      original_residue_count=old_structure.num_residues, name_changes=changes,
                      valid_counts=mask[:,:2].sum(0).tolist(), differences=differences,
                      protocols=protocols, replica_mapping=mappings, trajectory_sources=trajs,
                      input_sources={str(original_pdb): sha(original_pdb), str(original_label): sha(original_label)},
                      corrected_hashes={"structure": sha(newpdb), "labels": sha(label)},
                      seconds=time.monotonic()-started)
    except Exception as e:
        report.update(status="failed", error=repr(e), traceback=traceback.format_exc(), seconds=time.monotonic()-started)
    auditpath.parent.mkdir(parents=True,exist_ok=True)
    auditpath.write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding="utf-8")
    return report

def records():
    rows = [json.loads(x) for x in (ORIG / "artifacts/evaluation_atlas_hchi_v1/oof/prediction_manifest.jsonl").read_text().splitlines()]
    return [x["record"] for x in rows]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--limit", type=int)
    ap.add_argument("--sources", default="atlas,dynarepo,itrace")
    args = ap.parse_args()
    chosen = [r for r in records() if r["dataset"] in args.sources.split(",")]
    if args.limit:
        chosen = chosen[:args.limit]
    started = time.monotonic()
    results = []
    with ProcessPoolExecutor(args.workers) as pool:
        for i, result in enumerate(pool.map(process, chosen),1):
            results.append(result)
            if i%25 == 0 or result["status"] != "complete" or i == len(chosen):
                print(f"{i}/{len(chosen)} {result['sample_id']} {result['status']} elapsed={time.monotonic()-started:.0f}s",flush=True)
    summary = {"selected": len(chosen), "sources": args.sources, "completed": sum(r["status"]=="complete" for r in results),
               "failed": [r for r in results if r["status"] != "complete"], "elapsed_seconds": time.monotonic()-started,
               "valid_counts": [sum(r.get("valid_counts",[0,0])[k] for r in results) for k in [0,1]],
               "code_sha256": sha(__file__), "old_files_overwritten": False,
               "version": "chemical_identity_corrected_v1", "dyna_protocol": "1000-frame disjoint windows; full frame stride1",
               "itrace_protocol": "all frames per replica; three-replica mean and sample variance; exact chain identity map"}
    (HERE / f"extraction_summary_{args.sources.replace(',','_')}.json").write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps({k:v for k,v in summary.items() if k != "failed"}),flush=True)

if __name__ == "__main__":
    main()
