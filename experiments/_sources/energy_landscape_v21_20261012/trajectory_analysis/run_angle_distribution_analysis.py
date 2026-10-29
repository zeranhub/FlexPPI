"""Read-only raw-trajectory extraction and frozen angle-distribution analysis.

The existing 83-record sensitivity cohort is reused without selection changes.
Only byte-identical selected XTC frames are copied to disposable task staging;
MDAnalysis never receives a raw source trajectory path. No model is loaded.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import sys
import tempfile
import time
import traceback

for name in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[name] = "1"
os.environ["CUDA_VISIBLE_DEVICES"] = ""
HERE = Path(__file__).resolve().parent
SENS = HERE.parent.parent / "label_sensitivity"
sys.path.insert(0, str(SENS))
import run_label_sensitivity as prior
import numpy as np
import pandas as pd

sha256 = prior.sha256
write_json = prior.write_json
OFFSETS = (0, 30, 60, 90)


def utc():
    return datetime.now(timezone.utc).isoformat()


def freeze(output):
    output.mkdir(parents=True, exist_ok=True)
    source_plan = SENS / "sensitivity_plan.json"
    source = json.loads(source_plan.read_text(encoding="utf-8"))
    assert source["eligible_records"] == len(source["records"]) == 83
    assert sum(len(r["trajectories"]) for r in source["records"]) == 203
    assert (source["common_start_ps"], source["common_end_ps"], source["common_step_ps"], source["common_frames"]) == (0, 100000, 200, 501)
    assert sha256(SENS / "run_label_sensitivity.py") == source["analysis_code_sha256"]
    paths = {"prior_sensitivity_code": SENS / "run_label_sensitivity.py",
             "chemical_and_chain_mapping": prior.legacy_path,
             "constants": prior.CORRECTED / "training_payload/src/flexppi/constants.py",
             "dihedral_and_entropy": prior.CORRECTED / "training_payload/src/flexppi/data/labels.py",
             "pdb_parser": prior.CORRECTED / "training_payload/src/flexppi/data/structures.py"}
    plan = dict(protocol="frozen83_chi_distribution_v1", created_utc=utc(),
        inherited_plan_path=str(source_plan), inherited_plan_sha256=sha256(source_plan),
        analysis_code_sha256=sha256(__file__), dependency_hashes={k: dict(path=str(p), sha256=sha256(p)) for k,p in paths.items()},
        records=83, trajectories=203, source_counts={s: sum(r["dataset"] == s for r in source["records"]) for s in ("atlas", "dynarepo", "itrace")},
        cohort=[dict(sample_id=r["sample_id"], dataset=r["dataset"], fold=r["fold"],
                     replicas=[t["replica_index"] for t in r["trajectories"]],
                     trajectories=[dict(topology=t["topology"], topology_sha256=t["topology_sha256"], trajectory=t["trajectory"],
                          trajectory_sha256=t["trajectory_sha256"], selected_indices=t["selected_indices"],
                          selected_times_ps=[t["times_ps"][i] for i in t["selected_indices"]]) for t in r["trajectories"]]) for r in source["records"]],
        fixed_grid=dict(start_ps=0, end_ps=100000, step_ps=200, frames=501, interpolation=False),
        raw_access="All raw files read-only; full SHA256 checked; original selected XTC frames copied by frozen byte offsets; canonical topology and MDAnalysis offset files only in disposable task staging",
        coordinate_alignment="All CA atoms of each source chain; float64 Kabsch to first sampled frame, centered at that chain origin; chemistry/exact-identity mapping reused; float64 aligned coordinates stored in angstrom",
        angles="Original audited dihedral function, float32 degrees, padded to 4 chis with NaN; each available defined quadruplet retained; original joint label requires all defined chis",
        entropy="Joint 3**q categorical Shannon entropy using natural logs (nats); original interval [-180,-60),[-60,60),[60,180] with +180 in final bin",
        offset_degrees=list(OFFSETS), offset_definition="For positive shift s, wrap(theta-s) into [-180,180), then apply 120-degree bins. Original s=0 retains legacy +180 endpoint. All full-joint shifts are equal-width angular-boundary sensitivity, not physical basin correction. s=60 gives chi1 sectors centered -60,+60,180 degrees with the trans circular seam merged into one bin; compatibility with conventional chi1 centers is assessed separately only for non-Pro residues",
        chi1_centers=dict(original=[-120,0,120], shifted60=[-60,60,180]),
        conventional_chi1_exclusions="Pro excluded because ring-constrained chi1 near +/-30 degrees cannot be interpreted as the same conventional three rotamer sectors; each Pro residue is explicitly flagged. Pro remains in all legacy joint-label comparisons. No special chi2/chi3 partition is silently substituted",
        boundary_neighborhood="Fraction of observations within circular distance <=10 degrees of any bin boundary, both chi1 and any defined chi",
        apparent_changes="Number of changes between adjacent saved 200-ps categorical observations; no recrossing filter, no real transition-rate/energy-barrier interpretation",
        halves="Disjoint index halves [0,250) and [250,501), corresponding 0-49.8 ns and 50-100 ns; Jensen-Shannon divergence (natural log, range 0 to ln2)",
        windows="Five disjoint approximately 20-ns blocks, frame indices [0,100),[100,200),[200,300),[300,400),[400,501); report arithmetic and occupancy-weighted mean window H versus full-window H",
        repeats="All pairwise within-record independent saved replica distributions; DynaRepo has one replica, repeat JSD missing rather than zero",
        ca_decomposition="Aligned CA total mean-square fluctuation = occupancy-weighted within-state fluctuation + fluctuation of state means; use chi1 and joint categorical states; identity audited, descriptive rather than energetic",
        prediction_pairing="Frozen original corrected Full OOF mean unchanged; explicit output indices rmsf0/chi_entropy1, source labels rmsf0/chi_entropy1/sasa2; original target/mask/identities asserted; rescore identical paired valid intersection for rep1 and available-replica arithmetic label mean",
        aggregation="Per-trajectory/residue complete tables; record macro first, then source/all summaries; primary replica1 maintains common replica count; available replica means supplementary; no outlier trimming",
        source_exclusions=source["timestamp_exclusions"], no_new_exclusions="Inherited seven short trajectories remain excluded; all eligible83 retained; extraction failures explicit and completion requires83/203")
    path = output / "analysis_plan.json"
    if path.exists():
        frozen = json.loads(path.read_text(encoding="utf-8"))
        assert frozen["analysis_code_sha256"] == sha256(__file__)
        assert frozen["inherited_plan_sha256"] == sha256(source_plan)
        assert frozen["cohort"] == plan["cohort"]
        plan = frozen
    else:
        write_json(path, plan)
    return source, plan, sha256(path)


def chain_align(coordinates, chain_ids):
    aligned = np.full(coordinates.shape, np.nan, dtype=np.float64)
    ids_all = np.asarray(chain_ids)
    for chain in dict.fromkeys(ids_all.tolist()):
        ids = np.flatnonzero(ids_all == chain)
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
        aligned[:, ids] = np.einsum("tni,tij->tnj", centered, rotation, optimize=True)
    return aligned


def selected_xtc(item, destination):
    original = Path(item["trajectory"])
    before = original.stat()
    assert before.st_size == item["trajectory_bytes"]
    assert sha256(original) == item["trajectory_sha256"]
    offsets = item["frame_offsets"]
    assert len(offsets) == len(item["times_ps"]) == item["frames"]
    assert len(item["selected_indices"]) == 501
    with original.open("rb") as source, destination.open("wb") as target:
        for index in item["selected_indices"]:
            start = offsets[index]
            stop = offsets[index+1] if index+1 < len(offsets) else before.st_size
            assert 0 <= start < stop <= before.st_size
            source.seek(start)
            block = source.read(stop-start)
            assert len(block) == stop-start
            target.write(block)
    after = original.stat()
    assert (before.st_size, before.st_mtime_ns) == (after.st_size, after.st_mtime_ns)
    return dict(raw_sha256=item["trajectory_sha256"], raw_bytes=before.st_size,
                selected_xtc_sha256=sha256(destination), selected_bytes=destination.stat().st_size,
                raw_size_and_mtime_unchanged=True)


def states_for(angles, shift):
    # Keeping the legacy float32 arithmetic at offset0 is necessary for exact
    # reproduction at its stated +180 endpoint. Shifted circular bins wrap.
    if shift == 0:
        return np.clip(np.floor((angles + np.float32(180)) / np.float32(120)).astype(np.int16), 0, 2)
    wrapped = (angles.astype(np.float64) - shift + 180) % 360 - 180
    return np.clip(np.floor((wrapped + 180) / 120).astype(np.int16), 0, 2)


def joint_codes(angles, shift):
    states = states_for(angles, shift)
    return np.sum(states * (3 ** np.arange(states.shape[1])), axis=1).astype(np.int32)


def distribution(codes, size):
    return np.bincount(codes, minlength=size).astype(np.float64) / len(codes)


def entropy(p):
    selected = p > 0
    return float(-np.sum(p[selected] * np.log(p[selected])))


def jsd(p, q):
    value=entropy((p+q)/2) - (entropy(p)+entropy(q))/2
    assert -1e-12 <= value <= np.log(2)+1e-12
    return max(0.0,value)


def extract_one(row, item, output, plan_sha, reference):
    import MDAnalysis as mda
    identity = dict(sample_id=row["sample_id"], dataset=row["dataset"], fold=row["fold"], replica_index=item["replica_index"])
    destination = output / "cache" / (row["sample_id"] + "_rep" + str(item["replica_index"]) + ".npz")
    audit_path = output / "trajectory_audit" / (destination.stem + ".json")
    if audit_path.exists():
        existing = json.loads(audit_path.read_text(encoding="utf-8"))
        if existing.get("status") == "complete" and existing["plan_sha256"] == plan_sha:
            assert sha256(destination) == existing["cache_sha256"]
            return existing
    started = time.monotonic()
    staging = output / "_staging"
    staging.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="angle_", dir=staging) as folder:
        folder = Path(folder)
        assert folder.resolve().is_relative_to(staging.resolve())
        top, trajectory = folder / "canonical.pdb", folder / "selected.xtc"
        assert sha256(item["topology"]) == item["topology_sha256"]
        chemistry = prior.legacy.canonical_pdb(item["topology"], top)
        byte_audit = selected_xtc(item, trajectory)
        replica = prior.parse_pdb(str(top))
        universe = mda.Universe(str(top), str(trajectory))
        try:
            assert len(universe.trajectory) == 501
            residues = [r for r in universe.residues if r.resname in prior.THREE_TO_ONE and "CA" in r.atoms.names]
            keys = [(str(r.atoms.chainIDs[0]).strip() or "_", str(r.resid)+str(r.icode).strip(), str(r.resname)) for r in residues]
            assert keys == list(zip(replica.chain_ids, replica.author_residue_ids, replica.residue_names))
            ca, quads, quad_locations, q_counts = [], [], [], np.zeros(len(residues), dtype=np.int8)
            angle_available = np.zeros((len(residues),4), dtype=bool)
            for residue_index, residue in enumerate(residues):
                lookup = {a.name:a.index for a in residue.atoms}
                assert len(lookup) == len(residue.atoms), "Duplicate atom names"
                ca.append(lookup["CA"])
                definitions = prior.CHI_ATOMS.get(residue.resname, [])
                q_counts[residue_index] = len(definitions)
                assert len(definitions) <= 4
                for chi_index, atoms in enumerate(definitions):
                    if all(a in lookup for a in atoms):
                        quads.append([lookup[a] for a in atoms])
                        quad_locations.append((residue_index, chi_index))
                        angle_available[residue_index,chi_index] = True
            coordinates = np.empty((501,len(residues),3), dtype=np.float32)
            angles = np.full((501,len(residues),4), np.nan, dtype=np.float32)
            quad_array = np.asarray(quads, dtype=np.int64)
            bad_bonds = 0
            frame_times = []
            for frame_number, frame in enumerate(universe.trajectory):
                expected_time = item["times_ps"][item["selected_indices"][frame_number]]
                assert abs(frame.time-expected_time) <= 0.05
                coordinates[frame_number] = frame.positions[ca]
                frame_times.append(frame.time)
                if len(quads):
                    points = frame.positions[quad_array]
                    bad_bonds += int((np.linalg.norm(np.diff(points, axis=1), axis=-1) > 5).sum())
                    values = prior._dihedral_degrees(points)
                    assert np.isfinite(values).all()
                    for column, (ri, qi) in enumerate(quad_locations):
                        angles[frame_number,ri,qi] = values[column]
            aligned = chain_align(coordinates, replica.chain_ids)
            rmsf = np.sqrt(np.mean(np.sum((aligned-aligned.mean(axis=0))**2,axis=-1),axis=0))
            expected_rmsf = prior.legacy.chain_fit(coordinates, replica.chain_ids)
            np.testing.assert_allclose(rmsf, expected_rmsf, rtol=1e-12, atol=1e-12, equal_nan=True)
            mapping, mapping_stats = prior.legacy.replica_map(reference, replica)
            mapped_angles = np.full((501,reference.num_residues,4),np.nan,dtype=np.float32)
            mapped_ca = np.full((501,reference.num_residues,3),np.nan,dtype=np.float64)
            mapped_q = np.zeros(reference.num_residues,dtype=np.int8)
            mapped_available = np.zeros((reference.num_residues,4),dtype=bool)
            source_indices = np.full(reference.num_residues,-1,dtype=np.int32)
            for target_index, source_index in mapping:
                mapped_angles[:,target_index] = angles[:,source_index]
                mapped_ca[:,target_index] = aligned[:,source_index]
                mapped_q[target_index] = q_counts[source_index]
                mapped_available[target_index] = angle_available[source_index]
                source_indices[target_index] = source_index
            ca_valid = np.isfinite(mapped_ca).all(axis=(0,2)) & reference.frame_mask & prior.target_mask(reference,row["target_chains"])
            chi_valid = np.array([q>0 and mapped_available[i,:q].all() for i,q in enumerate(mapped_q)])
            chi_valid &= reference.frame_mask & prior.target_mask(reference,row["target_chains"])
            h = np.full(reference.num_residues,np.nan)
            for ri in np.flatnonzero(chi_valid):
                q = mapped_q[ri]
                h[ri] = prior._entropy(joint_codes(mapped_angles[:,ri,:q],0))
            mapped_rmsf = np.sqrt(np.mean(np.sum((mapped_ca-mapped_ca.mean(axis=0))**2,axis=-1),axis=0))
            destination.parent.mkdir(parents=True, exist_ok=True)
            np.savez_compressed(destination, chi_degrees=mapped_angles, aligned_ca_angstrom=mapped_ca,
                times_ps=np.asarray(frame_times), selected_original_indices=np.asarray(item["selected_indices"]),
                chi_count=mapped_q, chi_atom_available=mapped_available, ca_valid=ca_valid, chi_valid=chi_valid,
                source_residue_indices=source_indices, chain=np.asarray(reference.chain_ids),
                author_residue_id=np.asarray(reference.author_residue_ids), amino_acid=np.asarray(list(reference.sequence)),
                residue_name=np.asarray(reference.residue_names), uniform_rmsf=mapped_rmsf, uniform_hchi=h,
                analysis_plan_sha256=np.asarray(plan_sha), raw_trajectory_sha256=np.asarray(item["trajectory_sha256"]),
                topology_sha256=np.asarray(item["topology_sha256"]), units=np.asarray(["degrees","angstrom","ps","nats"]))
            report = identity | dict(status="complete", plan_sha256=plan_sha, cache_path=str(destination),
                cache_sha256=sha256(destination), mapping=mapping_stats, residue_count=reference.num_residues,
                mapped_residues=len(mapping), ca_valid=int(ca_valid.sum()), chi_valid=int(chi_valid.sum()),
                frames=501, first_ps=frame_times[0], last_ps=frame_times[-1], bad_chi_bonds_gt5A=bad_bonds,
                chemistry_changes=chemistry, byte_audit=byte_audit, elapsed_seconds=time.monotonic()-started)
            write_json(audit_path, report)
            return report
        finally:
            universe.trajectory.close()


def process_record(arguments):
    row, output_str, plan_sha = arguments
    output = Path(output_str)
    started = time.monotonic()
    report = dict(sample_id=row["sample_id"],dataset=row["dataset"],fold=row["fold"],plan_sha256=plan_sha)
    try:
        for name in ("structure","old_label","prediction"):
            assert sha256(row[name]) == row[name+"_sha256"]
        reference = prior.parse_pdb(row["structure"])
        audits = [extract_one(row,item,output,plan_sha,reference) for item in row["trajectories"]]
        samples=[]
        for audit in audits:
            with np.load(audit["cache_path"],allow_pickle=False) as cache:
                sample=np.full((reference.num_residues,3),np.nan)
                sample[:,0],sample[:,1]=cache["uniform_rmsf"],cache["uniform_hchi"]
                samples.append(sample)
        reproducibility=[]
        old_report=json.loads((SENS/"record_audit"/(row["sample_id"]+".json")).read_text(encoding="utf-8"))
        for variant,selected in (("uniform_rep1",samples[:1]),("uniform_available_replicas",samples)):
            values,variances,mask=prior._summarise_samples(selected,reference.num_residues)
            mask &= reference.frame_mask[:,None] & prior.target_mask(reference,row["target_chains"])[:,None]
            prior_label=old_report["labels"][variant]
            assert sha256(prior_label["path"])==prior_label["sha256"]
            with np.load(prior_label["path"],allow_pickle=False) as old:
                assert np.array_equal(mask[:,:2],old["mask"][:,:2])
                np.testing.assert_allclose(values[:,:2],old["values"][:,:2],rtol=1e-6,atol=1e-6)
                maximum=float(np.nanmax(np.abs(values[:,:2]-old["values"][:,:2])))
            reproducibility.append(dict(variant=variant, max_absolute_label_difference=maximum,mask_exact=True, prior_label_sha256=prior_label["sha256"]))
        report.update(status="complete",trajectories=audits,reproduction=reproducibility,elapsed_seconds=time.monotonic()-started)
    except Exception as error:
        report.update(status="failed",error=repr(error),traceback=traceback.format_exc(),elapsed_seconds=time.monotonic()-started)
    write_json(output/"record_audit"/(row["sample_id"]+".json"),report)
    return report


def boundary_fraction(angles, shift):
    distances = np.abs(((angles.astype(np.float64)[...,None] - (np.asarray([-180,-60,60])+shift) + 180)%360)-180)
    near=np.min(distances,axis=-1)<=10
    return float(near[:,0].mean()),float(near.any(axis=1).mean()),float(near.mean())


def variance_decomposition(ca,codes):
    center=ca.mean(axis=0)
    total=float(np.mean(np.sum((ca-center)**2,axis=1)))
    within,between=0.0,0.0
    for state in np.unique(codes):
        part=ca[codes==state]
        weight=len(part)/len(ca)
        mean=part.mean(axis=0)
        within+=weight*float(np.mean(np.sum((part-mean)**2,axis=1)))
        between+=weight*float(np.sum((mean-center)**2))
    assert abs(total-within-between) <= 1e-8*max(total,1)
    return total,within,between,between/total if total>0 else float("nan")


def analyze_record(arguments):
    row,output_str,plan_sha=arguments
    output=Path(output_str)
    started=time.monotonic()
    identity=dict(sample_id=row["sample_id"],dataset=row["dataset"],fold=row["fold"])
    try:
        report=json.loads((output/"record_audit"/(row["sample_id"]+".json")).read_text(encoding="utf-8"))
        assert report["status"]=="complete" and report["plan_sha256"]==plan_sha
        caches=[]
        for audit in report["trajectories"]:
            assert sha256(audit["cache_path"])==audit["cache_sha256"]
            with np.load(audit["cache_path"],allow_pickle=False) as z:
                caches.append({k:z[k] for k in z.files})
        with np.load(row["old_label"],allow_pickle=False) as z:
            old_target,old_mask=z["values"],z["mask"].astype(bool)
        with np.load(row["prediction"],allow_pickle=False) as z:
            mean,valid,target=z["mean"].astype(float),z["valid"].astype(bool),z["target"]
            assert mean.shape==(len(caches[0]["chi_count"]),2) and target.shape[1]==3
            np.testing.assert_array_equal(target,old_target)
            np.testing.assert_array_equal(valid,old_mask)
            for cache_key in ("chain","author_residue_id","amino_acid"):
                np.testing.assert_array_equal(z[cache_key],caches[0][cache_key])
        detailed,repairs,record_rows,prediction_rows,mask_rows=[],[],[],[],[]
        all_h={shift:[] for shift in OFFSETS}
        # Per-trajectory distribution probabilities retained for repeated JSD.
        distributions={}
        for replica_index,cache in enumerate(caches,start=1):
            angles=cache["chi_degrees"]
            q_counts=cache["chi_count"]
            eligible=cache["chi_valid"] & valid[:,1] & old_mask[:,1] & np.isfinite(mean[:,1])
            for shift in OFFSETS:
                h_array=np.full(len(q_counts),np.nan)
                replica_rows=[]
                for ri in np.flatnonzero(eligible):
                    q=int(q_counts[ri])
                    a=angles[:,ri,:q]
                    assert np.isfinite(a).all() and len(a)==501
                    joint=joint_codes(a,shift)
                    chi1=states_for(a[:,0:1],shift)[:,0]
                    p=distribution(joint,3**q)
                    p1=distribution(chi1,3)
                    half_js=jsd(distribution(joint[:250],3**q),distribution(joint[250:],3**q))
                    h=entropy(p)
                    h1=entropy(p1)
                    h_array[ri]=h
                    boundaries=boundary_fraction(a,shift)
                    windows=[joint[start:stop] for start,stop in ((0,100),(100,200),(200,300),(300,400),(400,501))]
                    window_h=np.asarray([entropy(distribution(c,3**q)) for c in windows])
                    weighted=float(np.average(window_h,weights=[len(c) for c in windows]))
                    assert weighted-h <= 1e-10
                    res_id=identity|dict(replica_index=replica_index,residue_index=int(ri),chain=cache["chain"][ri],
                        author_residue_id=cache["author_residue_id"][ri],amino_acid=cache["amino_acid"][ri],
                        residue_name=cache["residue_name"][ri],chi_count=q,shift_degrees=shift)
                    result=res_id|dict(h_joint=h,h_chi1=h1,chi1_p0=float(p1[0]),chi1_p1=float(p1[1]),chi1_p2=float(p1[2]),
                        occupied_joint_states=int((p>0).sum()),joint_max_occupancy=float(p.max()),chi1_max_occupancy=float(p1.max()),
                        apparent_joint_changes=int(np.count_nonzero(np.diff(joint))),apparent_chi1_changes=int(np.count_nonzero(np.diff(chi1))),
                        chi1_boundary_10deg_fraction=boundaries[0],any_chi_boundary_10deg_fraction=boundaries[1],
                        all_chi_boundary_10deg_fraction=boundaries[2],half_joint_jsd=half_js,
                        half_chi1_jsd=jsd(distribution(chi1[:250],3),distribution(chi1[250:],3)),
                        window_h_arithmetic_mean=float(window_h.mean()),window_h_weighted_mean=weighted,
                        full_minus_window_h_mean=h-float(window_h.mean()),full_minus_window_h_weighted=h-weighted,
                        old_oof_hchi_prediction=float(mean[ri,1]),old_hchi_label=float(old_target[ri,1]),
                        abs_oof_error_original=float(abs(mean[ri,1]-old_target[ri,1])),
                        abs_oof_error_uniform_offset=float(abs(mean[ri,1]-h)),
                        aligned_rmsf=float(cache["uniform_rmsf"][ri]))
                    conventional_eligible=cache["residue_name"][ri]!="PRO"
                    result.update(chi1_conventional_eligible=int(conventional_eligible),
                        chi1_conventional_exclusion_reason="Pro: ring-constrained chi1; conventional three-sector physical interpretation excluded" if not conventional_eligible else "",
                        conventional_chi1_p_gminus=float(p1[0]) if conventional_eligible and shift==60 else float("nan"),
                        conventional_chi1_p_gplus=float(p1[1]) if conventional_eligible and shift==60 else float("nan"),
                        conventional_chi1_p_trans=float(p1[2]) if conventional_eligible and shift==60 else float("nan"))
                    if cache["ca_valid"][ri]:
                        ca=cache["aligned_ca_angstrom"][:,ri]
                        total,within,between,fraction=variance_decomposition(ca,joint)
                        result.update(ca_total_variance_A2=total,ca_within_joint_variance_A2=within,ca_between_joint_variance_A2=between,
                                      ca_between_joint_fraction=fraction)
                        _,within1,between1,fraction1=variance_decomposition(ca,chi1)
                        result.update(ca_within_chi1_variance_A2=within1,ca_between_chi1_variance_A2=between1,ca_between_chi1_fraction=fraction1)
                    detailed.append(result)
                    replica_rows.append(result)
                    distributions[(replica_index,shift,ri)]=(p,p1)
                all_h[shift].append(h_array)
                values=pd.DataFrame(replica_rows)
                numerical=[c for c in values.select_dtypes(include="number") if c not in ("fold","replica_index","residue_index","chi_count","shift_degrees")]
                record_rows.append(identity|dict(replica_index=replica_index,shift_degrees=shift,residues=len(replica_rows),
                    **{c:float(values[c].mean()) for c in numerical}))
        for first in range(1,len(caches)+1):
            for second in range(first+1,len(caches)+1):
                for shift in OFFSETS:
                    for ri in np.flatnonzero(caches[first-1]["chi_valid"] & caches[second-1]["chi_valid"] & valid[:,1]):
                        if (first,shift,ri) not in distributions or (second,shift,ri) not in distributions:
                            continue
                        p,p1=distributions[(first,shift,ri)]
                        q,q1=distributions[(second,shift,ri)]
                        assert p.shape==q.shape
                        repairs.append(identity|dict(replica_first=first,replica_second=second,residue_index=int(ri),shift_degrees=shift,
                            chain=caches[0]["chain"][ri],author_residue_id=caches[0]["author_residue_id"][ri],
                            joint_repeat_jsd=jsd(p,q),chi1_repeat_jsd=jsd(p1,q1)))
        for variant,indices in (("uniform_rep1",[0]),("uniform_available_replicas",list(range(len(caches))))):
            for shift in OFFSETS:
                new_h=np.nanmean(np.asarray([all_h[shift][i] for i in indices]),axis=0)
                new_rmsf=np.nanmean(np.asarray([caches[i]["uniform_rmsf"] for i in indices]),axis=0)
                for mi,metric,new in ((0,"rmsf",new_rmsf),(1,"chi_entropy",new_h)):
                    new_mask=np.logical_and.reduce([caches[i]["ca_valid" if mi==0 else "chi_valid"] for i in indices])
                    common=new_mask & valid[:,mi] & old_mask[:,mi] & np.isfinite(new) & np.isfinite(mean[:,mi])
                    mask_rows.append(identity|dict(variant=variant,shift_degrees=shift,metric=metric,old_valid=int(valid[:,mi].sum()),
                        new_valid=int(new_mask.sum()),common=int(common.sum()),old_only=int((valid[:,mi]&~new_mask).sum()),new_only=int((new_mask&~valid[:,mi]).sum())))
                    before=prior.statistics(old_target[common,mi],mean[common,mi])
                    after=prior.statistics(new[common],mean[common,mi])
                    fields=identity|dict(variant=variant,shift_degrees=shift,metric=metric,count=int(common.sum()),
                        label_mae_shift=float(np.mean(abs(new[common]-old_target[common,mi]))),
                        label_bias_shift=float(np.mean(new[common]-old_target[common,mi])))
                    for key in ("mae","rmse","spearman","pearson","bias"):
                        fields.update({"old_"+key:before[key],"new_"+key:after[key],"new_minus_old_"+key:after[key]-before[key]})
                    prediction_rows.append(fields)
        tables=dict(per_residue_trajectory=detailed,repeat_jsd_per_residue=repairs,
                    per_record_trajectory=record_rows,prediction_effect_per_record=prediction_rows,valid_mask_audit=mask_rows)
        parts=output/"table_parts"/row["sample_id"]
        parts.mkdir(parents=True,exist_ok=True)
        paths={}
        for name,rows in tables.items():
            path=parts/(name+".csv.gz")
            pd.DataFrame(rows).to_csv(path,index=False,compression="gzip")
            paths[name]=dict(path=str(path),sha256=sha256(path),rows=len(rows))
        result=identity|dict(status="complete",plan_sha256=plan_sha,tables=paths,elapsed_seconds=time.monotonic()-started)
    except Exception as error:
        result=identity|dict(status="failed",plan_sha256=plan_sha,error=repr(error),traceback=traceback.format_exc(),elapsed_seconds=time.monotonic()-started)
    write_json(output/"distribution_audit"/(row["sample_id"]+".json"),result)
    return result


def run_pool(function,source,output,plan_sha,workers,phase):
    jobs=[(row,str(output),plan_sha) for row in source["records"]]
    results=[]
    with ProcessPoolExecutor(max_workers=workers) as pool:
        future_map={pool.submit(function,args):args[0]["sample_id"] for args in jobs}
        for index,future in enumerate(as_completed(future_map),start=1):
            result=future.result()
            results.append(result)
            print(phase,index,len(jobs),result["sample_id"],result["status"],round(result["elapsed_seconds"],2),flush=True)
            if result["status"]!="complete":
                print(result.get("traceback",result),flush=True)
    write_json(output/(phase+"_audit.json"),dict(status="complete" if all(r["status"]=="complete" for r in results) else "partial",
        plan_sha256=plan_sha,records=len(results),complete=sum(r["status"]=="complete" for r in results),
        failures=[r for r in results if r["status"]!="complete"],finished_utc=utc()))
    assert len(results)==83 and all(r["status"]=="complete" for r in results), phase+" incomplete; no exclusions or replacements allowed"


def summarize(source,output,plan_sha):
    reports=[json.loads((output/"distribution_audit"/(r["sample_id"]+".json")).read_text(encoding="utf-8")) for r in source["records"]]
    assert all(r["status"]=="complete" and r["plan_sha256"]==plan_sha for r in reports)
    hashes={}
    for name in ("per_residue_trajectory","repeat_jsd_per_residue","per_record_trajectory","prediction_effect_per_record","valid_mask_audit"):
        frames=[]
        for report in reports:
            table=report["tables"][name]
            assert sha256(table["path"])==table["sha256"]
            if table["rows"]:
                frames.append(pd.read_csv(table["path"]))
        frame=pd.concat(frames,ignore_index=True)
        destination=output/(name+(".csv.gz" if name in ("per_residue_trajectory","repeat_jsd_per_residue") else ".csv"))
        frame.to_csv(destination,index=False,compression="gzip" if destination.suffix==".gz" else None)
        hashes[name]=dict(path=str(destination),sha256=sha256(destination),rows=len(frame))
        if name=="per_record_trajectory":
            primary=frame[frame.replica_index==1].copy().assign(replica_mode="replica1")
            numerical=[c for c in frame.select_dtypes(include="number") if c not in ("fold","replica_index","shift_degrees")]
            available=frame.groupby(["sample_id","dataset","fold","shift_degrees"],as_index=False)[numerical].mean().assign(replica_mode="available_replica_mean")
            record_macro=pd.concat([primary.drop(columns="replica_index"),available],ignore_index=True)
            record_macro.to_csv(output/"distribution_record_macro.csv",index=False)
            all_frame=record_macro.copy().assign(dataset="all")
            by_source=pd.concat([record_macro,all_frame],ignore_index=True).groupby(["replica_mode","dataset","shift_degrees"])[numerical].agg(["mean","std","count"])
            by_source.columns=[a+"_"+b for a,b in by_source.columns]
            by_source.reset_index().to_csv(output/"distribution_source_summary.csv",index=False)
        elif name=="prediction_effect_per_record":
            numerical=[c for c in frame.select_dtypes(include="number") if c not in ("fold","shift_degrees")]
            by_source=pd.concat([frame,frame.copy().assign(dataset="all")],ignore_index=True).groupby(["variant","dataset","metric","shift_degrees"])[numerical].agg(["mean","std","count"])
            by_source.columns=[a+"_"+b for a,b in by_source.columns]
            by_source.reset_index().to_csv(output/"prediction_effect_source_summary.csv",index=False)
        elif name=="repeat_jsd_per_residue":
            repeat_macro=frame.groupby(["sample_id","dataset","fold","shift_degrees"])[["joint_repeat_jsd","chi1_repeat_jsd"]].mean().reset_index()
            repeat_macro.to_csv(output/"repeat_jsd_record_macro.csv",index=False)
            pd.concat([repeat_macro,repeat_macro.copy().assign(dataset="all")]).groupby(["dataset","shift_degrees"])[["joint_repeat_jsd","chi1_repeat_jsd"]].agg(["mean","std","count"]).to_csv(output/"repeat_jsd_source_summary.csv")
        elif name=="valid_mask_audit":
            assert (frame.old_only==0).all() and (frame.new_only==0).all()
            frame.groupby(["variant","dataset","metric","shift_degrees"])[["old_valid","new_valid","common","old_only","new_only"]].sum().to_csv(output/"valid_mask_source_audit.csv")
    extraction=[json.loads((output/"record_audit"/(r["sample_id"]+".json")).read_text(encoding="utf-8")) for r in source["records"]]
    trajectories=[t for r in extraction for t in r["trajectories"]]
    assert len(trajectories)==203 and all(t["frames"]==501 for t in trajectories)
    assert all(t["byte_audit"]["raw_size_and_mtime_unchanged"] for t in trajectories)
    summary_files=[p for p in output.glob("*.csv") if p.name not in [Path(v["path"]).name for v in hashes.values()]]
    hashes.update({p.stem:dict(path=str(p),sha256=sha256(p)) for p in summary_files})
    write_json(output/"completion_audit.json",dict(status="complete",completed_utc=utc(),records=83,trajectories=203,
        plan_sha256=plan_sha,analysis_code_sha256=sha256(__file__),total_selected_frames=203*501,
        ca_and_chi_cache_sha256={Path(t["cache_path"]).name:t["cache_sha256"] for t in trajectories},
        original_uniform_label_reproduction_max_absolute=max(x["max_absolute_label_difference"] for r in extraction for x in r["reproduction"]),
        original_uniform_masks_exact=True,all_raw_hashes_checked=True,all_raw_mtime_and_size_unchanged=True,
        all_source_prediction_and_label_hashes_checked=True,bad_chi_bonds_gt5A=sum(t["bad_chi_bonds_gt5A"] for t in trajectories),output_tables=hashes,
        interpretation="Descriptive finite-trajectory occupancy, angular-boundary sensitivity, apparent adjacent-frame changes, and OOF rescore only. No estimated true kinetic transition rates, free-energy barriers or Boltzmann inversion."))


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",type=Path,default=HERE)
    parser.add_argument("--phase",choices=("freeze","smoke","extract","analyze","summarize","all"),default="all")
    parser.add_argument("--workers",type=int,default=4)
    args=parser.parse_args()
    output=args.output.resolve()
    source,plan,plan_sha=freeze(output)
    print(json.dumps(dict(plan_sha256=plan_sha,code_sha256=sha256(__file__),records=83,trajectories=203,phase=args.phase)),flush=True)
    if args.phase=="smoke":
        results=[]
        for dataset in ("atlas","dynarepo","itrace"):
            row=next(r for r in source["records"] if r["dataset"]==dataset)
            arguments=(row,str(output),plan_sha)
            result=process_record(arguments)
            print("smoke_extract",row["sample_id"],result["status"],flush=True)
            assert result["status"]=="complete", result
            result=analyze_record(arguments)
            print("smoke_analyze",row["sample_id"],result["status"],flush=True)
            assert result["status"]=="complete", result
            results.append(result)
        write_json(output/"smoke_audit.json",dict(status="passed",plan_sha256=plan_sha,records=results,completed_utc=utc()))
    if args.phase in ("extract","all"):
        run_pool(process_record,source,output,plan_sha,args.workers,"extract")
    if args.phase in ("analyze","all"):
        run_pool(analyze_record,source,output,plan_sha,args.workers,"analyze")
    if args.phase in ("summarize","all"):
        summarize(source,output,plan_sha)


if __name__=="__main__":
    main()
