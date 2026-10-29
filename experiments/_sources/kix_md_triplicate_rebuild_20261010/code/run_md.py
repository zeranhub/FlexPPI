"""Auditable, restartable four-state KIX MD in an isolated OpenMM environment.

Run `prepare`, then `smoke`, then `queue`.  The queue owns only one GPU context
at a time.  A completed segment JSON is the transaction/commit marker; partial
files are moved aside before retrying from the preceding binary checkpoint.
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import math
import os
from pathlib import Path
import random
import shutil
import subprocess
import sys
import time

import numpy as np
import openmm as mm
from openmm import app, unit

STATES = {"apo": "B", "myb": "AB", "mll": "BC", "both": "ABC"}
PARAMS = dict(protocol="KIX_C36_NVT_v1", temperature_k=298.0, friction_per_ps=5.0,
              timestep_ps=0.002, cutoff_nm=1.2, switch_nm=1.0, salt_molar=0.1,
              edge_clearance_nm=1.0, equilibration_ns=10.0, production_ns=100.0,
              repeats=3, segment_ns=1.0, protein_frame_ps=1.0,
              restraint_equation="U=k*periodicdistance^2 (no factor 1/2)",
              equilibration_restraint_k_kj_mol_nm2=4184.0,
              minimization_restraint_k_kj_mol_nm2=16736.0,
              forcefield_files=["charmm36.xml", "charmm36/water.xml"],
              constraints="HBonds; rigidWater=True; tolerance=1e-6",
              thermostat="OpenMM LangevinIntegrator", pressure_control="none (NVT)",
              termini="charged, not capped", histidines="HID (all five KIX His)",
              minimizer="OpenMM L-BFGS, two stages, max 2000 iterations each")


def now():
    return dt.datetime.now(dt.timezone.utc).isoformat()


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    path = Path(path)
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")
    os.replace(temp, path)


def write_pdb(path, topology, positions):
    with Path(path).open("w") as stream:
        app.PDBFile.writeFile(topology, positions, stream, keepIds=True)


def atom_key(atom):
    res = atom.residue
    return (res.chain.id, res.id, res.insertionCode.strip(), atom.name)


def gpu_platform():
    platform = mm.Platform.getPlatformByName("CUDA")
    # CUDA_VISIBLE_DEVICES maps the selected physical GPU to local index 0.
    return platform, {"DeviceIndex": "0", "Precision": "mixed", "UseBlockingSync": "true"}


def integrator(seed):
    integ = mm.LangevinIntegrator(PARAMS["temperature_k"] * unit.kelvin,
                                 PARAMS["friction_per_ps"] / unit.picosecond,
                                 PARAMS["timestep_ps"] * unit.picoseconds)
    integ.setConstraintTolerance(1e-6)
    integ.setRandomNumberSeed(seed)
    return integ


def prepare(root):
    from pdbfixer import PDBFixer
    prepared = root / "prepared"
    prepared.mkdir(exist_ok=True)
    if (prepared / "READY.json").exists():
        print("Preparation already committed; refusing to replace inputs.", flush=True)
        return
    if any((prepared / state).exists() for state in STATES):
        raise RuntimeError("Uncommitted preparation exists: inspect before retry, do not overwrite.")
    source = root / "inputs" / "both.pdb"
    original = app.PDBFile(str(source))
    old = {atom_key(a): np.array(original.positions[a.index].value_in_unit(unit.nanometer))
           for a in original.topology.atoms()}
    fixer = PDBFixer(filename=str(source))
    fixer.platform = mm.Platform.getPlatformByName("CPU")
    fixer.findMissingResidues()
    if fixer.missingResidues:
        raise ValueError(f"Unexpected missing residues: {fixer.missingResidues}")
    fixer.findMissingAtoms()
    fixer.addMissingAtoms(seed=20260914)
    repaired = app.Modeller(fixer.topology, fixer.positions)
    new = {atom_key(a): np.array(repaired.positions[a.index].value_in_unit(unit.nanometer))
           for a in repaired.topology.atoms()}
    assert set(old).issubset(new), "Input atom lost during repair"
    error = max(float(np.linalg.norm(new[key] - value)) for key, value in old.items())
    assert error < 1e-6, f"Input heavy atoms moved during repair ({error} nm)"
    additions = sorted(set(new) - set(old))
    expected = {("C", "841", "", x) for x in ("CG", "OD1", "OD2")}
    expected.add(("C", "857", "", "OXT"))
    assert set(additions) == expected, additions
    write_pdb(prepared / "both_repaired_heavy.pdb", repaired.topology, repaired.positions)
    write_json(prepared / "repair_audit.json", dict(source_sha256=sha(source),
               existing_atom_max_displacement_nm=error, added_atoms=additions,
               note="Repair full ternary first, derive other states from identical heavy coordinates."))
    all_meta = {}
    for index, (state, chains) in enumerate(STATES.items()):
        folder = prepared / state
        folder.mkdir()
        model = app.Modeller(repaired.topology, repaired.positions)
        model.delete([c for c in model.topology.chains() if c.id not in chains])
        ff = app.ForceField(*PARAMS["forcefield_files"])
        variants = ["HID" if r.name == "HIS" else None for r in model.topology.residues()]
        random.seed(20260914 + index)
        applied_variants = model.addHydrogens(ff, pH=7.0, variants=variants,
                                            platform=mm.Platform.getPlatformByName("CPU"))
        # Atom repair is common; hydrogen building may relax H only. Keep both snapshots.
        write_pdb(folder / "protein.pdb", model.topology, model.positions)
        protein_atoms = list(model.topology.atoms())
        count = len(protein_atoms)
        heavy = [a.index for a in protein_atoms if a.element != app.element.hydrogen]
        xyz = np.asarray(model.positions.value_in_unit(unit.nanometer))
        side = float(np.ptp(xyz, axis=0).max() + 2 * PARAMS["edge_clearance_nm"])
        assert side > 2 * PARAMS["cutoff_nm"]
        # Explicit cube, centered bounding box: >= 1 nm on EACH side initially.
        xyz += side / 2 - (xyz.max(axis=0) + xyz.min(axis=0)) / 2
        model.positions = xyz * unit.nanometer
        model.addSolvent(ff, model="tip3p", boxSize=mm.Vec3(side, side, side) * unit.nanometer,
                         positiveIon="Na+", negativeIon="Cl-", ionicStrength=0.1 * unit.molar,
                         neutralize=True)
        system = ff.createSystem(model.topology, nonbondedMethod=app.PME,
                                 nonbondedCutoff=1.2 * unit.nanometer,
                                 switchDistance=1.0 * unit.nanometer,
                                 constraints=app.HBonds, rigidWater=True,
                                 ewaldErrorTolerance=0.0005)
        nb = next(f for f in system.getForces() if isinstance(f, mm.NonbondedForce))
        charge = sum(nb.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge)
                     for i in range(system.getNumParticles()))
        protein_charge = sum(nb.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge)
                             for i in range(count))
        assert abs(charge) < 1e-4, f"Non-neutral solvated system: {charge}"
        forces = []
        for force in system.getForces():
            info = dict(type=type(force).__name__, name=force.getName())
            for method in ["getUseSwitchingFunction", "getSwitchingDistance", "getCutoffDistance",
                           "getUseDispersionCorrection", "getUseLongRangeCorrection"]:
                if hasattr(force, method):
                    info[method] = str(getattr(force, method)())
            forces.append(info)
        restraint = mm.CustomExternalForce("restraint_k*periodicdistance(x,y,z,x0,y0,z0)^2")
        restraint.setName("Protein heavy atom positional restraint")
        restraint.addGlobalParameter("restraint_k", 0.0)
        for name in ["x0", "y0", "z0"]:
            restraint.addPerParticleParameter(name)
        positions = model.positions.value_in_unit(unit.nanometer)
        for atom_index in heavy:
            restraint.addParticle(atom_index, positions[atom_index])
        system.addForce(restraint)
        # Audit bonds: no covalent bonds between protein chains or to solvent.
        for a, b in model.topology.bonds():
            assert a.residue.chain == b.residue.chain, f"Unexpected inter-chain bond: {a} {b}"
        # A bounding-box lower bound on all periodic-image distances.
        span = np.ptp(np.asarray(positions[:count])[heavy], axis=0)
        periodic_gap_lower_bound = float(np.min(side - span))
        assert periodic_gap_lower_bound >= 2.0 - 1e-6
        write_pdb(folder / "solvated.pdb", model.topology, model.positions)
        # protein.pdb now uses the exact solvated coordinate frame / atom ordering.
        protein_top = app.Modeller(model.topology, model.positions)
        protein_top.delete(list(protein_top.topology.atoms())[count:])
        write_pdb(folder / "protein.pdb", protein_top.topology, protein_top.positions)
        (folder / "system.xml").write_text(mm.XmlSerializer.serialize(system))
        residues = [dict(index=r.index, chain=r.chain.id, author_resid=r.id, name=r.name,
                         variant=applied_variants[r.index]) for r in protein_top.topology.residues()]
        atoms = [dict(index=a.index, chain=a.residue.chain.id, author_resid=a.residue.id,
                      residue_name=a.residue.name, atom_name=a.name, element=a.element.symbol)
                 for a in protein_top.topology.atoms()]
        solvent_counts = {}
        for residue in list(model.topology.residues())[len(residues):]:
            solvent_counts[residue.name] = solvent_counts.get(residue.name, 0) + 1
        meta = dict(state=state, chains=chains, protein_atom_count=count,
                    protein_heavy_atom_count=len(heavy), total_atoms=system.getNumParticles(),
                    protein_charge_e=protein_charge, net_charge_e=charge, box_side_nm=side,
                    initial_periodic_heavy_gap_lower_bound_nm=periodic_gap_lower_bound,
                    solvent_counts=solvent_counts, residues=residues, protein_atoms=atoms,
                    nonbonded_forces=forces, system_sha256=sha(folder / "system.xml"),
                    protein_pdb_sha256=sha(folder / "protein.pdb"))
        write_json(folder / "metadata.json", meta)
        all_meta[state] = {k: v for k, v in meta.items() if k not in ["residues", "protein_atoms"]}
        print(json.dumps(all_meta[state]), flush=True)
    write_json(root / "protocol.json", PARAMS)
    data_dir = Path(app.__file__).parent / "data"
    ff_files = [data_dir / "charmm36.xml"] + sorted((data_dir / "charmm36").glob("*.xml"))
    provenance = dict(created_utc=now(), openmm_version=mm.__version__, parameters=PARAMS,
                      forcefield_sha256={str(p.relative_to(data_dir)): sha(p) for p in ff_files},
                      inputs_sha256={p.name: sha(p) for p in (root / "inputs").glob("*") if p.is_file()},
                      builder_sha256=sha(__file__), states=all_meta)
    write_json(prepared / "READY.json", provenance)


def simulation(root, state, seed):
    folder = root / "prepared" / state
    meta = json.loads((folder / "metadata.json").read_text())
    assert sha(folder / "system.xml") == meta["system_sha256"]
    assert sha(folder / "protein.pdb") == meta["protein_pdb_sha256"]
    pdb = app.PDBFile(str(folder / "solvated.pdb"))
    system = mm.XmlSerializer.deserialize((folder / "system.xml").read_text())
    platform, properties = gpu_platform()
    sim = app.Simulation(pdb.topology, system, integrator(seed), platform, properties)
    sim.context.setPositions(pdb.positions)
    return sim, meta


def thermodynamics(sim):
    state = sim.context.getState(getEnergy=True)
    potential = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    kinetic = state.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole)
    dof = 3 * sim.system.getNumParticles() - sim.system.getNumConstraints() - 3
    temp = 2 * kinetic / (dof * 0.008314462618)
    if not (math.isfinite(potential) and math.isfinite(kinetic) and 20 < temp < 1000):
        raise RuntimeError(f"Invalid thermodynamics: U={potential}, K={kinetic}, T={temp}")
    return dict(potential_kj_mol=potential, kinetic_kj_mol=kinetic, temperature_k=temp,
                time_ps=state.getTime().value_in_unit(unit.picosecond), step=sim.currentStep)


def minimize(sim, seed):
    # OpenMM L-BFGS is a documented departure from original steepest descent.
    sim.context.setParameter("restraint_k", PARAMS["minimization_restraint_k_kj_mol_nm2"])
    sim.minimizeEnergy(tolerance=10 * unit.kilojoule_per_mole / unit.nanometer, maxIterations=2000)
    sim.context.setParameter("restraint_k", PARAMS["equilibration_restraint_k_kj_mol_nm2"])
    sim.minimizeEnergy(tolerance=10 * unit.kilojoule_per_mole / unit.nanometer, maxIterations=2000)
    sim.context.setVelocitiesToTemperature(298 * unit.kelvin, seed)


def smoke(root):
    if not (root / "prepared" / "READY.json").exists():
        raise RuntimeError("Prepare all systems first")
    folder = root / "smoke"
    folder.mkdir(exist_ok=True)
    results = {}
    for index, state in enumerate(STATES):
        print(f"{now()} SMOKE {state}", flush=True)
        sim, meta = simulation(root, state, 20265000 + index)
        minimize(sim, 20265000 + index)
        sim.reporters.append(app.DCDReporter(str(folder / f"{state}_protein.dcd"), 500,
                             enforcePeriodicBox=False,
                             atomSubset=list(range(meta["protein_atom_count"]))))
        start = time.monotonic()
        samples = []
        # 20 ps restrained + 20 ps unrestrained. This is never production data.
        for block in range(20):
            if block == 10:
                sim.context.setParameter("restraint_k", 0.0)
            sim.step(1000)
            samples.append(thermodynamics(sim))
        elapsed = time.monotonic() - start
        close_reporters(sim)
        from MDAnalysis.coordinates.DCD import DCDReader
        from analyze_job import parse_topology, geometry_extrema
        protein_topology = parse_topology(root / "prepared" / state / "protein.pdb")
        with DCDReader(str(folder / f"{state}_protein.dcd")) as reader:
            assert len(reader) == 40 and reader.n_atoms == meta["protein_atom_count"]
            assert abs(reader.dt - 1.0) < 1e-5
            smoke_coordinates = np.stack([frame.positions.copy() for frame in reader])
            bond_min, bond_max = geometry_extrema(smoke_coordinates, protein_topology["bond_pairs"])
        # Exercise exact binary-checkpoint load before trusting long jobs.
        saved_step = sim.currentStep
        saved_energy = thermodynamics(sim)["potential_kj_mol"]
        saved = sim.context.createCheckpoint()
        sim.step(2)
        sim.context.loadCheckpoint(saved)
        assert sim.currentStep == saved_step
        restored_energy = thermodynamics(sim)["potential_kj_mol"]
        assert abs(restored_energy - saved_energy) < 0.05
        result = dict(state=state, samples=samples, wall_seconds=elapsed,
                      approximate_ns_day=0.040 / elapsed * 86400,
                      mean_temperature_last_20ps=float(np.mean([s["temperature_k"] for s in samples[10:]])),
                      total_atoms=meta["total_atoms"], platform=sim.context.getPlatform().getName(),
                      protein_dcd_round_trip="passed (40 frames, 1 ps, exact atom count)",
                      checkpoint_restore="passed (step and potential energy)",
                      protein_covalent_distance_range_angstrom=[bond_min, bond_max])
        if not 265 < result["mean_temperature_last_20ps"] < 330:
            raise RuntimeError(f"Smoke mean temperature outside acceptance band: {result}")
        write_pdb(folder / f"{state}_last.pdb", sim.topology,
                  sim.context.getState(getPositions=True).getPositions())
        write_json(folder / f"{state}.json", result)
        results[state] = result
        print(json.dumps({k: v for k, v in result.items() if k != "samples"}), flush=True)
        del sim
    write_json(folder / "PASSED.json", dict(finished_utc=now(), results=results,
               warning="40 ps technical smoke per state only; not evidence of equilibration/convergence."))


def checkpoint(sim, path):
    temp = Path(str(path) + ".tmp")
    sim.saveCheckpoint(str(temp))
    os.replace(temp, path)


def close_reporters(sim):
    for reporter in sim.reporters:
        out = getattr(reporter, "_out", None)
        if out is not None and hasattr(out, "close"):
            out.close()
    sim.reporters.clear()


def run_job(root, state, repeat):
    if json.loads((root / "protocol.json").read_text()) != PARAMS:
        raise RuntimeError("Protocol changed after preparation")
    seed = 20260914 + repeat * 100 + list(STATES).index(state)
    folder = root / "runs" / state / f"rep{repeat:02d}"
    folder.mkdir(parents=True, exist_ok=True)
    if (folder / "DONE.json").exists():
        return
    sim, meta = simulation(root, state, seed)
    write_json(folder / "job.json", dict(state=state, repeat=repeat, velocity_and_thermostat_seed=seed,
               system_sha256=meta["system_sha256"], code_sha256=sha(__file__), parameters=PARAMS))
    committed = []
    for phase in ("equilibration", "production"):
        phase_dir = folder / phase
        phase_dir.mkdir(exist_ok=True)
        files = sorted(phase_dir.glob("segment_*.json"))
        for expected, marker in enumerate(files):
            assert marker.stem == f"segment_{expected:04d}", "Checkpoint gap"
            record = json.loads(marker.read_text())
            ck = phase_dir / record["checkpoint"]
            assert ck.exists() and sha(ck) == record["checkpoint_sha256"]
            if phase == "production":
                assert (phase_dir / record["trajectory"]).stat().st_size == record["trajectory_bytes"]
            committed.append(ck)
    if committed:
        sim.loadCheckpoint(str(committed[-1]))
        print(f"{now()} RESUME {state} rep{repeat:02d} step={sim.currentStep}", flush=True)
    elif (folder / "minimized.chk").exists():
        sim.loadCheckpoint(str(folder / "minimized.chk"))
    else:
        minimize(sim, seed)
        checkpoint(sim, folder / "minimized.chk")
        write_json(folder / "minimized.json", thermodynamics(sim))
    start_wall = time.monotonic()
    start_step = sim.currentStep
    for phase, segment_count in [("equilibration", 10), ("production", 100)]:
        phase_dir = folder / phase
        sim.context.setParameter("restraint_k", 4184.0 if phase == "equilibration" else 0.0)
        for index in range(segment_count):
            stem = f"segment_{index:04d}"
            marker = phase_dir / f"{stem}.json"
            if marker.exists():
                continue
            # Never append uncertain trajectories: preserve any uncommitted files.
            leftovers = list(phase_dir.glob(f"{stem}.*"))
            if leftovers:
                recovery = folder / "interrupted" / str(time.time_ns())
                recovery.mkdir(parents=True)
                for p in leftovers:
                    shutil.move(str(p), str(recovery / p.name))
            steps_before = sim.currentStep
            expected_step = ((0 if phase == "equilibration" else 10) + index) * 500000
            assert steps_before == expected_step, (phase, index, steps_before, expected_step)
            dcd_path = phase_dir / f"{stem}.dcd"
            csv_path = phase_dir / f"{stem}.csv"
            if phase == "production":
                sim.reporters.append(app.DCDReporter(str(dcd_path), 500, enforcePeriodicBox=False,
                                                     atomSubset=list(range(meta["protein_atom_count"]))))
            sim.reporters.append(app.StateDataReporter(csv_path.open("w", buffering=1), 5000, step=True, time=True,
                                                       potentialEnergy=True, kineticEnergy=True,
                                                       temperature=True, volume=True, density=True,
                                                       speed=True, separator=","))
            block_start = time.monotonic()
            try:
                sim.step(500000)
                info = thermodynamics(sim)
            finally:
                close_reporters(sim)
            elapsed = time.monotonic() - block_start
            ck = phase_dir / f"{stem}.chk"
            checkpoint(sim, ck)
            info.update(phase=phase, segment=index, checkpoint=ck.name, checkpoint_sha256=sha(ck),
                        finished_utc=now(), wall_seconds=elapsed, speed_ns_day=86400 / elapsed)
            if phase == "production":
                info.update(trajectory=dcd_path.name, trajectory_bytes=dcd_path.stat().st_size,
                            expected_frames=1000, frame_interval_ps=1.0)
            # The JSON is written LAST. Its atomic replacement commits the segment.
            write_json(marker, info)
            current = dict(state=state, repeat=repeat, phase=phase, completed_segment=index + 1,
                           production_ns_completed=(index + 1 if phase == "production" else 0),
                           status="running", updated_utc=now(), latest=info,
                           elapsed_seconds_this_session=time.monotonic() - start_wall,
                           steps_this_session=sim.currentStep - start_step)
            write_json(root / "progress.json", current)
            print(json.dumps(current), flush=True)
            if (root / "STOP_AFTER_SEGMENT").exists():
                print("Graceful stop requested; latest committed checkpoint retained.", flush=True)
                raise SystemExit(75)
    sim.saveState(str(folder / "final_state.xml"))
    write_pdb(folder / "final_solvated.pdb", sim.topology,
              sim.context.getState(getPositions=True).getPositions())
    write_json(folder / "DONE.json", dict(state=state, repeat=repeat, finished_utc=now(),
               production_ns=100, equilibration_ns=10, final_step=sim.currentStep))
    del sim


def queue(root):
    import fcntl
    lock = (root / "queue.lock").open("a+")
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    if not (root / "smoke" / "PASSED.json").exists():
        raise RuntimeError("All four CUDA smoke tests must pass before queue launch")
    write_json(root / "queue.json", dict(pid=os.getpid(), started_utc=now(),
               cuda_visible_devices=os.environ.get("CUDA_VISIBLE_DEVICES"),
               order=[f"{state}/rep{r:02d}" for r in range(1, 4) for state in STATES],
               worker_count=1, code_sha256=sha(__file__)))
    for repeat in range(1, 4):
        for state in STATES:
            if (root / "STOP_AFTER_SEGMENT").exists():
                raise SystemExit(75)
            write_json(root / "progress.json", dict(status="starting_or_resuming", state=state,
                       repeat=repeat, updated_utc=now()))
            run_job(root, state, repeat)
            analysis_script = root / "code" / "analyze_job.py"
            analysis_complete = root / "runs" / state / f"rep{repeat:02d}" / "analysis" / "summary.json"
            if analysis_script.exists() and not analysis_complete.exists():
                result = subprocess.run([sys.executable, str(analysis_script),
                                         "--protein-pdb", str(root / "prepared" / state / "protein.pdb"),
                                         "--job-dir", str(root / "runs" / state / f"rep{repeat:02d}")],
                                        check=False)
                if result.returncode:
                    print(f"Analysis failed for {state}/{repeat}; MD data intact. Exit {result.returncode}", flush=True)
    write_json(root / "progress.json", dict(status="all_completed", finished_utc=now(),
               jobs=12, production_ns=1200, total_including_equilibration_ns=1320))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=["prepare", "smoke", "queue", "job"])
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--state", choices=STATES)
    parser.add_argument("--repeat", type=int)
    args = parser.parse_args()
    root = args.root.resolve()
    if args.command == "job":
        assert args.state and args.repeat in range(1, 4)
        run_job(root, args.state, args.repeat)
    else:
        globals()[args.command](root)


if __name__ == "__main__":
    main()
