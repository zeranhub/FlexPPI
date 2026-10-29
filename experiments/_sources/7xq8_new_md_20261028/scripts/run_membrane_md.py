"""Staged OpenMM membrane equilibration and one 100 ns reference trajectory.

No thermodynamic convergence claim is made by these numerical stability gates.
Restart creates new trajectory segments; manifest records physical start times.
"""
import argparse
import csv
import datetime
import json
import math
import os
from pathlib import Path
import time

import numpy as np
import openmm as mm
from openmm import app, unit


class PhysicalTimeDCDReporter(app.DCDReporter):
    """Write physical start time even across timestep changes and restarts."""
    def report(self, simulation, state):
        if self._dcd is None:
            topology = simulation.topology
            if self._atomSubset is not None:
                topology = app.Topology()
                topology.setPeriodicBoxVectors(simulation.topology.getPeriodicBoxVectors())
                chain=topology.addChain()
                residue=topology.addResidue('',chain)
                atoms=list(simulation.topology.atoms())
                for i in self._atomSubset:
                    topology.addAtom(atoms[i].name,atoms[i].element,residue)
            dt=simulation.integrator.getStepSize()
            first_step=round(state.getTime()/dt)
            self._dcd=app.DCDFile(self._out,topology,dt,first_step,self._reportInterval,self._append)
        super().report(simulation,state)


def atomic_json(path, data):
    tmp = path.with_suffix('.tmp')
    tmp.write_text(json.dumps(data, indent=2))
    tmp.replace(path)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--run-dir', type=Path, required=True)
    parser.add_argument('--gpu', default='1')
    parser.add_argument('--production-ns', type=float, default=100)
    args = parser.parse_args()
    root = args.run_dir
    out = root/'replica_1'
    out.mkdir(exist_ok=True)
    if (out/'complete.json').exists():
        print('Already complete', flush=True)
        return
    top = app.PDBxFile(str(root/'membrane_topology.pdbx')).topology
    initial = mm.XmlSerializer.deserialize((root/'membrane_state.xml').read_text())
    system = mm.XmlSerializer.deserialize((root/'membrane_system.xml').read_text())
    protein = json.loads((root/'protein_atom_indices.json').read_text())
    protein_set = set(protein)
    bb = mm.CustomExternalForce('0.5*kbb*periodicdistance(x,y,z,x0,y0,z0)^2')
    sc = mm.CustomExternalForce('0.5*ksc*periodicdistance(x,y,z,x0,y0,z0)^2')
    lipid = mm.CustomExternalForce('0.5*klipid*periodicdistance(x,y,z,x0,y0,z0)^2')
    for force,name in [(bb,'kbb'),(sc,'ksc'),(lipid,'klipid')]:
        force.addGlobalParameter(name,0.0)
        for c in ['x0','y0','z0']: force.addPerParticleParameter(c)
    xyz = initial.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
    for a in top.atoms():
        if a.element == app.element.hydrogen: continue
        if a.index in protein_set:
            force = bb if a.name in ['N','CA','C','O','OXT'] else sc
            force.addParticle(a.index,xyz[a.index])
        elif a.residue.name=='POPC' and a.name=='P':
            lipid.addParticle(a.index,xyz[a.index])
    for f in [bb,sc,lipid]:system.addForce(f)
    barostat = mm.MonteCarloMembraneBarostat(1*unit.bar,0*unit.bar*unit.nanometer,
                    310*unit.kelvin, mm.MonteCarloMembraneBarostat.XYIsotropic,
                    mm.MonteCarloMembraneBarostat.ZFree,25)
    barostat.setRandomNumberSeed(20261101)
    barostat.setFrequency(0)
    system.addForce(barostat)
    integrator = mm.LangevinMiddleIntegrator(310*unit.kelvin,1/unit.picosecond,0.001*unit.picoseconds)
    integrator.setConstraintTolerance(1e-6)
    integrator.setRandomNumberSeed(20261030)
    sim = app.Simulation(top,system,integrator,mm.Platform.getPlatformByName('CUDA'),
                         {'DeviceIndex':args.gpu,'Precision':'mixed','UseBlockingSync':'true'})
    # stage name, duration ps, timestep ps, backbone/sidechain/lipid k, pressure
    stages = [
        ('NVT_1',125,0.001,4000,2000,1000,False),
        ('NVT_2',125,0.001,2000,1000,700,False),
        ('NVT_3',125,0.001,1000,500,400,False),
        ('NPT_1',500,0.002,500,200,200,True),
        ('NPT_2',500,0.002,200,50,40,True),
        ('NPT_3',500,0.002,50,0,0,True),
        ('unrestrained_equilibration',10000,0.002,0,0,0,True),
        ('production',args.production_ns*1000,0.002,0,0,0,True),
    ]
    progress_path = out/'progress.json'
    progress = json.loads(progress_path.read_text()) if progress_path.exists() else {'stage_index':0,'stage_step':0,'generation':0,'segments':[]}
    if progress_path.exists():
        sim.loadCheckpoint(str(out/progress['checkpoint']))
        print('Resuming',progress['stage_index'],progress['stage_step'],flush=True)
    else:
        sim.context.setState(initial)
        for n,v in [('kbb',4000),('ksc',2000),('klipid',1000)]:sim.context.setParameter(n,v)
        print('Minimizing full membrane with position restraints',flush=True)
        sim.minimizeEnergy(tolerance=10*unit.kilojoule_per_mole/unit.nanometer,maxIterations=10000)
        sim.context.setVelocitiesToTemperature(310*unit.kelvin,20261030)
        sim.saveState(str(out/'minimized_state.xml'))
    base_box = np.array(initial.getPeriodicBoxVectors().value_in_unit(unit.nanometer))
    initial_volume = float(np.linalg.det(base_box))
    constraints = system.getNumConstraints()
    dof = 3*system.getNumParticles()-constraints-3
    gas_constant = unit.MOLAR_GAS_CONSTANT_R.value_in_unit(unit.kilojoule_per_mole/unit.kelvin)
    # The checkpoint and progress commit are paired by alternating slot names.
    def save(i,done):
        progress['generation'] += 1
        slot = f"checkpoint_{progress['generation']%2}.chk"
        temp = out/(slot+'.tmp')
        sim.saveCheckpoint(str(temp))
        temp.replace(out/slot)
        progress.update(stage_index=i,stage_step=done,checkpoint=slot,
                        total_step=sim.currentStep,physical_time_ps=sim.context.getState().getTime().value_in_unit(unit.picoseconds),
                        updated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
        atomic_json(progress_path,progress)
    atoms=list(top.atoms())
    sg_bonds=[(a.index,b.index) for a,b in top.bonds() if a.name==b.name=='SG']
    assert len(sg_bonds)==23, 'Expected all 23 S-S bonds in topology'
    phos=[a.index for a in atoms if a.residue.name=='POPC' and a.name=='P']
    assert len(phos)>100
    progress['generation']+=1
    attempt=progress['generation']
    with (out/f'numeric_checks_{attempt:04d}.csv').open('w',newline='') as audit_file:
        writer=csv.writer(audit_file)
        writer.writerow(['stage','step','time_ps','potential_kj_mol','temperature_K','volume_nm3','max_SS_nm','area_nm2'])
        for i,stage in enumerate(stages):
            if i<progress['stage_index']:continue
            name,duration,dt,kbb,ksc,klipid,pressure=stage
            if name=='production':
                # Gate on the last 5 ns of unrestrained equilibration. This is
                # a basic density/area stability screen, not proof of convergence.
                points={}
                for path in sorted(out.glob('numeric_checks_*.csv')):
                    with path.open(newline='') as handle:
                        for row in csv.DictReader(handle):
                            if row['stage']=='unrestrained_equilibration':
                                points[float(row['time_ps'])]=row
                times=sorted(points)
                tail=[points[t] for t in times if times and t>times[-1]-5000]
                assert len(tail)>=80, 'Insufficient unrestrained equilibration samples'
                ts=np.array([float(r['temperature_K']) for r in tail])
                vs=np.array([float(r['volume_nm3']) for r in tail])
                areas=np.array([float(r['area_nm2']) for r in tail])
                area_drift=abs(areas[-25:].mean()/areas[:25].mean()-1)
                volume_drift=abs(vs[-25:].mean()/vs[:25].mean()-1)
                gate={'temperature_mean_K':float(ts.mean()),'volume_cv':float(vs.std()/vs.mean()),
                      'area_cv':float(areas.std()/areas.mean()),'area_drift_fraction':float(area_drift),
                      'volume_drift_fraction':float(volume_drift),'window_ps':5000}
                gate['passed']=bool(300<ts.mean()<320 and vs.std()/vs.mean()<.03 and
                       areas.std()/areas.mean()<.03 and area_drift<.05 and volume_drift<.05)
                atomic_json(out/'production_entry_gate.json',gate)
                assert gate['passed'], 'Equilibration stability gate failed; production was not started'
            done=progress['stage_step'] if i==progress['stage_index'] else 0
            target=round(duration/dt)
            if done>=target:
                save(i+1,0)
                continue
            integrator.setStepSize(dt*unit.picoseconds)
            for n,v in [('kbb',kbb),('ksc',ksc),('klipid',klipid)]:sim.context.setParameter(n,v)
            barostat.setFrequency(25 if pressure else 0)
            # Frequency belongs to the Force object; updating the Context requires
            # reinitialization. Preserve positions, velocities, box, time and RNG.
            sim.context.reinitialize(preserveState=True)
            part=f'{name}_part_{attempt:04d}'
            start_ps=sim.context.getState().getTime().value_in_unit(unit.picoseconds)
            progress['segments'].append({'stage':name,'stem':part,'start_stage_step':done,'start_time_ps':start_ps,
                    'dt_ps':dt,'protein_interval_steps':round(10/dt),'full_interval_steps':round(100/dt)})
            sim.reporters=[app.StateDataReporter(str(out/(part+'_energy.csv')),round(10/dt),step=True,time=True,
                           potentialEnergy=True,kineticEnergy=True,totalEnergy=True,temperature=True,volume=True,speed=True),
                           PhysicalTimeDCDReporter(str(out/(part+'_protein.dcd')),round(10/dt),atomSubset=protein,enforcePeriodicBox=False)]
            if name=='production':
                sim.reporters.append(PhysicalTimeDCDReporter(str(out/(part+'_all.dcd')),round(100/dt),enforcePeriodicBox=False))
            print('STAGE',name,'remaining_ps',(target-done)*dt,'at',time.asctime(),flush=True)
            save(i,done)
            start=time.monotonic()
            stage_start=done
            while done<target:
                block=min(round(50/dt),target-done)
                sim.step(block)
                done+=block
                st=sim.context.getState(getEnergy=True,getPositions=True)
                pos=st.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
                energy=st.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
                kinetic=st.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole)
                temp=2*kinetic/(dof*gas_constant)
                box=np.array(st.getPeriodicBoxVectors().value_in_unit(unit.nanometer))
                volume=float(np.linalg.det(box))
                ss=max(np.linalg.norm(pos[a]-pos[b]) for a,b in sg_bonds)
                writer.writerow([name,sim.currentStep,st.getTime().value_in_unit(unit.picoseconds),energy,temp,volume,ss,box[0,0]*box[1,1]])
                audit_file.flush()
                if not np.isfinite(pos).all() or not math.isfinite(energy):raise RuntimeError('Nonfinite numerical state')
                if not 200<temp<450:raise RuntimeError(f'Temperature outside numerical safety range: {temp}')
                if not .55<volume/initial_volume<1.6:raise RuntimeError('Large volume drift')
                if not .17<ss<.28:raise RuntimeError(f'Abnormal S-S geometry {ss}')
                save(i,done)
                ns_day=(done-stage_start)*dt/1000/(time.monotonic()-start)*86400
                print(name,f'{done*dt/1000:.3f}/{duration/1000:.3f} ns',f'{ns_day:.2f} ns/day',flush=True)
            sim.reporters=[]
            sim.saveState(str(out/(name+'_final_state.xml')))
            save(i+1,0)
    atomic_json(out/'complete.json',{'production_ns':args.production_ns,'replicas':1,'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
                 'status':'trajectory_complete_requires_scientific_analysis','progress':progress})
    print('PRODUCTION_COMPLETE',flush=True)


if __name__=='__main__':
    main()
