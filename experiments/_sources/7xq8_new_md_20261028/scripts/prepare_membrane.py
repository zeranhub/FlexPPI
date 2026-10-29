"""Prepare the resolved six-chain 7XQ8 BCR only; no antigen or missing tails.

Sources: 7XQ8 SSBOND, VRC01 3NGB antibody bonds, UniProt P01871 C14-LC bond.
This is an explicitly unglycosylated POPC reference model, not a paper replica.
"""
import argparse
import collections
import hashlib
import io
import json
import os
from pathlib import Path
import time
import inspect
import textwrap

import numpy as np
import openmm as mm
from openmm import app, unit
from pdbfixer import PDBFixer


def dump(path, obj):
    path.write_text(json.dumps(obj, indent=2), encoding='utf-8')


def key(a):
    r = a.residue
    return (r.chain.id, r.id, r.insertionCode.strip(), a.name)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--run-dir', type=Path, required=True)
    parser.add_argument('--stage', choices=['protein', 'membrane'], required=True)
    parser.add_argument('--gpu', default='1')
    parser.add_argument('--platform', choices=['CUDA','CPU'], default='CUDA')
    args = parser.parse_args()
    root = args.run_dir
    os.environ['OPENMM_CPU_THREADS'] = '8'
    cuda = mm.Platform.getPlatformByName(args.platform)
    if args.platform == 'CUDA':
        cuda.setPropertyDefaultValue('DeviceIndex', args.gpu)
        cuda.setPropertyDefaultValue('Precision', 'mixed')
    ff = app.ForceField('charmm36_2024.xml', 'charmm36_2024/water.xml')
    print(time.asctime(), 'starting', args.stage, flush=True)
    if args.stage == 'protein':
        text = (root/'7xq8_original.pdb').read_text()
        clean = '\n'.join(l for l in text.splitlines() if l.startswith(('ATOM  ', 'TER   ')))+'\nEND\n'
        original = app.PDBFile(io.StringIO(clean))
        original_atoms = {key(a): np.array(original.positions[a.index].value_in_unit(unit.nanometer)) for a in original.topology.atoms()}
        original_res = [(r.chain.id, r.id, r.insertionCode, r.name) for r in original.topology.residues()]
        assert len(original_res) == 1889
        fixer = PDBFixer(pdbfile=io.StringIO(clean), platform=cuda)
        fixer.findMissingResidues()
        fixer.missingResidues = {}  # no signal peptides/tags/unresolved tails
        fixer.findNonstandardResidues()
        assert not fixer.nonstandardResidues, fixer.nonstandardResidues
        fixer.findMissingAtoms()
        missing = [{'chain': r.chain.id, 'resid': r.id, 'icode':r.insertionCode,
                    'atoms':[a.name for a in atoms]} for r, atoms in fixer.missingAtoms.items()]
        dump(root/'missing_heavy_atoms.json', missing)
        print('rebuilding heavy atoms', len(missing), flush=True)
        fixer.addMissingAtoms(seed=20261029)
        repaired_res = [(r.chain.id, r.id, r.insertionCode, r.name) for r in fixer.topology.residues()]
        assert repaired_res == original_res, 'residue mapping changed'
        existing_shift = max(np.linalg.norm(np.array(fixer.positions[a.index].value_in_unit(unit.nanometer))-original_atoms[key(a)])
                             for a in fixer.topology.atoms() if key(a) in original_atoms)
        assert existing_shift < 1e-6, existing_shift
        # Explicitly replace all automatically detected S-S bonds.
        pairs = []
        for l in text.splitlines():
            if l.startswith('SSBOND'):
                pairs.append(((l[15], l[17:21].strip(), l[21].strip()),
                              (l[29], l[31:35].strip(), l[35].strip()), '7XQ8_SSBOND'))
        for c in ['C', 'v']:
            for i,j in [(32,102),(42,108),(151,211)]:
                pairs.append(((c,str(i),''),(c,str(j),''),'VRC01_3NGB_and_IgM_homology'))
        for c in ['L','R']:
            for i,j in [(23,88),(136,196)]:
                pairs.append(((c,str(i),''),(c,str(j),''),'VRC01_3NGB'))
        for c,l in [('C','L'),('v','R')]:
            pairs.append(((c,'137',''),(l,'216',''),'UniProt_P01871_C14_light_chain_and_sequence_alignment'))
        assert len(pairs) == 23
        top = app.Topology()
        atom_map = {}
        for c in fixer.topology.chains():
            nc = top.addChain(c.id)
            for r in c.residues():
                nr = top.addResidue(r.name, nc, r.id, r.insertionCode)
                for a in r.atoms():
                    atom_map[a] = top.addAtom(a.name, a.element, nr, a.id)
        for a,b in fixer.topology.bonds():
            if not (a.name == b.name == 'SG'):
                top.addBond(atom_map[a], atom_map[b])
        sg = {(a.residue.chain.id,a.residue.id,a.residue.insertionCode.strip()):a for a in top.atoms() if a.name=='SG'}
        for a,b,source in pairs:
            top.addBond(sg[a],sg[b])
        dump(root/'disulfide_evidence.json', pairs)
        modeller = app.Modeller(top, fixer.positions)
        print('adding hydrogen atoms with CHARMM36m', flush=True)
        modeller.addHydrogens(ff, pH=7.0, platform=cuda)
        # Restrained local relaxation lets new rotamers and terminal light-chain
        # residues form the evidenced S-S bonds without replacing cryoEM backbone.
        system = ff.createSystem(modeller.topology, nonbondedMethod=app.CutoffNonPeriodic,
                                 nonbondedCutoff=1.2*unit.nanometer, constraints=app.HBonds)
        rest = mm.CustomExternalForce('0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)')
        for n in ['k','x0','y0','z0']:
            rest.addPerParticleParameter(n)
        xyz = modeller.positions.value_in_unit(unit.nanometer)
        for a in modeller.topology.atoms():
            if a.element == app.element.hydrogen:
                continue
            r = a.residue
            k = 4000.0 if a.name in ['N','CA','C','O'] else 1000.0
            if key(a) not in original_atoms: k = 0.0
            if r.name == 'CYS': k = min(k, 300.0)
            if r.chain.id in ['L','R'] and int(r.id)>=211: k = 50.0
            rest.addParticle(a.index, [k,*xyz[a.index]])
        system.addForce(rest)
        integ = mm.VerletIntegrator(0.001*unit.picoseconds)
        ctx = mm.Context(system, integ, cuda)
        ctx.setPositions(modeller.positions)
        print('local S-S geometry minimization', flush=True)
        mm.LocalEnergyMinimizer.minimize(ctx, 10, 10000)
        state = ctx.getState(getPositions=True,getEnergy=True,getForces=True)
        xyz2 = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
        sg2 = {(a.residue.chain.id,a.residue.id,a.residue.insertionCode.strip()):a for a in modeller.topology.atoms() if a.name=='SG'}
        distances = [{'a':a,'b':b,'source':s,'distance_nm':float(np.linalg.norm(xyz2[sg2[a].index]-xyz2[sg2[b].index]))} for a,b,s in pairs]
        backbone_shift = [{'key':key(a),'shift_nm':float(np.linalg.norm(xyz2[a.index]-original_atoms[key(a)]))}
                          for a in modeller.topology.atoms() if a.name=='CA' and key(a) in original_atoms]
        audit = {'openmm': mm.__version__, 'residue_count':len(original_res), 'chains':[c.id for c in modeller.topology.chains()],
                 'all_atoms':modeller.topology.getNumAtoms(), 'original_atom_shift_before_relaxation_nm':float(existing_shift),
                 'disulfides':distances,'CA_shifts':backbone_shift,'potential_kj_mol':state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole),
                 'termini':'Uncapped charged NH3+/COO- at all six resolved-chain termini; missing tails excluded.',
                 'glycans':'All NAG removed; unglycosylated exploratory POPC reference.'}
        dump(root/'protein_preparation_audit.json', audit)
        assert np.isfinite(xyz2).all()
        assert all(0.185 < x['distance_nm'] < 0.235 for x in distances), 'bad S-S geometry'
        assert max(x['shift_nm'] for x in backbone_shift)<0.6, 'large local backbone shift'
        general = [x['shift_nm'] for x in backbone_shift if not(x['key'][0] in ['L','R'] and int(x['key'][1])>=211)]
        assert max(general)<0.3, 'large backbone shift'
        axis = json.loads((root/'7xq8_TM_axis_proposal.json').read_text())
        rot = np.array(axis['rotation_rows'])
        center = np.array(axis['center_A'])/10
        oriented = (xyz2-center)@rot.T
        with (root/'protein_prepared.pdb').open('w') as f:
            app.PDBFile.writeFile(modeller.topology,oriented*unit.nanometer,f,keepIds=True)
        dump(root/'residue_mapping.json', [{'index':r.index,'chain':r.chain.id,'author_id':r.id,'insertion_code':r.insertionCode,'name':r.name} for r in modeller.topology.residues()])
        print('PROTEIN_READY', audit['all_atoms'], 'atoms; max CA shift',max(x['shift_nm'] for x in backbone_shift),flush=True)
    else:
        pdb = app.PDBFile(str(root/'protein_prepared.pdb'))
        assert pdb.topology.getNumResidues()==1889
        assert sum(a.name==b.name=='SG' for a,b in pdb.topology.bonds())==23, 'S-S lost on PDB roundtrip'
        modeller = app.Modeller(pdb.topology,pdb.positions)
        # OpenMM's half-size InflateGRO-like intermediate produced NaNs for this
        # large BCR. Use the same tiling/solvation/ion code with full-size cavity
        # deletion and a more conservative local membrane relaxation. This is a
        # project-local method; the installed OpenMM source is never modified.
        source=textwrap.dedent(inspect.getsource(app.Modeller.addMembrane))
        edits={
            'Vec3(0.5*p[0], 0.5*p[1], p[2])':'Vec3(p[0], p[1], p[2])',
            'overlapCutoff = 0.22':'overlapCutoff = 0.30',
            'LocalEnergyMinimizer.minimize(context, 10.0, 30)':'LocalEnergyMinimizer.minimize(context, 10.0, 2000)',
            'LangevinIntegrator(10.0, 50.0, 0.001)':'LangevinIntegrator(10.0, 50.0, 0.0005)',
        }
        for old,new in edits.items():
            assert source.count(old)==1, 'Unexpected OpenMM API source: '+old
            source=source.replace(old,new)
        (root/'addMembrane_fullsize_implementation.py').write_text(source)
        namespace=dict(app.Modeller.addMembrane.__globals__)
        exec(compile(source,str(root/'addMembrane_fullsize_implementation.py'),'exec'),namespace)
        app.Modeller.addMembrane=namespace['addMembrane']
        dump(root/'membrane_insertion_method.json',{'method':'full-size cavity deletion with 0.30nm overlap cutoff',
             'reason':'Default half-size inflation produced particle coordinate NaN',
             'source_openmm_version':mm.__version__,'source_sha256':hashlib.sha256(source.encode()).hexdigest(),
             'relaxation_timestep_ps':0.0005,'initial_minimize_max_iterations':2000,
             'installed_library_modified':False})
        print('adding POPC and 0.15M NaCl/water, 1.5nm minimum padding',flush=True)
        modeller.addMembrane(ff,lipidType='POPC',membraneCenterZ=0*unit.nanometer,
                             minimumPadding=1.5*unit.nanometer,ionicStrength=0.15*unit.molar,
                             neutralize=True,platform=cuda)
        positions_nm = np.asarray(modeller.positions.value_in_unit(unit.nanometer))
        phosphorus_z = [positions_nm[a.index,2] for a in modeller.topology.atoms()
                        if a.residue.name=='POPC' and a.name=='P']
        lower = float(np.median([z for z in phosphorus_z if z<0]))
        upper = float(np.median([z for z in phosphorus_z if z>0]))
        tm=[]
        for c,lo,hi in [('C',577,594),('v',577,594),('A',144,165),('B',161,181)]:
            z=[float(positions_nm[a.index,2]) for a in modeller.topology.atoms()
               if a.residue.chain.id==c and a.name=='CA' and lo<=int(a.residue.id)<=hi]
            tm.append({'chain':c,'author_range':[lo,hi],'z_min_nm':min(z),'z_max_nm':max(z),
                       'CA_between_phosphate_planes':sum(lower<v<upper for v in z)})
        geometry={'phosphate_plane_z_nm':[lower,upper],'TM_spans':tm,
                  'orientation_source':'four-TM consensus geometric normal; not experimental OPM annotation'}
        dump(root/'membrane_geometry_audit.json',geometry)
        assert 2.5<upper-lower<5.5, 'unexpected initial membrane thickness'
        assert all(t['z_min_nm']<0<t['z_max_nm'] and t['CA_between_phosphate_planes']>=8 for t in tm), 'TM insertion geometry failed'
        with (root/'membrane_topology.pdbx').open('w') as f:
            app.PDBxFile.writeFile(modeller.topology,modeller.positions,f,keepIds=True)
        print('building final PME system',modeller.topology.getNumAtoms(),flush=True)
        system = ff.createSystem(modeller.topology, nonbondedMethod=app.PME,
                                 nonbondedCutoff=1.2*unit.nanometer,switchDistance=1.0*unit.nanometer,
                                 constraints=app.HBonds,rigidWater=True,ewaldErrorTolerance=0.0005)
        (root/'membrane_system.xml').write_text(mm.XmlSerializer.serialize(system))
        integ = mm.VerletIntegrator(0.001*unit.picoseconds)
        ctx = mm.Context(system,integ,cuda)
        ctx.setPositions(modeller.positions)
        state = ctx.getState(getPositions=True,getEnergy=True)
        assert np.isfinite(state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))
        (root/'membrane_state.xml').write_text(mm.XmlSerializer.serialize(state))
        protein_indices=[a.index for a in modeller.topology.atoms() if a.residue.chain.index<6]
        dump(root/'protein_atom_indices.json',protein_indices)
        counts=dict(collections.Counter(r.name for r in modeller.topology.residues()))
        dump(root/'membrane_preparation_audit.json',{'atoms':system.getNumParticles(),'residue_counts':counts,
             'box_nm':np.array(state.getPeriodicBoxVectors().value_in_unit(unit.nanometer)).tolist(),
             'protein_atoms':len(protein_indices),'temperature_K':310,'NaCl_M':0.15,'pressure_bar':1,
             'forcefields':['charmm36_2024.xml','charmm36_2024/water.xml'],
             'scope':'7XQ8 resolved six chains, 1889 residues, unglycosylated POPC; single 100ns trajectory requested.'})
        print('MEMBRANE_READY',counts,flush=True)


if __name__=='__main__':
    main()
