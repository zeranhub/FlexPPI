"""Check newly rebuilt systems and all retained heavy coordinates before MD."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def atoms(path):
    result = {}
    for line in path.read_text().splitlines():
        if line.startswith('ATOM') and line[76:78].strip() not in ('H','D'):
            key=(line[21],line[22:27].strip(),line[12:16].strip())
            if key in result:
                raise ValueError('Duplicate atom key')
            result[key]=np.array([float(line[30:38]),float(line[38:46]),float(line[46:54])])
    return result


def audit(root):
    ready=json.loads((root/'prepared/READY.json').read_text())
    previous=json.loads((root/'historical_audit/prepared_READY.json').read_text())
    if ready['forcefield_sha256'] != previous['forcefield_sha256']:
        raise ValueError('Force-field files differ from historical OpenMM 8.4 batch')
    if ready['parameters'] != json.loads((root/'protocol.json').read_text()):
        raise ValueError('New prepared parameters differ from this batch protocol')
    results={}
    for state in ['apo','myb','mll','both']:
        source_path=root/'inputs'/f'{state}.pdb'
        source=atoms(source_path)
        folder=root/'prepared'/state
        prepared=atoms(folder/'protein.pdb')
        metadata=json.loads((folder/'metadata.json').read_text())
        if sha(source_path) != previous['inputs_sha256'][f'{state}.pdb']:
            raise ValueError('Input construct differs from frozen historical input')
        if not set(source).issubset(prepared):
            raise ValueError('A source heavy atom was lost')
        changes=np.stack([prepared[key]-value for key,value in source.items()])
        translation=changes.mean(axis=0)
        errors=np.linalg.norm(changes-translation,axis=1)
        if errors.max()>=0.002:
            raise ValueError('Retained heavy coordinates changed beyond PDB rounding and translation')
        if sha(folder/'protein.pdb')!=metadata['protein_pdb_sha256'] or sha(folder/'system.xml')!=metadata['system_sha256']:
            raise ValueError('Prepared system/protein checksum mismatch')
        if not (folder/'solvated.pdb').is_file() or abs(metadata['net_charge_e'])>=1e-4:
            raise ValueError('Missing solvated system or nonneutral system')
        histidines=[r for r in metadata['residues'] if r['chain']=='B' and r['name']=='HIS']
        if {int(r['author_resid']) for r in histidines}!={592,594,602,605,651} or not all(r['variant']=='HID' for r in histidines):
            raise ValueError('KIX histidine variant differs from prespecified HID')
        if metadata['initial_periodic_heavy_gap_lower_bound_nm']<2:
            raise ValueError('Periodic heavy-atom image gap too small')
        nonbonded=[f for f in metadata['nonbonded_forces'] if f['type'] in ('NonbondedForce','CustomNonbondedForce')]
        if len(nonbonded)!=2 or not all(f['getUseSwitchingFunction']=='True' and f['getCutoffDistance']=='1.2 nm' and f['getSwitchingDistance']=='1.0 nm' for f in nonbonded):
            raise ValueError('Nonbonded cutoff/switch protocol differs')
        results[state]=dict(passed=True,source_heavy_atoms=len(source),added_heavy_atoms=sorted(set(prepared)-set(source)),
            max_coordinate_deviation_after_translation_angstrom=float(errors.max()),translation_angstrom=translation.tolist(),
            protein_sha256=sha(folder/'protein.pdb'),system_sha256=sha(folder/'system.xml'),
            solvated_pdb_sha256=sha(folder/'solvated.pdb'),total_atoms=metadata['total_atoms'])
    report=dict(status='all_four_rebuilt_inputs_passed',batch='KIX_rebuilt_triplicate_20261010',
        rebuilt_batch=True,old_rep01_included=False,forcefield_files_match_historical=True,states=results,
        script_sha256=sha(Path(__file__)))
    (root/'prepared/input_audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,required=True)
    audit(parser.parse_args().root.resolve(strict=True))
