from pathlib import Path
import os,json
import numpy as np
ENV=Path('D:/Anaconda/envs/pvb_cpu');os.environ['PATH']=str(ENV/'Library/bin')+os.pathsep+os.environ.get('PATH','')
DLL_HANDLES=[os.add_dll_directory(str(p)) for p in [ENV/'Library/bin',ENV/'DLLs'] if p.exists()]
import pymol2
ROOT=Path(__file__).resolve().parent;P=ROOT/'pymol'
m=json.loads((P/'render_manifest.json').read_text(encoding='utf-8'))
cameras={case:m['assets'][case+'_overview']['view'][:9] for case in ['cam','ub']}
(P/'inputs/camera_rotations.json').write_text(json.dumps(cameras,indent=2),encoding='utf-8')
report={}
with pymol2.PyMOL() as pm:
    c=pm.cmd
    for case in ['cam','ub']:
        c.reinitialize();c.load(str(P/m['inputs'][case]['pdb']),'complex');c.remove('complex and (hydro or not polymer.protein)')
        original=c.get_model('complex');orig={(a.chain,a.resi,a.resn,a.name):np.array(a.coord) for a in original.atom}
        assert len(orig)==m['inputs'][case]['protein_heavy_atoms']
        for key,entry in m['assets'].items():
            if not key.startswith(case+'_'):continue
            c.reinitialize();c.load(str(P/entry['pse']))
            atoms=c.get_model('complex').atom;keys={(a.chain,a.resi,a.resn,a.name) for a in atoms};assert keys==set(orig)
            diff=max(float(np.max(abs(np.array(a.coord)-orig[(a.chain,a.resi,a.resn,a.name)]))) for a in atoms)
            assert diff==0
            sticks=c.get_model('complex and rep sticks').atom
            if key.endswith('overview'):assert len(sticks)==0
            else:
                expected={(ch,resi,resn) for ch,resi,resn in entry['whole_residues']}
                whole={(a.chain,a.resi,a.resn,a.name) for a in atoms if (a.chain,a.resi,a.resn) in expected}
                actual={(a.chain,a.resi,a.resn,a.name) for a in sticks};assert actual==whole
                assert len(sticks)==entry['stick_atom_count']
                site='41' if 'Q41' in key else '62'
                sc=np.array([a.coord for a in atoms if a.chain=='A' and a.resi==site and a.name not in {'N','CA','C','O'}])
                expected_distance={(a.chain,a.resi,a.resn) for a in atoms if np.linalg.norm(sc-np.array(a.coord),axis=1).min()<=5.0}
                assert expected_distance==expected
            report[key]=dict(actual_original_coordinates_identical=True,maximum_coordinate_difference_A=diff,protein_heavy_atoms=len(atoms),actual_stick_atoms=len(sticks),expected_stick_atoms=entry['stick_atom_count'],all_complete_residue_atoms_included=True)
            assert c.count_atoms('complex and rep surface')==0
    v1=m['assets']['ub_Q41_local']['view'];v2=m['assets']['ub_Q62_local']['view']
    assert v1[:12]==v2[:12]
report['two_local_views_share_rotation_and_optical_scale']=True
(P/'saved_session_geometry_audit.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(report,ensure_ascii=False))
