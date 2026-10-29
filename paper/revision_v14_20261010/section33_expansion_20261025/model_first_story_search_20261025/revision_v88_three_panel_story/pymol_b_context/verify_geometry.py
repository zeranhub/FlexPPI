from pathlib import Path
import os,json
import numpy as np
ENV=Path('D:/Anaconda/envs/pvb_cpu');os.environ['PATH']=str(ENV/'Library/bin')+os.pathsep+os.environ.get('PATH','')
DLL_HANDLES=[os.add_dll_directory(str(p)) for p in [ENV/'Library/bin',ENV/'DLLs'] if p.exists()]
import pymol2
OUT=Path(__file__).resolve().parent
m=json.loads((OUT/'render_manifest.json').read_text());report={}
with pymol2.PyMOL() as pm:
 c=pm.cmd
 for name,source in [('cam_bound','1MXE_CaM_native4_147_CaMKI294_318.pdb'),('autoinhibited','1A06_rigidly_aligned_to_1MXE.pdb')]:
  c.reinitialize();c.load(str(OUT/(name+'.pse')));actual=c.get_coords(name).copy();view=np.asarray(c.get_view());c.load(str(OUT/'inputs'/source),'source');expected=c.get_coords('source')
  chain='B' if name=='cam_bound' else 'A';reg=f'{name} and chain {chain} and resi 294-318'
  entry={'PSE_atoms_match_frozen_input_exactly':bool(np.array_equal(actual,expected)),'max_coordinate_displacement_A':float(np.abs(actual-expected).max()),'PSE_count':len(actual),'same_overview_rotation':bool(np.array_equal(view[:9],np.asarray(m['camera_rotation']))),'has_M316':c.count_atoms(f'{reg} and resi 316')>0,'residues_317_318_atom_count':c.count_atoms(f'{reg} and resi 317-318')}
  if name=='cam_bound':
   anchor=c.get_coords(f'{name} and chain B and resi 316 and not name N+CA+C+O');contacts=[]
   for n in [32,51,52,55,63,71]:
    a=c.get_coords(f'{name} and chain A and resi {n}');d=float(np.linalg.norm(anchor[:,None,:]-a[None,:,:],axis=2).min());contacts.append({'CaM_native_position':n,'min_M316_sidechain_to_residue_heavy_atom_distance_A':d})
   entry['six_contacts']=contacts;entry['all_six_within_4p5A']=all(x['min_M316_sidechain_to_residue_heavy_atom_distance_A']<=4.5 for x in contacts)
  report[name]=entry
 report['both_PSE_rotations_equal']=m['same_camera_rotation'];report['both_PSE_optical_scale_equal']=m['same_camera_scale'];report['matched_scale_crops_same_dimensions']=m['assets']['cam_bound']['matched_scale_trim']['width']==m['assets']['autoinhibited']['matched_scale_trim']['width'] and m['assets']['cam_bound']['matched_scale_trim']['height']==m['assets']['autoinhibited']['matched_scale_trim']['height']
(OUT/'independent_geometry_audit.json').write_text(json.dumps(report,indent=2),encoding='utf-8');print(json.dumps(report,indent=2))