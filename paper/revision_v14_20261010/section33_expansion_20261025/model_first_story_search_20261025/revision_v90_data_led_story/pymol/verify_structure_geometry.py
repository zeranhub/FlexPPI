from pathlib import Path
import json,importlib.util,numpy as np
from PIL import Image
P=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('render_base',P/'render_mechanism_structures.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
r=json.loads((P/'render_manifest.json').read_text());audit={}
with m.pymol2.PyMOL() as p:
 c=p.cmd
 for name,a in r['assets'].items():
  c.reinitialize();c.load(str(P/(name+'.pse')))
  obj,ch,src=('kinase','A','1A06_rigidly_aligned_to_1MXE.pdb') if name.startswith('1A06') else ('bound','B','1MXE_CaM_native4_147_CaMKI294_318.pdb')
  actual=c.get_coords(obj);c.load(str(P/'inputs'/src),'reference');reference=c.get_coords('reference')
  extra={}
  selections=({'CaM_L116':f'{obj} and chain A and resi 116 and name CD2','CaM_E120':f'{obj} and chain A and resi 120 and name OE2','CaM_M124':f'{obj} and chain A and resi 124 and name SD'} if obj=='bound' else {'Y184':f'{obj} and resi 184 and name OH','P216':f'{obj} and resi 216 and name CG'})
  for k,s in selections.items():
   pt=m.project(c,s,a['width'],a['height']);a['landmarks'][k]=pt;a['trim_landmarks'][k]=[pt[0]-a['crop_box'][0],pt[1]-a['crop_box'][1]]
  coords_equal=bool(np.array_equal(actual,reference));rotation=np.array(c.get_view())[:9]
  alpha=np.asarray(Image.open(P/a['file']).convert('RGBA'))[:,:,3]
  audit[name]={'atom_count':len(actual),'input_atom_count':len(reference),'coordinates_exactly_match_frozen_input':coords_equal,'max_abs_coordinate_change_A':float(abs(actual-reference).max()),'rotation_matches_shared_camera':bool(np.allclose(rotation,m.VIEW[:9],atol=1e-7)),'transparent_background':bool((alpha==0).any()),'has_all_four_complete_key_residues':all(c.count_atoms(f'{obj} and chain {ch} and resi {s}')==n for s,n in [(296,9),(298,11),(303,14),(316,8)])}
  assert coords_equal and audit[name]['rotation_matches_shared_camera'] and audit[name]['has_all_four_complete_key_residues']
(P/'render_manifest.json').write_text(json.dumps(r,indent=2),encoding='utf-8')
(P/'independent_geometry_audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
print(json.dumps(audit,indent=2))