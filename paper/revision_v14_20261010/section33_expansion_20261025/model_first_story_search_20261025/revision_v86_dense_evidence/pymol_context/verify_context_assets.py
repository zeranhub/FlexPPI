from pathlib import Path
import os,json,hashlib,numpy as np
ENV=Path('D:/Anaconda/envs/pvb_cpu');os.environ['PATH']=str(ENV/'Library/bin')+os.pathsep+os.environ.get('PATH','')
if hasattr(os,'add_dll_directory'):DLL_HANDLES=[os.add_dll_directory(str(p)) for p in [ENV/'Library/bin',ENV/'DLLs'] if p.exists()]
import pymol2
from PIL import Image
OUT=Path(__file__).resolve().parent
m=json.loads((OUT/'render_manifest.json').read_text(encoding='utf-8'))
def native(path):
 d={}
 for s in path.read_text().splitlines():
  if s.startswith('ATOM  '):d[(s[21],s[22:26].strip(),s[12:16].strip())]={'xyz':np.array([float(s[k:k+8]) for k in [30,38,46]]),'b':float(s[60:66])}
 return d
sources={'autoinhibited':native(OUT/'inputs/1A06_observed_native.pdb'),'CaM_bound':native(OUT/'inputs/1MXE_model_input_native.pdb')}
qa={'assets':{},'visual_review':'passed','visual_scope':'Final trimmed 1A06 and 1MXE inspected as actual rendered images; common fragment and key residues visible, no synthetic residues or depth clipping. Whole 1A06 kinase kept in frame.','numbering_audit':json.loads((OUT/'native_numbering_audit.json').read_text(encoding='utf-8'))}
with pymol2.PyMOL() as pm:
 c=pm.cmd
 for obj,e in m['assets'].items():
  c.reinitialize();c.load(str(OUT/(obj+'.pse')));atoms=c.get_model(obj).atom;s=sources[obj];keys=[(a.chain,a.resi,a.name) for a in atoms];assert set(keys)==set(s)
  pxyz=np.array([a.coord for a in atoms]);sxyz=np.array([s[k]['xyz'] for k in keys]);pairdiff=float(np.max(np.abs(np.linalg.norm(pxyz[:,None,:]-pxyz[None,:,:],axis=-1)-np.linalg.norm(sxyz[:,None,:]-sxyz[None,:,:],axis=-1))))
  bdiff=max(abs(a.b-s[k]['b']) for a,k in zip(atoms,keys));camdiff=float(np.max(np.abs(np.array(c.get_view())-np.array(e['view']))));im=Image.open(OUT/e['trim']['file']).convert('RGBA')
  checks={'atom_count':len(atoms),'native_atom_keys_match':True,'max_pairwise_distance_change_from_original_A':pairdiff,'max_native_B_factor_float_rounding_difference':bdiff,'max_camera_difference':camdiff,'trim_dimensions':list(im.size),'transparency_range':list(im.getchannel('A').getextrema()),'PNG_sha256':hashlib.sha256((OUT/e['file']).read_bytes()).hexdigest(),'trim_PNG_sha256':hashlib.sha256((OUT/e['trim']['file']).read_bytes()).hexdigest(),'PSE_sha256':hashlib.sha256((OUT/(obj+'.pse')).read_bytes()).hexdigest()}
  assert pairdiff<5e-5 and bdiff<1e-4 and camdiff<1e-5 and checks['transparency_range']==[0,255]
  qa['assets'][obj]=checks
(OUT/'structure_QA.json').write_text(json.dumps(qa,indent=2),encoding='utf-8');print(json.dumps(qa['assets'],indent=2))