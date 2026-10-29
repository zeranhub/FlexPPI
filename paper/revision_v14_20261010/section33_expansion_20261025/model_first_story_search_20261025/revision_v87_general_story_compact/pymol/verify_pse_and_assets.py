from pathlib import Path
import os,json,hashlib
import numpy as np
ENV=Path('D:/Anaconda/envs/pvb_cpu');os.environ['PATH']=str(ENV/'Library/bin')+os.pathsep+os.environ.get('PATH','')
if hasattr(os,'add_dll_directory'):DLL_HANDLES=[os.add_dll_directory(str(p)) for p in [ENV/'Library/bin',ENV/'DLLs'] if p.exists()]
import pymol2
from PIL import Image,ImageDraw
OUT=Path(__file__).resolve().parent
m=json.loads((OUT/'render_manifest.json').read_text(encoding='utf-8'))
source={}
for s in (OUT/'inputs/1MXE_CaM_native4_147_CaMKI294_318.pdb').read_text().splitlines():
 if s.startswith('ATOM'):source[(s[21],s[22:26].strip(),s[12:16].strip())]=np.array([float(s[k:k+8]) for k in [30,38,46]])
report={'assets':{},'visual_review':'pending','native_coordinate_source_heavy_atoms':len(source)}
with pymol2.PyMOL() as pm:
 c=pm.cmd
 for name,e in m['assets'].items():
  c.reinitialize();c.load(str(OUT/(name+'.pse')))
  obj='complex' if name=='overview' else name
  atoms=c.get_model(obj).atom;diff=max(float(np.max(np.abs(np.asarray(a.coord)-source[(a.chain,a.resi,a.name)]))) for a in atoms)
  cam=float(np.max(np.abs(np.asarray(c.get_view())-np.asarray(e['view']))))
  im=Image.open(OUT/e['file']).convert('RGBA')
  check={'atom_count_in_editable_display_object':len(atoms),'source_atom_keys_match':len(atoms)==len(source),'max_PyMOL_float_rounding_difference_A':diff,'max_saved_camera_difference':cam,'RGBA':im.mode=='RGBA','pixel_dimensions':list(im.size),'alpha_extrema':list(im.getchannel('A').getextrema()),'pse_file_bytes':(OUT/(name+'.pse')).stat().st_size,'png_sha256':hashlib.sha256((OUT/e['file']).read_bytes()).hexdigest(),'pse_sha256':hashlib.sha256((OUT/(name+'.pse')).read_bytes()).hexdigest()}
  assert check['source_atom_keys_match'] and diff<1e-4 and cam<1e-5 and check['RGBA'] and check['alpha_extrema']==[0,255]
  report['assets'][name]=check
report['rotation_identical_overview_and_locals']=m['camera_rotation_consistency'];report['optical_scale_identical_both_locals']=m['local_scale_consistency']
(OUT/'pse_and_assets_QA.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
canvas=Image.new('RGB',(1800,900),'white');draw=ImageDraw.Draw(canvas)
for i,key in enumerate(['overview','local_303','local_316']):
 e=m['assets'][key];path=e.get('trim',{}).get('file',e['file']);im=Image.open(OUT/path).convert('RGBA');im.thumbnail((580,850));x=i*600+(600-im.width)//2;y=35+(850-im.height)//2
 canvas.paste(im,(x,y),im);draw.text((i*600+10,10),key,fill='black')
canvas.save(OUT/'QA_structure_contact_sheet.png');print(json.dumps(report,indent=2))