from pathlib import Path
import os,json,hashlib
import numpy as np
from PIL import Image
OUT=Path(__file__).resolve().parent
ENV=Path('D:/Anaconda/envs/pvb_cpu')
os.environ['PATH']=str(ENV/'Library/bin')+os.pathsep+os.environ.get('PATH','')
DLL_HANDLES=[os.add_dll_directory(str(p)) for p in [ENV/'Library/bin',ENV/'DLLs'] if p.exists()]
import pymol2
m=json.loads((OUT/'render_manifest.json').read_text())
a=m['assets']['local'];t=a['trim'];im=Image.open(OUT/t['file']).convert('RGBA')
with pymol2.PyMOL() as pm:
 c=pm.cmd;c.load(str(OUT/'M316_contacts_local.pse'));saved=c.get_coords('autoinhibited').copy();rotation=list(c.get_view()[:9]);c.load(str(OUT/'inputs/1A06_rigidly_aligned_to_1MXE.pdb'),'reference');reference=c.get_coords('reference')
 qa={'PSE_atom_count':int(len(saved)),'PSE_coordinates_equal_input':bool(np.array_equal(saved,reference)),
     'PSE_camera_rotation_matches_manifest':rotation==a['view'][:9],
     'trim_size_matches_manifest':list(im.size)==[t['width'],t['height']],
     'all_landmarks_inside_trim':all(0<=p[0]<im.width and 0<=p[1]<im.height for p in t['landmarks'].values()),
     'no_rendering_canvas_edge_clipping':a['alpha_bbox'][0]>0 and a['alpha_bbox'][1]>0 and a['alpha_bbox'][2]<a['width'] and a['alpha_bbox'][3]<a['height'],
     'transparent_background':im.getchannel('A').getextrema()[0]==0,
     'PNG_sha256':hashlib.sha256((OUT/t['file']).read_bytes()).hexdigest()}
qa['passed']=all(v for k,v in qa.items() if isinstance(v,bool))
(OUT/'asset_QA.json').write_text(json.dumps(qa,indent=2),encoding='utf-8');print(json.dumps(qa,indent=2))