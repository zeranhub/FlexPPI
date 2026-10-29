"""Real-coordinate MEL5/alpha24beta17 molecular context; no structure generation."""
from pathlib import Path
import json, hashlib, os, sys
_pymol_dll_dir = Path(r'D:\Anaconda\envs\pvb_cpu\Library\bin')
_dll = os.add_dll_directory(str(_pymol_dll_dir)) if hasattr(os, 'add_dll_directory') and _pymol_dll_dir.is_dir() else None
import numpy as np
import pymol
pymol.finish_launching(['pymol', '-cq'])
from pymol import cmd

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT
sys.path.insert(0,str(ROOT.parent.parent))
from publication_style import BLUE, WARM, MD, LIGHT, TEAL
rgb=lambda color:[int(color[i:i+2],16)/255 for i in (1,3,5)]
cmd.load(str(SOURCE/'4jff_model_input.pdb'),'eng')
cmd.hide('everything','all')
cmd.bg_color('white')
cmd.set('orthoscopic',1)
cmd.set('ray_shadows',0)
cmd.set('specular',0.15)
cmd.set('antialias',2)
cmd.set('cartoon_fancy_helices',1)
cmd.set('cartoon_smooth_loops',1)
cmd.set('ray_opaque_background',0)
cmd.set_color('mhc',rgb(TEAL))
cmd.set_color('alpha',rgb(LIGHT))
cmd.set_color('beta',rgb(BLUE))
cmd.set_color('pep',rgb(WARM))
cmd.set_color('q72',rgb(TEAL))
for sel,col in [('chain A','mhc'),('chain D','alpha'),('chain E','beta'),('chain C','pep')]:
    cmd.color(col,sel)
cmd.show('cartoon','eng and ((chain A and resi 1-180) or (chain D+E and resi 1-110))')
cmd.show('sticks','eng and chain C and not hydro')
cmd.set('stick_radius',.2,'chain C')
cmd.show('spheres','eng and chain A and resi 72 and name CA')
cmd.set('sphere_scale',.35)
cmd.color('q72','chain A and resi 72')
xyz=lambda s:np.asarray(cmd.get_coords(s),float)
pep=xyz('chain C and name CA')
x=pep[-1]-pep[0]; x/=np.linalg.norm(x)
z=xyz('chain D+E and resi 1-110 and name CA').mean(0)-pep.mean(0)
z-=x*np.dot(z,x); z/=np.linalg.norm(z)
y=np.cross(z,x)
zview=.7*z-.714*y; zview/=np.linalg.norm(zview)
yview=np.cross(zview,x)
rot=np.stack([x,yview,zview])
cmd.set_view(tuple(rot.T.flatten())+(0,0,-150)+tuple(pep.mean(0)+z*10)+(100,200,0))
cmd.zoom('eng and ((chain A and resi 1-180) or (chain D+E and resi 1-110) or chain C)',buffer=1,complete=1)
overview=list(cmd.get_view())
cmd.turn('y',40)
overview=list(cmd.get_view())
cmd.png(str(ROOT/'4jff_overview.png'),1800,1500,dpi=450,ray=1)

cmd.hide('everything','all')
cmd.show('cartoon','eng and ((chain A and resi 69-75) or (chain E and resi 51-56))')
cmd.set('cartoon_transparency',.55)
cmd.show('sticks','eng and ((chain A and resi 72) or (chain E and resi 54-55)) and not hydro')
cmd.set('stick_radius',.18)
cmd.color('q72','chain A and resi 72')
cmd.color('beta','chain E')
cmd.color('oxygen','(chain A and resi 72 or chain E and resi 54-55) and elem O')
cmd.color('nitrogen','(chain A and resi 72 or chain E and resi 54-55) and elem N')
q=xyz('chain A and resi 72 and name CA')[0]
f=xyz('chain E and resi 54 and name CA')[0]
g=xyz('chain E and resi 55 and name CA')[0]
x=g-f; x/=np.linalg.norm(x)
y=(f+g)/2-q; y-=x*np.dot(y,x); y/=np.linalg.norm(y)
z=np.cross(x,y)
# Tilt slightly so the amide bonds are visible without flat projection.
yview=.91*y+.414*z; yview/=np.linalg.norm(yview)
zview=np.cross(x,yview)
rot=np.stack([x,yview,zview])
center=(q+f+g)/3
cmd.set_view(tuple(rot.T.flatten())+(0,0,-50)+tuple(center)+(20,80,0))
cmd.zoom('((chain A and resi 72) or (chain E and resi 54-55)) and not hydro',buffer=1.2,complete=1)
cmd.distance('d1','chain A and resi 72 and name OE1','chain E and resi 54 and name N')
cmd.distance('d2','chain A and resi 72 and name NE2','chain E and resi 55 and name O')
cmd.hide('labels','all')
cmd.set('dash_color','gray40');cmd.set('dash_width',2);cmd.set('dash_gap',.24);cmd.set('dash_length',.2)
local=list(cmd.get_view())
cmd.png(str(ROOT/'4jff_Q72_F53_G54.png'),2200,1500,dpi=450,ray=1)
cmd.save(str(ROOT/'MEL5_affinity_maturation.pse'))
names={'Q72_OE1':'chain A and resi 72 and name OE1','Q72_NE2':'chain A and resi 72 and name NE2','beta_F53_N':'chain E and resi 54 and name N','beta_G54_O':'chain E and resi 55 and name O'}
meta={'status':'passed','input_file':str(SOURCE/'4jff_model_input.pdb'),'input_SHA256':hashlib.sha256((SOURCE/'4jff_model_input.pdb').read_bytes()).hexdigest(),'common_numbering':{'HLA_Q72':'input A72','beta_Phe53':'input E54','beta_Gly54':'input E55'},'overview_view':overview,'local_view':local,'atom_coordinates':{name:xyz(sel)[0].tolist() for name,sel in names.items()},'distances_A':{'Q72_OE1_to_beta_Phe53_N':float(cmd.get_distance(names['Q72_OE1'],names['beta_F53_N'])),'Q72_NE2_to_beta_Gly54_O':float(cmd.get_distance(names['Q72_NE2'],names['beta_G54_O']))},'colors':{'HLA':TEAL,'TCR_alpha':LIGHT,'TCR_beta':BLUE,'peptide':WARM,'Q72':TEAL,'oxygen':'PyMOL conventional oxygen','nitrogen':'PyMOL conventional nitrogen'},'data_integrity':'Real existing MD input coordinates, no conformational editing. All-heavy structure shown for interpretation only; model uses N/CA/C.'}
(ROOT/'molecular_view.json').write_text(json.dumps(meta,indent=2),encoding='utf8')
cmd.quit()
