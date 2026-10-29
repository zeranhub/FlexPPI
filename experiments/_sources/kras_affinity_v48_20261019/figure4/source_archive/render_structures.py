"""Real KRAS-complex coordinates, MHC-superposed, with matched camera and scale."""
from pathlib import Path
import json
import numpy as np
import pymol
pymol.finish_launching(['pymol', '-cq'])
from pymol import cmd
OUT=Path(__file__).resolve().parent
SOURCE=OUT.parent.parent/'itrace_case_search_only_20261017'/'replica_recheck'
for pdb in ['7ow5','7ow6']:
    cmd.load(str(SOURCE/pdb/f'{pdb}_model_input.pdb'),'s'+pdb)
align=cmd.align('s7ow6 and chain A and resi 1-180 and name CA', 's7ow5 and chain A and resi 1-180 and name CA',cycles=0)
xyz=lambda sel:np.array(cmd.get_coords(sel),float)
pep=xyz('s7ow5 and chain C and name CA')
x=pep[-1]-pep[0];x/=np.linalg.norm(x)
z=xyz('s7ow5 and chain D+E and resi 1-110 and name CA').mean(0)-pep.mean(0)
z-=x*np.dot(z,x);z/=np.linalg.norm(z)
y=np.cross(z,x);zview=.6*z-.8*y;yview=np.cross(zview,x)
rot=np.stack([x,yview,zview])
cmd.hide('everything','all');cmd.bg_color('white')
cmd.set('orthoscopic',1);cmd.set('ray_shadows',0);cmd.set('specular',0.18)
cmd.set('antialias',2);cmd.set('cartoon_fancy_helices',1);cmd.set('cartoon_smooth_loops',1)
cmd.set_color('mhcgray',[.76,.77,.79]);cmd.set_color('tcrgray',[.58,.61,.64])
cmd.set_color('WT',[.8118,.5490,.1529]);cmd.set_color('G12D',[.4627,.1647,.5137])
for pdb,col in [('7ow5','WT'),('7ow6','G12D')]:
    obj='s'+pdb
    cmd.show('cartoon',f'{obj} and chain A and resi 1-180')
    cmd.color('mhcgray',f'{obj} and chain A')
    cmd.show('cartoon',f'{obj} and chain D+E and resi 1-110')
    cmd.color('tcrgray',f'{obj} and chain D+E')
    cmd.show('sticks',f'{obj} and chain C and not hydro')
    cmd.set('stick_radius',.19,f'{obj} and chain C')
    cmd.color(col,f'{obj} and chain C')
cmd.set_view(tuple(rot.T.flatten())+(0,0,-150)+tuple(pep.mean(0)+z*10)+(100,200,0))
cmd.zoom('(s7ow5 or s7ow6) and ((chain A and resi 1-180) or (chain D+E and resi 1-110) or chain C)',buffer=2,complete=1)
view=list(cmd.get_view());cmd.viewport(1500,1100)
for pdb in ['7ow5','7ow6']:
    cmd.disable('all');cmd.enable('s'+pdb);cmd.set_view(view)
    cmd.png(str(OUT/f'{pdb}_binding.png'),1500,1100,dpi=450,ray=1)
cmd.enable('all');cmd.hide('everything','all')
cmd.show('sticks','chain C and not hydro')
cmd.set('stick_radius',.13,'s7ow5 and chain C')
cmd.set('stick_radius',.17,'s7ow6 and chain C')
cmd.show('spheres','chain C and resi 6 and name CA')
cmd.set('sphere_scale',.22,'chain C and resi 6 and name CA')
cmd.set_view(tuple(rot.T.flatten())+(0,0,-150)+tuple(pep.mean(0))+(100,200,0))
cmd.zoom('chain C',buffer=2,complete=1)
peptideview=list(cmd.get_view())
cmd.png(str(OUT/'peptide_superposition.png'),1700,900,dpi=450,ray=1)
label_xyz={}
for name,sel in [('G6','s7ow5 and chain C and resi 6 and name CA'),('D6','s7ow6 and chain C and resi 6 and name CA')]:
    label_xyz[name]=xyz(sel).tolist()
cmd.save(str(OUT/'KRAS_WT_G12D_aligned.pse'))
audit={'alignment':'7OW6 to 7OW5: chain A C-alpha residues 1-180, no outlier rejection',
       'MHC_fit_RMSD_A':float(align[0]),'MHC_fit_n':int(align[1]),'shared_complex_camera':view,
       'peptide_camera':peptideview,'source_PDBs':[str(SOURCE/p/f'{p}_model_input.pdb') for p in ['7ow5','7ow6']],
       'state_colors':{'WT':'#CF8C27','G12D':'#762A83'},'CA6_coordinates':label_xyz,
       'transform_policy':'rigid MHC alignment only; no deformation, rescaling, smoothing or deletion of peptide atoms'}
(OUT/'molecular_view.json').write_text(json.dumps(audit,indent=2),encoding='utf8')
cmd.quit()
