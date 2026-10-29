"""Matched-camera molecular views, using real input coordinates in PyMOL."""
from pathlib import Path
import json, numpy as np
import pymol
pymol.finish_launching(['pymol','-cq'])
from pymol import cmd
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'figure6';OUT.mkdir(exist_ok=True)
for pdb in ['2ak4','4jry']:cmd.load(str(ROOT/'sb27_sb47'/f'{pdb}_model_input.pdb'),'s'+pdb)
cmd.align('s4jry and chain A and resi 1-180 and name CA','s2ak4 and chain A and resi 1-180 and name CA',cycles=0)
xyz=lambda sel:np.array(cmd.get_coords(sel),float)
pep=xyz('s2ak4 and chain C and name CA')
x=pep[-1]-pep[0];x/=np.linalg.norm(x)
z=xyz('s2ak4 and chain D+E and resi 1-110 and name CA').mean(0)-pep.mean(0);z-=x*np.dot(z,x);z/=np.linalg.norm(z)
y=np.cross(z,x);zview=.6*z-.8*y;yview=np.cross(zview,x)
rot=np.stack([x,yview,zview])
cmd.hide('everything','all')
cmd.bg_color('white')
cmd.set('orthoscopic',1);cmd.set('ray_shadows',0);cmd.set('specular',0.2)
cmd.set('antialias',2);cmd.set('cartoon_fancy_helices',1);cmd.set('cartoon_smooth_loops',1)
cmd.set_color('mhcgray',[.75,.77,.79]);cmd.set_color('tcra',[.30,.56,.71]);cmd.set_color('tcrb',[.43,.68,.57]);cmd.set_color('peptideorange',[.91,.58,.19]);cmd.set_color('q7red',[.70,.09,.17])
for pdb in ['2ak4','4jry']:
    obj='s'+pdb
    cmd.show('cartoon',f'{obj} and chain A and resi 1-180')
    cmd.color('mhcgray',f'{obj} and chain A')
    for ch,col in [('D','tcra'),('E','tcrb')]:
        cmd.show('cartoon',f'{obj} and chain {ch} and resi 1-110')
        cmd.color(col,f'{obj} and chain {ch}')
    cmd.show('sticks',f'{obj} and chain C and not hydro')
    cmd.set('stick_radius',.16,f'{obj} and chain C')
    cmd.color('peptideorange',f'{obj} and chain C')
    cmd.color('q7red',f'{obj} and chain C and resi 7')
cmd.set_view(tuple(rot.T.flatten())+(0,0,-150)+tuple(pep.mean(0)+z*10)+(100,200,0))
cmd.zoom('(s2ak4 or s4jry) and ((chain A and resi 1-180) or (chain D+E and resi 1-110) or chain C)',buffer=2,complete=1)
view=list(cmd.get_view()); cmd.viewport(1600,1100)
# Identical view and scale for both independently shown complexes.
for pdb in ['2ak4','4jry']:
    cmd.disable('all');cmd.enable('s'+pdb);cmd.set_view(view)
    cmd.png(str(OUT/f'{pdb}_binding.png'),1600,1100,dpi=450,ray=1)
    cmd.scene('SB27' if pdb=='2ak4' else 'SB47','store')
cmd.save(str(OUT/'SB27_SB47_binding_modes.pse'))
(OUT/'molecular_view.json').write_text(json.dumps({'shared_view':view,'MHC_alignment':'CA residues1–180, no outlier rejection','colors':{'MHC':'gray','TCRalpha':'blue','TCRbeta':'green','peptide':'orange','Q7':'red'}},indent=2))
cmd.quit()
