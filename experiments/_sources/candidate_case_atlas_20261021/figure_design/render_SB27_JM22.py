"""True model-input coordinates, matching cameras after full MHC C-alpha alignment."""
import os
os.add_dll_directory(r'D:\Anaconda\envs\pvb_cpu\Library\bin')
from pathlib import Path
import json, hashlib
import numpy as np
import pymol
pymol.finish_launching(['pymol','-cq'])
from pymol import cmd
ROOT = Path(__file__).resolve().parent
EXP = ROOT.parent.parent
CASES = [
    ('SB27_SB47', ['2ak4','4jry'], EXP/'tcr_case_v40_20261017'/'sb27_sb47', 7),
    ('JM22_F5L', ['1oga','5hhm'], EXP/'tcr_expanded_screen_v41_20261017'/'jm22_f5l_validation', 5),
]
for name, pdbs, inputdir, focal in CASES:
    out=ROOT/name
    out.mkdir(exist_ok=True)
    cmd.reinitialize()
    paths=[inputdir/f'{pdb}_model_input.pdb' for pdb in pdbs]
    for i,path in enumerate(paths):cmd.load(str(path),f's{i}')
    alignment=cmd.align('s1 and chain A and resi 1-180 and name CA',
                        's0 and chain A and resi 1-180 and name CA',cycles=0)
    xyz=lambda sel:np.asarray(cmd.get_coords(sel),float)
    pep=xyz('s0 and chain C and name CA')
    pepb=xyz('s1 and chain C and name CA')
    assert pep.shape==pepb.shape
    peptide_rmsd=float(np.sqrt(np.mean(np.sum((pep-pepb)**2,axis=1))))
    x=pep[-1]-pep[0];x/=np.linalg.norm(x)
    z=xyz('s0 and chain D+E and resi 1-110 and name CA').mean(0)-pep.mean(0)
    z-=x*np.dot(z,x);z/=np.linalg.norm(z)
    y=np.cross(z,x);zview=.6*z-.8*y;yview=np.cross(zview,x)
    rot=np.stack([x,yview,zview])
    cmd.hide('everything','all')
    cmd.bg_color('white')
    for option,value in [('orthoscopic',1),('ray_shadows',0),('specular',.2),
                         ('antialias',2),('cartoon_fancy_helices',1),('cartoon_smooth_loops',1)]:cmd.set(option,value)
    cmd.set_color('mhcgray',[.72,.77,.78])
    cmd.set_color('tcrgray',[.57,.59,.61])
    cmd.set_color('first',[207/255,140/255,39/255])
    cmd.set_color('second',[118/255,42/255,131/255])
    cmd.set_color('focal',[41/255,123/255,130/255])
    for i,col in enumerate(['first','second']):
        obj=f's{i}'
        cmd.show('cartoon',f'{obj} and chain A and resi 1-180')
        cmd.color('mhcgray',f'{obj} and chain A')
        cmd.show('cartoon',f'{obj} and chain D+E and resi 1-110')
        cmd.color('tcrgray',f'{obj} and chain D+E')
        cmd.show('sticks',f'{obj} and chain C and not hydro')
        cmd.set('stick_radius',.18,f'{obj} and chain C')
        cmd.color(col,f'{obj} and chain C')
        cmd.show('spheres',f'{obj} and chain C and resi {focal} and name CA')
        cmd.set('sphere_scale',.35,f'{obj} and chain C and resi {focal}')
        cmd.color('focal',f'{obj} and chain C and resi {focal} and name CA')
    cmd.set_view(tuple(rot.T.flatten())+(0,0,-150)+tuple(pep.mean(0)+z*10)+(100,200,0))
    cmd.zoom('(s0 or s1) and ((chain A and resi 1-180) or (chain D+E and resi 1-110) or chain C)',buffer=2,complete=1)
    view=list(cmd.get_view())
    for i,pdb in enumerate(pdbs):
        cmd.disable('all');cmd.enable(f's{i}');cmd.set_view(view)
        cmd.png(str(out/f'{pdb}_binding.png'),1600,1100,dpi=450,ray=1)
        cmd.scene(pdb,'store')
    cmd.enable('all')
    cmd.hide('cartoon','all')
    cmd.hide('spheres','all')
    cmd.show('spheres',f'chain C and resi {focal} and name CA')
    cmd.zoom('chain C',buffer=1.6,complete=1)
    cmd.png(str(out/'peptide_superposition.png'),1500,800,dpi=450,ray=1)
    cmd.save(str(out/f'{name}_model_input_views.pse'))
    metadata={'status':'passed','source':'Canonical preprocessed model-input PDBs, not a redrawn cartoon',
              'MHC_fit':'All matching chain A C-alpha residues1-180; cycles0 means no outlier rejection',
              'MHC_alignment_RMSD_A':float(alignment[0]),'MHC_alignment_atoms':int(alignment[1]),
              'peptide_CA_RMSD_after_MHC_alignment_A':peptide_rmsd,'shared_overview_view':view,
              'peptide_only_view':list(cmd.get_view()),'focal_position':focal,
              'input_paths':[str(p) for p in paths],
              'input_SHA256':[hashlib.sha256(p.read_bytes()).hexdigest() for p in paths],
              'no_independent_peptide_fit':True}
    (out/'molecular_view.json').write_text(json.dumps(metadata,indent=2),encoding='utf8')
    print(json.dumps({'case':name,'peptide_RMSD_A':peptide_rmsd,'status':'rendered'}),flush=True)
cmd.quit()
