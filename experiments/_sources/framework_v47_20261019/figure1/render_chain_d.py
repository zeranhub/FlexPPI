"""Re-render the true complete D chain with automatically assigned ribbons.

No atom coordinates, residue identities or prediction values are changed.
"""
from pathlib import Path
import os, json, hashlib
_dll=os.add_dll_directory(r'D:\Anaconda\envs\pvb_cpu\Library\bin')
import numpy as np
from PIL import Image
import pymol
from pymol import cmd

OUT=Path(__file__).resolve().parent
SOURCE=OUT/'source_archive'
SOURCE.mkdir(exist_ok=True)
OLD=Path(r'C:\Users\DELL\.codex\worktrees\4666\FlexPPI\paper\revision_20261002\figures\main_revision_v13_20261008')
import shutil
for relative in ['source_data/itrace_4g9f.pdb','source_data/local_frame_candidates.json','assets/projection_receipt.json','assets/whole_chain_receipt.json']:
    shutil.copy2(OLD/relative,SOURCE/Path(relative).name)
receipt=json.loads((SOURCE/'whole_chain_receipt.json').read_text(encoding='utf8'))
projection=json.loads((SOURCE/'projection_receipt.json').read_text(encoding='utf8'))
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
assert sha(SOURCE/'itrace_4g9f.pdb')==receipt['source_sha256']
pymol.finish_launching(['pymol','-cq'])
cmd.load(str(SOURCE/'itrace_4g9f.pdb'),'assembly')
cmd.hide('everything')
cmd.select('chain_d','assembly and polymer.protein and chain D')
oldcoords=cmd.get_coords('chain_d').copy()
cmd.dss('assembly and chain D')
model=cmd.get_model('chain_d and name CA')
ss={k:sum(a.ss==k for a in model.atom) for k in ['H','S','L','']}
assert len(model.atom)==204
cmd.show('cartoon','chain_d')
cmd.cartoon('automatic','chain_d')
cmd.set_color('protein_blue',[.62,.73,.79]);cmd.color('protein_blue','chain_d')
cmd.set_color('local_gold',[.84,.64,.32]);cmd.color('local_gold','chain_d and resi 166')
cmd.set('cartoon_fancy_helices',1)
cmd.set('cartoon_flat_sheets',1)
cmd.set('cartoon_smooth_loops',1)
cmd.set('cartoon_highlight_color',-1)
cmd.set('cartoon_loop_radius',.17)
cmd.set('cartoon_tube_radius',.17)
cmd.set('cartoon_rect_length',1.1)
cmd.set('cartoon_rect_width',.19)
cmd.bg_color('white');cmd.set('ray_opaque_background',0)
cmd.set('orthoscopic',1);cmd.set('antialias',2)
cmd.set('ray_shadows',0);cmd.set('depth_cue',0)
cmd.set('ambient',.52);cmd.set('direct',.40)
cmd.set('specular',.10);cmd.set('shininess',20)
cmd.set('ray_trace_mode',0)
cmd.set('ray_trace_gain',0)
view=receipt['view'];w,h=receipt['raw_image_size']
cmd.set_view(view);cmd.viewport(w,h)
cmd.png(str(OUT/'chain_d_ribbon_raw.png'),width=w,height=h,dpi=600,ray=1)
assert np.array_equal(oldcoords,cmd.get_coords('chain_d'))
cmd.save(str(OUT/'chain_d_ribbon.pse'))
raw=Image.open(OUT/'chain_d_ribbon_raw.png').convert('RGBA')
bb=raw.getchannel('A').getbbox()
crop=[max(0,bb[0]-45),max(0,bb[1]-45),min(w,bb[2]+45),min(h,bb[3]+45)]
raw.crop(crop).save(OUT/'chain_d_ribbon.png',dpi=(600,600))
new=dict(source_pdb_sha256=receipt['source_sha256'],atom_coordinates_changed=False,
         displayed_chain='D',displayed_residues=204,secondary_structure_assignment='PyMOL dss on actual assembly chain D',
         secondary_structure_counts=ss,cartoon='automatic (normal beta arrows and helices, not forced tube)',
         camera_view=view,raw_size=[w,h],crop=crop,image_size=[crop[2]-crop[0],crop[3]-crop[1]],
         pdb_sha256=sha(SOURCE/'itrace_4g9f.pdb'),image_sha256=sha(OUT/'chain_d_ribbon.png'),
         alpha_bbox_in_raw=list(bb),pymol_version=cmd.get_version()[0],
         pixels={k:[v[0]+receipt['crop'][0]-crop[0],v[1]+receipt['crop'][1]-crop[1]] for k,v in projection['pixels'].items()},
         camera_nominal_projection_px_per_A=projection['camera_nominal_projection_px_per_A'],
         camera_matrix=projection['PyMOL_rotation_convention']['camera_matrix'])
(OUT/'chain_d_ribbon_receipt.json').write_text(json.dumps(new,indent=2),encoding='utf8')
print(json.dumps(new),flush=True)
cmd.quit()
