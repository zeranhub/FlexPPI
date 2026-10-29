"""Render only existing KIX/ubiquitin molecular coordinates; no generated atoms."""
from pathlib import Path
import json,hashlib,os
_dll=os.add_dll_directory(r'D:\Anaconda\envs\pvb_cpu\Library\bin')
import numpy as np
import pymol
pymol.finish_launching(['pymol','-cq'])
from pymol import cmd
BASE=Path(r'E:\Project\FlexPPI\experiments\_sources')
SRC=BASE/'reader_revision_v39_20261016/english_figures/fig6'
OUT=BASE/'candidate_case_atlas_20261021/evidence/non_tcr_figures'
OUT.mkdir(exist_ok=True,parents=True)
renders=OUT/'structure_renders';renders.mkdir(exist_ok=True)
audit=[]
cmd.load(str(SRC/'kix/KIX_case_validation.pse'))
for scene,stem,size in [('KIX_context','kix_context',(1800,1500)),('KIX_L12_MD','kix_l12_md',(1800,1600)),('KIX_L12_FlexPPI_fold1','kix_l12_model',(1800,1600))]:
    assert scene in cmd.get_scene_list(),(scene,cmd.get_scene_list())
    cmd.scene(scene,'recall',animate=0)
    cmd.hide('labels','all')
    cmd.bg_color('white');cmd.set('ray_opaque_background',0);cmd.set('ray_shadows',0)
    path=renders/(stem+'.png')
    cmd.png(str(path),width=size[0],height=size[1],dpi=600,ray=1)
    audit.append(dict(stem=stem,pse=str(SRC/'kix/KIX_case_validation.pse'),scene=scene,labels='removed; vector labels in final composition',render=str(path),pixel_size=list(size)))
cmd.reinitialize()
ub=SRC/'ubiquitin/inputs/A003R_input_structure.pdb'
cmd.load(str(ub),'ub')
before=np.asarray(cmd.get_coords('ub'),dtype=float)
cmd.hide('everything','all');cmd.show('cartoon','ub')
cmd.bg_color('white');cmd.set('ray_opaque_background',0);cmd.set('ray_shadows',0)
cmd.set('orthoscopic',1);cmd.set('specular',.15);cmd.set('antialias',2)
cmd.set_color('ubgray',[.68,.71,.73]);cmd.color('ubgray','ub and chain A')
cmd.set_color('ubagreen',[.106,.471,.216]);cmd.color('ubagreen','ub and chain B')
for pos,col in [(2,'#762A83'),(41,'#F1A340')]:
    cmd.set_color('site'+str(pos),[int(col[i:i+2],16)/255 for i in [1,3,5]])
    sel=f'ub and chain A and resi {pos}'
    assert cmd.count_atoms(sel+' and name CA')==1
    cmd.show('sticks',sel+' and not hydro');cmd.set('stick_radius',.18,sel)
    cmd.show('spheres',sel+' and name CA');cmd.set('sphere_scale',.35,sel+' and name CA')
    cmd.color('site'+str(pos),sel)
view=json.loads((SRC/'ubiquitin/views.json').read_text(encoding='utf-8'))['01_complex']
cmd.set_view(view);cmd.zoom('ub',buffer=3,complete=1)
cmd.png(str(renders/'ubiquitin_context.png'),width=1800,height=1500,dpi=600,ray=1)
after=np.asarray(cmd.get_coords('ub'),dtype=float)
assert np.array_equal(before,after)
cmd.save(str(renders/'ubiquitin_Q2_Q41_input.pse'))
audit.append(dict(stem='ubiquitin_context',source=str(ub),source_sha256=hashlib.sha256(ub.read_bytes()).hexdigest(),chains={'A':cmd.count_atoms('chain A and name CA'),'B':cmd.count_atoms('chain B and name CA')},coordinates_unchanged=True,highlighted_sites=[2,41]))
(OUT/'structure_render_audit.json').write_text(json.dumps(dict(status='passed',renders=audit),ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(dict(status='passed',output=str(renders)),ensure_ascii=False))
cmd.quit()
