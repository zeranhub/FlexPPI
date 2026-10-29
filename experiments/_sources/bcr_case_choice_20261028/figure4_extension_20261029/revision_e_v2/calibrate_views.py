from pathlib import Path
import os,json
import numpy as np
from PIL import Image
_dll=os.add_dll_directory(r'D:\Anaconda\envs\pvb_cpu\Library\bin')
import pymol
pymol.finish_launching(['pymol','-cq'])
from pymol import cmd
ROOT=Path(__file__).resolve().parent
data=json.loads((ROOT/'candidate_views.json').read_text())
cmd.load(data['source'],'fab');original=np.array(cmd.get_coords('fab'))
W,H=data['size'];results={}
for name,angle in [('R96',-30),('R101',0)]:
    rec=next(r for r in data['views'] if r['name']==name and r['angle']==angle)
    cmd.delete('probe*');cmd.delete('holdout');cmd.disable('fab')
    view=rec['view'];center=np.array(view[12:15])
    probes=[center+q for q in [np.array([0.,0.,0.]),np.array([2.,0.,0.]),np.array([0.,2.,0.]),np.array([-2.,-2.,2.])]]
    for i,(p,c) in enumerate(zip(probes,['red','green','blue','yellow'])):
        cmd.pseudoatom(f'probe{i}',pos=p.tolist());cmd.show('spheres',f'probe{i}');cmd.set('sphere_scale',.12,f'probe{i}');cmd.color(c,f'probe{i}')
    cmd.set('ambient',1);cmd.set('direct',0);cmd.set('specular',0);cmd.set('ray_shadows',0);cmd.set('ray_opaque_background',0);cmd.set('orthoscopic',1);cmd.set('depth_cue',0)
    cmd.viewport(W,H);cmd.set_view(view);cmd.png(str(ROOT/f'{name}_probes.png'),W,H,ray=1)
    arr=np.array(Image.open(ROOT/f'{name}_probes.png').convert('RGBA'));r,g,b,a=[arr[:,:,k].astype(int) for k in range(4)]
    masks=[(r>180)&(g<70)&(b<70)&(a>150),(g>150)&(r<70)&(b<70)&(a>150),(b>180)&(r<70)&(g<70)&(a>150),(r>180)&(g>180)&(b<70)&(a>150)]
    screen=[]
    for mask in masks:
        yy,xx=np.nonzero(mask);assert len(xx);screen.append([float(xx.mean()),float(yy.mean())])
    mapping=np.linalg.solve(np.column_stack([probes,np.ones(4)]),np.array(screen))
    cmd.delete('probe*');hold=center+np.array([.7,-.9,.3]);cmd.pseudoatom('holdout',pos=hold.tolist());cmd.show('spheres','holdout');cmd.set('sphere_scale',.12,'holdout');cmd.color('magenta','holdout');cmd.set_view(view)
    cmd.png(str(ROOT/f'{name}_holdout.png'),W,H,ray=1)
    arr=np.array(Image.open(ROOT/f'{name}_holdout.png').convert('RGBA'));r,g,b,a=[arr[:,:,k].astype(int) for k in range(4)]
    yy,xx=np.nonzero((r>180)&(b>180)&(g<70)&(a>150));err=float(np.linalg.norm(np.array([xx.mean(),yy.mean()])-np.r_[hold,1.]@mapping));assert err<1.
    anchors={k:(np.r_[p,1.]@mapping).tolist() for k,p in rec['anchors_camera'].items()}
    results[name]=dict(**rec,anchors_pixels=anchors,projection_error_pixels=err,size=[W,H])
(ROOT/'final_projection.json').write_text(json.dumps(results,indent=2),encoding='utf-8')
print(json.dumps(results,indent=2),flush=True)
cmd.quit()
