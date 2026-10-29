"""Real peptide-bond context and accurate axes for the zoom inset."""
from pathlib import Path
import hashlib,json
import numpy as np
from PIL import Image
import pymol
from pymol import cmd

HERE=Path(__file__).resolve().parent;OUT=HERE/'assets'
SOURCE=HERE/'source_data/itrace_4g9f.pdb'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
pymol.finish_launching(['pymol','-cq'])
cmd.load(str(SOURCE),'assembly');cmd.remove('hydro');cmd.hide('everything')
cmd.select('shown','assembly and chain D and resi 165-167 and name N+CA+C+O')
cmd.select('target','assembly and chain D and resi 166')
cmd.show('sticks','shown');cmd.show('spheres','shown')
cmd.set_color('carbon_context',[.63,.76,.80]);cmd.color('carbon_context','shown')
cmd.set_color('target_gold',[.88,.64,.28]);cmd.color('target_gold','target')
cmd.set_color('nitrogen',[.34,.51,.73]);cmd.color('nitrogen','shown and elem N')
cmd.set_color('oxygen',[.76,.40,.43]);cmd.color('oxygen','shown and elem O')
cmd.set('stick_radius',.13);cmd.set('sphere_scale',.21)
cmd.set('sphere_scale',.29,'target and name CA')
cmd.set('sphere_quality',4);cmd.set('stick_quality',18);cmd.set('valence',0)
cmd.set('stick_ball',0);cmd.bg_color('white');cmd.set('ray_opaque_background',0)
cmd.set('orthoscopic',1);cmd.set('antialias',2);cmd.set('ray_shadows',0)
cmd.set('depth_cue',0);cmd.set('ambient',.52);cmd.set('specular',.10);cmd.set('ray_trace_mode',0)
v9=json.loads((HERE/'source_data/v9_local_projection.json').read_text())
cmd.set_view(v9['fixed_view'])
ca=np.array(cmd.get_atom_coords('target and name CA'))
n=np.array(cmd.get_atom_coords('target and name N'));c=np.array(cmd.get_atom_coords('target and name C'))
x=c-ca;x/=np.linalg.norm(x);y=n-ca-(n-ca).dot(x)*x;y/=np.linalg.norm(y);z=np.cross(x,y)
R=np.stack([x,y,z],axis=1)
outward=np.array([1.,-1.,1.]);outward/=np.linalg.norm(outward)
up=np.array([0.,0.,1.])-outward[2]*outward;up/=np.linalg.norm(up)
right=np.cross(up,outward)
camera=np.stack([right,up,outward])@R.T
new_view=list(cmd.get_view());new_view[:9]=camera.T.flatten().tolist();cmd.set_view(new_view)
cmd.zoom('shown',buffer=2.4)
cmd.viewport(1100,1050)
view=list(cmd.get_view())
cmd.png(str(OUT/'local_inset_raw.png'),width=1100,height=1050,dpi=600,ray=1)
cmd.save(str(OUT/'local_inset.pse'))
raw=Image.open(OUT/'local_inset_raw.png').convert('RGBA')
points={'origin':ca,'N':n,'C':c,'x':ca+3.6*x,'y':ca+3.6*y,'z':ca+3.6*z}
pixels={};cmd.hide('everything')
for key,coords in points.items():
    cmd.pseudoatom('probe',pos=coords.tolist(),vdw=.08);cmd.show('spheres','probe')
    cmd.color('red','probe');cmd.set_view(view)
    path=OUT/f'inset_probe_{key}.png'
    cmd.png(str(path),width=1100,height=1050,dpi=600,ray=1)
    bbox=Image.open(path).convert('RGBA').getchannel('A').getbbox();assert bbox is not None
    pixels[key]=[(bbox[0]+bbox[2])/2,(bbox[1]+bbox[3])/2]
    cmd.delete('probe')
bb=raw.getchannel('A').getbbox()
all_x=[p[0] for p in pixels.values()]+[bb[0],bb[2]]
all_y=[p[1] for p in pixels.values()]+[bb[1],bb[3]]
crop=[max(0,int(min(all_x))-75),max(0,int(min(all_y))-75),min(1100,int(max(all_x))+75),min(1050,int(max(all_y))+75)]
raw.crop(crop).save(OUT/'local_inset.png',dpi=(600,600))
pixels={k:[p[0]-crop[0],p[1]-crop[1]] for k,p in pixels.items()}
receipt=dict(source_sha256=sha(SOURCE),coordinates_modified=False,artificial_bonds_added=False,
             displayed_intervals=[['D',165,167]],axis_length_A=3.6,
             frame=R.tolist(),points={k:v.tolist() for k,v in points.items()},pixels=pixels,
             image_size=[crop[2]-crop[0],crop[3]-crop[1]],crop=crop,view=view,displayed_atoms=cmd.count_atoms('shown'),
             image_sha256=sha(OUT/'local_inset.png'),pymol_version=cmd.get_version()[0])
(OUT/'inset_receipt.json').write_text(json.dumps(receipt,indent=2),encoding='utf-8')
print('Inset rendered',flush=True);cmd.quit()
