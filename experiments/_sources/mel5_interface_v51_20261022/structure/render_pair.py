"""Real-coordinate two-state local views. A single common orthographic camera."""
from pathlib import Path
import json,os,math,hashlib
import numpy as np
from PIL import Image
_dll=os.add_dll_directory(r'D:\Anaconda\envs\pvb_cpu\Library\bin')
import pymol
pymol.finish_launching(['pymol','-cq'])
from pymol import cmd
ROOT=Path(__file__).resolve().parent
ROOT.mkdir(parents=True,exist_ok=True)
PAR=ROOT.parent
OLD=PAR.parent/'mel5_joint_case_v49_20261019/figure4/source_archive'
WT=OLD/'3hg1_model_input.pdb'
ENG=OLD/'4jff_model_input.pdb'
cmd.load(str(WT),'wt');cmd.load(str(ENG),'eng')
fit=json.loads((PAR/'analysis/geometry_audit.json').read_text(encoding='utf8'))['HLA_CA_alignment']
cmd.load_coords(np.asarray(cmd.get_coords('wt'),float)@np.array(fit['rotation_row_vector'])+np.array(fit['translation']),'wt')
cmd.hide('everything','all');cmd.bg_color('white')
for setting,value in [('orthoscopic',1),('ray_shadows',0),('specular',.12),('antialias',2),('cartoon_fancy_helices',1),('cartoon_smooth_loops',1),('ray_opaque_background',0),('depth_cue',0),('ray_trace_fog',0)]: cmd.set(setting,value)
cmd.set_color('mhc',[.16,.48,.51]);cmd.set_color('beta',[.4627,.1647,.5137]);cmd.set_color('q72',[.16,.48,.51]);cmd.set_color('WT_beta',[.8118,.5490,.1529])
for obj, b52,b53,b54 in [('wt',52,53,54),('eng',53,54,55)]:
    cmd.color('mhc',f'{obj} and chain A')
    cmd.color('WT_beta' if obj=='wt' else 'beta',f'{obj} and chain E')
    cmd.show('cartoon',f'{obj} and ((chain A and resi 68-77) or (chain E and resi {b52-3}-{b54+2}))')
    cmd.cartoon('tube',f'{obj} and ((chain A and resi 68-77) or (chain E and resi {b52-3}-{b54+2}))')
    cmd.set('cartoon_tube_radius',.20)
    cmd.set('cartoon_transparency',.60,f'{obj}')
    sticks=f'{obj} and ((chain A and resi 72) or (chain E and resi {b52}-{b54})) and not hydro'
    cmd.show('sticks',sticks);cmd.set('stick_radius',.18,sticks)
    # Mutated beta52 is visible, while the key beta53/beta54 pair remains visually dominant.
    cmd.set('stick_transparency',.35,f'{obj} and chain E and resi {b52}')
    cmd.color('q72',f'{obj} and chain A and resi 72')
    cmd.color('oxygen',f'{sticks} and elem O');cmd.color('nitrogen',f'{sticks} and elem N')
    cmd.distance(f'{obj}_d1',f'{obj} and chain A and resi 72 and name OE1',f'{obj} and chain E and resi {b53} and name N')
    cmd.distance(f'{obj}_d2',f'{obj} and chain A and resi 72 and name NE2',f'{obj} and chain E and resi {b54} and name O')
cmd.hide('labels','all');cmd.set('dash_color','gray45');cmd.set('dash_width',2.6);cmd.set('dash_gap',.22);cmd.set('dash_length',.19)

def xyz(sel):
    coords=cmd.get_coords(sel)
    return np.zeros((0,3)) if coords is None else np.asarray(coords,float)
q=xyz('eng and chain A and resi 72 and name CA')[0]
f=xyz('eng and chain E and resi 54 and name CA')[0]
g=xyz('eng and chain E and resi 55 and name CA')[0]
x=g-f;x/=np.linalg.norm(x)
y=(f+g)/2-q;y-=x*np.dot(y,x);y/=np.linalg.norm(y)
z=np.cross(x,y)
yview=.91*y+.414*z;yview/=np.linalg.norm(yview)
zview=np.cross(x,yview)
rot=np.stack([x,yview,zview])
center=(q+f+g)/3
cmd.viewport(1800,1500)
cmd.set_view(tuple(rot.T.flatten())+(0,0,-50)+tuple(center)+(20,80,0))
# Only this union selection determines scale and framing for both states.
union='((wt and chain E and resi 52-54) or (eng and chain E and resi 53-55) or (chain A and resi 72)) and not hydro'
cmd.zoom(union,buffer=2.0,complete=1)
cmd.set('orthoscopic',1)
view=list(cmd.get_view())
meta={'status':'generated_pending_visual_QA','sources':{'WT_original':str(WT),'mature':str(ENG),'alignment_audit':str(PAR/'analysis/geometry_audit.json')},'view':view,'image_size':[1800,1500],'camera':'Single orthographic camera, fitted on union of both states, then unchanged. WT is transformed in memory at full float precision; the separately written aligned PDB is not the numerical source for the rendering.', 'numbering':'Common MEL5 numbering; mature input beta = common+1.', 'structures':{}}
for obj,b53,b54 in [('wt',53,54),('eng',54,55)]:
    cmd.disable('all');cmd.enable(obj);cmd.enable(obj+'_d1');cmd.enable(obj+'_d2')
    cmd.set_view(view)
    name='MEL5_local_common_view.png' if obj=='wt' else 'alpha24beta17_local_common_view.png'
    cmd.png(str(ROOT/name),1800,1500,dpi=450,ray=1)
    sels={'Q72_CA':f'{obj} and chain A and resi 72 and name CA','Q72_OE1':f'{obj} and chain A and resi 72 and name OE1','Q72_NE2':f'{obj} and chain A and resi 72 and name NE2','beta52_CA':f'{obj} and chain E and resi {b53-1} and name CA','beta53_CA':f'{obj} and chain E and resi {b53} and name CA','beta53_N':f'{obj} and chain E and resi {b53} and name N','beta54_CA':f'{obj} and chain E and resi {b54} and name CA','beta54_O':f'{obj} and chain E and resi {b54} and name O'}
    coords={k:xyz(s)[0].tolist() for k,s in sels.items()}
    for label,resi in [('beta52_sidechain_center',b53-1),('beta53_sidechain_center',b53)]:
        sc=xyz(f'{obj} and chain E and resi {resi} and not hydro and not name N+CA+C+O+OXT')
        coords[label]=(sc.mean(0) if len(sc) else xyz(f'{obj} and chain E and resi {resi} and name CA')[0]).tolist()
    meta['structures'][obj]={'PNG':str(ROOT/name),'atom_coordinates':coords,'distances_A':{'Q72_OE1_beta53_N':float(cmd.get_distance(sels['Q72_OE1'],sels['beta53_N'])),'Q72_NE2_beta54_O':float(cmd.get_distance(sels['Q72_NE2'],sels['beta54_O']))}}
cmd.save(str(ROOT/'MEL5_dual_state_common_view.pse'))

# Empirical projection calibration: same view and canvas, three isolated small
# colored spheres. Their image centers validate the screen mapping rather than
# relying on an undocumented matrix transpose convention.
cmd.disable('all')
probe_points=[np.array(center),np.array(center)+np.array([3.,0.,0.]),np.array(center)+np.array([0.,3.,0.]),np.array(center)+np.array([0.,0.,3.])]
probe_colors=['red','green','blue','yellow']
for i,(p,c) in enumerate(zip(probe_points,probe_colors)):
    cmd.pseudoatom(f'probe{i}',pos=p.tolist());cmd.show('spheres',f'probe{i}');cmd.set('sphere_scale',.15,f'probe{i}');cmd.color(c,f'probe{i}')
cmd.set('ambient',1);cmd.set('direct',0);cmd.set('specular',0);cmd.set_view(view)
cmd.png(str(ROOT/'projection_calibration.png'),1800,1500,dpi=450,ray=1)
im=np.array(Image.open(ROOT/'projection_calibration.png').convert('RGBA'))
r,g,b,a=[im[:,:,k].astype(int) for k in range(4)]
masks=[(r>180)&(g<70)&(b<70)&(a>150),(g>150)&(r<70)&(b<70)&(a>150),(b>180)&(r<70)&(g<70)&(a>150),(r>180)&(g>180)&(b<70)&(a>150)]
centers=[]
for mask in masks:
    yy,xx=np.nonzero(mask)
    if len(xx)==0: raise RuntimeError('Projection probe not visible')
    centers.append([float(xx.mean()),float(yy.mean())])
design=np.column_stack([np.array(probe_points),np.ones(4)])
mapping=np.linalg.solve(design,np.array(centers))
# Holdout coordinate, independent fifth magenta sphere, same camera.
cmd.disable('all');holdout=np.mean(probe_points,axis=0)+np.array([1.0,-1.2,.7])
cmd.pseudoatom('holdout',pos=holdout.tolist());cmd.show('spheres','holdout');cmd.set('sphere_scale',.15,'holdout');cmd.color('magenta','holdout');cmd.set_view(view)
cmd.png(str(ROOT/'projection_holdout.png'),1800,1500,dpi=450,ray=1)
hi=np.array(Image.open(ROOT/'projection_holdout.png').convert('RGBA'));hr,hg,hb,ha=[hi[:,:,k].astype(int) for k in range(4)]
yy,xx=np.nonzero((hr>180)&(hb>180)&(hg<70)&(ha>150));observed=np.array([xx.mean(),yy.mean()]);predicted=np.r_[holdout,1.]@mapping
error=float(np.linalg.norm(observed-predicted))
if error>2.0: raise RuntimeError(f'Projection holdout error {error}')
for state in meta['structures'].values():
    state['atom_screen_pixels_top_left']={k:(np.r_[p,1.]@mapping).tolist() for k,p in state['atom_coordinates'].items()}
    state['PNG_sha256']=hashlib.sha256(Path(state['PNG']).read_bytes()).hexdigest()
    png_image=Image.open(state['PNG']).convert('RGBA')
    state['nontransparent_pixel_bbox']=list(png_image.getchannel('A').getbbox())
    state['alpha_range']=list(png_image.getchannel('A').getextrema())
    state['beta52_sidechain_center_note']='WT Gly52 has no heavy sidechain; its sidechain_center label anchor is CA. Mature Pro52 and beta53 use the mean of sidechain heavy atoms.'
meta['screen_projection']={'mapping_xyz1_to_xy':mapping.tolist(),'calibration_probe_world':np.array(probe_points).tolist(),'calibration_probe_screen':centers,'holdout_world':holdout.tolist(),'holdout_predicted':predicted.tolist(),'holdout_observed':observed.tolist(),'holdout_error_pixels':error,'note':'Pixel coordinates measured at original 1800x1500 canvas; origin upper-left.'}
(ROOT/'dual_state_molecular_view.json').write_text(json.dumps(meta,indent=2),encoding='utf8')
print(json.dumps({'status':'done','distances':{k:v['distances_A'] for k,v in meta['structures'].items()},'projection_holdout_error_pixels':error},indent=2))
cmd.quit()
