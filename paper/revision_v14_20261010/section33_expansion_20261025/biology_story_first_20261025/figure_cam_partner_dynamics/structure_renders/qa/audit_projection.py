"""Independently validate the manifest's 3-D to PNG projection with CA markers."""
from pathlib import Path
import os
import json
import numpy as np
from PIL import Image

ENV=Path('D:/Anaconda/envs/pvb_cpu')
os.environ['PATH']=str(ENV/'Library/bin')+os.pathsep+os.environ.get('PATH','')
if hasattr(os,'add_dll_directory'):
    DLL_HANDLES=[os.add_dll_directory(str(p)) for p in [ENV/'Library/bin',ENV/'DLLs'] if p.exists()]
import pymol2

OUT=Path(__file__).resolve().parent
ROOT=OUT.parent
manifest=json.loads((ROOT/'render_manifest.json').read_text())
asset=manifest['assets']['overview_nNOS']
width,height=asset['width'],asset['height']
view=np.array(asset['view'])
colors={'L18':(1.,0.,0.),'L39':(0.,1.,0.),'L105':(0.,0.,1.),'M109':(1.,0.,1.),'M144':(0.,1.,1.)}
rows=[]
with pymol2.PyMOL() as pm:
    c=pm.cmd
    c.load(str(ROOT/'cam_partner_dynamics.pse'))
    c.disable('all');c.enable('nNOS');c.hide('everything','all')
    for setting,value in dict(orthoscopic=1,ray_opaque_background=1,ray_shadows=0,
            depth_cue=0,ambient=1,direct=0,specular=0,antialias=2,max_threads=2,
            sphere_scale=.35).items(): c.set(setting,value)
    c.bg_color('white')
    for label,color in colors.items():
        selection=f'nNOS and chain A and resi {label[1:]} and name CA'
        assert c.count_atoms(selection)==1
        c.set_color('qa_'+label,list(color));c.color('qa_'+label,selection)
        c.show('spheres',selection)
    c.set_view(view.tolist())
    actual_view=list(c.get_view());fov=float(c.get('field_of_view'))
    c.png(str(OUT/'overview_nNOS_CA_projection_markers.png'),width=width,height=height,dpi=300,ray=1,quiet=1)
    pixels=np.asarray(Image.open(OUT/'overview_nNOS_CA_projection_markers.png').convert('RGB'))/255.
    halfh=-view[11]*np.tan(np.deg2rad(fov/2));halfw=halfh*width/height
    for label,color in colors.items():
        xyz=c.get_coords(f'nNOS and chain A and resi {label[1:]} and name CA',1)[0]
        qT=view[:9].reshape(3,3).T@(xyz-view[12:15])+view[9:12]
        qN=view[:9].reshape(3,3)@(xyz-view[12:15])+view[9:12]
        def project(q):return np.array([(.5+q[0]/(2*halfw))*width,(.5-q[1]/(2*halfh))*height])
        pT=project(qT);pN=project(qN)
        # Flat colored sphere pixels are isolated by proximity to the unique RGB.
        distance=np.linalg.norm(pixels-np.array(color)[None,None,:],axis=2)
        y,x=np.nonzero(distance<.25)
        assert len(x)>10,label
        centroid=np.array([x.mean()+.5,y.mean()+.5])
        stored=np.array(asset['landmarks'][label])
        rows.append(dict(label=label,xyz=xyz.tolist(),manifest_xy=stored.tolist(),
                         formula_transpose_xy=pT.tolist(),formula_no_transpose_xy=pN.tolist(),
                         rendered_marker_centroid_xy=centroid.tolist(),n_marker_pixels=len(x),
                         marker_bbox=[int(x.min()),int(y.min()),int(x.max()),int(y.max())],
                         manifest_error_px=float(np.linalg.norm(stored-centroid)),
                         transpose_error_px=float(np.linalg.norm(pT-centroid)),
                         no_transpose_error_px=float(np.linalg.norm(pN-centroid))))
payload=dict(asset='overview_nNOS',width=width,height=height,field_of_view=fov,
             method='Restore final PSE, keep nNOS only, replace molecular representations with five unique-color CA spheres, ray-render at the manifest dimensions, compare independent RGB marker centroids.',
             actual_view=actual_view,view_max_delta=float(np.max(np.abs(np.array(actual_view)-view))),
             pixel_convention='Detected centroids are pixel centers, +0.5 on integer PIL indices.',
             root_formula='view[:9].reshape(3,3).T @ (xyz-view[12:15])+view[9:12]; orthographic scale from -view[11]*tan(fov/2)',
             rows=rows,passed=all(r['manifest_error_px']<3 for r in rows))
(OUT/'projection_audit.json').write_text(json.dumps(payload,indent=2),encoding='utf-8')
print(json.dumps(payload,indent=2))
