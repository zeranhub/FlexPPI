"""Two native M316 contact contexts with the exact 1MXE overview camera rotation.
Input structures are frozen observed heavy atoms. No new fitting or atom edits.
1A06 is static literature context, without model prediction or MD.
"""
from pathlib import Path
import os, json, hashlib
import numpy as np
ENV=Path('D:/Anaconda/envs/pvb_cpu')
os.environ['PATH']=str(ENV/'Library/bin')+os.pathsep+os.environ.get('PATH','')
DLL_HANDLES=[os.add_dll_directory(str(p)) for p in [ENV/'Library/bin',ENV/'DLLs'] if p.exists()]
import pymol2
from PIL import Image
OUT=Path(__file__).resolve().parent
BASE=json.loads((OUT/'inputs/1MXE_camera_manifest.json').read_text())
ROTATION=BASE['assets']['overview']['view'][:9]
COLORS={'background':'#8A9299','partner':'#8B72A5','key':'#3A8D8C'}
INPUTS={'cam_bound':OUT/'inputs/1MXE_CaM_native4_147_CaMKI294_318.pdb','autoinhibited':OUT/'inputs/1A06_rigidly_aligned_to_1MXE.pdb'}
W,H=1600,1600

def project_xyz(c,xyz):
    v=np.asarray(c.get_view());q=v[:9].reshape(3,3).T@(np.asarray(xyz)-v[12:15])+v[9:12]
    hh=-v[11]*np.tan(np.deg2rad(float(c.get('field_of_view'))/2));hw=hh*W/H
    return [float((.5+q[0]/(2*hw))*W),float((.5-q[1]/(2*hh))*H)]

def selectpoint(c,sel):return c.get_coords(sel).mean(0)

def init(c):
    for k,v in dict(ray_opaque_background=0,opaque_background=0,ray_trace_mode=0,ray_shadows=0,antialias=3,orthoscopic=1,depth_cue=0,specular=.12,shininess=35,ambient=.65,direct=.35,reflect=.12,light_count=2,cartoon_fancy_helices=1,cartoon_fancy_sheets=1,cartoon_flat_sheets=1,cartoon_smooth_loops=1,cartoon_sampling=20,cartoon_loop_radius=.23,cartoon_tube_radius=.23,cartoon_oval_length=1.25,cartoon_oval_width=.28,max_threads=4,stick_radius=.20,valence=0).items():c.set(k,v)
    for k,h in COLORS.items():c.set_color(k,[int(h[i:i+2],16)/255 for i in [1,3,5]])
    c.bg_color('white');c.viewport(W,H)

def style(c,name):
    obj=name;c.load(str(INPUTS[name]),obj);c.hide('everything');c.dss(obj)
    peptide_chain='B' if name=='cam_bound' else 'A'
    reg=f'{obj} and chain {peptide_chain} and resi 312-318'
    anchor=f'{obj} and chain {peptide_chain} and resi 316'
    c.show('cartoon',reg);c.color('partner',reg);c.set('cartoon_color','partner',reg)
    if name=='cam_bound':
        neighbors=f'{obj} and chain A and resi 32+51+52+55+63+71'
        local=f'byres ({obj} and chain A and polymer.protein within 7.0 of ({anchor}))'
        c.show('cartoon',local);c.color('background',local);c.set('cartoon_color','background',local);c.set('cartoon_transparency',.58,local)
        marks={f'{aa}{n}':f'{obj} and chain A and resi {n} and name '+('SD' if aa=='M' else 'CG' if aa=='L' else 'CG1' if aa=='I' else 'CB') for aa,n in [('L',32),('M',51),('I',52),('V',55),('I',63),('M',71)]}
    else:
        neighbors=f'{obj} and chain A and resi 45+99'
        local=neighbors
        marks={'L45':f'{obj} and chain A and resi 45 and name CD1','S99':f'{obj} and chain A and resi 99 and name OG'}
    c.show('sticks',neighbors+' and not name N+C+O');c.color('background',neighbors)
    c.set('stick_radius',.17,neighbors)
    c.show('sticks',anchor+' and not name N+C+O');c.color('key',anchor+' and not name N+C+O');c.set('stick_radius',.27,anchor)
    marks['M316']=anchor+' and name CE';marks['M316_CA']=anchor+' and name CA';marks['M316_SD']=anchor+' and name SD'
    visible=reg+' or '+local+' or '+neighbors
    c.set_view(BASE['assets']['overview']['view']);c.zoom(visible,buffer=2.0,complete=1)
    v=list(c.get_view());v[:9]=ROTATION;v[15]=1.;v[16]=-v[11]+180.;c.set_view(v)
    return visible,marks

def main():
    manifest={'scope':'M316 static contact contexts from observed heavy atoms. Both objects share all nine original 1MXE overview rotation values. The pre-aligned 1A06 input undergoes no further coordinate changes.',
      'limitations':['1A06 is a static autoinhibited kinase structure; no model prediction or MD here.','1MXE contains a CaM-bound regulatory peptide, not a full activated kinase.','1A06 residues317-318 are unobserved and absent.','Dashed distances are heavy-atom separations, without a hydrogen-bond assignment.'],
      'palette':COLORS,'camera_rotation':ROTATION,'assets':{},'geometry_audit':{},'contacts':[]}
    # Determine one optical scale, in the same camera orientation, for both contexts.
    views={}
    with pymol2.PyMOL() as pm:
        c=pm.cmd;init(c)
        for name in INPUTS:
            c.delete('all');visible,marks=style(c,name);views[name]=list(c.get_view())
    common_depth=max(-v[11] for v in views.values())
    for name,path in INPUTS.items():
        with pymol2.PyMOL() as pm:
            c=pm.cmd;init(c);visible,marks=style(c,name)
            before=c.get_coords(name).copy();v=list(c.get_view());v[:9]=ROTATION;v[9:12]=[0.,0.,-common_depth];v[15]=1.;v[16]=common_depth+180.;c.set_view(v)
            landmarks={key:project_xyz(c,selectpoint(c,sel)) for key,sel in marks.items()}
            if name=='autoinhibited':
                for partner,leftatom,rightatom in [(45,'SD','CD1'),(99,'CE','OG')]:
                    left=f'{name} and chain A and resi 316 and name {leftatom}';right=f'{name} and chain A and resi {partner} and name {rightatom}'
                    dist=c.distance(f'distance_{partner}',left,right,mode=0);c.hide('labels',f'distance_{partner}')
                    for k,value in [('dash_color','background'),('dash_width',3.6),('dash_gap',.15),('dash_length',.18)]:c.set(k,value,f'distance_{partner}')
                    midpoint=(selectpoint(c,left)+selectpoint(c,right))/2
                    landmarks[f'distance_{partner}_midpoint']=project_xyz(c,midpoint)
                    manifest['contacts'].append({'M316_atom':leftatom,'partner_residue':partner,'partner_atom':rightatom,'distance_A':float(dist)})
            c.png(str(OUT/(name+'.png')),width=W,height=H,dpi=600,ray=1,quiet=1)
            c.save(str(OUT/(name+'.pse')))
            im=Image.open(OUT/(name+'.png')).convert('RGBA');bbox=list(im.getchannel('A').getbbox());pad=36
            crop=[max(0,bbox[0]-pad),max(0,bbox[1]-pad),min(W,bbox[2]+pad),min(H,bbox[3]+pad)]
            im.crop(crop).save(OUT/(name+'_trim.png'))
            manifest['assets'][name]={'file':name+'.png','width':W,'height':H,'view':list(c.get_view()),'alpha_bbox':bbox,'landmarks':landmarks,
             'trim':{'file':name+'_trim.png','crop_box':crop,'width':crop[2]-crop[0],'height':crop[3]-crop[1],'landmarks':{key:[p[0]-crop[0],p[1]-crop[1]] for key,p in landmarks.items()}}}
            after=c.get_coords(name)
            manifest['geometry_audit'][name]={'input_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'atom_count':len(before),'coordinate_change_max_A':float(np.abs(before-after).max()),'coordinates_exactly_unchanged':bool(np.array_equal(before,after)),'rotation_exactly_inherited':bool(np.array_equal(np.asarray(c.get_view()[:9]),np.asarray(ROTATION))),'visible_atoms':c.count_atoms(visible)}
            print('RENDERED',name,flush=True)
    manifest['shared_optical_scale_depth_A']=common_depth
    manifest['same_camera_rotation']=manifest['assets']['cam_bound']['view'][:9]==manifest['assets']['autoinhibited']['view'][:9]
    manifest['same_camera_scale']=manifest['assets']['cam_bound']['view'][9:12]==manifest['assets']['autoinhibited']['view'][9:12]
    crop_w=max(a['trim']['width'] for a in manifest['assets'].values());crop_h=max(a['trim']['height'] for a in manifest['assets'].values())
    for name,a in manifest['assets'].items():
        im=Image.open(OUT/a['file']).convert('RGBA');b=a['alpha_bbox'];x=round((b[0]+b[2]-crop_w)/2);y=round((b[1]+b[3]-crop_h)/2)
        crop=[x,y,x+crop_w,y+crop_h];fn=name+'_matched_scale_trim.png';im.crop(crop).save(OUT/fn)
        a['matched_scale_trim']={'file':fn,'width':crop_w,'height':crop_h,'crop_box':crop,'landmarks':{key:[p[0]-x,p[1]-y] for key,p in a['landmarks'].items()}}
    (OUT/'render_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    print('COMPLETE',flush=True)
if __name__=='__main__':main()