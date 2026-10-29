"""Exact CaM inputs in a shared, auditable molecular frame.

The complete IQCG complex is rigidly fitted to NHE1 using 24 C-alpha atoms
in CaM helices I and IV (mature residues 7-19 and 65-75). These two helices
form the more stable shared N-lobe frame. One additional common rigid
rotation defines the viewing direction. Neither lobe is transformed alone.
No atom, loop, calcium ion, or experimental conformation is reconstructed.
"""
from pathlib import Path
import os, json, hashlib, argparse
import numpy as np

OUT=Path(__file__).resolve().parent
ENV=Path('D:/Anaconda/envs/pvb_cpu')
os.environ['PATH']=str(ENV/'Library/bin')+os.pathsep+os.environ.get('PATH','')
if hasattr(os,'add_dll_directory'):
    DLL_HANDLES=[os.add_dll_directory(str(p)) for p in [ENV/'Library/bin',ENV/'DLLs'] if p.exists()]
import pymol2
from PIL import Image,ImageDraw,ImageFont

COLORS={'N_lobe':'#2166AC','C_lobe':'#85ABD1','linker':'#A6AEB4',
        'partner':'#B87936','key':'#3A8D8C','sulfur':'#B18F31',
        'oxygen':'#A83D3D','nitrogen':'#447BA7','context':'#A6AEB4'}
FIT_RES=list(range(7,20))+list(range(65,76))

def rgb(h):
    return [int(h[i:i+2],16)/255 for i in (1,3,5)]

def project(c,sel,w,h):
    xyz=c.get_coords(sel).mean(0);v=np.asarray(c.get_view())
    q=v[:9].reshape(3,3).T@(xyz-v[12:15])+v[9:12]
    halfh=-v[11]*np.tan(np.deg2rad(float(c.get('field_of_view'))/2))
    halfw=halfh*w/h
    return [float((.5+q[0]/(2*halfw))*w),float((.5-q[1]/(2*halfh))*h)]

def style(c,obj):
    c.hide('everything',obj);c.show('cartoon',obj+' and polymer.protein')
    for key,selection in [('N_lobe','chain A and resi 4-75'),
                          ('linker','chain A and resi 76-81'),
                          ('C_lobe','chain A and resi 82-148'),
                          ('partner','chain P')]:
        c.color(key,obj+' and '+selection);c.set('cartoon_color',key,obj+' and '+selection)
    for site in [72,109]:
        focus=f'{obj} and chain A and resi {site}'
        c.show('sticks',focus+' and not name N+C+O')
        c.color('key',focus+' and elem C');c.color('sulfur',focus+' and elem S')
        c.set('stick_radius',.26,focus)

def fit(a,b):
    ac=a.mean(0);bc=b.mean(0)
    u,s,v=np.linalg.svd((a-ac).T@(b-bc))
    r=u@np.diag([1.,1.,np.linalg.det(u@v)])@v
    t=bc-ac@r
    d=np.linalg.norm(a@r+t-b,axis=1)
    return r,t,float(np.sqrt(np.mean(d*d)))

def coords_for(c,obj,ids):
    return np.asarray([c.get_coords(f'{obj} and chain A and resi {i} and name CA')[0] for i in ids],float)

def identity_view(c):
    v=list(c.get_view());v[:9]=[1.,0.,0.,0.,1.,0.,0.,0.,1.];c.set_view(v)

def render_pair(c,base,basis,angle,prefix,w,h):
    theta=np.deg2rad(angle)
    around_x=np.array([[1,0,0],[0,np.cos(theta),-np.sin(theta)],[0,np.sin(theta),np.cos(theta)]],float)
    rotation=basis@around_x
    for label in base:
        c.load_coords(base[label]@rotation,label)
    c.viewport(w,h);identity_view(c);c.zoom('NHE1 or IQCG',buffer=1.6,complete=1)
    # Fit the projected union, rather than a loose 3-D bounding sphere. Both
    # states keep the same camera and scale; no molecule is enlarged alone.
    union=np.concatenate([c.get_coords(label) for label in base])
    lower=union.min(0);upper=union.max(0);center=(lower+upper)/2
    half_height=max((upper[0]-lower[0])/2/(w/h),(upper[1]-lower[1])/2)*1.075
    distance=half_height/np.tan(np.deg2rad(float(c.get('field_of_view'))/2))
    common=list(c.get_view());common[9:12]=[0.,0.,-float(distance)];common[12:15]=center.tolist()
    common[15]=float(distance-(upper[2]-center[2])-10.)
    common[16]=float(distance-(lower[2]-center[2])+10.)
    c.set_view(common)
    entries={}
    for label in ['NHE1','IQCG']:
        c.disable('all');c.enable(label);c.set_view(common)
        filename=f'{prefix}_{label}.png'
        c.png(str(OUT/filename),width=w,height=h,dpi=450,ray=1,quiet=1)
        im=Image.open(OUT/filename).convert('RGBA')
        entries[label]={'file':filename,'width':w,'height':h,'view':list(c.get_view()),
                        'alpha_bbox':list(im.getchannel('A').getbbox()),
                        'landmarks':{f'M{i}':project(c,f'{label} and chain A and resi {i} and name CE',w,h) for i in [72,109]}}
    boxes=[v['alpha_bbox'] for v in entries.values()]
    # Identical crop for both states preserves angstrom-to-pixel scale and position.
    pad=max(12,round(w*.015))
    bb=[max(0,min(b[0] for b in boxes)-pad),max(0,min(b[1] for b in boxes)-pad),
        min(w,max(b[2] for b in boxes)+pad),min(h,max(b[3] for b in boxes)+pad)]
    for label,entry in entries.items():
        filename=f'{prefix}_{label}_trim.png'
        Image.open(OUT/entry['file']).crop(bb).save(OUT/filename)
        entry['paired_trim']={'file':filename,'crop_box':bb,'width':bb[2]-bb[0],'height':bb[3]-bb[1],
          'landmarks':{k:[p[0]-bb[0],p[1]-bb[1]] for k,p in entry['landmarks'].items()}}
    print('RENDERED',prefix,flush=True)
    return entries,rotation

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--preview',action='store_true');parser.add_argument('--angle',type=float,default=240.)
    args=parser.parse_args()
    manifest={'scope':'Exact model/simulation-input atomic snapshots, not MD averages. No rebuilt atoms or calcium.',
              'palette':COLORS,'fit_residues_mature':FIT_RES,'fit_description':'Complete IQCG complex rigidly fitted to NHE1 CaM helices I and IV C-alpha atoms (7-19,65-75).'}
    with pymol2.PyMOL() as pm:
        c=pm.cmd;c.bg_color('white')
        for k,v in dict(ray_opaque_background=0,opaque_background=0,ray_trace_mode=0,ray_shadows=0,
                        antialias=3,orthoscopic=1,depth_cue=0,specular=.12,shininess=35,
                        ambient=.65,direct=.35,reflect=.12,light_count=2,
                        cartoon_fancy_helices=1,cartoon_fancy_sheets=1,cartoon_flat_sheets=1,
                        cartoon_smooth_loops=1,cartoon_sampling=20,cartoon_loop_radius=.24,
                        cartoon_tube_radius=.24,cartoon_oval_length=1.25,cartoon_oval_width=.28,
                        max_threads=4,stick_radius=.17,valence=0).items():c.set(k,v)
        for k,col in COLORS.items():c.set_color(k,rgb(col))
        for label in ['NHE1','IQCG']:
            c.load(str(OUT/'inputs'/f'{label}_mature_CaM_and_partner.pdb'),label)
            c.remove(label+' and hydro');c.dss(label);style(c,label)
        original={label:np.array(c.get_coords(label),float) for label in ['NHE1','IQCG']}
        a=coords_for(c,'IQCG',FIT_RES);b=coords_for(c,'NHE1',FIT_RES)
        r,t,rms=fit(a,b);base={'NHE1':original['NHE1'],'IQCG':original['IQCG']@r+t}
        n=coords_for(c,'NHE1',list(range(4,76))).mean(0)
        cc=coords_for(c,'NHE1',list(range(82,149))).mean(0)
        x=(cc-n)/np.linalg.norm(cc-n)
        partner=c.get_coords('NHE1 and chain P and name CA')
        residual=(partner-partner.mean(0))-(partner-partner.mean(0))@x[:,None]*x
        _,_,vh=np.linalg.svd(residual,full_matrices=False)
        y=vh[0];y-=np.dot(y,x)*x;y/=np.linalg.norm(y)
        z=np.cross(x,y);basis=np.stack([x,y,z],axis=1)
        center=n
        base={label:xyz-center for label,xyz in base.items()}
        manifest['fit_RMSD_A']=rms;manifest['fit_rotation_row_vector']=r.tolist();manifest['fit_translation_A']=t.tolist()
        manifest['input_sha256']={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (OUT/'inputs').glob('*.pdb')}
        manifest['heavy_atom_count']={label:len(xyz) for label,xyz in original.items()}
        manifest['coordinate_audit']={}
        for label in original:
            # All pairwise distances across each complete complex are invariant.
            # A dense stride gives broad coverage without a large quadratic array.
            s0=original[label][::7];s1=base[label][::7]
            d0=np.linalg.norm(s0[:,None]-s0[None,:],axis=2)
            d1=np.linalg.norm(s1[:,None]-s1[None,:],axis=2)
            manifest['coordinate_audit'][label]={'max_pairwise_distance_error_A':float(np.max(np.abs(d0-d1))),
                                                 'complete_complex_rigid_transform':True}
            assert np.max(np.abs(d0-d1))<1e-6
        if args.preview:
            previews={}
            for angle in [0,60,120,180,240,300]:
                entries,_=render_pair(c,base,basis,angle,f'preview_{angle}',900,940)
                previews[str(angle)]=entries
            canvas=Image.new('RGB',(1200,3*610),'white');draw=ImageDraw.Draw(canvas)
            try:font=ImageFont.truetype('C:/Windows/Fonts/arial.ttf',24)
            except OSError:font=ImageFont.load_default()
            for j,angle in enumerate([0,60,120,180,240,300]):
                x0=(j%2)*600;y0=(j//2)*610
                draw.text((x0+8,y0+5),f'Common view: {angle} deg',fill='black',font=font)
                for k,label in enumerate(['NHE1','IQCG']):
                    im=Image.open(OUT/previews[str(angle)][label]['paired_trim']['file']).convert('RGBA')
                    im.thumbnail((290,520));xx=x0+k*300+(300-im.width)//2;yy=y0+65
                    canvas.paste(im,(xx,yy),im)
                    draw.text((x0+k*300+65,y0+35),label,fill='black',font=font)
            canvas.save(OUT/'shared_view_contact_sheet.png')
            manifest['previews']=previews
            (OUT/'preview_manifest.json').write_text(json.dumps(manifest,indent=2))
        else:
            entries,rotation=render_pair(c,base,basis,args.angle,'whole',2600,2718)
            manifest['assets']={'whole_'+label:entry for label,entry in entries.items()}
            manifest['view_angle_degrees']=args.angle;manifest['shared_display_rotation_row_vector']=rotation.tolist()
            manifest['common_reference_origin_A']=center.tolist()
            manifest['projection']='Identical orthographic camera, scale, canvas and union crop; CaM helix-I/IV reference coordinates coincide.'
            for label in original:
                observed=np.asarray(c.get_coords(label),float)
                expected=base[label]@rotation
                error=float(np.max(np.linalg.norm(observed-expected,axis=1)))
                manifest['coordinate_audit'][label]['display_coordinate_max_rounding_error_A']=error
                manifest['coordinate_audit'][label]['heavy_atom_count_preserved']=(len(observed)==len(original[label]))
                assert error<1e-4 and len(observed)==len(original[label])
            v=entries['NHE1']['view']
            half_height=-v[11]*np.tan(np.deg2rad(float(c.get('field_of_view'))/2))
            manifest['shared_pixels_per_A']=2718/(2*half_height)
            assert np.array_equal(entries['NHE1']['view'],entries['IQCG']['view'])
            assert entries['NHE1']['paired_trim']['crop_box']==entries['IQCG']['paired_trim']['crop_box']
            c.disable('all');c.enable('NHE1');c.scene('NHE1_shared_view','store')
            c.disable('all');c.enable('IQCG');c.scene('IQCG_shared_view','store')
            c.disable('all');c.enable('NHE1')
            c.save(str(OUT/'CaM_shared_reference_v82.pse'))
            manifest['pymol_version']=str(c.get_version()[0])
            (OUT/'render_manifest.json').write_text(json.dumps(manifest,indent=2))
    print('COMPLETE',str(OUT),flush=True)

if __name__=='__main__':main()
