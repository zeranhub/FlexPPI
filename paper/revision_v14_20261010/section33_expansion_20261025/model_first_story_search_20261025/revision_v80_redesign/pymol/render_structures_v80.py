"""Continuous-cartoon PyMOL assets from the exact model/simulation inputs.

No atom or residue is rebuilt. Rigid-body fitting changes object display only.
Local views retain each complete chain; no spatial selection is used for cartoons.
Scientific image rendering: inherit the actual previous PyMOL asset, not a bitmap.
"""
from pathlib import Path
import os, json, hashlib, shutil
import numpy as np

ENV = Path('D:/Anaconda/envs/pvb_cpu')
os.environ['PATH'] = str(ENV/'Library/bin') + os.pathsep + os.environ.get('PATH','')
if hasattr(os, 'add_dll_directory'):
    DLL_HANDLES = [os.add_dll_directory(str(p)) for p in [ENV/'Library/bin',ENV/'DLLs'] if p.exists()]
import pymol2
from PIL import Image

OUT = Path(__file__).resolve().parent
OLD = OUT.parent.parent/'figure_cam_model_story'
COLORS = {'N_lobe':'#2166AC', 'C_lobe':'#85ABD1', 'linker':'#A6AEB4',
          'partner':'#B87936', 'key':'#3A8D8C', 'sulfur':'#B18F31',
          'oxygen':'#A83D3D', 'nitrogen':'#447BA7', 'context':'#A6AEB4'}

def rgb(h):
    return [int(h[i:i+2],16)/255 for i in (1,3,5)]

def project(c, sel, w, h):
    xyz = c.get_coords(sel).mean(0)
    v = np.asarray(c.get_view())
    q = v[:9].reshape(3,3).T @ (xyz-v[12:15]) + v[9:12]
    halfh = -v[11]*np.tan(np.deg2rad(float(c.get('field_of_view'))/2))
    halfw = halfh*w/h
    return [float((.5+q[0]/(2*halfw))*w),float((.5-q[1]/(2*halfh))*h)]

def style(c, obj, local_site=None):
    c.hide('everything',obj)
    c.show('cartoon',obj+' and polymer.protein')
    # Representation-specific cartoon colours prevent an emphasized side chain
    # from creating a misleading green band in the helix ribbon itself.
    for key, selection in [('N_lobe','chain A and resi 4-75'),
                           ('linker','chain A and resi 76-81'),
                           ('C_lobe','chain A and resi 82-148'),
                           ('partner','chain P')]:
        c.color(key,obj+' and '+selection)
        c.set('cartoon_color',key,obj+' and '+selection)
    sites = [72,109] if local_site is None else [local_site]
    for site in sites:
        focus = f'{obj} and chain A and resi {site}'
        c.show('sticks',focus+' and not name N+C+O')
        c.color('key',focus+' and elem C')
        c.color('sulfur',focus+' and elem S')
        c.set('stick_radius', .26 if local_site is None else .21,focus)
    if local_site is not None:
        focus = f'{obj} and chain A and resi {local_site}'
        # Neighbour selection applies ONLY to sticks. All peptide cartoons stay
        # continuous, including the entire partner chain and the second CaM lobe.
        neighbours = (f'{obj} and byres (({obj} and not ({focus})) '
                      f'within 4.3 of ({focus} and not name N+CA+C+O))')
        c.show('sticks',neighbours+' and not name N+C+O')
        c.set('stick_radius',.11,neighbours)
        c.set('stick_transparency',.12,neighbours)
        # Clear local lobe, complete faded opposite lobe, no sphere-cut cartoons.
        opposite = '82-148' if local_site < 76 else '4-75'
        c.set('cartoon_transparency',.72,obj+' and chain A and resi '+opposite)
        c.set('cartoon_transparency',.4,obj+' and chain A and resi 76-81')
        c.set('cartoon_transparency',.12,obj+' and chain P')
    c.set('cartoon_transparency',0,obj+' and chain A and resi '+('4-75' if local_site==72 else '82-148')) if local_site else None

def render(c,name,w,h,landmarks):
    c.png(str(OUT/(name+'.png')),width=w,height=h,dpi=450,ray=1,quiet=1)
    im=Image.open(OUT/(name+'.png')).convert('RGBA')
    alpha=im.getchannel('A')
    bbox=alpha.getbbox()
    entry={'file':name+'.png','width':w,'height':h,'view':list(c.get_view()),
           'field_of_view':float(c.get('field_of_view')), 'alpha_bbox':list(bbox),
           'landmarks':{k:project(c,s,w,h) for k,s in landmarks.items()}}
    c.scene(name,'store')
    print('RENDERED',name,flush=True)
    return entry

def paired_trim(assets,names):
    boxes=[assets[n]['alpha_bbox'] for n in names]
    w,h=assets[names[0]]['width'],assets[names[0]]['height']
    bb=[max(0,min(b[0] for b in boxes)-35),max(0,min(b[1] for b in boxes)-35),
        min(w,max(b[2] for b in boxes)+35),min(h,max(b[3] for b in boxes)+35)]
    for n in names:
        im=Image.open(OUT/(n+'.png'));name=n+'_trim.png';im.crop(bb).save(OUT/name)
        assets[n]['paired_trim']={'file':name,'crop_box':bb,'width':bb[2]-bb[0],'height':bb[3]-bb[1],
           'landmarks':{k:[v[0]-bb[0],v[1]-bb[1]] for k,v in assets[n]['landmarks'].items()}}

def main():
    OUT.mkdir(parents=True,exist_ok=True);(OUT/'inputs').mkdir(exist_ok=True)
    manifest={'palette':COLORS,'assets':{},'scope':'Exact simulation/model-input atomic structures; no rebuilt atoms; no calcium added; static snapshots are not MD conformations.'}
    for label in ['NHE1','IQCG']:
        name=label+'_mature_CaM_and_partner.pdb'
        if not (OUT/'inputs'/name).exists(): shutil.copy2(OLD/'inputs'/name,OUT/'inputs'/name)
    with pymol2.PyMOL() as pm:
        c=pm.cmd;c.bg_color('white')
        for k,v in dict(ray_opaque_background=0,opaque_background=0,ray_trace_mode=0,ray_shadows=0,
                        antialias=3,orthoscopic=1,depth_cue=0,specular=.12,shininess=35,
                        ambient=.65,direct=.35,reflect=.12,light_count=2,
                        cartoon_fancy_helices=1,cartoon_fancy_sheets=1,cartoon_flat_sheets=1,
                        cartoon_smooth_loops=1,cartoon_sampling=20,cartoon_loop_radius=.24,
                        cartoon_tube_radius=.24,cartoon_oval_length=1.25,cartoon_oval_width=.28,
                        max_threads=4,stick_radius=.17,valence=0).items():
            c.set(k,v)
        for n,col in COLORS.items():c.set_color(n,rgb(col))
        for label in ['NHE1','IQCG']:
            path=OUT/'inputs'/(label+'_mature_CaM_and_partner.pdb')
            c.load(str(path),label);c.remove(label+' and hydro');c.dss(label);style(c,label)
        ali=c.align('IQCG and chain A and resi 4-75 and name CA','NHE1 and chain A and resi 4-75 and name CA',cycles=0)
        manifest['global_N_lobe_alignment_RMSD_A']=float(ali[0])
        c.viewport(2200,2000)
        c.orient('NHE1');c.turn('x',20);c.turn('y',20);c.turn('z',-16)
        c.zoom('NHE1 or IQCG',buffer=1.6,complete=1)
        global_view=c.get_view()
        for label in ['NHE1','IQCG']:
            c.disable('all');c.enable(label);c.set_view(global_view)
            marks={f'M{k}':f'{label} and chain A and resi {k} and name CE' for k in [72,109]}
            marks.update(N_lobe=label+' and chain A and resi 4-75 and name CA',
                         C_lobe=label+' and chain A and resi 82-148 and name CA',
                         partner=label+' and chain P and name CA')
            manifest['assets']['whole_'+label]=render(c,'whole_'+label,2200,2000,marks)
        paired_trim(manifest['assets'],['whole_NHE1','whole_IQCG'])
        # Same complete-chain representations in local comparison. Both states
        # share the local backbone fit and camera; every atom remains available.
        for site, domain in [(72,'4-81'),(109,'76-148')]:
            objs={}
            for label in ['NHE1','IQCG']:
                obj=label+'_M'+str(site);objs[label]=obj;c.create(obj,label);style(c,obj,site)
            left,right=objs['NHE1'],objs['IQCG']
            ali=c.align(f'{right} and chain A and resi {site-5}-{site+5} and name N+CA+C',
                        f'{left} and chain A and resi {site-5}-{site+5} and name N+CA+C',cycles=0)
            manifest[f'local_M{site}_backbone_fit_RMSD_A']=float(ali[0])
            c.viewport(2200,1800)
            c.orient(f'{left} and chain A and resi {domain}')
            c.turn('y',35);c.turn('x',18);c.turn('z',-12)
            # Fit COMPLETE chains in the viewport; local differences are read
            # from opacity plus sticks, never from disconnected ribbon pieces.
            c.zoom(f'{left} or {right}',buffer=1.1,complete=1)
            localview=c.get_view()
            for label,obj in objs.items():
                c.disable('all');c.enable(obj);c.set_view(localview)
                focus=f'{obj} and chain A and resi {site}'
                marks={'M'+str(site):focus+' and name CE','M'+str(site)+'_CA':focus+' and name CA',
                       'partner':obj+' and chain P and name CA'}
                name=f'M{site}_{label}'
                manifest['assets'][name]=render(c,name,2200,1800,marks)
            paired_trim(manifest['assets'],[f'M{site}_NHE1',f'M{site}_IQCG'])
        c.disable('all');c.enable('NHE1');c.set_view(global_view)
        c.save(str(OUT/'CaM_continuous_context_v80.pse'))
        manifest['pymol_version']=str(c.get_version()[0])
        manifest['input_sha256']={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (OUT/'inputs').glob('*.pdb')}
        manifest['representation_audit']={'cartoons':'complete polymer chains, no spatial cut',
            'neighbour_sticks':'within 4.3 A of focused side-chain heavy atoms',
            'hydrogens':'hidden/removed only in display objects',
            'coordinates':'input PDB coordinates unmodified; rigid fitting only on PyMOL object copies',
            'local_camera':'same view and scale for each matched pair; whole chains retained in viewport',
            'ribbon_colour':'representation specific; key-residue sticks do not recolour ribbon bands'}
    (OUT/'render_manifest.json').write_text(json.dumps(manifest,indent=2))
    print('COMPLETE',str(OUT),flush=True)

if __name__=='__main__':main()
