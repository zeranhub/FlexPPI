"""Reproducible 1A06 M316/L45/S99 local rendering; inherited v86 camera rotation.
No modelling, coordinate fitting, atom generation or coordinate transformation.
"""
from pathlib import Path
import os, json, hashlib
import numpy as np
ENV = Path('D:/Anaconda/envs/pvb_cpu')
os.environ['PATH'] = str(ENV/'Library/bin') + os.pathsep + os.environ.get('PATH','')
DLL_HANDLES = [os.add_dll_directory(str(p)) for p in [ENV/'Library/bin', ENV/'DLLs'] if p.exists()]
import pymol2
from PIL import Image
OUT = Path(__file__).resolve().parent
INPUT = OUT/'inputs/1A06_rigidly_aligned_to_1MXE.pdb'
BASE = json.loads((OUT/'inputs/v86_render_manifest.json').read_text(encoding='utf-8'))
COLORS = BASE['palette']

def project(cmd, selection, width, height):
    xyz = cmd.get_coords(selection).mean(0)
    v = np.asarray(cmd.get_view())
    q = v[:9].reshape(3,3).T @ (xyz-v[12:15]) + v[9:12]
    hh = -v[11]*np.tan(np.deg2rad(float(cmd.get('field_of_view'))/2))
    hw = hh*width/height
    return [float((.5+q[0]/(2*hw))*width), float((.5-q[1]/(2*hh))*height)]

def main():
    manifest = {'input': str(INPUT.relative_to(OUT)), 'input_sha256': hashlib.sha256(INPUT.read_bytes()).hexdigest(),
                'provenance': '1A06 experimentally observed atoms, rigidly aligned in v86 to 1MXE; unchanged here',
                'limitations': ['1A06 is a static structural context, with no model prediction or MD.',
                                'Residues317-318 have no observed atoms and are not generated.',
                                'Dashed segments denote minimum heavy-atom separations, not assigned hydrogen bonds.'],
                'palette': COLORS, 'assets': {}}
    with pymol2.PyMOL() as pm:
        c = pm.cmd
        c.load(str(INPUT), 'autoinhibited')
        xyz_before = c.get_coords('autoinhibited').copy()
        for key,val in dict(ray_opaque_background=0, opaque_background=0, ray_trace_mode=0, ray_shadows=0,
                           antialias=3, orthoscopic=1, depth_cue=0, specular=.12, shininess=35, ambient=.65,
                           direct=.35, reflect=.12, light_count=2, cartoon_fancy_helices=1, cartoon_fancy_sheets=1,
                           cartoon_flat_sheets=1, cartoon_smooth_loops=1, cartoon_sampling=20,
                           cartoon_loop_radius=.24, cartoon_tube_radius=.24, cartoon_oval_length=1.25,
                           cartoon_oval_width=.28, max_threads=4, stick_radius=.24, valence=0).items():
            c.set(key,val)
        for key,h in COLORS.items(): c.set_color(key,[int(h[k:k+2],16)/255 for k in [1,3,5]])
        c.bg_color('white'); c.hide('everything'); c.dss('autoinhibited')
        reg = 'autoinhibited and chain A and resi 294-316'
        local_reg = 'autoinhibited and chain A and resi 312-316'
        local_kinase = 'autoinhibited and chain A and resi 42-48+97-100'
        c.show('cartoon',local_reg)
        c.color('background',local_kinase); c.set('cartoon_color','background',local_kinase)
        c.color('regulatory',local_reg); c.set('cartoon_color','regulatory',local_reg)
        c.set('cartoon_transparency',.25,local_kinase)
        c.show('sticks','autoinhibited and chain A and resi 45+99+316 and not name N+C+O')
        c.color('background','autoinhibited and chain A and resi 45+99')
        c.color('key','autoinhibited and chain A and resi 316')
        c.set('stick_radius',.29,'autoinhibited and chain A and resi 316')
        # Distance guides are actual closest heavy-atom pairs, without an interaction assignment.
        pairs = [(45,'SD','CD1'), (99,'CE','OG')]
        facts = []
        for residue,anchor_atom,partner_atom in pairs:
            left = f'autoinhibited and chain A and resi 316 and name {anchor_atom}'
            right = f'autoinhibited and chain A and resi {residue} and name {partner_atom}'
            dist = c.distance(f'distance_{residue}',left,right,mode=0)
            c.hide('labels',f'distance_{residue}')
            c.set('dash_color','background',f'distance_{residue}'); c.set('dash_width',4.4,f'distance_{residue}')
            c.set('dash_gap',.18,f'distance_{residue}'); c.set('dash_length',.2,f'distance_{residue}')
            facts.append({'M316_atom':anchor_atom,'partner_residue':residue,'partner_atom':partner_atom,'distance_A':float(dist)})
        width,height = 2400,2100
        c.viewport(width,height)
        # Zoom changes camera translation, depth and clipping only. Preserve all nine rotation elements.
        c.zoom(local_reg+' or (autoinhibited and chain A and resi 45+99 and not name N+C+O)',buffer=2,complete=1)
        view = list(c.get_view()); view[:9] = BASE['camera_common_rotation']; view[15] = 1; view[16] = -view[11]+180; c.set_view(view)
        landmarks = {key: project(c,sel,width,height) for key,sel in {
            'M316':'autoinhibited and chain A and resi 316 and name CE',
            'M316_CA':'autoinhibited and chain A and resi 316 and name CA',
            'L45':'autoinhibited and chain A and resi 45 and name CD1',
            'S99':'autoinhibited and chain A and resi 99 and name OG',
            'M316_SD':'autoinhibited and chain A and resi 316 and name SD',
        }.items()}
        c.png(str(OUT/'M316_contacts_local.png'),width=width,height=height,dpi=600,ray=1,quiet=1)
        c.save(str(OUT/'M316_contacts_local.pse'))
        im = Image.open(OUT/'M316_contacts_local.png').convert('RGBA')
        bbox = list(im.getchannel('A').getbbox()); pad=32
        crop = [max(0,bbox[0]-pad),max(0,bbox[1]-pad),min(width,bbox[2]+pad),min(height,bbox[3]+pad)]
        im.crop(crop).save(OUT/'M316_contacts_local_trim.png')
        manifest['assets']['local'] = {'file':'M316_contacts_local.png','width':width,'height':height,'view':list(c.get_view()),'alpha_bbox':bbox,'landmarks':landmarks,
            'trim': {'file':'M316_contacts_local_trim.png','crop_box':crop,'width':crop[2]-crop[0],'height':crop[3]-crop[1],
                     'landmarks':{key:[xy[0]-crop[0],xy[1]-crop[1]] for key,xy in landmarks.items()}}}
        manifest['contacts'] = facts
        manifest['geometry_audit'] = {'input_atom_count':len(xyz_before),'output_atom_count':c.count_atoms('autoinhibited'),
            'max_coordinate_change_A':float(np.abs(xyz_before-c.get_coords('autoinhibited')).max()),
            'coordinates_exactly_unchanged':bool(np.array_equal(xyz_before,c.get_coords('autoinhibited'))),
            'view_rotation_inherited_exactly':bool(np.array_equal(np.asarray(c.get_view()[:9]),np.asarray(BASE['camera_common_rotation']))),
            'no_generated_atoms':c.count_atoms('autoinhibited')==len(xyz_before)}
        manifest['pymol_version']=str(c.get_version()[0])
    (OUT/'render_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    print(json.dumps(manifest,indent=2))
if __name__=='__main__': main()