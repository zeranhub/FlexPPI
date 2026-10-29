"""Ray-trace deposited structures; no invented conformers or motion arrows.

Clean, sequence-mapped first-model structures are prepared separately. Structure
snapshots locate the NMR probes; their conformations do not encode dynamics.
"""
from pathlib import Path
import json
import os
import sys
import numpy as np

ENV = Path('D:/Anaconda/envs/pvb_cpu')
os.environ['PATH'] = str(ENV/'Library/bin')+os.pathsep+os.environ.get('PATH', '')
if hasattr(os, 'add_dll_directory'):
    DLL_HANDLES = [os.add_dll_directory(str(p)) for p in [ENV/'Library/bin', ENV/'DLLs'] if p.exists()]
import pymol2

ROOT = Path(__file__).resolve().parent
OUT = ROOT/'structure_renders'
SPECS = [('CaMKKalpha', '1CKK', '#762A83'), ('eNOS', '1NIW', '#328C8D'),
         ('nNOS', '2O60', '#2166AC'), ('smMLCK', '1CDL', '#B87936')]
KEYS = [18, 39, 105, 109, 144]


def rgb(h):
    return [int(h[i:i+2], 16)/255 for i in [1, 3, 5]]


def setup(cmd):
    cmd.bg_color('white')
    for k, v in dict(ray_opaque_background=0, ray_trace_mode=0, ray_trace_gain=.025,
                     ray_shadows=0, antialias=2, orthoscopic=1, depth_cue=0,
                     specular=.16, ambient=.58, direct=.42, cartoon_fancy_helices=1,
                     cartoon_fancy_sheets=1, cartoon_flat_sheets=1,
                     cartoon_smooth_loops=1, cartoon_sampling=16, max_threads=4,
                     stick_radius=.18, sphere_scale=.22, valence=0).items():
        cmd.set(k, v)
    cmd.set_color('cam_gray', rgb('#CBD0D4'))
    cmd.set_color('probe_carbon', rgb('#AA2838'))
    cmd.set_color('calcium_gray', rgb('#839096'))


def style(cmd, obj, color):
    cmd.set_color('partner_'+obj, rgb(color))
    cmd.hide('everything', obj)
    cmd.show('cartoon', obj+' and polymer.protein')
    cmd.color('cam_gray', obj+' and chain A')
    cmd.color('partner_'+obj, obj+' and chain P')
    cmd.set('cartoon_color', 'cam_gray', obj+' and chain A')
    cmd.set('cartoon_color', 'partner_'+obj, obj+' and chain P')
    cmd.show('spheres', obj+' and inorganic')
    cmd.color('calcium_gray', obj+' and inorganic')
    selected = obj+' and chain A and resi '+'+'.join(map(str, KEYS))
    cmd.show('sticks', selected+' and not name N+C+O')
    cmd.color('probe_carbon', selected+' and elem C')
    cmd.color('blue', selected+' and elem N')
    cmd.color('red', selected+' and elem O')
    cmd.color('yellow', selected+' and elem S')
    cmd.show('spheres', selected+' and (name CD1+CD2+CE)')


def project(cmd, selection, width, height):
    pts = cmd.get_coords(selection, state=1)
    xyz = pts.mean(axis=0)
    view = np.asarray(cmd.get_view())
    q = view[:9].reshape(3, 3).T @ (xyz-view[12:15])+view[9:12]
    halfh = -view[11]*np.tan(np.deg2rad(float(cmd.get('field_of_view'))/2))
    halfw = halfh*width/height
    return [float((.5+q[0]/(2*halfw))*width), float((.5-q[1]/(2*halfh))*height)]


def render(cmd, filename, width, height, selections):
    cmd.png(str(OUT/filename), width=width, height=height, dpi=300, ray=1, quiet=1)
    return dict(file=filename, width=width, height=height, view=list(cmd.get_view()),
                landmarks={key:project(cmd, sel, width, height) for key, sel in selections.items()})


def main():
    OUT.mkdir(exist_ok=True)
    manifest = dict(description='Deposited first models; mapped mature-CaM numbering; no dynamic conformer overlay',
                    common_view=True, assets={}, pymol_version=None)
    with pymol2.PyMOL() as pm:
        c = pm.cmd
        setup(c)
        for label, pid, color in SPECS:
            c.load(str(ROOT/f'structure_data/{pid}_CaM_A_peptide_P_clean.pdb'), label)
            c.dss(label)
            if pid == '1NIW':
                assert c.count_atoms(label+' and resn MSE and polymer.protein') == c.count_atoms(label+' and resn MSE')
        for label, pid, color in SPECS:
            if label != 'nNOS':
                align = c.align(f'{label} and chain A and name CA', 'nNOS and chain A and name CA', cycles=0, transform=1)
                manifest.setdefault('whole_CaM_alignment_rmsd_A', {})[label] = float(align[0])
            style(c, label, color)
        c.orient('nNOS and chain A')
        c.turn('x', 12)
        c.turn('y', 25)
        c.turn('z', -8)
        c.zoom('all and polymer.protein', buffer=3.5, complete=1)
        overview = c.get_view()
        for label, pid, color in SPECS:
            c.disable('all'); c.enable(label); c.set_view(overview)
            sels = {f'{"L" if k in [18,39,105] else "M"}{k}':f'{label} and chain A and resi {k} and name CA' for k in KEYS}
            sels.update(N_lobe=f'{label} and chain A and resi 10-65 and name CA',
                        C_lobe=f'{label} and chain A and resi 90-140 and name CA',
                        peptide=f'{label} and chain P and name CA')
            manifest['assets']['overview_'+label] = render(c, 'overview_'+label+'.png', 1400, 1300, sels)
            c.scene('whole_'+label, 'store')
            print('rendered overview', label, flush=True)
        # One genuine local environment for each of the three cleanest probes.
        for site in [39, 105, 144]:
            name = ('L' if site != 144 else 'M')+str(site)
            # Locally fit the CaM backbone, so each local pair is viewed in a
            # shared local frame rather than confusing lobe closure with rotamers.
            c.align(f'smMLCK and chain A and resi {site-6}-{site+6} and name N+CA+C',
                    f'nNOS and chain A and resi {site-6}-{site+6} and name N+CA+C',
                    cycles=0, transform=1)
            c.disable('all');c.enable('nNOS')
            c.hide('spheres', 'all')
            focus = f'nNOS and chain A and resi {site}'
            c.orient(focus+' or (nNOS and chain P within 8 of ('+focus+'))')
            c.turn('y', 25 if site == 39 else -25)
            c.turn('x', -10)
            c.zoom(focus, buffer=7.0, complete=1)
            view = c.get_view()
            for label in ['nNOS', 'smMLCK']:
                c.disable('all');c.enable(label)
                c.hide('everything', label)
                c.show('cartoon', label+' and polymer.protein')
                c.set('cartoon_transparency', .58, label+' and chain A')
                c.show('sticks', f'{label} and chain A and resi {site} and not name N+C+O')
                c.show('sticks', f'{label} and chain P and byres (({label} and chain A and resi {site}) around 4.5)')
                methyl = 'CD2' if site == 39 else 'CD1+CD2' if site == 105 else 'CE'
                c.show('spheres', f'{label} and chain A and resi {site} and name {methyl}')
                c.set('sphere_scale', .32, f'{label} and chain A and resi {site}')
                c.set_view(view)
                sels = {name:f'{label} and chain A and resi {site} and name CA',
                        'methyl':f'{label} and chain A and resi {site} and name '+('CD2' if site == 39 else 'CD1' if site == 105 else 'CE')}
                manifest['assets'][f'local_{name}_{label}'] = render(c, f'local_{name}_{label}.png', 1000, 1100, sels)
                c.scene(f'local_{name}_{label}', 'store')
                print('rendered local', name, label, flush=True)
        c.save(str(OUT/'cam_partner_dynamics.pse'))
        manifest['pymol_version'] = str(c.get_version()[0])
    (OUT/'render_manifest.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')


if __name__ == '__main__':
    main()
