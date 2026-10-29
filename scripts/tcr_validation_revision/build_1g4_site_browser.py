"""Interactive, self-contained PyMOL site browser; not a dynamics illustration.

All protein atoms are unmodified deposited 2BNR coordinates. Annotation CGOs
and pseudoatoms are separate objects. Colors identify regions, not performance.
Run with D:/Anaconda/envs/pvb_cpu/python.exe.
"""
from pathlib import Path
import csv
import hashlib
import json
import os

import numpy as np

ENV = Path('D:/Anaconda/envs/pvb_cpu')
os.environ['PATH'] = str(ENV / 'Library/bin') + os.pathsep + os.environ.get('PATH', '')
if hasattr(os, 'add_dll_directory'):
    DLL_HANDLES = [os.add_dll_directory(str(p)) for p in [ENV / 'Library/bin', ENV / 'DLLs'] if p.exists()]
import pymol2
from pymol.cgo import BEGIN, END, LINES, COLOR, VERTEX, LINEWIDTH

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / 'artifacts/tcr_pmhc_reassessment_20260915/sources/2BNR.cif'
OUT = ROOT / 'output/tcr_prediction_reference_revision_20260915/pymol_1g4_sites_20260916'
OBJECTS = {'D': 'TCR_alpha_D', 'E': 'TCR_beta_E', 'A': 'HLA_A2_A', 'B': 'beta2m_B', 'C': 'NYESO1_peptide_C'}
SITES = [
    dict(name='site_01_alpha_96_99', chain='D', ids=[96, 97, 98, 99], sequence='SGGS', color='#D9792F', title='alpha 96-99', subtitle='D:96-99 | SGGS', scene='04_alpha_96_99', side=-1, dy=6),
    dict(name='site_02_alpha_198_202', chain='D', ids=list(range(198, 203)), sequence='EDTFF', color='#B53656', title='alpha 198-202', subtitle='D:198-202 | EDTFF', scene='05_alpha_198_202', side=-1, dy=-4),
    dict(name='site_03_beta_230_233', chain='E', ids=list(range(230, 234)), sequence='TQIV', color='#3E78B4', title='beta 230-233', subtitle='E:230-233 | TQIV', scene='06_beta_230_233', side=1, dy=-8),
    dict(name='site_04_beta_29_30', chain='E', ids=[29, 30], sequence='EY', color='#408C42', title='beta 29-30', subtitle='E:29-30 | EY', scene='07_beta_29_30', side=1, dy=12),
    dict(name='site_05_alpha_51', chain='D', ids=[51], sequence='Q', color='#8056A2', title='alpha Q51', subtitle='D:51', scene='08_alpha_Q51', side=-1, dy=15),
    dict(name='site_06_beta_100_101', chain='E', ids=[100, 101], sequence='EL', color='#8C693D', title='beta 100-101', subtitle='E:100-101 | EL', scene='09_beta_100_101', side=1, dy=-2),
    dict(name='site_07_alpha_Y100', chain='D', ids=[100], sequence='Y', color='#008A91', title='alpha Y100', subtitle='D:100 | Y', scene='10_alpha_Y100', side=-1, dy=-5),
]
AA = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}


def rgb(hex_color):
    return [int(hex_color[i:i+2], 16) / 255 for i in (1, 3, 5)]


def camera(cmd, selection, basis, pad=1.15, width=1600, height=1250):
    xyz = cmd.get_coords(selection)
    projected = xyz @ basis.T
    center_world = (projected.min(0) + projected.max(0)) / 2 @ basis
    span = np.ptp(projected, axis=0)
    half_height = max(span[1] / 2, span[0] * height / width / 2) * pad + 2
    distance = half_height / np.tan(np.deg2rad(10))
    clip = max(span[2] / 2 + 40, 40)
    cmd.viewport(width, height)
    cmd.set_view([*basis.T.reshape(-1), 0, 0, -distance, *center_world, distance-clip, distance+clip, 20])


def style(cmd):
    cmd.bg_color('white')
    for key, value in dict(orthoscopic=1, depth_cue=0, ray_shadows=0, antialias=2,
            ray_opaque_background=1, ray_trace_mode=0, ray_trace_gain=.03,
            ambient=.55, direct=.45, specular=.14, max_threads=2,
            cartoon_fancy_helices=1, cartoon_fancy_sheets=1, cartoon_flat_sheets=1,
            cartoon_smooth_loops=1, cartoon_loop_radius=.24, stick_radius=.19,
            sphere_scale=.35, label_font_id=7, label_size=18, label_color='black',
            label_outline_color='white', label_bg_color='white', label_bg_transparency=.12,
            label_padding=(.18,.18,0), seq_view=0, scene_buttons=1).items():
        cmd.set(key, value)


def base_display(cmd):
    cmd.enable('TCR'); cmd.enable('Partners')
    cmd.disable('Overview_labels'); cmd.disable('Overview_guides')
    cmd.disable('Detail_labels_*'); cmd.disable('Detail_guides_*')
    cmd.hide('everything', 'protein_atoms')
    cmd.show('cartoon', 'protein_atoms')
    cmd.set('cartoon_transparency', 0, 'protein_atoms')
    cmd.set('cartoon_transparency', .38, 'HLA_A2_A or beta2m_B')
    cmd.hide('cartoon', 'NYESO1_peptide_C')
    cmd.show('sticks', 'NYESO1_peptide_C or marked_sites')
    cmd.show('spheres', 'marked_sites and name CA')
    cmd.hide('labels', 'protein_atoms')
    cmd.deselect()


def detail_labels(cmd, name, items, context, basis):
    """Camera-aware, nonoverlapping labels, independent of protein atoms."""
    pcoords=cmd.get_coords(context) @ basis.T
    assignments=[]
    for i,(selection,label,color) in enumerate(items):
        center=cmd.get_coords(selection).mean(0)
        point=center @ basis.T
        side=-1 if i%2==0 else 1
        assignments.append(dict(center=center,point=point,side=side,label=label,color=color))
    for side in [-1,1]:
        prev=None
        for item in sorted([r for r in assignments if r['side']==side],key=lambda r:-r['point'][1]):
            desired=item['point'][1]
            item['y']=desired if prev is None else min(desired,prev-3.5)
            prev=item['y']
    guides=[LINEWIDTH,1.,BEGIN,LINES]
    obj='Detail_labels_'+name; lines='Detail_guides_'+name
    for i,item in enumerate(assignments):
        side=item['side'];x=pcoords[:,0].min()-3 if side<0 else pcoords[:,0].max()+3
        end=np.array([x,item['y'],pcoords[:,2].max()+5]) @ basis
        line_end=end-side*1.8*basis[0]
        cmd.pseudoatom(obj,pos=end.tolist(),name=f'L{i+1}',resi=str(i+1),label=item['label'])
        guides += [COLOR,*rgb(item['color']),VERTEX,*item['center'],VERTEX,*line_end]
    guides += [END];cmd.load_cgo(guides,lines)
    cmd.hide('everything',obj);cmd.show('labels',obj);cmd.set('label_size',18,obj)
    cmd.group('Annotations',obj+' '+lines)
    camera(cmd,context+' or '+obj,basis,pad=1.30)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    manifest = dict(source=str(SOURCE), source_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        pdb='2BNR', structure='Deposited holo crystal structure; no CD3; not an MD frame or motion',
        palette_semantics='Region identity only; not binding energy or prediction accuracy',
        chain_roles=OBJECTS, sites=SITES, scenes=[])
    with pymol2.PyMOL() as pm:
        cmd = pm.cmd; style(cmd); cmd.set('cif_use_auth', 1)
        cmd.load(str(SOURCE), 'source_2BNR')
        source_atoms = cmd.get_model('source_2BNR and polymer.protein', state=1).atom
        source_map = {(a.chain,a.resi,a.name,a.alt): np.array(a.coord) for a in source_atoms}
        assert len(source_map) == len(source_atoms) == 6625
        for ch, name in OBJECTS.items():
            cmd.create(name, f'source_2BNR and polymer.protein and chain {ch}', 1, 1)
        cmd.delete('source_2BNR')
        cmd.group('TCR', 'TCR_alpha_D TCR_beta_E')
        cmd.group('Partners', 'HLA_A2_A beta2m_B NYESO1_peptide_C')
        cmd.select('protein_atoms', 'TCR_alpha_D or TCR_beta_E or HLA_A2_A or beta2m_B or NYESO1_peptide_C')
        cmd.dss('protein_atoms')
        for color, value, sel in [('tcr_alpha_base','#A8B8C2','TCR_alpha_D'),('tcr_beta_base','#C4CCC3','TCR_beta_E'),
                ('mhc_base','#CAD0D7','HLA_A2_A'),('b2m_base','#DEDFE2','beta2m_B'),('antigen_gold','#D9B326','NYESO1_peptide_C')]:
            cmd.set_color(color, rgb(value)); cmd.color(color, sel)
        for i, site in enumerate(SITES):
            expr = f"{OBJECTS[site['chain']]} and resi {'+'.join(map(str,site['ids']))}"
            cmd.select(site['name'], expr)
            residues = sorted(cmd.get_model(expr+' and name CA').atom, key=lambda a:int(a.resi))
            assert ''.join(AA[a.resn] for a in residues) == site['sequence']
            assert [int(a.resi) for a in residues] == site['ids']
            color_name = f'region_color_{i+1}'
            cmd.set_color(color_name, rgb(site['color'])); cmd.color(color_name, expr)
            for atom in residues:
                label = f"{AA[atom.resn]}{atom.resi} ({'alpha' if site['chain']=='D' else 'beta'}; {site['chain']})"
                atomsel = expr+f' and resi {atom.resi} and name CA'
                cmd.label(atomsel, repr(label))
                rows.append(dict(region=site['name'], chain=atom.chain, author_id=atom.resi,
                    amino_acid=AA[atom.resn], label=label, color=site['color'], scene=site['scene']))
        cmd.select('marked_sites', ' or '.join(s['name'] for s in SITES))
        assert cmd.count_atoms('marked_sites and name CA') == 19
        cmd.select('pmhc_contact_4p5A', 'byres ((TCR_alpha_D or TCR_beta_E) within 4.5 of (HLA_A2_A or beta2m_B or NYESO1_peptide_C))')
        assert cmd.count_atoms('pmhc_contact_4p5A and name CA') == 28
        # Camera: variable domains above constant domains, alpha to left of beta.
        centroid = lambda sel: cmd.get_coords(sel).mean(0)
        ac = centroid('TCR_alpha_D and resi 120-190 and name CA')
        bc = centroid('TCR_beta_E and resi 125-225 and name CA')
        v = centroid('(TCR_alpha_D or TCR_beta_E) and resi 5-105 and name CA')
        yy = v - (ac+bc)/2; yy /= np.linalg.norm(yy)
        xx = bc-ac; xx -= np.dot(xx,yy)*yy; xx /= np.linalg.norm(xx)
        basis = np.array([xx,yy,np.cross(xx,yy)])
        pcoords = cmd.get_coords('TCR_alpha_D or TCR_beta_E') @ basis.T
        leaders = [LINEWIDTH, 1.1, BEGIN, LINES]
        label_ys={}
        for side in [-1,1]:
            previous=None
            ordered=sorted([s for s in SITES if s['side']==side],
                key=lambda s:-(centroid(s['name']+' and name CA') @ basis.T)[1]-s['dy'])
            for site in ordered:
                desired=(centroid(site['name']+' and name CA') @ basis.T)[1]+site['dy']
                label_ys[site['name']]=desired if previous is None else min(desired,previous-11)
                previous=label_ys[site['name']]
        for i, site in enumerate(SITES):
            center = centroid(site['name']+' and name CA')
            p = center @ basis.T
            edge = pcoords[:,0].min()-16 if site['side']<0 else pcoords[:,0].max()+16
            end = np.array([edge,label_ys[site['name']],pcoords[:,2].max()+8]) @ basis
            cmd.pseudoatom('Overview_labels', pos=end.tolist(), name=f'L{i+1}', resi=str(i+1),
                label=site['title']+'\n'+site['subtitle'], color=f'region_color_{i+1}')
            line_end=end-site['side']*7*basis[0]
            leaders += [COLOR,*rgb(site['color']),VERTEX,*center,VERTEX,*line_end]
        leaders += [END]
        cmd.load_cgo(leaders, 'Overview_guides')
        cmd.hide('everything', 'Overview_labels'); cmd.show('labels', 'Overview_labels')
        cmd.group('Annotations', 'Overview_labels Overview_guides')
        cmd.set('label_size', 20, 'Overview_labels')
        def store(name, description, preview=False):
            cmd.deselect(); cmd.scene(name, 'store', message=description)
            manifest['scenes'].append(dict(name=name, description=description, view=list(cmd.get_view())))
            if preview:
                cmd.png(str(OUT/(name+'.png')), width=1600, height=1250, dpi=200, ray=1)
            print('Stored',name,flush=True)
        base_display(cmd); cmd.enable('Annotations'); cmd.enable('Overview_labels'); cmd.enable('Overview_guides')
        camera(cmd,'protein_atoms or Overview_labels',basis,pad=1.30)
        store('01_overview','All six requested regions plus alpha Y100, on the deposited holo 2BNR structure. Label leaders are annotations, not molecular bonds.',True)
        base_display(cmd); cmd.disable('Partners'); cmd.enable('Overview_labels'); cmd.enable('Overview_guides')
        camera(cmd,'TCR_alpha_D or TCR_beta_E or Overview_labels',basis,pad=1.30)
        store('02_TCR_only','Same coordinates and colors, pMHC hidden for clearer inspection of the two TCR chains.',True)
        base_display(cmd)
        view_selection = '(TCR_alpha_D or TCR_beta_E) and resi 20-105'
        cmd.hide('cartoon', 'protein_atoms and not (byres (NYESO1_peptide_C around 18))')
        cmd.show('cartoon', 'protein_atoms and byres (NYESO1_peptide_C around 18)')
        cmd.hide('sticks', 'site_02_alpha_198_202 or site_03_beta_230_233')
        cmd.hide('spheres', 'site_02_alpha_198_202 or site_03_beta_230_233')
        face='NYESO1_peptide_C or site_01_alpha_96_99 or site_04_beta_29_30 or site_05_alpha_51 or site_06_beta_100_101 or site_07_alpha_Y100'
        face_items=[(s['name']+' and name CA',s['title'],s['color']) for s in SITES if s['ids'][0]<150]
        detail_labels(cmd,'face',face_items,face,basis)
        store('03_binding_face','Binding-face overview. Positions are author numbered; color is only a region identifier.',True)
        for site in SITES:
            base_display(cmd)
            context = f"byres (protein_atoms within 7 of {site['name']})"
            cmd.hide('everything','protein_atoms'); cmd.show('cartoon',context)
            cmd.show('sticks',context); cmd.set('stick_radius',.12,'protein_atoms'); cmd.set('stick_radius',.23,site['name'])
            cmd.show('spheres',site['name']+' and name CA')
            # Nearby peptide/MHC and other TCR-chain residues stay visible.
            items=[]
            for atom in sorted(cmd.get_model(site['name']+' and name CA').atom,key=lambda a:int(a.resi)):
                items.append((site['name']+f' and resi {atom.resi} and name CA',
                    f"{AA[atom.resn]}{atom.resi} ({'alpha' if atom.chain=='D' else 'beta'}; {atom.chain})",site['color']))
            detail_labels(cmd,site['scene'],items,context,basis)
            store(site['scene'],f"{site['title']} ({site['subtitle']}). Surrounding deposited residues within 7 A shown; no hydrogen bonds or energetic hotspots are inferred.",site['name'] in ['site_01_alpha_96_99','site_05_alpha_51'])
        cmd.scene('01_overview','recall',animate=0)
        cmd.set('scene_buttons',1)
        cmd.deselect()
        cmd.save(str(OUT/'1G4_2BNR_marked_sites.pse'))
        cmd.save(str(OUT/'2BNR_protein_coordinates.pdb'),'protein_atoms',state=1)
        final_atoms=cmd.get_model('protein_atoms').atom
        assert len(final_atoms)==len(source_atoms)
        assert all(np.array_equal(np.array(a.coord),source_map[a.chain,a.resi,a.name,a.alt]) for a in final_atoms)
        manifest['protein_atoms']=len(final_atoms); manifest['marked_residues']=len(rows)
        manifest['protein_coordinates_unchanged']=True
        manifest['pymol_version']=list(cmd.get_version())
        manifest['session_scene_names']=cmd.get_scene_list()
    with (OUT/'residue_color_scene_index.csv').open('w',encoding='utf-8-sig',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    # Reopen in a fresh PyMOL instance: portable session, no source reloading.
    with pymol2.PyMOL() as check:
        c=check.cmd;c.load(str(OUT/'1G4_2BNR_marked_sites.pse'))
        assert c.count_atoms('protein_atoms')==6625
        assert c.count_atoms('marked_sites and name CA')==19
        assert len(c.get_scene_list())==10
        for site in SITES:
            assert c.count_atoms(site['name']+' and name CA')==len(site['ids'])
            c.scene(site['scene'],'recall',animate=0)
        manifest['fresh_session_reload_passed']=True
    (OUT/'session_audit.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
    print('Saved and independently reloaded',OUT/'1G4_2BNR_marked_sites.pse',flush=True)


if __name__=='__main__':
    main()
