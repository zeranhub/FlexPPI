"""Authentic PDZ3 apo structure and local alpha3 contacts, PyMOL ray tracing."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from render_mechanistic_structures import pymol2

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'artifacts/pdz3_alpha3_focus_20260914'
OUT=ROOT/'output/flexppi_pdz3_alpha3_focus_20260914/structure_assets'

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    specs={s['id']:s for s in json.loads((BASE/'states.json').read_text())}
    df=pd.read_csv(BASE/'residue_summary.csv')
    sub=df[(df.background=='apo1bfe')&(df.kind=='deletion')&(df.context=='d7')&(df.mutant=='WT')&df.hchi_valid]
    with pymol2.PyMOL() as pm:
        c=pm.cmd;c.load(specs['apo1bfe_full_WT']['path'],'pdz3')
        c.hide('everything');c.show('cartoon');c.dss();c.bg_color('white')
        for k,v in dict(ray_opaque_background=1,ray_trace_mode=1,ray_trace_gain=.05,ray_shadows=0,antialias=2,
            orthoscopic=1,depth_cue=0,specular=.12,ambient=.6,direct=.4,cartoon_fancy_helices=1,
            cartoon_fancy_sheets=1,cartoon_flat_sheets=1,cartoon_smooth_loops=1,cartoon_sampling=14,max_threads=4,
            label_color='black',label_size=20,label_outline_color='white',label_font_id=7).items():c.set(k,v)
        for name,color in [('nochi',[.73,.75,.77]),('blue',[33/255,102/255,172/255]),('zero',[.97]*3),
                           ('red',[178/255,24/255,43/255]),('alpha',[.95,.64,.25]),('carbon',[.50,.53,.56])]:c.set_color(name,color)
        c.color('nochi');c.alter('all','b=-999')
        for r in sub.itertuples():c.alter(f'resi {r.author_number}',f'b={r.hchi_mean}')
        c.spectrum('b','blue zero red','b > -998',minimum=-.55,maximum=.55)
        c.color('alpha','resi 396-402')
        c.select('anchor','resi 334+337+358+397+400')
        c.show('sticks','anchor and not name N+C+O');c.set('stick_radius',.14)
        # Orient by the helix-facing core, then retain identical orientation in the zoom.
        c.orient('resi 318-323+330-338+355-358+394-402');c.turn('y',20);c.turn('z',-20)
        c.zoom('pdz3',buffer=3,complete=1)
        views={}
        views['whole']=list(c.get_view());c.scene('whole','store')
        c.png(str(OUT/'pdz3_whole.png'),width=1700,height=1300,dpi=300,ray=1)
        c.hide('sticks');c.hide('cartoon','not (resi 324-340+354-359+392-402)')
        c.show('sticks','resi 334+337+358+396+397+399+400 and not name N+C+O')
        c.set('cartoon_transparency',.18)
        c.zoom('resi 334+337+358+396+397+399+400',buffer=3,complete=1)
        views['local']=list(c.get_view());c.scene('local','store')
        c.png(str(OUT/'pdz3_local.png'),width=1700,height=1300,dpi=300,ray=1)
        # Separate chemically interpretable atom pair views, no dynamic edges.
        c.hide('everything');c.show('sticks','resi 334+337+397+399+400')
        c.color('carbon','elem C');c.color('alpha','resi 397+399+400 and elem C');c.color('blue','elem N');c.color('red','elem O')
        c.distance('static_pair','resi 334 and name OE1','resi 399 and name NH2')
        c.set('dash_color','carbon');c.set('dash_width',2);c.hide('labels','static_pair')
        c.zoom('resi 334+337+397+399+400',buffer=4,complete=1)
        views['contacts']=list(c.get_view());c.scene('contacts','store')
        c.png(str(OUT/'pdz3_contacts.png'),width=1700,height=1300,dpi=300,ray=1)
        c.save(str(OUT/'pdz3_alpha3.pse'))
        atoms={}
        for n in [323,328,334,337,355,358,393,396,397,399,400,401]:
            atoms[str(n)]=list(c.get_atom_coords(f'pdz3 and resi {n} and name CA'))
        record=dict(source=specs['apo1bfe_full_WT'],metric='delta Hchi d7-WT, five-model mean',
            alpha3='396-402 shown in orange as reference; absent from deleted input',scale=[-.55,.55],
            saturated_residues=int((sub.hchi_mean.abs()>.55).sum()),views=views,label_ca_coordinates=atoms,
            raster_size=[1700,1300],local_view='same orientation as whole',
            chemical_contact='E334 OE1-R399 NH2; static crystallographic distance only')
        (OUT/'render_provenance.json').write_text(json.dumps(record,indent=2),encoding='utf-8')
        print('Rendered whole / local / static contacts')

if __name__=='__main__':main()
