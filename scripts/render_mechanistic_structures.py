"""Ray-trace authentic input structures with model deltas; no generated geometry.

Run with the existing PyMOL conda environment. Deposited NMR model1 is used only
as a spatial reference; quantitative KIX panels include all20 models.
"""
from pathlib import Path
import os
import sys
import json
import csv
from collections import defaultdict
import numpy as np

ENV=Path('D:/Anaconda/envs/pvb_cpu')
os.environ['PATH']=str(ENV/'Library/bin')+os.pathsep+os.environ.get('PATH','')
if hasattr(os,'add_dll_directory'):
    DLL_HANDLES=[os.add_dll_directory(str(p)) for p in [ENV/'Library/bin',ENV/'DLLs'] if p.exists()]
import pymol2

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'artifacts/mechanistic_case_studies_20260913'
OUT=ROOT/'output/flexppi_mechanistic_case_studies_20260913'
ASSETS=OUT/'structure_assets'
JOBS=[
 dict(name='pdz3',state='pdz_wt_free',contrast='pdz_deletion_free',metric='chi_entropy',targets='A',focus='chain A and resi 330-336+389-402',labels=[['A','332','D332'],['A','397','removed 396-402']],accent='chain A and resi 396-402'),
 dict(name='tcr',state='tcr_bound_wtseq',contrast='tcr_apo_holo',metric='rmsf',targets='DE',focus='chain D and resi 116-206 or chain E and resi 218-231',labels=[['D','150','C alpha'],['E','226','beta 218-231']]),
 dict(name='tcr_context',state='tcr_bound_wtseq',contrast='tcr_context',metric='rmsf',targets='DE',focus='chain D and resi 116-206 or chain E and resi 218-231',labels=[['D','150','C alpha'],['E','226','beta 218-231']]),
 dict(name='thrombin',state='th_tm456',contrast='th_tm4_increment',metric='rmsf',targets='AM',focus='chain M and resi 60-64',labels=[['M','60A','60s loop'],['I','360','TM4']],accent='chain I and resi 345-387'),
 dict(name='kix',state='kix_m01_both',contrast='kix_m01_mll_effect',metric='rmsf',targets='B',focus='chain B and resi 610-625',labels=[['B','617','L12 loop'],['A','303','cMyb'],['C','849','MLL']]),
 dict(name='giv',state='giv_bound',contrast='giv_binding',metric='rmsf',targets='A',focus='chain A and resi 201-216',labels=[['A','211','W211 / SwII'],['C','1681','GIV']]),
 dict(name='bap1',state='bap_wt_bound',contrast='bap_mutation_free',metric='chi_entropy',targets='L',focus='chain L and resi 15+45-53+168+170+228',labels=[['L','49','L49'],['L','168','F168'],['L','170','F170']]),
]

def main():
    ASSETS.mkdir(parents=True,exist_ok=True)
    specs={s['id']:s for s in json.loads((BASE/'states.json').read_text())}
    groups=defaultdict(lambda:defaultdict(list))
    with (BASE/'residue_predictions.csv').open() as f:
        for r in csv.DictReader(f):
            for metric in ['rmsf','chi_entropy']:
                if r['delta_'+metric] and not(metric=='chi_entropy' and r['aa'] in 'AGP'):
                    groups[(r['contrast'],metric)][(r['chain'],r['author_id'])].append(float(r['delta_'+metric]))
    records=[]
    for job in JOBS:
        values={k:float(np.mean(v)) for k,v in groups[(job['contrast'],job['metric'])].items()}
        limit=1.0 if job['metric']=='rmsf' else .2
        with pymol2.PyMOL() as pm:
            cmd=pm.cmd
            cmd.load(specs[job['state']]['path'],'protein')
            cmd.hide('everything');cmd.show('cartoon','protein');cmd.dss('protein')
            cmd.bg_color('white')
            for k,v in dict(ray_opaque_background=1,ray_trace_mode=1,ray_trace_gain=.08,ray_shadows=0,antialias=2,
                orthoscopic=1,depth_cue=0,specular=.13,ambient=.55,direct=.45,cartoon_fancy_helices=1,
                cartoon_fancy_sheets=1,cartoon_flat_sheets=1,cartoon_smooth_loops=1,cartoon_sampling=14,max_threads=4,
                label_color='black',label_size=18,label_outline_color='white',label_font_id=7).items():
                cmd.set(k,v)
            cmd.set_color('missing',[.69,.72,.75]);cmd.color('missing','all')
            target=' or '.join('chain '+c for c in job['targets'])
            cmd.select('target','protein and ('+target+')')
            cmd.set_color('partner',[.80,.72,.57]);cmd.color('partner','protein and not target')
            if job['name']=='kix':
                cmd.set_color('myb',[.56,.49,.67]);cmd.color('myb','chain A')
            cmd.alter('protein','b=-999')
            for (ch,author),v in values.items():
                cmd.alter(f'protein and chain {ch} and resi {author}',f'b={v:.10f}')
            cmd.set_color('delta_blue',[33/255,102/255,172/255])
            cmd.set_color('delta_zero',[247/255]*3)
            cmd.set_color('delta_red',[178/255,24/255,43/255])
            cmd.spectrum('b','delta_blue delta_zero delta_red','target and b > -998',minimum=-limit,maximum=limit)
            if job.get('accent'):
                cmd.set_color('accent',[.95,.64,.25]);cmd.color('accent',job['accent'])
            if job['name'].startswith('tcr'):
                sites=json.loads((BASE/'tcr/selected_nmr_site_mapping.json').read_text())
                ids=[s['pdb_author_id'] for s in sites if s.get('pdb_author_id')]
                cmd.create('nmr_markers','protein and chain E and name CA and resi '+'+'.join(ids))
                cmd.hide('everything','nmr_markers')
                cmd.show('spheres','nmr_markers');cmd.set('sphere_scale',.65,'nmr_markers')
                cmd.set_color('nmr_purple',[118/255,42/255,131/255]);cmd.color('nmr_purple','nmr_markers')
            cmd.select('focus','protein and ('+job['focus']+')')
            for n,(ch,r,label) in enumerate(job['labels']):
                sel=f'protein and chain {ch} and resi {r} and name CA'
                if cmd.count_atoms(sel):
                    cmd.pseudoatom('tag'+str(n),pos=cmd.get_atom_coords(sel),label=label)
                    cmd.set('label_position',[1.5,1.5,2],'tag'+str(n))
            # Avoid occluding the protein with raster text; vector panel titles
            # and component keys carry labels in the composed figure.
            cmd.hide('labels','all')
            cmd.orient('target');cmd.zoom('protein',buffer=4,complete=1)
            fullview=cmd.get_view()
            cmd.scene('whole','store')
            full=ASSETS/(job['name']+'_whole.png')
            cmd.png(str(full),width=1500,height=1100,dpi=300,ray=1)
            # Local context: actual atoms within6 A of literature-defined region.
            cmd.hide('labels','all')
            if job['name'].startswith('tcr'):
                cmd.hide('everything','nmr_markers')
            cmd.select('local_context','protein and (byres (focus around 6) or focus)')
            cmd.hide('cartoon','protein and not local_context')
            cmd.show('sticks','focus and not name N+C+O')
            cmd.set('stick_radius',.12)
            cmd.zoom('local_context',buffer=3,complete=1)
            localview=cmd.get_view()
            local=ASSETS/(job['name']+'_local.png')
            cmd.png(str(local),width=1500,height=1100,dpi=300,ray=1)
            cmd.scene('local','store')
            cmd.save(str(ASSETS/(job['name']+'.pse')))
            # Explicit PML allows camera and presentation to be reproduced.
            pml=f"load {specs[job['state']]['path'].replace(chr(92),'/')}, protein\n# Full styling: scripts/render_mechanistic_structures.py\nset_view ({','.join(map(str,fullview))})\n"
            (ASSETS/(job['name']+'.pml')).write_text(pml,encoding='utf-8')
            record=dict(job,scale=[-limit,limit],colored_residues=len(values),clipped_residues=sum(abs(v)>limit for v in values.values()),
                        source=specs[job['state']],full_view=list(fullview),local_view=list(localview),
                        whole=str(full),local=str(local),coordinate_ensemble=False)
            records.append(record)
            print('Rendered',job['name'],flush=True)
    (ASSETS/'render_provenance.json').write_text(json.dumps(records,indent=2),encoding='utf-8')

if __name__=='__main__':
    main()
