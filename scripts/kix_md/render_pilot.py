# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# Protein cartoons/local sticks: no compatible production asset; real PyMOL
# geometry, inheriting the existing audited render_mechanistic_structures.py.
"""Render fixed 2AGH model-1 geometry, colored by actual pilot/model observables."""
from pathlib import Path
import os
import csv
import json
import numpy as np

ENV = Path('D:/Anaconda/envs/pvb_cpu')
os.environ['PATH'] = str(ENV/'Library/bin') + os.pathsep + os.environ.get('PATH', '')
DLLS = [os.add_dll_directory(str(p)) for p in [ENV/'Library/bin', ENV/'DLLs'] if p.exists()]
import pymol2

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT/'artifacts/kix_md_20260914/pilot_results_20260916'
OUT = ROOT/'output/flexppi_kix_md_pilot_20260916/structure_assets'
PALETTES = {'rank': ['#F7FBFF','#6BAED6','#08306B'], 'delta': ['#2166AC','#F7F7F7','#B2182B']}

def colors(value, low, high, palette):
    t = np.clip((value-low)/(high-low), 0, 1)*2
    k = min(int(t), 1)
    rgb = [np.array([int(c[i:i+2],16)/255 for i in [1,3,5]]) for c in palette]
    return (rgb[k]*(1-(t-k))+rgb[k+1]*(t-k)).tolist()

def rows(name):
    with (BASE/name).open(encoding='utf-8',newline='') as f:
        return list(csv.DictReader(f))

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    ranks = [r for r in rows('relative_residues.csv') if r['state']=='both' and r['metric']=='rmsf' and r['normalization']=='percentile' and r['model_fold']=='ensemble_mean']
    delta = [r for r in rows('relative_deltas.csv') if r['metric']=='rmsf' and r['model_fold']=='ensemble_mean']
    jobs = [dict(name='both_rank', rows=ranks, field='value', low=0, high=100, palette='rank'),
            dict(name='both_delta', rows=[r for r in delta if r['contrast']=='both_minus_apo' and r['normalization']=='raw'],field='delta',low=-4.2,high=4.2,palette='delta'),
            dict(name='partner_rank_shift',rows=[r for r in delta if r['contrast']=='mll_minus_myb' and r['normalization']=='percentile'],field='delta',low=-60,high=60,palette='delta')]
    manifest=[]
    for job in jobs:
        assert len(job['rows'])==87
        for source in ['md','model']:
            with pymol2.PyMOL() as pm:
                c=pm.cmd
                c.load(str(BASE.parent/'remote_audit/both_protein.pdb'),'protein')
                c.hide('everything');c.show('cartoon','protein');c.dss('protein');c.bg_color('white')
                for k,v in dict(ray_opaque_background=1,ray_trace_mode=1,ray_trace_gain=.08,ray_shadows=0,
                                antialias=2,orthoscopic=1,depth_cue=0,specular=.13,ambient=.55,direct=.45,
                                cartoon_fancy_helices=1,cartoon_fancy_sheets=1,cartoon_flat_sheets=1,
                                cartoon_smooth_loops=1,cartoon_sampling=14,max_threads=4).items():
                    c.set(k,v)
                c.set_color('partnerA',[.65,.60,.71]);c.color('partnerA','chain A')
                c.set_color('partnerC',[.82,.72,.53]);c.color('partnerC','chain C')
                for row in job['rows']:
                    value=float(row[source+'_'+job['field']]); rid=int(row['author_residue_number'])
                    name='res'+str(rid)
                    c.set_color(name, colors(value,job['low'],job['high'],PALETTES[job['palette']]))
                    c.color(name,f'chain B and resi {rid}')
                c.select('focus','chain B and resi 614-621')
                c.show('sticks','focus and not hydro and not name N+C+O')
                c.set('stick_radius',.11)
                c.orient('chain B');c.zoom('protein',buffer=3,complete=1)
                # Same deterministic camera for MD and model; colors do not change geometry.
                whole_view=list(c.get_view())
                c.scene('whole','store')
                c.png(str(OUT/f'{job["name"]}_{source}_whole.png'),width=1350,height=1100,dpi=300,ray=1)
                c.select('local_context','protein and (byres (focus around 6) or focus)')
                c.hide('cartoon','protein and not local_context')
                c.zoom('local_context',buffer=2.5,complete=1)
                local_view=list(c.get_view())
                c.scene('L12','store')
                c.png(str(OUT/f'{job["name"]}_{source}_local.png'),width=1350,height=1050,dpi=300,ray=1)
                c.save(str(OUT/f'{job["name"]}_{source}.pse'))
                manifest.append(dict(name=job['name'],source=source,scale=[job['low'],job['high']],
                                     colored_residues=87,clipped_residues=sum(not(job['low']<=float(r[source+'_'+job['field']])<=job['high']) for r in job['rows']),
                                     whole_view=whole_view,local_view=local_view,reference_geometry='prepared 2AGH model1 ternary; NOT an MD snapshot',
                                     focus='B614-621, all atoms within 6 A by residue', palette=PALETTES[job['palette']]))
                print('Rendered',job['name'],source,flush=True)
    (OUT/'render_provenance.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')

if __name__=='__main__':
    main()
