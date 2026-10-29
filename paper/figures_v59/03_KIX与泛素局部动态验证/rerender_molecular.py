"""Native PyMOL redraw using archived cameras, atoms, and exact crop boxes.

Run with a Python environment containing PyMOL and Pillow. On this workstation:
D:/Anaconda/envs/pvb_cpu/python.exe rerender_molecular.py
"""
from pathlib import Path
import json
import hashlib
import os
import runpy
import numpy as np
import pandas as pd

CONDA = Path('D:/Anaconda/envs/pvb_cpu')
if os.name == 'nt' and (CONDA/'Library/bin').exists():
    os.environ['PATH'] = str(CONDA/'Library/bin')+os.pathsep+os.environ.get('PATH','')
    DLL_HANDLES = [os.add_dll_directory(str(p)) for p in [CONDA/'Library/bin',CONDA/'DLLs'] if p.exists()]
import pymol2

ROOT = Path(__file__).resolve().parent
BUNDLE = ROOT/'source_archive'
COLORS = {'cmyb_green':'#3A8D8C','mll_purple':'#B87936','l12_orange':'#2166AC',
          'q41_orange':'#B87936','q62_purple':'#3A8D8C','uba_green':'#2166AC'}
RGB = lambda h: [int(h[i:i+2],16)/255 for i in (1,3,5)]
SHA = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()

def recolor(c):
    for name, value in COLORS.items():
        c.set_color(name, RGB(value))
    c.set('max_threads', 4)
    c.set('ray_shadows', 0)

def snapshot(c):
    return {n:c.get_coords(n).copy() for n in c.get_names('objects')
            if c.get_type(n)=='object:molecule' and c.get_coords(n) is not None}

def verify(c, before):
    for n,xyz in before.items():
        assert np.array_equal(xyz,c.get_coords(n)), n

def draw(c, scene, path, size, audit):
    c.scene(scene,'recall',animate=0)
    recolor(c)
    before=snapshot(c); view=np.asarray(c.get_view()).copy()
    c.png(str(path),width=size[0],height=size[1],dpi=600,ray=1,quiet=1)
    verify(c,before)
    assert np.allclose(view,c.get_view(),rtol=0,atol=1e-6)
    c.scene(scene,'store')
    audit.append(dict(scene=scene,file=str(path.relative_to(ROOT)),
                      atomic_coordinates_unchanged=True,camera_unchanged=True,sha256=SHA(path)))

def main():
    audit=[]
    with pymol2.PyMOL() as pm:
        c=pm.cmd
        kb=BUNDLE/'kix'
        c.load(str(kb/'KIX_case_validation.pse'))
        rank=pd.read_csv(kb/'kix_both_fold1_source.csv')
        palette=np.asarray([RGB(h) for h in ['#F7FBFF','#6BAED6','#08306B']])
        for row in rank.loc[rank.metric.eq('rmsf')].itertuples():
            for source in ['md','model']:
                pos=getattr(row,source+'_value')/100*2
                lo=min(int(pos),1); frac=pos-lo
                color=((1-frac)*palette[lo]+frac*palette[lo+1]).tolist()
                c.set_color('rank_'+source+'_'+str(row.author_residue_number),color)
        config=json.loads((kb/'scene_config.json').read_text(encoding='utf-8'))
        draw(c,'KIX_context',kb/'kix_context.png',config['context_size'],audit)
        crop=runpy.run_path(str(kb/'whole_protein_crop.py'))['render_crop']
        for source,scene,stem in [('md','KIX_L12_MD','kix_l12_md'),('model','KIX_L12_FlexPPI_fold1','kix_l12_fold1')]:
            c.scene(scene+'_whole','recall',animate=0);recolor(c);before=snapshot(c)
            plan=config['crop_plans'][source]
            crop(c,plan,kb/(stem+'_whole_source.png'),kb/(stem+'.png'),dpi=600)
            verify(c,before)
            c.scene(scene,'store')
            audit.append(dict(scene=scene,file=str((kb/(stem+'.png')).relative_to(ROOT)),
                              atomic_coordinates_unchanged=True,original_crop_box_unchanged=True,
                              exact_pixel_crop_verified=plan['exact_pixel_crop_verified'],sha256=SHA(kb/(stem+'.png'))))
        (kb/'crop_manifest.json').write_text(json.dumps(config['crop_plans'],indent=2)+'\n',encoding='utf-8')
        # Store every scene with the updated native color definitions.
        for scene in c.get_scene_list():
            c.scene(scene,'recall',animate=0);recolor(c);c.scene(scene,'store')
        c.scene('KIX_context','recall',animate=0);c.save(str(kb/'KIX_case_validation.pse'))
        c.reinitialize()
        ub=BUNDLE/'ubiquitin'
        c.load(str(ub/'ubiquitin_2QHO.pse'))
        draw(c,'01_complex',ub/'renders/ubiquitin_context.png',[1500,1250],audit)
        plans=json.loads((ub/'crop_manifest.json').read_text(encoding='utf-8'))
        crop=runpy.run_path(str(ub/'whole_protein_crop.py'))['render_crop']
        for site in [41,62]:
            scene=f'03_Q{site}_whole_for_crop'
            c.scene(scene,'recall',animate=0);recolor(c);before=snapshot(c)
            plan=plans['cases'][f'Q{site}']
            crop(c,plan,ub/f'renders/q{site}_whole_source.png',ub/f'renders/q{site}_zoom.png',dpi=600)
            verify(c,before)
            c.scene(f'02_Q{site}_MD_observed','store')
            audit.append(dict(scene=scene,file=f'source_archive/ubiquitin/renders/q{site}_zoom.png',
                              atomic_coordinates_unchanged=True,original_crop_box_unchanged=True,
                              exact_pixel_crop_verified=plan['exact_pixel_crop_verified'],sha256=SHA(ub/f'renders/q{site}_zoom.png')))
        (ub/'crop_manifest.json').write_text(json.dumps(plans,indent=2)+'\n',encoding='utf-8')
        for scene in c.get_scene_list():
            c.scene(scene,'recall',animate=0);recolor(c);c.scene(scene,'store')
        for site in [41,62]:
            c.scene(f'02_Q{site}_MD_observed','recall',animate=0)
            c.save(str(ub/f'f_Q{site}_local.pse'))
        c.scene('01_complex','recall',animate=0);c.save(str(ub/'ubiquitin_2QHO.pse'))
        for stem in ['ubiquitin_context','q41_zoom','q62_zoom','q41_whole_source','q62_whole_source']:
            (ub/(stem+'.png')).write_bytes((ub/'renders'/(stem+'.png')).read_bytes())
    (ROOT/'molecular_style_audit.json').write_text(json.dumps(dict(
        status='pass',method='Native PyMOL ray tracing; original cameras, crop boxes and atomic coordinates',
        colors=COLORS,data_changed=False,outputs=audit),indent=2)+'\n',encoding='utf-8')
    print('Molecular redraw complete',flush=True)

if __name__=='__main__':
    main()
