"""Only recolor and ray-trace existing verified PyMOL scenes; no new frames."""
import os
from pathlib import Path
ENV = Path(r'D:\Anaconda\envs\pvb_cpu')
os.environ['PATH'] = str(ENV/'Library/bin') + os.pathsep + os.environ.get('PATH','')
_dll_handles = [os.add_dll_directory(str(p)) for p in [ENV/'Library/bin', ENV/'DLLs'] if p.is_dir()]
import json
import hashlib
import numpy as np
import pymol2

SOURCE = Path(r'E:\Project\FlexPPI\experiments\_case_validation_sources\ubiquitin')
OUT = Path(__file__).resolve().parent / 'structure_renders'
OUT.mkdir(parents=True, exist_ok=True)
VIEWS = json.loads((SOURCE/'views.json').read_text(encoding='utf-8'))

def recolor(cmd,overview=False):
    for name, hexcode in [('ubiquitin_blue','#2166AC'),('uba_grey','#A6AEB4'),('q41_orange','#B87936'),('q62_purple','#3A8D8C')]:
        cmd.set_color(name, [int(hexcode[i:i+2],16)/255 for i in (1,3,5)])
    if overview:
        cmd.color('ubiquitin_blue','complex and chain A')
        cmd.color('uba_grey','complex and chain B')
        cmd.color('q41_orange','complex and chain A and resi 41')
        cmd.color('q62_purple','complex and chain A and resi 62')

with pymol2.PyMOL() as session:
    cmd = session.cmd
    cmd.load(str(SOURCE/'ubiquitin_2QHO.pse'))
    source_coords = {name:np.array(cmd.get_coords(name),copy=True) for name in cmd.get_names('objects') if cmd.count_atoms(name)}
    audit = {'source_pse':str(SOURCE/'ubiquitin_2QHO.pse'),'source_pse_sha256':hashlib.sha256((SOURCE/'ubiquitin_2QHO.pse').read_bytes()).hexdigest(),'render_rule':'Existing PSE objects, MD representatives, geometry and cameras unchanged; color definitions only.','images':{}}
    cmd.set('max_threads', 4)
    for scene, name, size in [('01_complex','ubiquitin_context',(1500,1250)),('02_Q41_MD_observed','q41_zoom',(1740,1450)),('02_Q62_MD_observed','q62_zoom',(1740,1450))]:
        cmd.scene(scene,'recall',animate=0)
        recolor(cmd,scene=='01_complex')
        cmd.set_view(VIEWS[scene])
        # Raster labels would shrink below 8 pt at print size; the compositor
        # supplies editable 8 pt vector labels at the same structural context.
        cmd.hide('labels','all')
        before = np.array(cmd.get_view())
        path = OUT/(name+'.png')
        cmd.png(str(path),width=size[0],height=size[1],dpi=600,ray=1,quiet=1)
        audit['images'][name]={'scene':scene,'pixels':list(size),'camera_max_difference':float(np.max(abs(before-np.array(VIEWS[scene])))),'path':str(path),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}
    for name, coords in source_coords.items():
        assert np.array_equal(coords,cmd.get_coords(name)), name
    audit['all_object_coordinates_unchanged'] = True
    cmd.scene('01_complex','recall',animate=0); recolor(cmd,True)
    cmd.hide('labels','all');cmd.scene('01_complex','store')
    cmd.save(str(OUT/'ubiquitin_Q41_Q62_recolored.pse'))
(OUT/'recolor_render_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(audit,ensure_ascii=False,indent=2))
