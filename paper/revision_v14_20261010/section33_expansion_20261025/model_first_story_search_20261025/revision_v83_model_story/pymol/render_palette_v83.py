"""Re-render saved CaM coordinates in a palette separate from state comparisons.

Load the two previous auditable PyMOL sessions unchanged. Reuse every camera,
canvas and crop from the source manifest. No fitting, coordinate changes or
representation additions are performed here. Only molecular colours change.
"""
from pathlib import Path
import os, json, copy, hashlib
import numpy as np

OUT = Path(__file__).resolve().parent
ENV = Path('D:/Anaconda/envs/pvb_cpu')
os.environ['PATH'] = str(ENV/'Library/bin') + os.pathsep + os.environ.get('PATH','')
if hasattr(os, 'add_dll_directory'):
    DLL_HANDLES = [os.add_dll_directory(str(p)) for p in [ENV/'Library/bin',ENV/'DLLs'] if p.exists()]
import pymol2
from PIL import Image

COLORS = {
    'N_lobe':'#8A9299', 'C_lobe':'#BAC1C7', 'linker':'#D9DDE0',
    'partner':'#8B72A5', 'key':'#3A8D8C', 'sulfur':'#B18F31',
    'oxygen':'#A83D3D', 'nitrogen':'#447BA7', 'context':'#A6AEB4',
}

def rgb(value):
    return [int(value[i:i+2],16)/255. for i in (1,3,5)]

def recolor(c):
    for name, value in COLORS.items():
        c.set_color(name, rgb(value))
    for obj in c.get_names('objects'):
        for color, selection in [
            ('N_lobe','chain A and resi 4-75'),
            ('linker','chain A and resi 76-81'),
            ('C_lobe','chain A and resi 82-148'),
            ('partner','chain P'),
        ]:
            s = f'{obj} and ({selection})'
            c.color(color,s+' and elem C')
            c.set('cartoon_color',color,s)
        c.color('sulfur',f'{obj} and elem S')
        c.color('oxygen',f'{obj} and elem O')
        c.color('nitrogen',f'{obj} and elem N')
        sites = [int(obj.split('_M')[1])] if '_M' in obj else [72,109]
        for site in sites:
            c.color('key',f'{obj} and chain A and resi {site} and elem C')

def main():
    source = json.loads((OUT/'inputs/source_render_manifest.json').read_text())
    manifest = copy.deepcopy(source)
    manifest['palette'] = COLORS
    manifest['palette_semantics'] = {
        'molecular_components':'CaM neutral grey, partner muted purple, highlighted Met teal; all states use the same component colours.',
        'state_comparison_colours':'No NHE1/IQCG state colours are used in molecular cartoons. State titles and numerical plots are composed separately.',
    }
    manifest['recolor_coordinate_audit'] = {}
    manifest['recolor_source_sha256'] = {p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (OUT/'inputs').glob('*.pse')}
    groups = [
        ('CaM_shared_reference_v82.pse','CaM_shared_reference_v83.pse',
         [('whole_'+state,state) for state in ['NHE1','IQCG']]),
        ('CaM_local_context_from_v80.pse','CaM_local_context_v83.pse',
         [(f'detail_focus_M{site}_{state}',f'{state}_M{site}') for site in [72,109] for state in ['NHE1','IQCG']]),
    ]
    with pymol2.PyMOL() as pm:
        c = pm.cmd
        for source_name, output_name, targets in groups:
            c.reinitialize()
            c.load(str(OUT/'inputs'/source_name))
            before = {obj:np.asarray(c.get_coords(obj),float).copy() for obj in c.get_names('objects')}
            for scene in c.get_scene_list():
                c.scene(scene,'clear')
            recolor(c)
            for key, obj in targets:
                entry = manifest['assets'][key]
                c.disable('all');c.enable(obj)
                c.viewport(entry['width'],entry['height'])
                c.set_view(entry['view'])
                assert np.allclose(c.get_view(),entry['view'],atol=1e-5)
                c.png(str(OUT/entry['file']),width=entry['width'],height=entry['height'],dpi=450,ray=1,quiet=1)
                im = Image.open(OUT/entry['file']).convert('RGBA')
                entry['alpha_bbox'] = list(im.getchannel('A').getbbox())
                entry['palette'] = 'v83_component_palette'
                if key.startswith('detail_'):
                    entry['focus_carbon_color'] = COLORS['key']
                trim = entry['paired_trim']
                im.crop(trim['crop_box']).save(OUT/trim['file'])
                c.scene(key,'store')
                print('RENDERED',key,flush=True)
            for obj, xyz in before.items():
                after = np.asarray(c.get_coords(obj),float)
                error = float(np.max(np.abs(after-xyz)))
                assert np.array_equal(xyz,after), (obj,error)
                manifest['recolor_coordinate_audit'][source_name+':'+obj] = {
                    'atom_count':len(xyz),'maximum_coordinate_change_A':error,'exactly_identical_coordinates':True,
                }
            c.save(str(OUT/output_name))
        manifest['pymol_version'] = str(c.get_version()[0])
    manifest['recolor_camera_audit'] = 'All six source views, viewports, paired crop boxes and landmark pixels retained verbatim.'
    (OUT/'render_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    print('COMPLETE',str(OUT),flush=True)

if __name__ == '__main__':
    main()

