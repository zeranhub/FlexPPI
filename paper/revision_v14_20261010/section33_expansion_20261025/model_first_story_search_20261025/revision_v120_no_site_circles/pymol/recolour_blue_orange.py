"""Native PyMOL recolouring only. Coordinates and camera must stay identical."""
from pathlib import Path
from hashlib import sha256
import json,os
import numpy as np
ENV=Path('D:/Anaconda/envs/pvb_cpu')
os.environ['PATH']=str(ENV/'Library/bin')+os.pathsep+os.environ.get('PATH','')
DLL_HANDLES=[os.add_dll_directory(str(p)) for p in [ENV/'Library/bin',ENV/'DLLs'] if p.exists()]
import pymol2
ROOT=Path(__file__).resolve().parent
manifest=json.loads((ROOT/'render_manifest.json').read_text(encoding='utf-8'))
BLUE='#2166AC';ORANGE='#B87936'
def rgb(h):return [int(h[i:i+2],16)/255 for i in [1,3,5]]
audit={}
for key,chain,sites,expected in [('cam_overview','B',[303,316],1344),('ub_overview','A',[41,62],935)]:
    with pymol2.PyMOL() as pm:
        c=pm.cmd;c.load(str(ROOT/(key+'.pse')))
        before=c.get_coords('all').copy();old_view=list(c.get_view())
        assert len(before)==expected
        c.set_color('target_chain',rgb(BLUE));c.set_color('key',rgb(ORANGE))
        c.color('target_chain',f'chain {chain}')
        c.color('key',f'chain {chain} and resi '+ '+'.join(map(str,sites)))
        c.set('max_threads',4);c.set_view(old_view);c.viewport(1600,1400)
        c.png(str(ROOT/(key+'.png')),width=1600,height=1400,dpi=800,ray=1,quiet=1)
        c.save(str(ROOT/(key+'.pse')))
        assert np.array_equal(before,c.get_coords('all'))
        assert np.array_equal(np.asarray(old_view),np.asarray(c.get_view()))
        z=manifest['assets'][key]
        z['sha256']=sha256((ROOT/(key+'.png')).read_bytes()).hexdigest()
        z['target_chain_colour']=BLUE
        audit[key]={'all_atom_coordinates_exactly_preserved':True,'all18_camera_values_exactly_preserved':True,
            'atoms':expected,'target_chain':chain,'target_colour':BLUE,'key_sites':sites,'key_colour':ORANGE}
manifest['colors']['target_chain']=BLUE
manifest['colors']['key']=ORANGE
manifest['recolouring_only_v118']=audit
(ROOT/'render_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(audit,ensure_ascii=False))
