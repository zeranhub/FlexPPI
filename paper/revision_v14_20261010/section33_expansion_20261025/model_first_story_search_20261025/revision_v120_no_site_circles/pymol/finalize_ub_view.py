from pathlib import Path
from hashlib import sha256
import json
import numpy as np
from PIL import Image
from ub_view_candidates import pymol2,fit,project
ROOT=Path(__file__).resolve().parent
manifest_path=ROOT/'render_manifest.json'
manifest=json.loads(manifest_path.read_text(encoding='utf-8'))
candidates=json.loads((ROOT/'view_candidates/candidates.json').read_text())
selection='orient_x30'
W,H=1600,1400
with pymol2.PyMOL() as pm:
    c=pm.cmd;c.load(str(ROOT/'ub_overview.pse'))
    original=c.get_coords('all').copy();assert len(original)==935
    c.set_view(candidates[selection]['view'])
    c.set('ray_opaque_background',0);c.set('opaque_background',0);c.set('max_threads',4)
    view=fit(c,W,H);c.viewport(W,H)
    c.png(str(ROOT/'ub_overview.png'),width=W,height=H,dpi=800,ray=1,quiet=1)
    c.save(str(ROOT/'ub_overview.pse'))
    assert np.array_equal(original,c.get_coords('all'))
    px=project(original,view,W,H)
    assert ((px[:,0]>=0)&(px[:,0]<W)&(px[:,1]>=0)&(px[:,1]<H)).all()
    landmarks={f'Q{n}':project(c.get_coords(f'chain A and resi {n} and name CA'),view,W,H)[0].tolist() for n in [41,62]}
    im=Image.open(ROOT/'ub_overview.png').convert('RGBA')
    bbox=im.getchannel('A').getbbox();pad=45
    crop=[max(0,bbox[0]-pad),max(0,bbox[1]-pad),min(W,bbox[2]+pad),min(H,bbox[3]+pad)]
    manifest['assets']['ub_overview'].update(png='ub_overview.png',pse='ub_overview.pse',width=W,height=H,
        view=list(c.get_view()),landmarks=landmarks,crop=crop,sha256=sha256((ROOT/'ub_overview.png').read_bytes()).hexdigest(),
        selected_atom_count=935,selected_atoms_in_frame=True,stick_atom_count=0,
        view_selection=selection,atomic_coordinates_unchanged=True,camera_candidates_reviewed=12)
manifest_path.write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(manifest['assets']['ub_overview'],ensure_ascii=False))
