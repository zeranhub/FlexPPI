"""Verify 2-D structure callouts against PyMOL-rendered CA landmark centers."""
from pathlib import Path
import json
import numpy as np
from PIL import Image
from render_mechanistic_structures import pymol2

OUT=Path(__file__).resolve().parents[1]/'output/flexppi_pdz3_alpha3_focus_20260914/structure_assets'
REC=json.loads((OUT/'render_provenance.json').read_text())

def main():
    checks=[]
    with pymol2.PyMOL() as pm:
        c=pm.cmd;c.load(str(OUT/'pdz3_alpha3.pse'))
        c.bg_color('white');c.set('sphere_scale',.15)
        for scene in ['whole','local']:
            c.scene(scene,'recall');c.set_view(REC['views'][scene])
            for r in [334,337,397]:
                c.hide('everything');c.show('spheres',f'pdz3 and resi {r} and name CA')
                c.color('magenta',f'pdz3 and resi {r} and name CA')
                path=OUT/f'qa_landmark_{scene}_{r}.png'
                c.png(str(path),width=1700,height=1300,ray=1)
                im=np.asarray(Image.open(path).convert('RGB'))
                ys,xs=np.where(np.any(im<240,axis=2))
                actual=np.array([(xs.min()+xs.max())/2,(ys.min()+ys.max())/2])
                view=np.array(REC['views'][scene]);rot=view[:9].reshape(3,3).T
                q=rot@(np.array(REC['label_ca_coordinates'][str(r)])-view[12:15])+view[9:12]
                halfh=-view[11]*np.tan(np.deg2rad(20/2));halfw=halfh*1700/1300
                expected=np.array([(.5+q[0]/(2*halfw))*1700,(.5-q[1]/(2*halfh))*1300])
                checks.append(dict(scene=scene,residue=r,rendered_pixel=actual.tolist(),
                    projected_pixel=expected.tolist(),error_px=float(np.linalg.norm(actual-expected))))
    (OUT/'projection_checks.json').write_text(json.dumps(checks,indent=2),encoding='utf-8')
    print(json.dumps(checks,indent=2))

if __name__=='__main__':main()
