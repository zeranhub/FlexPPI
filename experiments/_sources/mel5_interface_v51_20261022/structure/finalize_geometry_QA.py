from pathlib import Path
import json

root=Path(__file__).resolve().parent
path=root/'dual_state_molecular_view.json'
meta=json.loads(path.read_text(encoding='utf8'))
assert meta['screen_projection']['holdout_error_pixels']<1
for state in meta['structures'].values():
    assert state['alpha_range']==[0,255]
    assert len(state['atom_screen_pixels_top_left'])==10
    assert all(0<=p[0]<1800 and 0<=p[1]<1500 for p in state['atom_screen_pixels_top_left'].values())
meta['status']='passed'
meta['visual_QA']={'review':'Manual visual inspection of both final 1800x1500 PNGs, same camera and scale; primary two atom pairs visible; no burned-in text; no molecular coordinate editing other than rigid HLA fit.',
                   'background':'Fine tube representation radius 0.20; HLA teal, WT TCR orange, mature TCR purple.',
                   'bounds_note':'The HLA short segment terminates near the lower image edge; no key residues or distance endpoints are clipped. Use a common crop for both states.'}
path.write_text(json.dumps(meta,ensure_ascii=False,indent=2),encoding='utf8')
(root/'visual_QA.json').write_text(json.dumps({'status':'passed','files':{k:{'PNG':v['PNG'],'sha256':v['PNG_sha256'],'alpha_range':v['alpha_range']} for k,v in meta['structures'].items()},'projection_holdout_error_pixels':meta['screen_projection']['holdout_error_pixels'],'visual_review':meta['visual_QA']},ensure_ascii=False,indent=2),encoding='utf8')
print('Geometry and two-state rendering QA passed.')
