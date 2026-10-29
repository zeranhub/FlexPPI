from pathlib import Path
import os, json, hashlib
_dll=os.add_dll_directory(r'D:\Anaconda\envs\pvb_cpu\Library\bin')
import pymol
pymol.finish_launching(['pymol','-cq'])
from pymol import cmd
root=Path(__file__).resolve().parent
reference=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010\section33_expansion_20261025\model_first_story_search_20261025\revision_v124_figure4_compact_palette\structure\MEL5_affinity_maturation.pse')
cmd.load(str(reference))
names=['orthoscopic','ray_shadows','specular','antialias','cartoon_fancy_helices','cartoon_smooth_loops','ray_opaque_background','ambient','direct','reflect','depth_cue','ray_trace_fog','ray_trace_mode','cartoon_loop_radius','cartoon_tube_radius','cartoon_flat_sheets','cartoon_sampling','stick_quality','cartoon_fancy_sheets','ray_trace_gain','two_sided_lighting','dash_width','dash_radius','dash_gap','dash_length']
settings={name:cmd.get(name,quiet=1) for name in names}
(root/'panel_a_render_style.json').write_text(json.dumps({'reference':str(reference),'sha256':hashlib.sha256(reference.read_bytes()).hexdigest(),'settings':settings},indent=2),encoding='utf-8')
print(json.dumps(settings))
cmd.quit()
