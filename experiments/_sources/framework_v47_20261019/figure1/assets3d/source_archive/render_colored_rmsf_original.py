"""Color the unchanged assembly by all 840 saved OOF RMSF predictions."""
from pathlib import Path
import csv
import hashlib
import json
import shutil

import numpy as np
import pymol
from pymol import cmd
from PIL import Image

HERE = Path(__file__).resolve().parent
OLD = HERE.parent / 'reference_rebuild_20261005' / 'flexppi_main'
DATA = HERE / 'source_data'
OUT = HERE / 'rmsf_structure'
DATA.mkdir(exist_ok=True)
OUT.mkdir(exist_ok=True)
for name in ('itrace_4g9f.pdb', 'itrace_4g9f_prediction.npz', 'itrace_4g9f_label.npz'):
    target = DATA / name
    if not target.exists():
        shutil.copy2(OLD / 'source_data' / name, target)


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


anchors = ['#1D5C93', '#57B9BB', '#F6D77B', '#DD8B50', '#B34342']
rgb = np.array([[int(h[k:k+2], 16)/255 for k in (1, 3, 5)] for h in anchors])


def color(value):
    assert 0 <= value <= 4
    k = min(int(value), 3)
    fraction = value-k
    return rgb[k]*(1-fraction)+rgb[k+1]*fraction


files = [DATA / 'itrace_4g9f.pdb', DATA / 'itrace_4g9f_prediction.npz', DATA / 'itrace_4g9f_label.npz']
hashes = {p.name: sha(p) for p in files}
old_manifest = DATA / 'render_source_manifest_v1.json'
if not old_manifest.exists():
    shutil.copy2(OLD / 'assets/render_source_manifest.json', old_manifest)
old_meta = json.loads(old_manifest.read_text(encoding='utf-8'))
assert hashes[files[0].name] == old_meta['pdb']['sha256']
assert hashes[files[1].name] == old_meta['prediction']['sha256']
assert hashes[files[2].name] == old_meta['label']['sha256']
pred = np.load(files[1], allow_pickle=False)
label = np.load(files[2], allow_pickle=False)
assert np.array_equal(pred['target'], label['values'], equal_nan=True)
assert np.array_equal(pred['valid'], label['mask'])
assert pred['mean'].shape == (840, 2)
assert pred['valid'][:, 0].sum() == 840
values = pred['mean'][:, 0].astype(float)
assert np.isfinite(values).all() and (values >= 0).all() and (values <= 4).all()
identities = list(zip(pred['chain'].tolist(), pred['author_residue_id'].tolist()))
assert len(set(identities)) == 840

pymol.finish_launching(['pymol', '-cq'])
cmd.load(str(files[0]), 'assembly')
cmd.remove('not polymer.protein')
cas = cmd.get_model('assembly and name CA').atom
assert len(cas) == 840
lookup = {(a.chain, a.resi): a for a in cas}
assert set(lookup) == set(identities)
cmd.hide('everything')
cmd.show('cartoon', 'assembly')
cmd.set('cartoon_fancy_helices', 1)
cmd.set('cartoon_smooth_loops', 1)
cmd.set('cartoon_sampling', 14)
cmd.set('cartoon_loop_radius', .25)
cmd.set('cartoon_transparency', 0)
cmd.set('ray_opaque_background', 0)
cmd.bg_color('white')
cmd.set('orthoscopic', 1)
cmd.set('antialias', 2)
cmd.set('ray_shadows', 0)
cmd.set('ambient', .43)
cmd.set('specular', .15)
cmd.set('ray_trace_mode', 0)
cmd.set('ray_trace_color', 'gray45')
cmd.set_view(old_meta['view'])
rows = []
for i, ((chain, residue), value) in enumerate(zip(identities, values)):
    name = 'rmsf_' + str(i)
    channels = color(value)
    cmd.set_color(name, channels.tolist())
    selected = f'assembly and chain {chain} and resi {residue}'
    assert cmd.count_atoms(selected + ' and name CA') == 1
    cmd.color(name, selected)
    atom = lookup[chain, residue]
    rows.append([i, chain, residue, str(pred['amino_acid'][i]), *atom.coord, value,
                 float(pred['target'][i, 0]), *channels.tolist()])
cmd.png(str(OUT / 'predicted_rmsf_full.png'), width=2200, height=1900, dpi=600, ray=1)
assert np.allclose(cmd.get_view(), old_meta['view'], atol=1e-6, rtol=0)
image = Image.open(OUT / 'predicted_rmsf_full.png').convert('RGBA')
bound = image.getchannel('A').getbbox()
assert bound
pad = 35
x0, y0, x1, y1 = bound
crop = (max(0, x0-pad), max(0, y0-pad), min(image.width, x1+pad), min(image.height, y1+pad))
image.crop(crop).save(OUT / 'predicted_rmsf_colored.png', dpi=(600, 600))
cmd.save(str(OUT / 'predicted_rmsf_colored.pse'))
with (OUT / 'all_840_rmsf_colors.csv').open('w', encoding='utf-8', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['global_index', 'chain', 'author_residue_id', 'amino_acid',
                     'ca_x_A', 'ca_y_A', 'ca_z_A', 'rmsf_prediction_A', 'rmsf_reference_A',
                     'red', 'green', 'blue'])
    writer.writerows(rows)
assert all(sha(p) == hashes[p.name] for p in files)
receipt = dict(sample_id='itrace_4g9f', residues=840, valid_RMSF=840,
               view=list(cmd.get_view()), palette_anchors_hex=anchors,
               color_range_A=[0, 4], palette_anchor_positions_A=[0, 1, 2, 3, 4],
               prediction_min_A=float(values.min()), prediction_max_A=float(values.max()),
               color_definition='Linear interpolation of RGB anchor colors at absolute RMSF values 0,1,2,3,4 Angstrom. Each residue retains the exact saved OOF prediction. No smoothing, normalization or clipping.',
               geometry_definition='Unchanged static protein structure with constant cartoon thickness; color encodes the saved predicted RMSF amplitude.',
               source_sha256=hashes, new_inference=0, new_simulation=0,
               crop_bbox=list(crop), image_size=list(image.crop(crop).size),
               outputs={p.name: sha(p) for p in OUT.iterdir() if p.is_file()},
               pymol_version=list(cmd.get_version()))
receipt['chains'] = {c: sum(k[0] == c for k in identities) for c in sorted(set(pred['chain']))}
(OUT / 'render_receipt.json').write_text(json.dumps(receipt, indent=2), encoding='utf-8')
print('Colored all 840 residues using exact saved RMSF values', flush=True)
cmd.quit()
