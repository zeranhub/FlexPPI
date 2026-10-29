"""Verify that the visual redesign preserves all archived source values and coordinates."""
from pathlib import Path
import hashlib
import json

import numpy as np

HERE = Path(__file__).resolve().parent
BASELINE = Path(r'E:/Project/FlexPPI/paper/figures_v58') / HERE.name
sha = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
checks = {}
for path in sorted(HERE.glob('source_*.csv')):
    previous = BASELINE / path.name
    assert previous.is_file(), path.name
    assert sha(path) == sha(previous), path.name
    checks[path.name] = {'sha256': sha(path), 'unchanged_from_v58': True}

current = np.load(HERE / 'source_heatmap_arrays.npz')
previous = np.load(BASELINE / 'source_heatmap_arrays.npz')
for key in current.files:
    assert np.array_equal(current[key], previous[key], equal_nan=True), key
    checks[key] = {'shape': list(current[key].shape), 'unchanged_from_v58': True}
assert current['CDR_delta_RMSF'].shape == (3, 45)
assert current['interface_delta_Hchi'].shape == (3, 30)
assert np.isfinite(current['interface_delta_Hchi']).sum(axis=1).tolist() == [23, 23, 22]
assert max(np.nanmax(np.abs(current[key])) for key in current.files) < 0.45

structure = json.loads((HERE / 'structure/molecular_view.json').read_text(encoding='utf8'))
old_structure = json.loads((BASELINE / 'structure/molecular_view.json').read_text(encoding='utf8'))
assert structure['input_SHA256'] == old_structure['input_SHA256']
for key in ['atom_coordinates', 'distances_A', 'common_numbering', 'overview_view', 'local_view']:
    assert structure[key] == old_structure[key], key
assert sha(HERE / 'structure/4jff_model_input.pdb') == structure['input_SHA256']
checks['molecular_coordinates_views_and_distances'] = 'unchanged_from_v58'

audit_path = HERE / 'figure_data_audit.json'
audit = json.loads(audit_path.read_text(encoding='utf8'))
assert audit['minimum_font_pt'] >= 6.5
assert all(item['inside_canvas'] for item in audit['visible_text'])
assert audit['dimensions_mm'] == [183, 175]
for ext in ['png', 'svg', 'pdf']:
    path = HERE / f'fig4_mel5_compact_joint.{ext}'
    assert path.is_file() and path.stat().st_size > 1000
    checks[path.name] = {'bytes': path.stat().st_size, 'sha256': sha(path)}
checks['all_text_inside_canvas'] = True
checks['heatmap_clipped_values'] = 0
checks['maximum_absolute_heatmap_value'] = float(max(np.nanmax(np.abs(current[key])) for key in current.files))
audit['source_integrity_checks'] = checks
audit['status'] = 'data_verified_pending_visual_QA'
audit_path.write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding='utf8')
print(json.dumps({'data_integrity': 'passed', 'checks': len(checks), 'clipped_values': 0}))
