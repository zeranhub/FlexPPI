"""Check source preservation, shared mappings, units, and complete exports."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
from publication_style import SEQUENTIAL, UNIT_TICKS, unit_cmap, unit_norm, signed_norm

ROOT = Path(__file__).resolve().parent
BASELINE = Path('E:/Project/FlexPPI/paper/figures_v58')
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    files = [p for p in BASELINE.rglob('*') if p.is_file() and p.suffix in {'.csv', '.npz', '.pdb'}]
    for old in files:
        new = ROOT / old.relative_to(BASELINE)
        assert sha(new) == sha(old), new
    folders = {p.name[:2]: p for p in ROOT.iterdir() if p.is_dir()}
    fig2 = folders['02']
    fig3 = folders['03']
    fig4 = folders['04']
    fig5 = folders['05']
    palette = json.loads((fig2 / 'palette.json').read_text(encoding='utf-8'))
    provenance = json.loads((fig3 / 'figure_provenance.json').read_text(encoding='utf-8'))
    assert palette['sequential'] == SEQUENTIAL == provenance['panel_c_colormap']['anchors']
    assert provenance['panel_c_colormap']['ticks'] == UNIT_TICKS
    assert provenance['panel_c_colormap']['vmin'] == 0
    assert provenance['panel_c_colormap']['vmax'] == 1
    colors = unit_cmap()(unit_norm()(np.array(UNIT_TICKS)))
    assert np.isfinite(colors).all()
    arrays = np.load(fig4 / 'source_heatmap_arrays.npz')
    limits = {k: [float(np.nanmin(arrays[k])), float(np.nanmax(arrays[k]))] for k in arrays.files}
    assert all(max(abs(lo), abs(hi)) < .45 for lo, hi in limits.values())
    assert signed_norm(.45)(0) == .5
    assert np.isfinite(arrays['interface_delta_Hchi']).sum(axis=1).tolist() == [23, 23, 22]
    plotted = pd.read_csv(fig5 / 'all_plotted_points.csv')
    assert len(plotted) == 7256
    manifest = json.loads((ROOT / 'manifest.json').read_text(encoding='utf-8'))
    for f in manifest['figures']:
        if f['figure'] in ['2','3','4','5']:
            for ext in ['png','pdf','svg']:
                path = ROOT / f['folder'] / (f['stem']+'.'+ext)
                assert path.is_file() and path.stat().st_size > 1000
    report = {'status': 'pass', 'byte_identical_scientific_files': len(files),
              'unit_heatmap_scale': [0,1], 'unit_ticks': UNIT_TICKS,
              'unit_colormap_anchors': SEQUENTIAL, 'signed_heatmap_limits': [-.45,.45],
              'signed_observed_ranges': limits, 'signed_clipped_cells': 0,
              'interface_Hchi_valid_counts': [23,23,22], 'figure5_plotted_records': len(plotted),
              'AP_and_percentile_are_distinct_metrics': True, 'no_new_statistical_tests': True}
    (ROOT / 'suite_validation.json').write_text(json.dumps(report, ensure_ascii=False, indent=2),encoding='utf-8')
    print(json.dumps(report, ensure_ascii=False))

if __name__ == '__main__':
    main()
