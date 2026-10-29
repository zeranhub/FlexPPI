"""Read-only scientific/format verification with a generated audit record."""
from pathlib import Path
import json
import hashlib
import re
import subprocess
import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
from analyze_predictions import BASE,OUT


def main():
    states=json.loads((BASE/'states.json').read_text())
    cache_files=[]
    for s in states:
        for fold in range(5):
            p=BASE/'prediction_cache'/(s['id']+f'_fold{fold}.npz')
            with np.load(p,allow_pickle=False) as z:
                assert z['prediction'].shape==(s['n_residues'],2) and z['mask'].all()
                assert np.isfinite(z['prediction']).all()
                provenance=json.loads(str(z['provenance']))
                assert provenance['input_sha256']==hashlib.sha256(Path(s['path']).read_bytes()).hexdigest()
            cache_files.append(str(p))
    d=pd.read_csv(OUT/'residue_predictions.csv',dtype={'author_id':str})
    for metric in ['rmsf','hchi']:
        np.testing.assert_allclose(d['comparison_'+metric]-d['reference_'+metric],d['delta_'+metric],atol=1e-6)
    counts={}
    for case,n in [('B4.2.3',421),('1G4-CD3',366)]:
        q=d[d.common_all & d.case.eq(case) & d.kind.eq('real_apo_holo')]
        size=q.groupby(['contrast','fold']).size()
        assert (size==n).all()
        counts[case]=n
    figs=json.loads((OUT/'figure_export_audit.json').read_text())['figures']
    exports=[]
    for f in figs:
        stem=f['file']
        for ext in ['pdf','png','svg']:
            path=OUT/(stem+'.'+ext)
            assert path.stat().st_size>10000
        svg=ET.parse(OUT/(stem+'.svg'))
        ntext=len(svg.findall('.//{http://www.w3.org/2000/svg}text'))
        assert ntext>=10,'SVG text was outlined or rasterized'
        fonts=subprocess.check_output(['pdffonts',str(OUT/(stem+'.pdf'))],text=True)
        rows=[line for line in fonts.splitlines()[2:] if line.strip()]
        assert rows and all('TrueType' in line and re.search(r'yes\s+yes\s+yes',line) for line in rows),fonts
        info=subprocess.check_output(['pdfinfo',str(OUT/(stem+'.pdf'))],text=True)
        a,b=map(float,re.search(r'Page size:\s+([0-9.]+) x ([0-9.]+)',info).groups())
        assert a*25.4/72 <= 183.05,'Export wider than the figure contract'
        assert max(a,b)*25.4/72 <= 247
        exports.append(dict(stem=stem,n_editable_svg_text=ntext,pdf_width_mm=a*25.4/72,pdf_height_mm=b*25.4/72,
            fonts=rows,sha256={ext:hashlib.sha256((OUT/(stem+'.'+ext)).read_bytes()).hexdigest() for ext in ['pdf','png','svg']}))
    render_checks=[]
    for stem,n in [('b423',421),('1g4',366)]:
        meta=json.loads((OUT/'structure_assets'/(stem+'_render_provenance.json')).read_text())
        assert meta['n_CA_for_rigid_fit']==n and not meta['local_refit']
        assert meta['delta_maps'][0]['n_residues']==n
        for a in meta['assets']:
            camera=a['camera']
            actual=np.array(camera['view'][:9]).reshape(3,3).T
            np.testing.assert_allclose(actual,camera['basis_rows'],atol=1e-6)
            for l in a['landmarks']:
                xyz=(np.array(l['xyz_A'])-camera['center_world']) @ actual.T
                scale=camera['height_px']/(2*camera['half_height_A'])
                expected=[camera['width_px']/2+xyz[0]*scale,camera['height_px']/2-xyz[1]*scale]
                np.testing.assert_allclose(expected,l['pixel_xy_top_left'],atol=.001)
        render_checks.append(dict(system=stem,shared_residues=n,projection_verified=True,saturation=meta['delta_maps'][0]['n_clipped']))
    validation=[]
    for fold in range(5):
        tests=json.loads((BASE/f'chunk_validation_fold{fold}.json').read_text())
        validation.append(dict(fold=fold,maximum_absolute_native_chunk_error=max(t['maximum_absolute_error'] for t in tests)))
    report=dict(status='passed',n_prediction_runs=len(cache_files),n_source_rows=len(d),primary_common_counts=counts,
        native_chunk_checks=validation,render_checks=render_checks,exports=exports,
        no_new_training=True,no_new_MD=True,warning='Automated correctness/format checks do not establish biological validity; see interpretive report and visual QA.')
    (OUT/'deliverable_validation.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps(dict(status='passed',predictions=len(cache_files),figures=len(exports),format_checks=exports,render_checks=render_checks),indent=2),flush=True)


if __name__=='__main__':
    main()
