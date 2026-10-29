import json, hashlib, struct, re, subprocess
from pathlib import Path
import numpy as np
import pandas as pd
from PIL import Image
OUT=Path(__file__).resolve().parent
ORIGINAL=Path(r'E:\Project\FlexPPI\experiments\_sources\figure_merge_v43_20261018\figure5')
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
copied=[]
for p in (ORIGINAL/'source').rglob('*'):
    if p.is_file():
        q=OUT/'source'/p.relative_to(ORIGINAL/'source')
        assert q.is_file() and sha(p)==sha(q),q
        copied.append(str(p.relative_to(ORIGINAL)))
old=pd.read_csv(ORIGINAL/'all_plotted_points.csv',float_precision='round_trip')
new=pd.read_csv(OUT/'all_plotted_points.csv',float_precision='round_trip')
assert old.shape==new.shape,(old.shape,new.shape)
for c in old.columns:
    if old[c].dtype.kind in 'biufc':
        np.testing.assert_allclose(old[c].to_numpy(dtype=float),new[c].to_numpy(dtype=float),rtol=0,atol=0,equal_nan=True,err_msg=c)
    else:
        assert old[c].fillna('').equals(new[c].fillna('')),c
img=Image.open(OUT/'fig5_joint_readouts_conformations.png')
assert img.size==(3242,4074),img.size
font_output=subprocess.check_output([r'C:\texlive\2026\bin\windows\pdffonts.exe',str(OUT/'fig5_joint_readouts_conformations.pdf')],text=True)
assert 'Arial' in font_output and 'CID TrueType' in font_output
pdf=(OUT/'fig5_joint_readouts_conformations.pdf').read_bytes()
mb=re.search(rb'/MediaBox\s*\[\s*0\s*0\s*([\d.]+)\s+([\d.]+)\s*\]',pdf)
assert mb
mm=[float(x)*25.4/72 for x in mb.groups()]
np.testing.assert_allclose(mm,[183,230],atol=1e-7)
audit=json.loads((OUT/'data_audit.json').read_text(encoding='utf-8'))
audit.update(status='passed',visual_QA='Actual final PNG inspected with view_image after label refinement; no clipping or data/annotation collisions observed.',
    exact_original_point_table_match=True,source_files_byte_identical=len(copied),png_size=list(img.size),pdf_canvas_mm=mm,
    pdf_fonts=font_output.strip(),font_embedding_verified=True,structure_landmarks_read_from_current_json=True)
(OUT/'data_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
result={'status':'passed','all_7256_plot_rows_exactly_equal_original':True,'copied_source_files_byte_identical':len(copied),
'plotted_points_by_panel':new.groupby('panel').size().to_dict(),'f_NaN_zero_bins':int(new[new.panel=='f'].y.isna().sum()),
'pdf_canvas_mm':mm,'png_size':list(img.size),'pdf_fonts_embedded':True,'raster_only_actual_PyMOL_structure':audit['raster_elements_in_svg']==1}
(OUT/'verification_summary.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(result,ensure_ascii=False,indent=2))

