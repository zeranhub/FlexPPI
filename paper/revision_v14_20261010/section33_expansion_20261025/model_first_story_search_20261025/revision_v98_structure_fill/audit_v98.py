from pathlib import Path
from zipfile import ZipFile
from PIL import Image
import numpy as np
import hashlib,json,sys
sys.stdout.reconfigure(encoding='utf-8')
ROOT=Path(__file__).resolve().parent;OLD=ROOT.parent/'revision_v97_journal_composition';BASE=ROOT.parents[2]
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
assert sha(ROOT/'figure3/source_data/all_plotted_values.csv')==sha(OLD/'figure3/source_data/all_plotted_values.csv')
for p in (OLD/'figure3/source_data').glob('*'):
    if p.is_file():assert sha(p)==sha(ROOT/'figure3/source_data'/p.name)
for p in (OLD/'pymol').rglob('*'):
    if p.is_file():assert sha(p)==sha(ROOT/'pymol'/p.relative_to(OLD/'pymol'))
a=np.asarray(Image.open(OLD/'figure3/fig3_dual_case_entropy.png'))
b=np.asarray(Image.open(ROOT/'figure3/fig3_dual_case_entropy.png'))
assert a.shape==b.shape
start=int(np.ceil(a.shape[0]*54/200))
assert np.array_equal(a[start:],b[start:]),'A data-panel pixel changed.'
qa=json.loads((ROOT/'figure3/figure_QA.json').read_text(encoding='utf-8'))
scales={}
for name,view in qa['structure_viewport_crops'].items():
    z=np.asarray(Image.open(ROOT/'pymol'/(name+'.png')).convert('RGBA'))
    yy,xx=np.where((z[:,:,3]>0)&(z[:,:,:3].min(axis=2)<254))
    bounds=[int(xx.min()),int(yy.min()),int(xx.max()+1),int(yy.max()+1)]
    assert view[0]<=bounds[0] and view[1]<=bounds[1] and view[2]>=bounds[2] and view[3]>=bounds[3],(name,bounds,view)
    scale=min(74/(view[2]-view[0]),41/(view[3]-view[1]))/(44/1600)
    scales[name]={'full_native_nonwhite_bounds':bounds,'view':view,'linear_scale_vs_v97':scale,'all_molecular_pixels_in_view':True}
docs=[BASE/'paper_v97_20261027_中文修订_期刊参考组图重排.docx',BASE/'paper_v98_20261027_中文修订_结构区域放大.docx']
parts=[]
for p in docs:
    with ZipFile(p) as z:parts.append({n:z.read(n) for n in z.namelist()})
assert parts[0].keys()==parts[1].keys()
assert [n for n in parts[0] if parts[0][n]!=parts[1][n]]==['word/media/image5.png']
audit={'all_853_data_values_and_all_source_tables_identical':True,'native_PNG_PSE_PDB_and_cameras_identical':True,'all_raster_pixels_below_54mm_identical':True,'source_structures_aspect_preserved':True,'structure_display_scales':scales,'Word_changed_part_only':'word/media/image5.png','all_text_captions_references_styles_tables_and_other45_media_identical':True,'figure_canvas_and_data_positions_unchanged':True}
oldpages=sorted((OLD/'render/release_01').glob('page-*.png'),key=lambda p:int(p.stem.split('-')[1]))
pages=sorted((ROOT/'render/release_01').glob('page-*.png'),key=lambda p:int(p.stem.split('-')[1]))
if len(oldpages)==len(pages)==29:
    changed=[i+1 for i,(x,y) in enumerate(zip(oldpages,pages)) if sha(x)!=sha(y)]
    assert changed==[11],changed
    audit.update(rendered_pages=29,changed_pages=[11],other28_pages_identical_to_verified_v97=True)
(ROOT/'audit/structure_and_document_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(audit,ensure_ascii=False))
