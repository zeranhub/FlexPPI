"""Validate native unchanged panels and frozen rendered curve geometry."""
from pathlib import Path
import json
import fitz
import numpy as np

ROOT=Path(__file__).resolve().parent
PT=72/25.4
W,H,TOP,ROW=183.,214.44,86.,65.
new=fitz.open(ROOT/'fig2_grouped.pdf')
src=fitz.open(ROOT/'source_archive'/'original_fig2_compact.pdf')
clip=fitz.Rect(0,0,W*PT,TOP*PT)
def pixels(page,clip=None):
    p=page.get_pixmap(dpi=160,clip=clip,alpha=False)
    return np.frombuffer(p.samples,dtype=np.uint8).reshape(p.height,p.width,3)
a=pixels(new[0],clip);b=pixels(src[0],clip)
mask=np.ones(a.shape[:2],dtype=bool)
assert np.array_equal(a[mask],b[mask])

bottom=fitz.open(ROOT/'source_archive'/'original_fig2_panel_d.pdf')
ref=fitz.open();p=ref.new_page(width=W*PT,height=H*PT)
p.show_pdf_page(fitz.Rect(0,(TOP+ROW)*PT,W*PT,H*PT),bottom,0)
a=pixels(new[0]);b=pixels(p)
mask=np.zeros(a.shape[:2],dtype=bool)
# Only the legend counts and isolated panel letter changed in the bottom group.
mask[int((TOP+ROW+12)*160/25.4):]=True
mask[int((TOP+ROW+5)*160/25.4):int((TOP+ROW+13)*160/25.4),
     int(7*160/25.4):int(14*160/25.4)]=False
assert np.array_equal(a[mask],b[mask])

oldrow=fitz.open(ROOT.parent/'figure2_relayout_20261024'/'fig2_second_row.pdf')
newrow=fitz.open(ROOT/'fig2_second_row.pdf')
oldpaths=oldrow[0].get_drawings();newpaths=newrow[0].get_drawings()
def points(path):
    return np.array([[v.x,v.y] for item in path['items'] for v in item[1:]
                     if isinstance(v,fitz.Point)])
curve_checks=[]
for panel,left in [('c',15),('d_PR',76),('d_ROC',137)]:
    rect=fitz.Rect(left*PT,14*PT,(left+42)*PT,54*PT)
    rect=rect+(-.01,-.01,.01,.01)
    selected=[q for q in oldpaths if q['type']=='s' and .69<q['width']<1.7
        and rect.contains(q['rect']) and q['rect'].width>40*PT]
    # Four internal and seven external curves, with all original vertices.
    assert len(selected)==(4 if panel=='c' else 7),(panel,len(selected))
    for q in selected:
        pp=points(q);pp[:,1]-=4*PT
        if panel=='d_PR':pp[:,0]+=5*PT
        candidates=[z for z in newpaths if z['type']=='s' and z['color']==q['color']
            and abs(z['width']-q['width'])<1e-6 and len(points(z))==len(pp)]
        assert any(np.allclose(points(z),pp,atol=2e-4,rtol=0) for z in candidates),(panel,q['rect'])
    curve_checks.append(dict(panel=panel,series=len(selected),geometry='unchanged except 4mm up; PR additionally 5mm right to group with ROC'))

q=json.loads((ROOT/'figure_QA.json').read_text(encoding='utf-8'))
q['status']='READY'
q['visual_review']='Complete grouped main figure inspected. b has its independent four-method key; c retains its internal legend; external PR/ROC share d, with PR closer to ROC; bottom is e with a centered concise group legend. No curve occlusion or clipped text.'
q['unchanged_panel_verification']={
    'top_a_b_including_independent_b_legend_pixel_identical_160dpi':True,
    'bottom_bars_axes_labels_excluding_repositioned_legend_and_panel_letter_pixel_identical_160dpi':True,
    'curve_native_geometry_comparison_against_v65':curve_checks,
    'raster_images_in_pdf':len(new[0].get_images(full=True)),
    'raster_images_in_svg':(ROOT/'fig2_grouped.svg').read_text(encoding='utf-8').count('<image')}
assert q['unchanged_panel_verification']['raster_images_in_pdf']==0
assert q['unchanged_panel_verification']['raster_images_in_svg']==0
(ROOT/'figure_QA.json').write_text(json.dumps(q,indent=2,ensure_ascii=False),encoding='utf-8')
print(json.dumps(q['unchanged_panel_verification'],indent=2))
