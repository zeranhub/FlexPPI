# Academic Figure Skill: native-vector layout refinement; no data recomputation.
# a/b/c/d/e label grid is Arial Bold 11pt black. All axes remain physically fixed.
from pathlib import Path
import hashlib
import json
import shutil
import fitz
import numpy as np

OUT=Path(__file__).resolve().parent
ARCH=OUT/'source_archive'
SOURCE=OUT.parent/'figure2_title_spacing_20261024'/'fig2_title_spacing.pdf'
ARCH.mkdir(exist_ok=True)
FROZEN=ARCH/'fig2_title_spacing.pdf'
if not FROZEN.exists():shutil.copy2(SOURCE,FROZEN)
PT=72/25.4
W,H=183.0,192.0
BAR_SHIFT=5.5
# All letters share this exact distance above their group's axis top.
LABEL_TO_AXIS=5.766
LETTER_BASELINES=dict(a=11-LABEL_TO_AXIS,b=11-LABEL_TO_AXIS,
    c=87-LABEL_TO_AXIS,d=87-LABEL_TO_AXIS,e=149.505-LABEL_TO_AXIS)
LETTER_X=dict(a=4.,b=94.,c=4.,d=65.,e=4.)

doc=fitz.open(FROZEN);page=doc[0]
streams={x[0]:doc.xref_stream(x[0]) for x in page.get_xobjects() if doc.xref_stream(x[0])}
original_streams=dict(streams)
changes=[]
for xref,content in streams.items():
    changed=content
    # Remove old individual letters, retaining their native surrounding stream.
    for letter in 'abcd':
        token=b'[ ('+letter.encode('utf-16-be')+b') ] TJ'
        if token in changed:
            changes.append(dict(xref=xref,operation='remove_old_letter',letter=letter,count=changed.count(token)))
            changed=changed.replace(token,b'[ () ] TJ')
    if b'/relayoutArialBold' in changed:
        assert changed.count(b'[<0048>]TJ')==1
        changed=changed.replace(b'[<0048>]TJ',b'[<>]TJ')
        changes.append(dict(xref=xref,operation='remove_old_e'))
    # Existing top titles descend 1 mm to equalize label-to-axis spacing.
    if b'387.4961860236 cm' in changed:
        for x in ['42.5196850394','297.6377952756']:
            old=f'1 0 -0 1 {x} 387.4961860236 cm'.encode()
            new=f'1 0 -0 1 {x} {387.4961860236-PT:.10f} cm'.encode()
            assert changed.count(old)==1
            changed=changed.replace(old,new)
        changes.append(dict(xref=xref,operation='top_title_down_1mm'))
    if len(content)<10000 and '1 interface'.encode('utf-16-be') in content:
        # After the 5.5 mm bottom-group shift, the 3-key legend and e letter
        # share a common baseline; only the legend's placement changes.
        old_legend_baseline=151.67
        legend_up=(old_legend_baseline-BAR_SHIFT)-LETTER_BASELINES['e']
        changed=(f'q\n1 0 0 1 0 {legend_up*PT:.10f} cm\n'.encode()+changed+b'\nQ\n')
        changes.append(dict(xref=xref,operation='legend_up',mm=legend_up))
    if changed!=content:doc.update_stream(xref,changed)

def rect_mm(values):return fitz.Rect(*(v*PT for v in values))
# Tight source rectangles do not overlap after placement, so their native white
# background cannot erase adjacent axis-label text.
PARTS=[
    dict(name='top_a_b',clip=[0,0,183,74.5],shift=[0,0]),
    dict(name='c',clip=[0,76.5,61.5,136.5],shift=[0,0]),
    dict(name='d_PR',clip=[70.5,76.5,125.5,136.5],shift=[-5,0]),
    dict(name='d_ROC',clip=[126.5,76.5,182,136.5],shift=[-3,0]),
    dict(name='e_bars',clip=[0,146,183,197.5],shift=[0,-BAR_SHIFT]),
]
out=fitz.open();p=out.new_page(width=W*PT,height=H*PT)
for part in PARTS:
    source=part['clip'];dx,dy=part['shift']
    target=[source[0]+dx,source[1]+dy,source[2]+dx,source[3]+dy]
    p.show_pdf_page(rect_mm(target),doc,0,clip=rect_mm(source))
    part['destination']=target
p.insert_font(fontname='figureArialBold',fontfile=r'C:\Windows\Fonts\arialbd.ttf')
for letter in 'abcde':
    p.insert_text((LETTER_X[letter]*PT,LETTER_BASELINES[letter]*PT),letter,
        fontname='figureArialBold',fontsize=11,color=(0,0,0))
out.subset_fonts()
out.save(OUT/'fig2_compact_aligned.pdf',garbage=4,deflate=True)
out.close();doc.close()
result=fitz.open(OUT/'fig2_compact_aligned.pdf');p=result[0]
p.get_pixmap(dpi=450,alpha=False).save(OUT/'fig2_compact_aligned.png')
p.get_pixmap(dpi=160,alpha=False).save(OUT/'fig2_compact_aligned_preview.png')
svg=p.get_svg_image(text_as_path=False)
(OUT/'fig2_compact_aligned.svg').write_text(svg,encoding='utf-8')
assert not p.get_images(full=True) and '<image' not in svg
ss=[s for b in p.get_text('dict')['blocks'] if 'lines' in b for l in b['lines'] for s in l['spans']]
labels={}
for letter in 'abcde':
    found=[s for s in ss if s['text']==letter]
    assert len(found)==1,(letter,len(found))
    s=found[0]
    assert abs(s['size']-11)<1e-5 and s['color']==0
    labels[letter]=dict(origin_mm=[v/PT for v in s['origin']],bbox_mm=[v/PT for v in s['bbox']],font=s['font'],size=s['size'],color=s['color'])
assert np.isclose(labels['a']['origin_mm'][0],labels['c']['origin_mm'][0])
assert np.isclose(labels['a']['origin_mm'][0],labels['e']['origin_mm'][0])
assert np.isclose(labels['a']['origin_mm'][1],labels['b']['origin_mm'][1])
assert np.isclose(labels['c']['origin_mm'][1],labels['d']['origin_mm'][1])

def unique(name):
    f=[s for s in ss if s['text']==name]
    assert len(f)==1,(name,len(f))
    return fitz.Rect(f[0]['bbox'])
assert not unique('Full vs MD').intersects(unique('MD Hχ (nats)'))
# Actual content-to-content gaps, including the unified third-row letter.
first_bottom=max(s['bbox'][3]/PT for s in ss if s['bbox'][1]/PT<75)
second_top=min(s['bbox'][1]/PT for s in ss if 75<s['bbox'][1]/PT<137)
second_bottom=max(s['bbox'][3]/PT for s in ss if 75<s['bbox'][1]/PT<137)
third_top=min(s['bbox'][1]/PT for s in ss if s['bbox'][1]/PT>137)
last_bottom=max(s['bbox'][3]/PT for s in ss)
audit=dict(status='numeric_passed_independent_visual_review_pending',width_mm=W,height_mm=H,
    source=str(SOURCE),source_sha256=hashlib.sha256(FROZEN.read_bytes()).hexdigest(),
    part_rectangles_mm=PARTS,labels=labels,label_to_axis_baseline_distance_mm=LABEL_TO_AXIS,
    measured_vertical_gaps_mm=dict(row1_to_row2=second_top-first_bottom,row2_to_row3=third_top-second_bottom,bottom_margin=H-last_bottom),
    row_content_bounds_mm=[[min(s['bbox'][1]/PT for s in ss),first_bottom],[second_top,second_bottom],[third_top,last_bottom]],
    axes_positions_mm=dict(a=[30,11,48,48],b=[105,11,71,48],c=[15,87,42,40],d_PR=[76,87,42,40],d_ROC=[134,87,42,40]),
    axes_horizontal_gaps_mm=dict(c_to_PR_before=24,c_to_PR_after=19,PR_to_ROC_before=14,PR_to_ROC_after=16),
    title_to_MD_title_horizontal_clearance_mm=(unique('MD Hχ (nats)').x0-unique('Full vs MD').x1)/PT,
    source_stream_changes=changes,no_data_or_statistics_recomputed=True,
    axes_sizes_fonts_colors_markers_and_data_paths_unchanged=True,raster_images_pdf=0,raster_images_svg=0)
(OUT/'layout_QA.json').write_text(json.dumps(audit,indent=2,ensure_ascii=False),encoding='utf-8')
(OUT/'source_manifest.json').write_text(json.dumps([dict(file=str(FROZEN.relative_to(OUT)),source=str(SOURCE),sha256=audit['source_sha256'])],indent=2),encoding='utf-8')
print(json.dumps({k:audit[k] for k in ['status','width_mm','height_mm','measured_vertical_gaps_mm','row_content_bounds_mm','labels']},indent=2))
