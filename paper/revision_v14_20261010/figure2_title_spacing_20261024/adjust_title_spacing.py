# Academic Figure Skill: cosmetic native-PDF adjustment; no data or plot redraw.
# a/b titles and panel letters move down 9 mm; the whole figure is then cropped
# by 9 mm at the top. Original vector paths, fonts, sizes, colors and data remain.
from pathlib import Path
import hashlib
import json
import shutil
import fitz
import numpy as np

OUT=Path(__file__).resolve().parent
ARCH=OUT/'source_archive'
SOURCE=OUT.parent/'figure2_grouped_layout_20261024'/'fig2_grouped.pdf'
PT=72/25.4
SHIFT_MM=9.0
SHIFT=SHIFT_MM*PT
ARCH.mkdir(exist_ok=True)
FROZEN=ARCH/'fig2_grouped.pdf'
if not FROZEN.exists():shutil.copy2(SOURCE,FROZEN)
source_hash=hashlib.sha256(FROZEN.read_bytes()).hexdigest()
doc=fitz.open(FROZEN)
page=doc[0]
matches=[x[0] for x in page.get_xobjects() if x[1]=='fullpage'
    and doc.xref_stream(x[0]) and b'413.1954970472 cm' in doc.xref_stream(x[0])]
assert len(matches)==1,matches
xref=matches[0]
original=doc.xref_stream(xref)
modified=original
replacements=[]
for x,y,label in [('11.3385826772','413.1954970472','a'),
    ('42.5196850394','413.0079970472','Full vs MD'),
    ('266.4566929134','413.1954970472','b'),
    ('297.6377952756','413.0079970472','Record AP distributions')]:
    old=f'1 0 -0 1 {x} {y} cm'.encode('ascii')
    new=f'1 0 -0 1 {x} {float(y)-SHIFT:.10f} cm'.encode('ascii')
    assert modified.count(old)==1,label
    modified=modified.replace(old,new)
    replacements.append(dict(label=label,old_matrix=old.decode(),new_matrix=new.decode()))
reverse=modified
for r in replacements:reverse=reverse.replace(r['new_matrix'].encode(),r['old_matrix'].encode())
assert reverse==original
# Every byte outside these four text-placement matrices is unchanged.
doc.update_stream(xref,modified)
width,height=page.rect.width,page.rect.height
clip=fitz.Rect(0,SHIFT,width,height)
out=fitz.open();p=out.new_page(width=width,height=height-SHIFT)
p.show_pdf_page(p.rect,doc,0,clip=clip)
out.save(OUT/'fig2_title_spacing.pdf',garbage=4,deflate=True)
out.close();doc.close()

result=fitz.open(OUT/'fig2_title_spacing.pdf');p=result[0]
p.get_pixmap(dpi=450,alpha=False).save(OUT/'fig2_title_spacing.png')
p.get_pixmap(dpi=160,alpha=False).save(OUT/'fig2_title_spacing_preview.png')
svg=p.get_svg_image(text_as_path=False)
(OUT/'fig2_title_spacing.svg').write_text(svg,encoding='utf-8')
assert not p.get_images(full=True) and '<image' not in svg

def spans(pg):
    return [s for b in pg.get_text('dict')['blocks'] if 'lines' in b
        for line in b['lines'] for s in line['spans']]
new_spans=spans(p)
source_doc=fitz.open(FROZEN);old_spans=spans(source_doc[0])
def span_one(ss,name):
    candidates=[s for s in ss if s['text']==name]
    assert len(candidates)==1,(name,len(candidates))
    return candidates[0]
new_a=fitz.Rect(span_one(new_spans,'Full vs MD')['bbox'])
new_b=fitz.Rect(span_one(new_spans,'Record AP distributions')['bbox'])
md=fitz.Rect(span_one(new_spans,'MD Hχ (nats)')['bbox'])
assert not new_a.intersects(md)
assert (md.x0-new_a.x1)/PT>12
old_a=fitz.Rect(span_one(old_spans,'Full vs MD')['bbox'])
old_b=fitz.Rect(span_one(old_spans,'Record AP distributions')['bbox'])
def axis_top(pg):
    candidates=[q for q in pg.get_drawings() if q['type']=='s'
        and q['width'] is not None and abs(q['width']-.6)<1e-5
        and abs(q['rect'].x0/PT-105)<.01 and q['rect'].width<.01
        and abs(q['rect'].height/PT-48)<.01]
    assert len(candidates)==1,len(candidates)
    return candidates[0]['rect'].y0
before_top=axis_top(source_doc[0]);after_top=axis_top(p)
assert abs((before_top-after_top)/PT-SHIFT_MM)<1e-4

# Verify every source XObject stream survives exactly, except the four title
# matrices intentionally edited above. No axis, data path or legend is redrawn.
expected_streams=[]
for item in source_doc[0].get_xobjects():
    content=source_doc.xref_stream(item[0])
    if content:expected_streams.append(modified if item[0]==xref else content)
new_streams=[result.xref_stream(item[0]) for item in p.get_xobjects()]
assert all(content in new_streams for content in expected_streams)
audit=dict(status='native_and_measurement_passed_visual_review_pending',
    source=str(SOURCE),source_sha256=source_hash,
    width_mm=p.rect.width/PT,height_mm=p.rect.height/PT,
    title_shift_down_mm=SHIFT_MM,top_crop_mm=SHIFT_MM,
    title_to_axes_gap_mm=dict(
        a_before=(before_top-old_a.y1)/PT,a_after=(after_top-new_a.y1)/PT,
        b_before=(before_top-old_b.y1)/PT,b_after=(after_top-new_b.y1)/PT),
    a_title_to_MD_axis_title_horizontal_clearance_mm=(md.x0-new_a.x1)/PT,
    title_boxes_nonoverlapping=True,
    original_vector_streams_preserved_except_four_text_matrices=True,
    data_or_statistics_recomputed=False,
    plot_dimensions_fonts_colors_and_legend_data_unchanged=True,
    other_panels_global_translation_up_mm=SHIFT_MM,
    raster_images_pdf=0,raster_images_svg=0,text_matrix_edits=replacements)
(OUT/'title_spacing_QA.json').write_text(json.dumps(audit,indent=2,ensure_ascii=False),encoding='utf-8')
(OUT/'source_manifest.json').write_text(json.dumps([dict(file=str(FROZEN.relative_to(OUT)),
    original_source=str(SOURCE),sha256=source_hash)],indent=2),encoding='utf-8')
print(json.dumps(audit,indent=2,ensure_ascii=False))
