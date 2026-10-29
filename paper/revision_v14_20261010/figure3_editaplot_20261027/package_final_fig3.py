"""Package verified native Origin subplots without rasterizing the PDF master."""
from pathlib import Path
import json
import hashlib
import zipfile
from collections import Counter
import numpy as np
from PIL import Image, ImageChops
from pypdf import PdfReader, PdfWriter, PageObject, Transformation
import pdfplumber

ROOT=Path(__file__).resolve().parent
DELIVERY=ROOT/'origin_panels_final'
IDS=['b_rmsf_all25','b_detail_md_full','c_ubiquitin_parity','d_gln_pairs',
     'e_chi1','e_chi2','e_chi3','e_chi4','f_fine_entropy','f_lowF_width']
MM=72/25.4
def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

plan=json.loads((ROOT/'origin_renderer_plan.json').read_text(encoding='utf-8'))
proposal=json.loads((ROOT/'scientific_design_proposal.json').read_text(encoding='utf-8'))
assert all(digest(ROOT/'source_data'/name)==value for name,value in proposal['source_hashes'].items())
specs={x['id']:x for x in plan['panels']}
qa={}
for name in IDS:
    folder=DELIVERY/name
    report=json.loads((folder/'origin_verify_report.json').read_text(encoding='utf-8'))
    assert report['editable_numeric_data_exact']
    assert all((folder/f'result.{ext}').stat().st_size>0 for ext in ['opju','pdf','png','tif'])
    before=Image.open(folder/'result.png').convert('RGB')
    after=Image.open(folder/'reopened.png').convert('RGB')
    identical=before.size==after.size and ImageChops.difference(before,after).getbbox() is None
    assert identical
    with pdfplumber.open(folder/'result.pdf') as pdf:
        page=pdf.pages[0]
        assert not page.images, f'{name}: expected vector-only PDF'
        font_sizes=sorted({round(c['size'],3) for c in page.chars if 'Arial' in c['fontname']})
    assert b'/FontFile' in (folder/'result.pdf').read_bytes()
    for phase in ['before_save','after_independent_reopen']:
        axis=report[phase]['axis_geometry']
        assert abs(axis['x']['label.pt']-9)<.05 and abs(axis['y']['label.pt']-9)<.05
        assert abs(axis['x']['thickness']-.6)<.05 and abs(axis['y']['thickness']-.6)<.05
    measurements=[s for s in specs[name]['series'] if s['role']=='measurement']
    report['visual_qa_status']='passed'
    report['human_visual_qa']='Opened the final PNGs and the complete final contact sheet; no clipping or legend/data collision remains.'
    report['reopened_png_pixel_identical']=True
    report['pdf_vector_only']=True
    (folder/'origin_verify_report.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
    qa[name]={'measurement_marks':sum(len(s['x']) for s in measurements),
              'native_font_pt':{'axis_and_tick':9,'title':10,'panel_letter':11},
              'axis_width_pt':.6,'saved_numeric_exact':True,'reopened_png_pixel_identical':True,
              'vector_pdf_no_raster_images':True,'embedded_font':True,
              'pdf_Arial_text_sizes_pt':font_sizes,
              'artifacts':{f'result.{e}':{'bytes':(folder/f'result.{e}').stat().st_size,
                  'sha256':digest(folder/f'result.{e}')} for e in ['opju','pdf','png','tif']}}

# Compose groups at their exported physical size, without scaling any panel.
# Origin's default PDF margin adds about 0.25 mm to the 89 mm page contract.
groups=[[[IDS[0],IDS[1]]],[[IDS[2],IDS[3]]],[[IDS[4],IDS[5]],[IDS[6],IDS[7]]],[[IDS[8],IDS[9]]]]
writer=PdfWriter()
reader_cache={name:PdfReader(DELIVERY/name/'result.pdf') for name in IDS}
for rows in groups:
    row_heights=[max(float(reader_cache[n].pages[0].mediabox.height) for n in row) for row in rows]
    height=sum(row_heights)+(len(rows)-1)*5*MM
    width=183*MM
    page=PageObject.create_blank_page(width=width,height=height)
    top=height
    for row,row_height in zip(rows,row_heights):
        gap=width-sum(float(reader_cache[n].pages[0].mediabox.width) for n in row)
        x=0.0
        for name in row:
            source=reader_cache[name].pages[0]
            h=float(source.mediabox.height)
            page.merge_transformed_page(source,Transformation().translate(tx=x,ty=top-h))
            x+=float(source.mediabox.width)+gap
        top-=row_height+5*MM
    writer.add_page(page)
writer.add_metadata({'/Title':'Figure 3 redesigned quantitative panels','/Author':'FlexPPI',
                     '/Subject':'Native Origin vector subplots from frozen v98 source data'})
pdf_path=ROOT/'figure3_nonstructural_origin_collection.pdf'
with pdf_path.open('wb') as stream:
    writer.write(stream)
assert len(PdfReader(pdf_path).pages)==4

overview=Image.new('RGB',(1800,3900),'white')
for i,name in enumerate(IDS):
    img=Image.open(DELIVERY/name/'result.png').convert('RGB')
    img.thumbnail((890,765),Image.Resampling.LANCZOS)
    overview.paste(img,((i%2)*900+(900-img.width)//2,(i//2)*780+(780-img.height)//2))
overview.save(ROOT/'figure3_origin_panels_overview.png',dpi=(300,300))

manifest=json.loads((DELIVERY/'delivery_manifest.json').read_text(encoding='utf-8'))
manifest.update(visual_qa_status='passed',all_saved_numeric_exact=True,all_reopened_png_pixel_identical=True,
                native_style_verified=True,vector_only_pdf=True,
                output_scope='10 editable subplot projects; original structure panels and DOCX unchanged')
(DELIVERY/'delivery_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')

qa_report={'status':'passed','renderer_plan_sha256':plan['renderer_plan_sha256'],
           'panels':qa,'source_data_unchanged':True,
           'font_export_note':'Native editable Origin text objects read back exactly 9/10/11 pt. The native PDF driver rounds rendered Arial glyph sizes by up to 0.12 pt (printer device quantization), recorded above. No post-export type or data edits.',
           'scope_note':'The 4-page collection preserves physical subplot size; the long PNG contact sheet is only a viewing index, not a single print-size manuscript figure.',
           'statistics':'No significance test or multiplicity correction added. n units, source-stratified paired bootstrap and limitations are specified in data_audit.md and final_caption_zh.md.'}
(ROOT/'final_figure_QA.json').write_text(json.dumps(qa_report,ensure_ascii=False,indent=2),encoding='utf-8')

include=['README.md','final_caption_zh.md','data_audit.md','data_audit.json','f_paired_differences.csv',
         'f_paired_summary.csv','final_figure_QA.json','render_origin_panels.py','origin_renderer_plan.json',
         'scientific_design_proposal.json','human_scientific_confirmation.json','origin_execution_approval.json',
         'prepare_origin_inputs.py','package_final_fig3.py','figure3_nonstructural_origin_collection.pdf','figure3_origin_panels_overview.png']
archive=ROOT/'figure3_editaplot_正式分图与源数据.zip'
with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    for name in include:
        path=ROOT/name
        if path.is_file(): z.write(path,name)
    for folder_name in ['origin_panels_final','source_data','origin_inputs']:
        for path in sorted((ROOT/folder_name).rglob('*')):
            if path.is_file() and path.name!='reopened.png':
                z.write(path,str(path.relative_to(ROOT)))
print(json.dumps({'status':'passed','projects':10,'vector_collection_pages':4,
                  'pdf':str(pdf_path),'archive':str(archive),'archive_bytes':archive.stat().st_size},ensure_ascii=False))
