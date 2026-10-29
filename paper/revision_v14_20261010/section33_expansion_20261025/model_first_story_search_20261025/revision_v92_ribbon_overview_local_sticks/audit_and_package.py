from pathlib import Path
import json,hashlib,shutil,zipfile,subprocess,importlib.util,inspect,sys
from dataclasses import asdict
import pandas as pd
import numpy as np
from PIL import Image,ImageChops
ROOT=Path(__file__).resolve().parent
OLD=ROOT.parent/'revision_v91_reference_driven_layout'
BASE=ROOT.parents[2]
F=ROOT/'figure3'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
source=F/'source_data'/'all25_dynamics_and_static_lobe_contacts.csv'
assert sha(source)==sha(OLD/'figure3'/'source_data'/source.name)
df=pd.read_csv(source);pl=pd.read_csv(F/'source_data'/'all_plotted_dynamics.csv');prev=pd.read_csv(F/'source_data'/'previous_v91_curves.csv')
assert pl.equals(prev) and len(pl)==92
errors=[]
for _,q in pl.iterrows():
    row=df[df.CaMKI_native_position==q.position].iloc[0]
    val=row[('pred' if q.model=='Full' else 'intra')+'_'+q.metric]
    errors.append(abs(float(val)-q.value))
assert max(errors)==0
im=Image.open(F/'fig3_ribbon_local_sticks.png').convert('RGB')
oldim=Image.open(OLD/'figure3'/'fig3_spatial_sequence_study.png').convert('RGB')
assert im.size==oldim.size==(3242,2515)
cut=round(65.4/25.4*450)
diff=ImageChops.difference(im.crop((0,cut,im.width,im.height)),oldim.crop((0,cut,im.width,im.height)))
assert diff.getbbox() is None,'Lower profiles differ from V91'
fonts=subprocess.check_output(['C:/texlive/2026/bin/windows/pdffonts.exe',str(F/'fig3_ribbon_local_sticks.pdf')],text=True,encoding='utf-8',errors='replace')
(F/'pdf_fonts.txt').write_text(fonts,encoding='utf-8')
assert 'ArialMT' in fonts and 'yes yes yes' in fonts
qa=json.loads((F/'figure_QA.json').read_text(encoding='utf-8'))
qa.update(curves_exactly_equal_V91=True,source_CSV_exactly_equal_V91=True,lower_panel_pixel_identical_V91=True,lower_panel_comparison_from_mm=65.4,unique_curve_values=92,maximum_source_value_difference=max(errors),fonts_embedded=True,pixels=list(im.size),manuscript_modified=False,new_inference_or_MD=False)
geometry=json.loads((ROOT/'pymol'/'independent_geometry_representation_color_audit.json').read_text(encoding='utf-8'))
assert geometry['assets']['overview']['stick_atom_count']==0
for key,n in [('zoom_N297',61),('zoom_N308',125)]:
    item=geometry['assets'][key]
    assert item['sticks_represented_atom_count']==item['expected_whole_residue_heavy_atom_count']==n
    assert item['all_expected_whole_residue_atoms_in_sticks_representation'] and item['no_sticks_outside_selected_whole_residues']
    assert item['loaded_source_coordinates_exactly_equal'] and item['surface_atom_count']==0
assert geometry['both_local_views_same_rotation'] and geometry['both_local_views_same_optical_distance'] and geometry['both_local_views_same_dimensions']
qa.update(independent_representation_coordinate_color_audit='passed',all_selected_local_residue_heavy_atoms_drawn_as_sticks=True,overview_actual_stick_atom_count=0,local_stick_atom_counts={'N297':61,'N308':125})
if sha(F/'fig3_ribbon_local_sticks.png')=='c6bc19418dd5e79f52108caffc98daff282bc038791ccd948869a1fddc1d1ecc':
    qa.update(visual_review='Root inspected final exact PNG; independent dynamics reviewer inspected same layout and native structure author inspected final 2800px assets',reviewed_png_sha256=sha(F/'fig3_ribbon_local_sticks.png'))
(F/'figure_QA.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8')
spec=importlib.util.spec_from_file_location('figure_qa',Path('C:/Users/DELL/.codex/skills/academic-figure-skill/scripts/qa_validator.py'))
mod=importlib.util.module_from_spec(spec);sys.modules[spec.name]=mod;spec.loader.exec_module(mod)
style=[]
for name,func in inspect.getmembers(mod,inspect.isfunction):
    if name.startswith('check_'):
        found=func((F/'compose_fig3_1MXE.py').read_text(encoding='utf-8'))
        style.extend(asdict(x) for x in (found if isinstance(found,list) else [found]))
(F/'automated_style_QA_raw.json').write_text(json.dumps(style,ensure_ascii=False,indent=2),encoding='utf-8')
shutil.copy2(OLD/'figure3'/'descriptive_statistics.json',F/'descriptive_statistics.json')
for ext in ['png','pdf','svg']:
    shutil.copy2(F/('fig3_ribbon_local_sticks.'+ext),BASE/('figure3_v92_20261026_整体丝带与局部棍状.'+ext))
archive=BASE/'figure3_v92_20261026_源文件与PyMOL.zip'
with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    for p in sorted(ROOT.rglob('*')):
        if p.is_file() and '__pycache__' not in p.parts:z.write(p,p.relative_to(ROOT.parent))
print(json.dumps({'curves_unchanged':True,'lower_panel_pixel_identical':True,'curve_values':92,'curve_value_max_difference':max(errors),'pixels':im.size,'zip_bytes':archive.stat().st_size,'PNG_sha256':sha(F/'fig3_ribbon_local_sticks.png')},ensure_ascii=False,indent=2))
