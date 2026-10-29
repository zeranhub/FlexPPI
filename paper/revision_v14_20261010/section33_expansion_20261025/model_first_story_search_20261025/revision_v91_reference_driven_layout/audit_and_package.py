from pathlib import Path
import json, hashlib, shutil, zipfile, subprocess, inspect, importlib.util, sys
from dataclasses import asdict
import numpy as np
import pandas as pd
from PIL import Image

ROOT=Path(__file__).resolve().parent
BASE=ROOT.parents[2]
FIG=ROOT/'figure3'
DATA=FIG/'source_data'
P=ROOT/'pymol_quantitative'
REF=ROOT.parent/'figure_reference_research_20261026'
EXPECTED_REVIEWED_PNG='f95455c91400a3b5eaa13a88ca9254c8329ef57a4549d2949f538d4d773950d6'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
src=pd.read_csv(DATA/'all25_dynamics_and_static_lobe_contacts.csv')
pl=pd.read_csv(DATA/'all_plotted_dynamics.csv')
assert len(src)==25 and src.valid_chi_entropy.sum()==21 and len(pl)==92
assert not pl.duplicated(['metric','position','model']).any()
errors=[]
for _,q in pl.iterrows():
    row=src[src.CaMKI_native_position==q.position].iloc[0]
    pref='pred' if q.model=='Full' else 'intra'
    expected=float(row[pref+'_'+q.metric])
    errors.append(abs(expected-float(q.value)))
assert max(errors)<1e-12
assert set(src.loc[~src.valid_chi_entropy,'CaMKI_native_position'])=={299,306,309,311}
manifest=json.loads((P/'render_manifest.json').read_text(encoding='utf-8'))
mapping=manifest['residue_mapping']
assert len(mapping)==25
n_maps=0
for q in mapping:
    row=src[src.CaMKI_native_position==q['native_residue']].iloc[0]
    assert abs(q['RMSF']-float(row.pred_rmsf))<1e-12
    n_maps+=1
    if q['Hchi_valid']:
        assert abs(q['Hchi']-float(row.pred_chi_entropy))<1e-12
        n_maps+=1
    else:assert q['Hchi'] is None
assert n_maps==46
catalog=json.loads((REF/'catalog_all60.json').read_text(encoding='utf-8'))
assert len(catalog)==60
n_papers=len({q['paper_title'].strip().casefold() for q in catalog})
assert n_papers==34
assert all(Path(q['local_path']).is_file() and q['visually_viewed'] for q in catalog)
spec=importlib.util.spec_from_file_location('figure_qa',Path('C:/Users/DELL/.codex/skills/academic-figure-skill/scripts/qa_validator.py'))
mod=importlib.util.module_from_spec(spec);sys.modules[spec.name]=mod;spec.loader.exec_module(mod)
code=(FIG/'compose_fig3_1MXE.py').read_text(encoding='utf-8')
findings=[]
for name,func in inspect.getmembers(mod,inspect.isfunction):
    if name.startswith('check_'):
        result=func(code)
        findings.extend(asdict(x) for x in (result if isinstance(result,list) else [result]))
(FIG/'automated_style_QA_raw.json').write_text(json.dumps(findings,ensure_ascii=False,indent=2),encoding='utf-8')
fonts=subprocess.check_output(['C:/texlive/2026/bin/windows/pdffonts.exe',str(FIG/'fig3_spatial_sequence_study.pdf')],text=True,encoding='utf-8',errors='replace')
(FIG/'pdf_fonts.txt').write_text(fonts,encoding='utf-8')
assert 'ArialMT' in fonts and 'yes yes yes' in fonts
im=Image.open(FIG/'fig3_spatial_sequence_study.png')
assert abs(im.width-183/25.4*450)<1.1 and abs(im.height-142/25.4*450)<1.1
qa=json.loads((FIG/'figure_QA.json').read_text(encoding='utf-8'))
assert sha(FIG/'fig3_spatial_sequence_study.png')==EXPECTED_REVIEWED_PNG
qa.update(dict(
    visual_review='Root and independent figure reviewers inspected this exact PNG',
    reviewed_png_sha256=EXPECTED_REVIEWED_PNG,
    independent_value_checks={'curve_values':92,'structure_values':46,'max_curve_abs_difference':max(errors),'undefined_Ala':[299,306,309,311]},
    b_environment='Actual CaM atoms within 5 A of Asn sidechain heavy atoms; selected local surface; not whole pocket or SASA quantitative display',
    final_b_asset='pymol_quantitative/b_environment_candidates/1MXE_Hchi_local_surface_focus.png',
    static_geometry_and_colors_independently_audited=True,
    raster_pixels=list(im.size),
    vector_pdf_fonts_embedded=True,
    manuscript_modified=False,
    new_inference_or_MD=False
))
(FIG/'figure_QA.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8')
stats={}
core=src.CaMKI_native_position.between(303,316)
for metric in ['rmsf','chi_entropy']:
    valid=src['valid_'+metric]
    stats[metric]={label:{'count':int(valid.sum()),'minimum':float(src.loc[valid,prefix+'_'+metric].min()),'maximum':float(src.loc[valid,prefix+'_'+metric].max()),'mean':float(src.loc[valid,prefix+'_'+metric].mean())} for label,prefix in [('Full','pred'),('Intra','intra'),('MD_reference','target')]}
stats['core_303_316_Full_RMSF']={'count':14,'min':float(src.loc[core,'pred_rmsf'].min()),'max':float(src.loc[core,'pred_rmsf'].max())}
stats['inferential_tests']='Not used; individual residue outputs are not independent biological replicates.'
(FIG/'descriptive_statistics.json').write_text(json.dumps(stats,ensure_ascii=False,indent=2),encoding='utf-8')
for ext in ['png','pdf','svg']:
    shutil.copy2(FIG/('fig3_spatial_sequence_study.'+ext),BASE/('figure3_v91_20261026_定量结构与完整序列.'+ext))
dest=BASE/'figure3_v91_20261026_源文件与PyMOL.zip'
with zipfile.ZipFile(dest,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    for f in sorted(ROOT.rglob('*')):
        if f.is_file() and '__pycache__' not in f.parts:
            z.write(f,f.relative_to(ROOT.parent))
print(json.dumps({'n_reference_groups':len(catalog),'n_papers':n_papers,'unique_curve_values':len(pl),'unique_structure_values':n_maps,'maximum_curve_value_difference':max(errors),'PNG_pixels':im.size,'source_zip_bytes':dest.stat().st_size,'raw_style_findings':[q for q in findings if not q['pass_']],'delivered_png_sha256':EXPECTED_REVIEWED_PNG},ensure_ascii=False,indent=2))
