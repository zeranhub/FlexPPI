from pathlib import Path
import shutil
ROOT=Path(__file__).resolve().parent
OLD=ROOT.parent/'revision_v97_journal_composition'
for name in ['figure3/source_data','pymol','render','audit','manuscript']:
    (ROOT/name).mkdir(parents=True,exist_ok=True)
for p in (OLD/'figure3/source_data').glob('*'):
    if p.is_file():shutil.copy2(p,ROOT/'figure3/source_data'/p.name)
for p in (OLD/'pymol').rglob('*'):
    if p.is_file():
        out=ROOT/'pymol'/p.relative_to(OLD/'pymol');out.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,out)
for name in ['render_docx_word_wrapper.py','export_word_readonly.ps1']:
    shutil.copy2(OLD/'render'/name,ROOT/'render'/name)
s=(OLD/'figure3/compose_fig3_v97.py').read_text(encoding='utf-8')
needle="z=manifest['assets'][key];a=ax(x,y,w,h);a.imshow(Image.open(R/z['png']).convert('RGBA'));a.axis('off')"
replace="""z=manifest['assets'][key];a=ax(x,y,w,h);a.imshow(Image.open(R/z['png']).convert('RGBA'));a.axis('off')
    # Display viewport only: keep original PNG/PSE and camera/coordinates intact.
    crop={'cam_overview':(200,110,1430,1315),'ub_overview':(90,205,1535,1195)}[key]
    a.set_xlim(crop[0],crop[2]);a.set_ylim(crop[3],crop[1]);a.set_aspect('equal')"""
assert needle in s;s=s.replace(needle,replace)
old="""txt(20,2,'CaM–CaMKI (1MXE)',FS)
txt(112,2,'Ubiquitin–UBA (2QHO)',FS)
molecule('cam_overview',17,7,44,38.5,{'W303':(.12,.61),'M316':(.99,.28)})
molecule('ub_overview',111,7,44,38.5,{'Q41':(.56,.94),'Q62':(.96,.09)})"""
new="""txt(45,2,'CaM–CaMKI (1MXE)',FS,ha='center')
txt(124,2,'Ubiquitin–UBA (2QHO)',FS,ha='center')
molecule('cam_overview',8,6,74,41,{'W303':(.14,.62),'M316':(1.12,.28)})
molecule('ub_overview',87,6,74,41,{'Q41':(.55,.95),'Q62':(.95,.12)})"""
assert old in s;s=s.replace(old,new)
needle="black_labels_and_colored_keylines=True,"
assert needle in s
s=s.replace(needle,"black_labels_and_colored_keylines=True,original_native_assets_unchanged=True,structure_viewport_crops={'cam_overview':[200,110,1430,1315],'ub_overview':[90,205,1535,1195]},structure_aspect_preserved=True,data_panel_positions_unchanged=True,")
(ROOT/'figure3/compose_fig3_v98.py').write_text(s,encoding='utf-8')
print('Prepared v98: tighten display viewports, enlarge native molecular content, preserve all data positions.')
