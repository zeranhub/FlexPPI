from pathlib import Path
import shutil
root=Path(__file__).resolve().parent;old=root.parent/'revision_v108_dual_axis_case_profiles';out=root/'figure3';out.mkdir(exist_ok=True)
shutil.copytree(old/'figure3/source_data',out/'source_data',dirs_exist_ok=True)
shutil.copytree(old/'pymol',root/'pymol',dirs_exist_ok=True)
shutil.copy2(old/'figure3/publication_style.py',out/'publication_style.py')
shutil.copy2(old/'plot_overlay.py',root/'plot_overlay.py')
source=(old/'figure3/compose_fig3_v108.py').read_text(encoding='utf-8')
source=source.replace('W,H=183.,230.','W,H=183.,238.')
source=source.replace('from case_profiles import draw_case,method_handles','from case_profiles import draw_case,method_handles,draw_sasa')
source=source.replace("panel('b','CaMKI dynamics and shielding'","panel('b','CaMKI dynamics and solvent exposure'")
start=source.index("a=ax(111,44,60,9);style_axis(a,'b_static_shielding')")
end=source.index('# c/d: complete sequence accuracy',start)
replacement="""a,sasa_values=draw_sasa(fig,(W,H),(111,44,60,9),'CaMKI')
values('b','static_sidechain_SASA_A2',sasa_values.residue_id,sasa_values.sidechain_SASA_A2,series='starting_structure')

"""
source=source[:start]+replacement+source[end:]
source=source.replace("(14,74,69,30),'Ubiquitin')","(14,74,69,22),'Ubiquitin',xlabel=False)")
needle="ordered=gln.sort_values('MD_Hchi_nats')"
source=source.replace(needle,"""a,sasa_values=draw_sasa(fig,(W,H),(14,100,69,9),'Ubiquitin')
values('c','static_sidechain_SASA_A2',sasa_values.residue_id,sasa_values.sidechain_SASA_A2,series='starting_structure')
"""+needle)
# Expand the second row by8mm, preserving relative positions within e/f.
replacements={
 "panel('e','Per-χ occupancy association',3,115)":"panel('e','Per-χ occupancy association',3,123)",
 '1-118/H':'1-126/H','a=ax(left,132,35.5,25)':'a=ax(left,140,35.5,25)',
 'txt(left,126.5,':'txt(left,134.5,','txt(3.1,144.5,':'txt(3.1,152.5,','txt(97,163,':'txt(97,171,',
 "panel('f','RMSF–Hχ relationship',3,173)":"panel('f','RMSF–Hχ relationship',3,181)",
 '(14,190,157,30)':'(14,198,157,30)',
 "{'b':117,'c':276,'d':18,'e':302}":"{'b':117,'c':348,'d':18,'e':302}",
 "'a_d_e_and_f_unchanged_from_v107':True":"'a_d_e_f_data_unchanged_from_v108':True,'b_c_static_SASA_tracks':True"
}
for a,b in replacements.items():
    assert a in source,a
    source=source.replace(a,b)
(out/'compose_fig3_v109.py').write_text(source,encoding='utf-8')
print('Prepared SASA tracks below both case profiles; b shielding replaced, c SASA added.')
