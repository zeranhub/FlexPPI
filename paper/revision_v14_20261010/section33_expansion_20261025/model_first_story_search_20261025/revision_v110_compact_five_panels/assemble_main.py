from pathlib import Path
import shutil
root=Path(__file__).resolve().parent;old=root.parent/'revision_v109_case_SASA_tracks';out=root/'figure3';out.mkdir(exist_ok=True)
for folder in ['figure3/source_data','pymol','source_inputs','sasa_audit']:
    shutil.copytree(old/folder,root/folder,dirs_exist_ok=True)
for name in ['case_profiles.py','plot_overlay.py','case_profiles_all97_with_masks.csv','case_profile_manifest.json',
             'case_SASA_all97_starting_complex.csv','standalone_all368_plotted_values.csv',
             'same_residue_MD_and_Full_all824988.npz','direct_relationship_statistics.json','cohort_all3063_records.csv']:
    shutil.copy2(old/name,root/name)
shutil.copy2(old/'figure3/publication_style.py',out/'publication_style.py')
source=(old/'figure3/compose_fig3_v109.py').read_text(encoding='utf-8')
source=source.replace('W,H=183.,238.','W,H=183.,181.')
source=source.replace("panel('d','Same side chain, different environments',98,62)","")
start=source.index("ordered=gln.sort_values('MD_Hchi_nats')")
end=source.index('# e: all four angle comparisons',start)
source=source[:start]+'''# d: compact single-axis RMSF/Hchi distribution, formerly panel f.
from plot_overlay import draw_overlay,load_data
panel('d','RMSF–Hχ relationship',98,62)
density_audit=draw_overlay(fig,(W,H),(111,80,60,31),xmax=10,gridsize=(50,35),marginal_mm=4.)
direct_arrays,direct_stats=load_data()

'''+source[end:]
start=source.index('# f: one-axis overlay')
end=source.index('df=pd.DataFrame(plotted)',start)
source=source[:start]+source[end:]
source=source.replace("{'b':117,'c':348,'d':18,'e':302}","{'b':117,'c':348,'e':302}")
source=source.replace("'a_d_e_f_data_unchanged_from_v108':True","'a_b_c_e_data_unchanged_from_v109':True,'old_d_removed':True,'old_f_relocated_to_d':True,'density_axis_box_mm':[111,80,60,31]")
source=source.replace("'f':'single-axis overlay of MD and Full RMSF-Hchi densities and marginals'","'d':'compact single-axis overlay of MD and Full RMSF-Hchi densities and marginals'")
source=source.replace("'common_f_positions':824988","'common_d_positions':824988").replace("'f_records':3063","'d_records':3063")
(out/'compose_fig3_v110.py').write_text(source,encoding='utf-8')
print('Removed old d; old f moved to compact right column as d; 5 panels,183x181mm.')
