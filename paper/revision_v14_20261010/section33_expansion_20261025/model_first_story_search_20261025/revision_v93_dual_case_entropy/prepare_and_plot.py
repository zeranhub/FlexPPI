from pathlib import Path
import shutil,subprocess,sys
ROOT=Path(__file__).resolve().parent;S=ROOT.parent
F=ROOT/'figure3';P=F/'source_data';P.mkdir(parents=True,exist_ok=True)
shutil.copy2(S/'revision_v90_data_led_story/figure3/source_data/all25_dynamics_and_static_lobe_contacts.csv',P/'cam_all25.csv')
for f in ['ubiquitin_72_position_profile_with_original_mask.csv','ubiquitin_all6_Gln_Hchi_SASA_plot.csv']:
    shutil.copy2(ROOT/'ubiquitin_audit'/f,P/f)
for f in ['retained_absolute_correlations_all498.csv','retained_prediction_groups_all498.csv','retained_panel_statistics.csv']:
    shutil.copy2(ROOT/'entropy_audit'/f,P/f)
source=(S/'revision_v90_data_led_story/figure3/compose_fig3_1MXE.py').read_text(encoding='utf-8')
head=source.split('# Inherited audited full-profile plotting;')[0]
assert 'save_cns_figure' in head and 'font.sans-serif' in head
(F/'compose_fig3_v93.py').write_text(head+(ROOT/'compose_main_body.py').read_text(encoding='utf-8'),encoding='utf-8')
subprocess.run([sys.executable,str(F/'compose_fig3_v93.py')],check=True)
