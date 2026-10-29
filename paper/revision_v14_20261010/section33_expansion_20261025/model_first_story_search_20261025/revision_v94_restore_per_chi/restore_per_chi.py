from pathlib import Path
import shutil,json,subprocess,sys
ROOT=Path(__file__).resolve().parent;OLD=ROOT.parent/'revision_v93_dual_case_entropy'
F=ROOT/'figure3';P=F/'source_data';P.mkdir(parents=True,exist_ok=True)
for source in (OLD/'figure3/source_data').iterdir():
    if source.is_file() and source.name!='all_plotted_values.csv':shutil.copy2(source,P/source.name)
for name in ['defined_angles_long.csv','per_chi_correlations.csv']:
    shutil.copy2(OLD/'entropy_audit/supplement_source'/name,P/name)
if not (ROOT/'pymol').exists():shutil.copytree(OLD/'pymol',ROOT/'pymol')
old=(OLD/'figure3/compose_fig3_v93.py').read_text(encoding='utf-8')
head=old.split('# Production inheritance:')[0]
body=(OLD/'compose_main_body.py').read_text(encoding='utf-8')
body=body.replace('W,H=183.,225.','W,H=183.,247.')
prefix=body.split('# g: absolute type-controlled correlations;')[0]
tail=body.split("df=pd.DataFrame(plotted);")[1]
tail="df=pd.DataFrame(plotted);"+tail
tail=tail.replace("{'b':100,'c':132,'e':12,'f':332,'g':332}","{'b':100,'c':132,'e':12,'f':302,'g':332}")
middle=(ROOT/'per_chi_panel_body.py').read_text(encoding='utf-8')
(F/'compose_fig3_v94.py').write_text(head+prefix+middle+'\n'+tail,encoding='utf-8')
subprocess.run([sys.executable,str(F/'compose_fig3_v94.py')],check=True)
