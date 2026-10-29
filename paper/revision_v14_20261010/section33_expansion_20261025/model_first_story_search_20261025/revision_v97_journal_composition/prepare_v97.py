from pathlib import Path
import shutil
ROOT=Path(__file__).resolve().parent
OLD=ROOT.parent/'revision_v96_visual_hierarchy'
for name in ['figure3/source_data','pymol','render','audit','manuscript']:
    (ROOT/name).mkdir(parents=True,exist_ok=True)
for p in (OLD/'figure3/source_data').glob('*'):
    if p.is_file():shutil.copy2(p,ROOT/'figure3/source_data'/p.name)
for p in (OLD/'pymol').rglob('*'):
    if p.is_file():
        out=ROOT/'pymol'/p.relative_to(OLD/'pymol');out.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,out)
for name in ['render_docx_word_wrapper.py','export_word_readonly.ps1']:
    shutil.copy2(OLD/'render'/name,ROOT/'render'/name)
old=(OLD/'figure3/compose_fig3_v96.py').read_text(encoding='utf-8')
baseline=old[old.index('# Academic Figure Skill Typography Baseline'):old.index('# Production inheritance:')]
body=(ROOT/'compose_body_v97.py').read_text(encoding='utf-8')
(ROOT/'figure3/compose_fig3_v97.py').write_text('''# Academic Figure Skill Asset Confirmation (verified against existing native production scripts)
# a -> original v96 PyMOL overview PNG/PSE -> native run, no image editing
# b/c -> v96 absolute residue profiles / LineTrend -> param inherit
# d/e/f -> v96 complete scatter and association plots -> param inherit
# All model and MD arrays are frozen; layout, typography, and annotations only.
'''+baseline+body,encoding='utf-8')
print('Prepared v97 from v96; all frozen tables and intact native overviews copied.')
