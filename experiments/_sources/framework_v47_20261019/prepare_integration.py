"""Reuse preservation/render machinery; the v46 graphic itself stays rejected."""
from pathlib import Path
import shutil
ROOT=Path(__file__).resolve().parent
PRIOR=ROOT.parent/'framework_v46_20261019'
for name in ['build_v46.py','audit_package.py','export_v46_word.ps1','render_v46.py','finalize_v46.py']:
    out=ROOT/name.replace('v46','v47')
    source=(PRIOR/name).read_text(encoding='utf8').replace('v46','v47')
    out.write_text(source,encoding='utf8')
(ROOT/'audit').mkdir(exist_ok=True)
shutil.copy2(PRIOR/'audit/prose_clarifications.json',ROOT/'audit/prose_clarifications.json')
contract=(PRIOR/'artifact.md').read_text(encoding='utf8')
contract+='\nVisual revision v47: the user rejected the flattened v46 graphic. Keep the original Figure 1 overall composition and substantive component hierarchy, and refine the real 3D protein/atom renders, component geometry, typography, colors and connectors. v46 remains a draft; v47 starts from the accepted v45 document.\n'
(ROOT/'artifact.md').write_text(contract,encoding='utf8')
print('Prepared v47 integration scripts, preserving accepted v45 source.')
