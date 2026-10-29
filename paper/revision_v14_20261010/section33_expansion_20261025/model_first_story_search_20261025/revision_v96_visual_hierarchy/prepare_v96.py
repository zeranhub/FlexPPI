"""Inherit frozen v95 values and native structures, then apply a coherent grid."""
from pathlib import Path
import shutil
ROOT=Path(__file__).resolve().parent
OLD=ROOT.parent/'revision_v95_grouped_structures'
for folder in ['figure3/source_data','pymol','render','manuscript','audit']:
    (ROOT/folder).mkdir(parents=True,exist_ok=True)
shutil.copytree(OLD/'figure3/source_data',ROOT/'figure3/source_data',dirs_exist_ok=True)
shutil.copytree(OLD/'pymol',ROOT/'pymol',dirs_exist_ok=True)
for name in ['render_docx_word_wrapper.py','export_word_readonly.ps1']:
    shutil.copy2(OLD/'render'/name,ROOT/'render'/name)
s=(OLD/'figure3/compose_fig3_v95.py').read_text(encoding='utf-8')
s=s[:s.index('# a/b: the two intact')]
s=s.replace('W,H=183.,247.','W,H=183.,229.')
s=s.replace("'font.size':9.0,'axes.labelsize':9.0,'xtick.labelsize':9.0,'ytick.labelsize':9.0,'legend.fontsize':9.0", "'font.size':9.2,'axes.labelsize':9.2,'xtick.labelsize':9.2,'ytick.labelsize':9.2,'legend.fontsize':9.2")
s=s.replace('def txt(x,y,s,size=9,','def txt(x,y,s,size=9.2,')
s=s.replace("txt(x,y,label,11.5", "txt(x,y,label,12")
s=s.replace("title,9.6", "title,10")
s=s.replace("color=TEAL,fontsize=9,", "color=TEAL,fontsize=9.2,")
s=s.replace('# (a,b) native structure overview', '# (a) native structure overview').replace('# (c,d) absolute residue profiles','# (b,c) absolute residue profiles').replace('# (e,f,g) complete-point scatter','# (d,e,f) complete-point scatter')
s+=(ROOT/'compose_fig3_v96_body.py').read_text(encoding='utf-8')
(ROOT/'figure3/compose_fig3_v96.py').write_text(s,encoding='utf-8')
print('Prepared coherent six-panel composition with unchanged source values.')
