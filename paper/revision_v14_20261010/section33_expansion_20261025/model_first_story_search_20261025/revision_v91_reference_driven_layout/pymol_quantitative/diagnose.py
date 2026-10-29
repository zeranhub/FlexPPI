from pathlib import Path
from PIL import Image
import numpy as np,json,os
root=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010\section33_expansion_20261025\model_first_story_search_20261025\revision_v91_reference_driven_layout\pymol_quantitative')
a=np.array(Image.open(root/'1MXE_Hchi_axis090.png'))
print('mode shape',a.shape,'alpha',np.unique(a[:,:,3],return_counts=True))
print('CaM pixel',a[200,470],a[550,600],a[390,500])
env=Path('D:/Anaconda/envs/pvb_cpu');os.environ['PATH']=str(env/'Library/bin')+os.pathsep+os.environ.get('PATH','');handles=[os.add_dll_directory(str(p)) for p in [env/'Library/bin',env/'DLLs'] if p.exists()]
import pymol2
with pymol2.PyMOL() as p:
 c=p.cmd;c.load(str(root/'1MXE_Hchi_axis090.pse'))
 print('cam cartoon_transparency',c.get('cartoon_transparency','cam'),'peptide',c.get('cartoon_transparency','peptide'))
 print('global',c.get('cartoon_transparency'),'opaque',c.get('ray_opaque_background'))