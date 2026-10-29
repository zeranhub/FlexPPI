from pathlib import Path
import os
ENV=Path('D:/Anaconda/envs/pvb_cpu');os.environ['PATH']=str(ENV/'Library/bin')+os.pathsep+os.environ.get('PATH','');H=[os.add_dll_directory(str(p)) for p in [ENV/'Library/bin',ENV/'DLLs'] if p.exists()]
import pymol2
pse=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010\section33_expansion_20261025\model_first_story_search_20261025\revision_v92_ribbon_overview_local_sticks\pymol\zoom_N297.pse')
with pymol2.PyMOL() as p:
 c=p.cmd;c.load(str(pse));ss=c.get_session()
 print('top',list(ss.keys()));print('names records',[(x[0],len(x)) for x in ss['names'] if x])
 n=c.count_atoms('bound');obj=[x for x in ss['names'] if x and x[0]=='bound'][0]
 def find(v,path):
  if isinstance(v,list):
   if len(v)==n:print('NATOMS LIST',path,'firsttypes',[type(a).__name__ for a in v[:3]],'firstlens',[len(a) if isinstance(a,(list,tuple)) else -1 for a in v[:3]])
   for i,x in enumerate(v):
    if isinstance(x,list):find(x,path+[i])
 find(obj,[])
 try:print('rep sticks selector',c.count_atoms('bound and rep sticks'))
 except Exception as e:print('rep selector unavailable',str(e))