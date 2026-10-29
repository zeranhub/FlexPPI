from pathlib import Path
import json,shutil
r=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010\section33_expansion_20261025\model_first_story_search_20261025\revision_v82_shared_view');f=r/'figure3';p=r/'pymol';a=f/'structure_assets'
old=json.loads((a/'render_manifest.json').read_text(encoding='utf-8'));new=json.loads((p/'render_manifest.json').read_text(encoding='utf-8'))
for k in ['whole_NHE1','whole_IQCG']:
    old['assets'][k]=new['assets'][k]; item=new['assets'][k];name=item.get('paired_trim',item)['file'];shutil.copy2(p/name,a/name)
old['global_display_alignment']={k:v for k,v in new.items() if k not in ['assets','palette']};old['local_display_alignment']='Unchanged local M72/M109 views: +/-5 residues, N/CA/C local backbone fitting, same camera/scale per pair.'
(a/'render_manifest.json').write_text(json.dumps(old,ensure_ascii=False,indent=2),encoding='utf-8')
for name in ['shared_view_body.py','compose_fig3.py']:
    path=f/name;s=path.read_text(encoding='utf-8');s=s.replace('x=1+i*47;text(x+22,10,','x=1+i*50;text(x+24.5,10,').replace("x,14,44,47)","x,14,49,47)").replace('text(x+22,63,','text(x+24.5,63,').replace('text(46,67.5,','text(50,67.5,')
    s=s.replace("panel('b','Local views and absolute readouts',98,2)","panel('b','Local views and absolute readouts',105,2)").replace('text(99,top,','text(106,top,').replace('98+j*21,top+10,20,18','106+j*18.5,top+10,17.5,18').replace('text(108+j*21,top+6,','text(114.75+j*18.5,top+6,').replace('text(141,y,label,6.4)','text(144,y,label,6.4)')
    path.write_text(s,encoding='utf-8')
cap=f/'caption_zh.md';cap.write_text(cap.read_text(encoding='utf-8').replace('两态输入结构，以N叶瓣','两态输入结构，显示时以N叶瓣'),encoding='utf-8')
print('Installed shared-reference PyMOL renders and widened whole-complex views.')