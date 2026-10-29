from pathlib import Path
import numpy as np,json,csv
S=Path(r'E:/Project/FlexPPI/paper/revision_v14_20261010/section33_expansion_20261025/model_first_story_search_20261025')
out=S/'revision_v89_mechanism_story/pymol_mechanism';out.mkdir(parents=True,exist_ok=True)
def read(path):
 a=[]
 for l in path.read_text().splitlines():
  if l[:6].strip()=='ATOM' and l[76:78].strip()!='H':
   a.append({'chain':l[21],'resi':int(l[22:26]),'resn':l[17:20].strip(),'name':l[12:16].strip(),'xyz':np.array([float(l[i:i+8]) for i in (30,38,46)])})
 return a
inp=S/'revision_v88_two_panel_story/source_analysis/M316_static_context/inputs'
r={};rows=[]
for code,fn,ch in [('1A06','1A06_rigidly_aligned_to_1MXE.pdb','A'),('1MXE','1MXE_CaM_native4_147_CaMKI294_318.pdb','B')]:
 aa=read(inp/fn);context=[a for a in aa if (a['chain']=='A' and a['resi']<294) if code=='1A06'] if code=='1A06' else [a for a in aa if a['chain']=='A']
 r[code]={}
 for site in [296,298,303,316]:
  atoms=[a for a in aa if a['chain']==ch and a['resi']==site]
  rr={};r[code][site]={'atom_names':[a['name'] for a in atoms]}
  for kind in ['all_heavy','sidechain']:
   at=atoms if kind=='all_heavy' else [a for a in atoms if a['name'] not in ['N','CA','C','O']]
   for b in context:
    dd=np.linalg.norm(np.array([a['xyz'] for a in at])-b['xyz'],axis=1);i=dd.argmin();d=float(dd[i]);key=(b['chain'],b['resi'],b['resn'])
    if key not in rr or d<rr[key][0]: rr[key]=(d,at[i]['name'],b['name'])
   ordered=sorted(rr.items(),key=lambda kv:kv[1][0]);r[code][site][kind]=[{'chain':k[0],'resi':k[1],'resn':k[2],'distance_A':round(v[0],4),'site_atom':v[1],'context_atom':v[2]} for k,v in ordered[:10]]
   for k,v in ordered:
    rows.append({'structure':code,'site':site,'scope':kind,'chain':k[0],'residue':k[1],'resname':k[2],'min_distance_A':v[0],'site_atom':v[1],'context_atom':v[2]})
   rr={}
(out/'K296_F298_W303_M316_contacts.json').write_text(json.dumps(r,indent=2),encoding='utf-8')
with (out/'all_context_min_distances.csv').open('w',newline='',encoding='utf-8') as f:
 w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
print(json.dumps(r,indent=2))