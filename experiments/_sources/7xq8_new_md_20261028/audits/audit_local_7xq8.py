from pathlib import Path
from collections import OrderedDict, Counter
import json, math

root=Path(__file__).parent
pdb=Path(r'E:\Project\FlexPPI\tmp\bcr\7xq8.pdb')
lines=pdb.read_text().splitlines()
chains=OrderedDict()
for line in lines:
    if line.startswith('ATOM  '):
        chain=line[21]; key=(int(line[22:26]),line[26])
        res=chains.setdefault(chain,OrderedDict()).setdefault(key,{'name':line[17:20],'atoms':{}})
        res['atoms'][line[12:16].strip()]=[float(line[30:38]),float(line[38:46]),float(line[46:54])]
missing={}
for line in lines:
    if line.startswith('REMARK 465') and len(line)>26:
        try:
            n=int(line[22:26]); c=line[19]; name=line[15:18]
        except ValueError: continue
        missing.setdefault(c,[]).append([n,line[26:27],name])

out={'source':str(pdb),'chains':{},'ssbond_records':[x for x in lines if x.startswith('SSBOND')], 'links':[x for x in lines if x.startswith('LINK')], 'hetero_residues':[]}
hetero=OrderedDict()
for line in lines:
    if line.startswith('HETATM'):
        key=(line[21],line[22:27].strip(),line[17:20]); hetero[key]=1
out['hetero_residues']=[list(k) for k in hetero]
for c,rr in chains.items():
    items=list(rr.items()); ca=[r['atoms']['CA'] for _,r in items]
    gaps=[]
    for (key1,r1),(key2,r2) in zip(items[:-1],items[1:]):
        d=math.dist(r1['atoms']['C'],r2['atoms']['N'])
        if d>2.: gaps.append({'left':list(key1),'right':list(key2),'C_N_A':d})
    cyst=[]
    for key,r in items:
        if r['name']=='CYS': cyst.append({'id':list(key),'atoms':list(r['atoms']),'CB':r['atoms'].get('CB'),'SG':r['atoms'].get('SG')})
    out['chains'][c]={'n_res':len(items),'first':list(items[0][0]),'last':list(items[-1][0]),'xyz_ca_min':[min(x[i] for x in ca) for i in range(3)], 'xyz_ca_max':[max(x[i] for x in ca) for i in range(3)], 'peptide_gaps':gaps,'remark465':missing.get(c,[]),'cysteines':cyst, 'residues_only_backbone_CB':sum(set(r['atoms'])<=set(['N','CA','C','O','CB']) and r['name'] not in ['GLY','ALA'] for _,r in items),'tail40':[{'id':list(k),'name':r['name'],'CA':r['atoms']['CA']} for k,r in items[-40:]]}
out['cysteine_CB_pairs_under_5A']=[]
allc=[(c,k,r) for c,rr in chains.items() for k,r in rr.items() if r['name']=='CYS']
for i,(c,k,r) in enumerate(allc):
    for c2,k2,r2 in allc[i+1:]:
        if 'CB' in r['atoms'] and 'CB' in r2['atoms']:
            d=math.dist(r['atoms']['CB'],r2['atoms']['CB'])
            if d<5:out['cysteine_CB_pairs_under_5A'].append({'a':[c,*k],'b':[c2,*k2],'cb_distance_A':d,'has_sg':('SG'in r['atoms'] and 'SG'in r2['atoms'])})
(root/'7xq8_local_structure_audit.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf-8')
for c,r in out['chains'].items():
    print(c, 'nres',r['n_res'],'first/last',r['first'],r['last'],'range',r['xyz_ca_min'],r['xyz_ca_max'],'gaps',r['peptide_gaps'],'backboneCB only',r['residues_only_backbone_CB'],'missing',len(r['remark465']))
print('SSBOND',len(out['ssbond_records']),'LINK',len(out['links']),'HETATM',out['hetero_residues'])
print(json.dumps(out['cysteine_CB_pairs_under_5A'],indent=2))
