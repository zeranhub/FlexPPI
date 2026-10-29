from pathlib import Path
from collections import OrderedDict
import urllib.request,json,sys,hashlib
ROOT=Path(__file__).resolve().parent;OUT=ROOT/'catchbond_validation';OUT.mkdir(exist_ok=True)
sys.path.insert(0,str(ROOT.parent/'label_sensitivity/runtime'))
from Bio.Align import PairwiseAligner
AA=dict(zip('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split(),'ARNDCQEGHILKMFPSTWYV'))
AA.update(HID='H',HIE='H',HIP='H',HSD='H',HSE='H',HSP='H',CYX='C',MSE='M')
aliases={'HID':'HIS','HIE':'HIS','HIP':'HIS','HSD':'HIS','HSE':'HIS','HSP':'HIS','CYX':'CYS','MSE':'MET'}
def parse(path):
    chains=OrderedDict()
    for l in path.read_text().splitlines():
        if l.startswith('ENDMDL'):break
        if l.startswith(('ATOM  ','HETATM')) and l[17:20].strip() in AA and l[16:17] in (' ','A'):
            ch=l[21];key=l[22:27].strip();chains.setdefault(ch,OrderedDict()).setdefault(key,[]).append(l)
    return {ch:[(key,AA[ls[0][17:20].strip()],ls) for key,ls in rr.items() if any(l[12:16].strip()=='CA' for l in ls)] for ch,rr in chains.items()}
al=PairwiseAligner(mode='global',match_score=2,mismatch_score=-2,open_gap_score=-6,extend_gap_score=-.5)
al.target_end_gap_score=0;al.query_end_gap_score=0
audit=[]
for cond,pdb in [('WT','9nmu'),('S32M','9nmx')]:
    raw=OUT/(pdb+'_experimental_full.pdb')
    if not raw.exists():raw.write_bytes(urllib.request.urlopen(f'https://files.rcsb.org/download/{pdb.upper()}.pdb',timeout=60).read())
    prepared=ROOT/'other_data_inventory'/f'catchbond_{cond}_prepared_protein.pdb'
    reference=parse(prepared);source=parse(raw)
    print(cond,[(k,len(v),''.join(r[1] for r in v)[:15]) for k,v in source.items()])
    cropped=[];mapping=[];claimed=set()
    for target,rr in reference.items():
        query=''.join(x[1] for x in rr);options=[]
        for ch,ss in source.items():
            if ch in claimed:continue
            seq=''.join(x[1] for x in ss)
            a=al.align(query,seq)[0];ij=[(int(i),int(j)) for i,j in zip(*a.indices) if i>=0 and j>=0]
            identical=sum(query[i]==seq[j] for i,j in ij);coverage=len(ij)/len(query);identity=identical/len(ij) if ij else 0
            options.append((identity*coverage,identity,coverage,ch,ij))
        score,identity,coverage,ch,ij=max(options)
        assert identity>=.99 and coverage>=.95,(cond,target,ch,identity,coverage)
        claimed.add(ch)
        for i,j in ij:
            key,aa,lines=source[ch][j];target_key=rr[i][0]
            for l in lines:
                if not l[12:16].strip().lstrip('0123456789').startswith(('H','D')):
                    rn=aliases.get(l[17:20],l[17:20]);new=l[:16]+' '+rn+l[20:21]+target+f'{int(target_key):4d}'+' '+l[27:]
                    cropped.append(new)
            mapping.append({'condition':cond,'target_chain':target,'target_residue_id':target_key,'source_chain':ch,'source_residue_id':key,'aa':aa,'identity':identity,'coverage':coverage})
        cropped.append('TER')
    path=OUT/f'catchbond_{cond}_experimental_matched.pdb';path.write_text('\n'.join(cropped)+'\nEND\n')
    import shutil
    dest=OUT/f'catchbond_{cond}_prepared.pdb';shutil.copy2(prepared,dest)
    audit.append({'condition':cond,'PDB':pdb,'source_sha256':hashlib.sha256(raw.read_bytes()).hexdigest(),'experimental_matched_pdb':str(path),'prepared_pdb':str(dest),'residue_mapping':mapping})
(OUT/'structure_input_audit.json').write_text(json.dumps(audit,indent=2),encoding='utf8')
