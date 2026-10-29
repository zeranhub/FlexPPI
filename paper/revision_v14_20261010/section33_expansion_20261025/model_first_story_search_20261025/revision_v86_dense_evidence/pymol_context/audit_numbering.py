from pathlib import Path
import json,hashlib
out=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010\section33_expansion_20261025\model_first_story_search_20261025\revision_v86_dense_evidence\pymol_context')
p=out/'inputs/1A06_raw_official.pdb';lines=p.read_text().splitlines();atoms=[]
for s in lines:
 if s.startswith('ATOM  '):atoms.append({'chain':s[21],'resn':s[17:20],'resi':int(s[22:26]),'icode':s[26],'name':s[12:16].strip(),'alt':s[16],'xyz':[float(s[k:k+8]) for k in [30,38,46]]})
ca={a['resi']:a for a in atoms if a['chain']=='A' and a['name']=='CA'}
seqmap={'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLU':'E','GLN':'Q','GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
obs={n:ca[n]['resn'] if n in ca else None for n in range(294,319)}
needed={n:['N','CA','C','O']+({'ASN':['CB','CG','OD1','ND2'],'TRP':['CB','CG','CD1','CD2','NE1','CE2','CE3','CZ2','CZ3','CH2'],'MET':['CB','CG','SD','CE']}.get(ca.get(n,{}).get('resn',''),[])) for n in [297,303,308,316]}
audit={'pdb_id':'1A06','download_url':'https://files.rcsb.org/download/1A06.pdb','source_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'observed_chain_A_residue_count':len(ca),'observed_range':[min(ca),max(ca)],'observed_residues_294_318':obs,'observed_294_318_sequence':''.join(seqmap.get(obs[n],'?') for n in range(294,319)),'expected_1MXE_peptide_sequence':'IKKNFAKSKWKQAFNATAVVRHMRK','missing_residues_in_294_318':[n for n in range(294,319) if n not in ca],'key_atoms_missing':{n:[name for name in names if not any(a['chain']=='A' and a['resi']==n and a['name']==name for a in atoms)] for n,names in needed.items()}}
(out/'native_numbering_audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8');print(json.dumps(audit,indent=2))