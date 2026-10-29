"""Prepare auditable CaM/peptide structural illustrations, without model inference.

Only the first deposited MODEL is retained. CaM residues are aligned by sequence
to the 148-residue mature CaM sequence in all four raw BMRB input audits. The
output CaM chain is A, numbered by that sequence; partner chain is P, retaining
its original author residue numbers. Only original protein atoms and Ca2+ ions
assigned to the original CaM chain are retained; no atom is modelled or rebuilt.
"""
from pathlib import Path
from collections import OrderedDict, Counter
import csv
import hashlib
import json
import math

OUT = Path(__file__).resolve().parent / 'structure_data'
RAW = Path(r'E:\Project\FlexPPI\artifacts\apo_holo_exploration_20260913\raw')
ROOT = Path(__file__).resolve().parent.parent
AA = dict(zip('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split(), 'ARNDCQEGHILKMFPSTWYV'))
AA['MSE']='M'  # Sequence-equivalent only; original Se chemistry is retained.
HEAVY = {
 'ALA':'N CA C O CB', 'ARG':'N CA C O CB CG CD NE CZ NH1 NH2',
 'ASN':'N CA C O CB CG OD1 ND2', 'ASP':'N CA C O CB CG OD1 OD2',
 'CYS':'N CA C O CB SG', 'GLN':'N CA C O CB CG CD OE1 NE2',
 'GLU':'N CA C O CB CG CD OE1 OE2', 'GLY':'N CA C O',
 'HIS':'N CA C O CB CG ND1 CD2 CE1 NE2', 'ILE':'N CA C O CB CG1 CG2 CD1',
 'LEU':'N CA C O CB CG CD1 CD2', 'LYS':'N CA C O CB CG CD CE NZ',
 'MET':'N CA C O CB CG SD CE', 'PHE':'N CA C O CB CG CD1 CD2 CE1 CE2 CZ',
 'PRO':'N CA C O CB CG CD', 'SER':'N CA C O CB OG',
 'THR':'N CA C O CB OG1 CG2', 'TRP':'N CA C O CB CG CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2',
 'TYR':'N CA C O CB CG CD1 CD2 CE1 CE2 CZ OH', 'VAL':'N CA C O CB CG1 CG2'}
HEAVY['MSE']='N CA C O CB CG SE CE'
CASES = [('CaMKK', '1CKK', 'A', 'B', 15183), ('eNOS','1NIW','A','B',15185),
         ('nNOS','2O60','A','B',15191), ('smMLCK','1CDL','A','E',15186)]
KEY = {18:'L',39:'L',105:'L',109:'M',144:'M'}

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def write_csv(name, rows):
    with (OUT / name).open('w', newline='', encoding='utf-8-sig') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)

def align(ref, observed):
    """Needleman-Wunsch global alignment; no author-number assumption."""
    n,m = len(ref), len(observed)
    score = [[0]*(m+1) for _ in range(n+1)]
    pointer = [[0]*(m+1) for _ in range(n+1)]
    for i in range(1,n+1): score[i][0],pointer[i][0] = -4*i,1
    for j in range(1,m+1): score[0][j],pointer[0][j] = -4*j,2
    for i in range(1,n+1):
        for j in range(1,m+1):
            opts = [score[i-1][j-1]+(3 if ref[i-1]==observed[j-1] else -5), score[i-1][j]-4, score[i][j-1]-4]
            best=max(opts); score[i][j]=best; pointer[i][j]=opts.index(best)
    i,j=n,m; pairs=[]; ar=[]; ao=[]
    while i or j:
        step=pointer[i][j]
        if step==0:
            pairs.append((i-1,j-1)); ar.append(ref[i-1]); ao.append(observed[j-1]); i-=1;j-=1
        elif step==1: ar.append(ref[i-1]);ao.append('-');i-=1
        else: ar.append('-');ao.append(observed[j-1]);j-=1
    return dict(reversed(pairs)), ''.join(reversed(ar)), ''.join(reversed(ao)),score[n][m]

def parse_first_model(path):
    text=path.read_text(); lines=text.splitlines()
    model_count=sum(line.startswith('MODEL ') for line in lines)
    started=False; atoms=[]; model_id=1
    for line in lines:
        if line.startswith('MODEL '):
            if started: break
            started=True;model_id=int(line[10:14]);continue
        if line.startswith('ENDMDL') and started: break
        if line[:6] not in ('ATOM  ','HETATM'): continue
        atom=dict(line=line, record=line[:6].strip(),name=line[12:16].strip(),altloc=line[16],resname=line[17:20].strip(),chain=line[21],
                  author_residue_id=int(line[22:26]),insertion_code=line[26].strip(),
                  x=float(line[30:38]),y=float(line[38:46]),z=float(line[46:54]),
                  occupancy=float(line[54:60]),element=line[76:78].strip())
        if not atom['element']:
            atom['element']='H' if atom['name'].lstrip('0123456789').startswith('H') else atom['name'][0]
        atoms.append(atom)
    return atoms, dict(model_count=model_count or 1, selected_model_id=model_id,
        expdta=' '.join(line[10:].strip() for line in lines if line.startswith('EXPDTA')),
        title=' '.join(line[10:].strip() for line in lines if line.startswith('TITLE ')),
        compnd=' '.join(line[10:].strip() for line in lines if line.startswith('COMPND')),
        dbref=[line.strip() for line in lines if line.startswith('DBREF ')],
        seqadv=[line.strip() for line in lines if line.startswith('SEQADV')])

def residues(atoms,chain):
    groups=OrderedDict()
    for a in atoms:
        if a['chain']==chain and a['resname'] in AA and (a['record']=='ATOM' or a['resname']=='MSE'):
            key=(a['author_residue_id'],a['insertion_code'],a['resname'])
            groups.setdefault(key,[]).append(a)
    return groups

def choose_altloc(group):
    # Resolve each named atom by occupancy, blank preferred, A preferred thereafter.
    by_name=OrderedDict()
    for a in group:
        old=by_name.get(a['name'])
        def rank(x): return (x['occupancy'],x['altloc']==' ',x['altloc']=='A')
        if old is None or rank(a)>rank(old): by_name[a['name']]=a
    return list(by_name.values())

def dist(a,b):
    return math.sqrt(sum((a[k]-b[k])**2 for k in ('x','y','z')))

def format_atom(a,serial,chain,resid,icode=''):
    line=a['line'].ljust(80)
    return f"{line[:6]}{serial:5d}{line[11:16]} {line[17:21]}{chain}{resid:4d}{(icode or ' '):1s}{line[27:]}"

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    nmr_audit=json.loads((ROOT/'raw_nmr_input_audit.json').read_text(encoding='utf-8'))
    seqs={x['label']:x['cam_sequence'] for x in nmr_audit['inputs']}
    assert len(set(seqs.values()))==1
    reference=next(iter(seqs.values())); assert len(reference)==148
    gate={r['label']:r for r in csv.DictReader((ROOT/'four_partners_peptide_match_gate.csv').open(encoding='utf-8-sig'))}
    allmap=[]; keyrows=[]; audit=[]
    for partner,pdb,cam_chain,pep_chain,bmrb in CASES:
        path=RAW/(pdb+'.pdb');atoms,meta=parse_first_model(path)
        cr=residues(atoms,cam_chain);pr=residues(atoms,pep_chain)
        cam_list=list(cr); pep_list=list(pr)
        observed=''.join(AA[k[2]] for k in cam_list)
        pep_seq=''.join(AA[k[2]] for k in pep_list)
        mapping,ar,ao,score=align(reference,observed)
        assert len(set(mapping.values()))==len(observed), f'{pdb}: unmapped structure residues'
        mismatches=[dict(bmrb=i+1,expected=reference[i],observed=observed[j]) for i,j in mapping.items() if reference[i]!=observed[j]]
        assert not mismatches, f'{pdb}: sequence mismatch {mismatches}'
        inv={j:i+1 for i,j in mapping.items()}
        peptide_ref=gate[partner]['peptide_sequence_bmrb']
        pm,pref,pobs,pscore=align(peptide_ref,pep_seq)
        peptide_mismatches=[dict(nmr_sequence_position=i+1,nmr_residue=peptide_ref[i],pdb_residue=pep_seq[j]) for i,j in pm.items() if peptide_ref[i]!=pep_seq[j]]
        # Peptide terminal differences are documented, never silently equated.
        clean_cam=[];clean_pep=[];sidechain_missing=[]
        for i in range(len(reference)):
            j=mapping.get(i)
            if j is None:
                row=dict(partner=partner,pdb_id=pdb,bmrb=bmrb,bmrb_residue_id=i+1,bmrb_amino_acid=reference[i],observed=False,original_chain=cam_chain,
                         original_author_residue_id='',original_insertion_code='',original_resname='',clean_chain='A',clean_residue_id=i+1,missing_heavy_atoms='',selected_atom_names='')
            else:
                key=cam_list[j]; selected=choose_altloc(cr[key]); names={a['name'] for a in selected if a['element']!='H'}
                missing=sorted(set(HEAVY[key[2]].split())-names)
                if missing: sidechain_missing.append(dict(bmrb_residue_id=i+1,resname=key[2],missing_heavy_atoms=missing))
                clean_cam.extend((a,i+1) for a in selected)
                row=dict(partner=partner,pdb_id=pdb,bmrb=bmrb,bmrb_residue_id=i+1,bmrb_amino_acid=reference[i],observed=True,original_chain=cam_chain,
                         original_author_residue_id=key[0],original_insertion_code=key[1],original_resname=key[2],clean_chain='A',clean_residue_id=i+1,
                         missing_heavy_atoms=';'.join(missing),selected_atom_names=';'.join(sorted(names)))
            allmap.append(row)
        for key,group in pr.items():
            clean_pep.extend((a,key[0],key[1]) for a in choose_altloc(group))
        pep_heavy=[a for a,_,_ in clean_pep if a['element']!='H']
        for pos,aa in KEY.items():
            assert reference[pos-1]==aa
            j=mapping.get(pos-1);assert j is not None
            key=cam_list[j]; selected=choose_altloc(cr[key]); heavy=[a for a in selected if a['element']!='H'];names={a['name'] for a in heavy}
            expected_methyl=['CD1','CD2'] if aa=='L' else ['CE']
            assert set(expected_methyl)<=names,f'{pdb} {pos}: missing methyl'
            sidechain=[a for a in heavy if a['name'] not in {'N','CA','C','O','OXT'}]
            closest=min(((dist(a,b),a,b) for a in heavy for b in pep_heavy),key=lambda v:v[0])
            side_closest=min(((dist(a,b),a,b) for a in sidechain for b in pep_heavy),key=lambda v:v[0])
            row=dict(partner=partner,pdb_id=pdb,bmrb=bmrb,bmrb_residue_id=pos,amino_acid=aa,original_chain=cam_chain,original_author_residue_id=key[0],
                     original_insertion_code=key[1],clean_chain='A',clean_residue_id=pos,methyl_atoms=';'.join(expected_methyl),methyl_atoms_complete=True,
                     residue_heavy_atoms_complete=set(HEAVY[key[2]].split())<=names,min_residue_to_peptide_heavy_distance_A=round(closest[0],4),
                     min_sidechain_to_peptide_heavy_distance_A=round(side_closest[0],4),closest_cam_atom=closest[1]['name'],
                     closest_peptide_chain=pep_chain,closest_peptide_author_residue_id=closest[2]['author_residue_id'],
                     closest_peptide_resname=closest[2]['resname'],closest_peptide_atom=closest[2]['name'],
                     CA_x=next(a['x'] for a in heavy if a['name']=='CA'),CA_y=next(a['y'] for a in heavy if a['name']=='CA'),CA_z=next(a['z'] for a in heavy if a['name']=='CA'))
            keyrows.append(row)
        calcium=[a for a in atoms if a['record']=='HETATM' and a['chain']==cam_chain and a['resname']=='CA' and a['element']=='CA']
        assert len(calcium)==4,f'{pdb}: expected four original calcium ions, found {len(calcium)}'
        pdbout=OUT/f'{pdb}_CaM_A_peptide_P_clean.pdb'
        lines=[f'HEADER    CAM PARTNER STRUCTURAL ILLUSTRATION         {pdb}',
               'REMARK   1 FIRST DEPOSITED MODEL ONLY. ORIGINAL COORDINATES; NO ATOM REBUILD.',
               'REMARK   1 CAM CHAIN A NUMBERED BY MATURE BMRB CAM SEQUENCE (ALA 1).',
               f'REMARK   1 PEPTIDE CHAIN P WAS {pep_chain}; ORIGINAL PEPTIDE AUTHOR NUMBERING.',
               'REMARK   1 CA2+ CHAIN A NUMBERED 501-504; WATER/OTHER MOLECULES REMOVED.',
               'REMARK   1 STRUCTURAL CONTEXT ONLY; NOT CLAIMED TO BE IDENTICAL NMR SAMPLE.']
        serial=0
        for a,pos in clean_cam:
            serial+=1;lines.append(format_atom(a,serial,'A',pos))
        lines.append('TER')
        for a,pos,icode in clean_pep:
            serial+=1;lines.append(format_atom(a,serial,'P',pos,icode))
        lines.append('TER')
        for k,a in enumerate(calcium,501):
            serial+=1;lines.append(format_atom(a,serial,'A',k))
        lines.append('END');pdbout.write_text('\n'.join(lines)+'\n',encoding='ascii')
        audit.append(dict(partner=partner,pdb_id=pdb,bmrb=bmrb,source_path=str(path),source_sha256=sha(path),output_path=str(pdbout),output_sha256=sha(pdbout),
                          original_cam_chain=cam_chain,original_peptide_chain=pep_chain,clean_cam_chain='A',clean_peptide_chain='P',
                          reference_sequence=reference,observed_cam_sequence=observed,alignment_reference=ar,alignment_observed=ao,alignment_score=score,
                          cam_residue_count=len(cr),cam_matched_identity=1.0,cam_missing_mature_positions=[i+1 for i in range(len(reference)) if i not in mapping],
                          cam_sequence_mismatches=mismatches,cam_missing_heavy_atom_details=sidechain_missing,
                          cam_modified_residues=[dict(original_author_residue_id=k[0],original_resname=k[2],bmrb_residue_id=inv[j],sequence_equivalent='M') for j,k in enumerate(cam_list) if k[2]=='MSE'],
                          peptide_observed_sequence=pep_seq,peptide_nmr_sequence=peptide_ref,peptide_alignment_reference=pref,peptide_alignment_observed=pobs,
                          peptide_sequence_mismatches=peptide_mismatches,peptide_mapped_count=len(pm),peptide_residue_count=len(pr),
                          retained_calcium_count=len(calcium),retained_calcium_original_author_numbers=[a['author_residue_id'] for a in calcium],
                          alternative_location_atoms_present=sum(a['altloc']!=' ' for a in atoms if a['chain'] in (cam_chain,pep_chain)),
                          clean_atom_count=serial,**meta))
    write_csv('cam_all_residue_mapping.csv',allmap)
    write_csv('cam_key_site_structure_context.csv',keyrows)
    payload=dict(method='Global full-sequence alignment to mature BMRB CaM; first model; no coordinate changes.',
                 disclaimer='Static structural illustrations use matched partner peptide structures, not the actual dynamic NMR samples. Static distances do not establish causes of S2 differences.',
                 notes=['1CDL COMPND peptide description conflicts with DBREF/sequence; DBREF MYLK_CHICK and peptide sequence establish smMLCK matching.',
                        '1NIW has nine selenomethionines (MSE). MSE is aligned to Met for sequence mapping only; the original MSE/Se atom identities and coordinates are retained, and this chemistry differs from the NMR CaM sample.',
                        '1MXE is excluded; no structures other than 1CKK,1NIW,2O60,1CDL are used.',
                        'NMR mature CaM sequence positions are retained; absent coordinates are left absent.',
                        'Alternative locations selected by maximum occupancy, then blank, then A. Hydrogen atoms are retained if deposited.'],cases=audit)
    (OUT/'structure_mapping_audit.json').write_text(json.dumps(payload,indent=2,ensure_ascii=False),encoding='utf-8')
    print(json.dumps([dict(partner=a['partner'],pdb=a['pdb_id'],n=a['cam_residue_count'],missing=a['cam_missing_mature_positions'],key_complete=True,peptide_mismatch=a['peptide_sequence_mismatches']) for a in audit],ensure_ascii=False,indent=2))

if __name__=='__main__': main()
