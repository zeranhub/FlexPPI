from pathlib import Path
import sys,json,csv,hashlib
sys.path.insert(0,r'E:\Project\FlexPPI\experiments\_sources\label_sensitivity\runtime')
sys.path.insert(0,r'E:\Project\FlexPPI\experiments\_sources\completed_results_20261012\plot_runtime')
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
P=Path(__file__).parent;R=P/'raw';out=[]
for pdb in ['4FTV','5E9D']:
    f=R/(pdb+'.cif');d=MMCIF2Dict(str(f))
    rec={'pdb_id':pdb,'source':str(f),'sha256':hashlib.sha256(f.read_bytes()).hexdigest(),'struct_title':d.get('_struct.title'),'entity_polymers':[],'chains':[]}
    for i,e in enumerate(d['_entity_poly.entity_id']):
        rec['entity_polymers'].append({'entity_id':e,'author_chains':d['_entity_poly.pdbx_strand_id'][i],'sequence':d['_entity_poly.pdbx_seq_one_letter_code_can'][i].replace('\n',''),'description':dict(zip(d['_entity.id'],d['_entity.pdbx_description'])).get(e)})
    groups={}
    rows=[]
    for i,ch in enumerate(d['_atom_site.auth_asym_id']):
        if d['_atom_site.pdbx_PDB_model_num'][i]!='1':continue
        if d['_atom_site.group_PDB'][i]!='ATOM':continue
        alt=d['_atom_site.label_alt_id'][i]
        if alt not in ['.','A','1']:continue
        key=(ch,d['_atom_site.auth_seq_id'][i],d['_atom_site.pdbx_PDB_ins_code'][i])
        if key not in groups:
            groups[key]={'pdb_id':pdb,'chain':ch,'author_residue_number':key[1],'insertion_code':'' if key[2] in ['?','.'] else key[2],'label_seq_id':d['_atom_site.label_seq_id'][i],'label_asym_id':d['_atom_site.label_asym_id'][i],'entity_id':d['_atom_site.label_entity_id'][i],'resname':d['_atom_site.label_comp_id'][i],'atoms':[]}
            rows.append(groups[key])
        groups[key]['atoms'].append(d['_atom_site.label_atom_id'][i])
    for ch in dict.fromkeys(r['chain'] for r in rows):
        x=[r for r in rows if r['chain']==ch]
        rec['chains'].append({'author_chain':ch,'label_chain':x[0]['label_asym_id'],'entity_id':x[0]['entity_id'],'n_observed_residues':len(x),'first_res':x[0],'last_res':x[-1]})
    out.append(rec)
    with (P/(pdb+'_residue_author_mapping.csv')).open('w',encoding='utf-8-sig',newline='') as h:
        w=csv.DictWriter(h,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
(P/'parent_structures.json').write_text(json.dumps(out,indent=2),encoding='utf8')
for r in out:
    print(r['pdb_id'],r['struct_title']);print('entities',r['entity_polymers']);print('chains',[(x['author_chain'],x['entity_id'],x['n_observed_residues']) for x in r['chains']])
