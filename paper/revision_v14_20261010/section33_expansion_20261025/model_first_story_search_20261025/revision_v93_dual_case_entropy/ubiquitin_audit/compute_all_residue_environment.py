"""Read-only SASA/contact analysis of the exact frozen A003R input structure."""
from pathlib import Path
import csv, json, math, hashlib
import pymol
from pymol import cmd

OUT=Path(__file__).resolve().parent
SOURCE=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010\section33_expansion_20261025\model_first_story_search_20261025\revision_v80_redesign\pymol\inputs\A003R_input_structure.pdb')
pymol.finish_launching(['pymol','-cq'])
cmd.load(str(SOURCE),'complex')
cmd.remove('hydro')
cmd.set('dot_solvent',1)
cmd.set('solvent_radius',1.4)
cmd.set('dot_density',4)
area=cmd.get_area('complex',load_b=1)
atoms=cmd.get_model('complex').atom
rows=[]
bb={'N','CA','C','O','OXT','OT1','OT2'}
groups={}
for a in atoms:
    groups.setdefault((a.chain,int(a.resi),a.resn),[]).append(a)
for (chain,resid,resn),group in groups.items():
    side=[a for a in group if a.name not in bb]
    other=[a for a in atoms if (a.chain,int(a.resi))!=(chain,resid)]
    row={'chain':chain,'input_residue_id':resid,'resname':resn,
         'whole_residue_SASA_A2':sum(a.b for a in group),
         'sidechain_SASA_A2':sum(a.b for a in side),'n_sidechain_heavy_atoms':len(side)}
    for d in [4.5,5.0,8.0]:
        hits={}
        pairs=0
        for a in side:
            for b in other:
                dist=math.dist(a.coord,b.coord)
                if dist<=d:
                    key=(b.chain,int(b.resi),b.resn)
                    hits[key]=min(hits.get(key,float('inf')),dist)
                    pairs+=1
        tag=str(d).replace('.','p')
        row[f'sc_residue_contacts_{tag}A_intra']=sum(k[0]==chain for k in hits)
        row[f'sc_residue_contacts_{tag}A_intra_nonlocal']=sum(k[0]==chain and abs(k[1]-resid)>1 for k in hits)
        row[f'sc_residue_contacts_{tag}A_partner']=sum(k[0]!=chain for k in hits)
        row[f'sc_atom_contact_pairs_{tag}A']=pairs
        if d==5.0:
            row['sc_contact_residue_list_5A']=';'.join(f'{k[0]}:{k[1]}:{k[2]}:{v:.3f}' for k,v in sorted(hits.items()))
    rows.append(row)
with (OUT/'all119_input_environment.csv').open('w',newline='',encoding='utf-8-sig') as f:
    w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
meta={'source_structure':str(SOURCE),'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
      'software':str(cmd.get_version()),'n_heavy_atoms':len(atoms),'n_residues':len(rows),
      'SASA_probe_A':1.4,'SASA_dot_density':4,'total_area_A2':area,
      'sidechain_definition':'all residue heavy atoms except N,CA,C,O,OXT,OT1,OT2',
      'contact_definition':'any sidechain heavy atom within cutoff of a different residue heavy atom; same-chain nonlocal additionally excludes +/- 1 sequence neighbor',
      'scope':'Exact single input structure; not MD averaged; no source mutation or model rerun'}
(OUT/'environment_provenance.json').write_text(json.dumps(meta,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(meta,ensure_ascii=False))
cmd.quit()
