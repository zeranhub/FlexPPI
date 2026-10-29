from pathlib import Path
import csv, json, re
import xml.etree.ElementTree as ET
import numpy as np
from Bio.PDB import PDBParser

BASE = Path(r'E:\Project\FlexPPI\experiments\_sources')
OUT = BASE/'opposite_affinity_rmsf_cases_20261022'/'spr_mechanism'
OLD = BASE/'local_antigen_screening_20261021'
DIRECT = BASE/'local_peptide_max_rmsf_drop_20261022'/'direct_affinity'

def txt(e):
    return re.sub(r'\s+', ' ', ''.join(e.itertext())).strip()

for name,path in [('piepenbrink2013',OUT/'piepenbrink2013_fullText.xml'),
                  ('pettmann2021', OLD/'experimental_labels'/'raw'/'pettmann2021_fullText.xml'),
                  ('aleksic2010',OLD/'experimental_labels'/'raw'/'aleksic2010_fullText.xml')]:
    root=ET.parse(path).getroot()
    paras=[]
    for e in root.iter():
        if e.tag in ['p','title','ref','fig','table-wrap']:
            s=txt(e)
            if len(s)>15 and re.search(r'G4A|L3A|L3Y|Thr98|Asp99|Arg65|Bmax|SPR|37.?C|temperature|C-terminal|6.9|4G4|2BNQ|2BNR|glycine at|position 4|3QH3',s,re.I):
                paras.append({'tag':e.tag,'id':e.get('id'), 'text':s})
    (OUT/f'{name}_selected_source_passages.json').write_text(json.dumps(paras,ensure_ascii=False,indent=2),encoding='utf-8')
    print(name, 'selected',len(paras))
    if name=='piepenbrink2013':
        for p in paras:
            if any(t in p['text'] for t in ['Thr98','position 4','glycine at']): print(p['id'],p['text'][:4500])

rows=list(csv.DictReader((DIRECT/'scores_and_original_study_KD_all_versions.csv').open(encoding='utf-8-sig')))
sel=[r for r in rows if r['primary_template']=='True' and r['SPR_variant']=='reported_hybrid' and r['mutation'] in ['G4A','L3A','L3Y']]
for r in sel:
    print('CASE',json.dumps({k:r[k] for k in ['panel_id','mutation','study_id','KD_uM','parent_KD_uM','error_value','error_type','n','delta_pKD_vs_parent','max_drop_positions','score_max_rmsf_drop_A','primary_locator']},ensure_ascii=False))

templates=json.loads((OLD/'structures'/'fixed_backbone_manifest.json').read_text(encoding='utf-8'))['templates']
parser=PDBParser(QUIET=True)
contacts=[]
for t in templates:
    if t['template_pdb'] not in ['1AO7','2BNQ']:continue
    m=parser.get_structure(t['template_pdb'],t['raw_pdb'])[0]
    roles=t['roles']; pc=m[roles['peptide']]
    print('RESIDS',t['template_pdb'],[(c.id,[(r.id,r.resname) for r in c if r.id[0]==' ' and r.id[1] in [65,98,99]]) for c in m if c.id in ['A','D']])
    for res in pc:
        if res.id[0]!=' ':continue
        atoms=[a for a in res if a.element not in ['H','D']]
        for role in ['hla_heavy','tcr_alpha','tcr_beta']:
            ch=m[roles[role]]
            for partner in ch:
                if partner.id[0]!=' ':continue
                pa=[a for a in partner if a.element not in ['H','D']]
                distances=[(float(np.linalg.norm(a.coord-b.coord)),a.name,b.name) for a in atoms for b in pa]
                d,an,bn=min(distances)
                if d<=4.5:
                    contacts.append({'template':t['template_pdb'],'peptide_chain':pc.id,'peptide_pos':res.id[1],'peptide_aa':res.resname,'partner_role':role,'partner_chain':ch.id,'partner_resid':str(partner.id[1])+partner.id[2].strip(),'partner_aa':partner.resname,'min_distance_A':d,'peptide_atom':an,'partner_atom':bn,'contact_definition':'minimum heavy-atom distance <=4.5 A; parent crystal structure only'})
    if t['template_pdb']=='1AO7':
        for rn in [98,99]:
            rr=m['D'][rn]
            ar=m['A'][65]
            for a in rr:
                if a.element not in ['N','O']:continue
                for b in ar:
                    if b.element not in ['N','O']:continue
                    d=float(np.linalg.norm(a.coord-b.coord))
                    if d<=3.5:print('POLAR_GEOMETRY',rn,rr.resname,a.name,'Arg65',b.name,d)

with (OUT/'parent_crystal_peptide_contacts.csv').open('w',encoding='utf-8-sig',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(contacts[0]));w.writeheader();w.writerows(contacts)
for r in contacts:
    if (r['template']=='1AO7' and r['peptide_pos']==4) or (r['template']=='2BNQ' and r['peptide_pos'] in [1,3]):print('CONTACT',r)
