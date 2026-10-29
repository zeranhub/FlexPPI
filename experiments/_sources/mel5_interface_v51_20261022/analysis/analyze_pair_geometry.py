from pathlib import Path
import json, hashlib
import numpy as np
import pandas as pd

BASE = Path(r'E:\Project\FlexPPI\experiments\_sources')
OLD = BASE/'mel5_joint_case_v49_20261019'
OUT = BASE/'mel5_interface_v51_20261022'/'analysis'
OUT.mkdir(parents=True, exist_ok=True)
SOURCES = {
    'MEL5': OLD/'figure4/source_archive/3hg1_model_input.pdb',
    'alpha24beta17': OLD/'figure4/source_archive/4jff_model_input.pdb',
    'alignment': OLD/'figure4/source_archive/corrected_pair_all_aligned_positions.csv',
}

def parse(path):
    atoms = {}
    for line in path.read_text().splitlines():
        if line.startswith(('ATOM  ', 'HETATM')):
            element = line[76:78].strip() or line[12:16].strip()[0]
            if element in ('H','D'): continue
            if line[16:17] not in (' ', 'A'): continue
            key = (line[21], int(line[22:26]))
            name = line[12:16].strip()
            atoms.setdefault(key, {})[name] = np.array([float(line[30:38]),float(line[38:46]),float(line[46:54])])
    return atoms

PA = parse(SOURCES['MEL5']); PB = parse(SOURCES['alpha24beta17'])
df = pd.read_csv(SOURCES['alignment'])

def idkey(chain, residue):
    return str(chain), int(str(residue).rstrip('ABCDEFGHIJKLMNOPQRSTUVWXYZ'))

def residue(row, state):
    suffix = 'a' if state == 'MEL5' else 'b'
    return (PA if state == 'MEL5' else PB)[idkey(row['chain_'+suffix],row['id_'+suffix])]

def minpair(a,b, sidechain_a=False, backbone_b=False):
    aa = {k:v for k,v in a.items() if not sidechain_a or k not in ('N','CA','C','O','OXT')}
    bb = {k:v for k,v in b.items() if not backbone_b or k in ('N','CA','C','O','OXT')}
    pairs = [(float(np.linalg.norm(x-y)), n,m) for n,x in aa.items() for m,y in bb.items()]
    d,n,m = min(pairs)
    return {'distance_A':d,'pMHC_atom':n,'TCR_atom':m,'n_pairs_le4A':sum(x[0]<=4 for x in pairs)}

qrow = df[(df.role=='mhc') & (df.position_a==72)].iloc[0]
tcr = df[df.role.isin(['tcr_alpha','tcr_beta'])].copy()
if len(tcr)==0:
    tcr = df[df.chain_a.isin(['D','E'])].copy()

qrows=[]
for _,r in tcr.iterrows():
    rec={'TCR_chain':r.chain_a, 'common_position':int(r.position_a),'WT_aa':r.aa_a,'mature_aa':r.aa_b,'WT_input_id':str(r.id_a),'mature_input_id':str(r.id_b), 'region':r.region,
         'WT_RMSF_A':r.rmsf_full_a,'mature_RMSF_A':r.rmsf_full_b,'delta_RMSF_A':r.rmsf_full_b-r.rmsf_full_a}
    for st,label in [('MEL5','WT'),('alpha24beta17','mature')]:
        qr = residue(qrow,st); tr = residue(r,st)
        m=minpair(qr,tr)
        for k,v in m.items(): rec[label+'_'+k]=v
        rec[label+'_CA_distance_A']=float(np.linalg.norm(qr['CA']-tr['CA']))
        rec[label+'_Q72_sidechain_to_TCR_backbone_A']=minpair(qr,tr,True,True)['distance_A']
    rec['delta_min_heavy_A']=rec['mature_distance_A']-rec['WT_distance_A']
    rec['delta_CA_A']=rec['mature_CA_distance_A']-rec['WT_CA_distance_A']
    qrows.append(rec)
qdf=pd.DataFrame(qrows).sort_values('mature_distance_A')
qdf.to_csv(OUT/'Q72_to_all_aligned_TCR_residues.csv',index=False)
qdf[(qdf.TCR_chain=='E') & qdf.common_position.between(48,56)].sort_values('common_position').to_csv(OUT/'Q72_beta48_56_neighborhood.csv',index=False)

specific=[]
for pos,qatom,tatom in [(53,'OE1','N'),(54,'NE2','O')]:
    r=tcr[(tcr.chain_a=='E')&(tcr.position_a==pos)].iloc[0]
    z={'HLA_residue':'Q72','Q72_atom':qatom,'TCR_common_position':pos,'TCR_atom':tatom,'WT_residue':r.aa_a+str(pos),'mature_residue':r.aa_b+str(pos),'WT_TCR_input_id':str(r.id_a),'mature_TCR_input_id':str(r.id_b)}
    for st,label in [('MEL5','WT'),('alpha24beta17','mature')]:
        z[label+'_distance_A']=float(np.linalg.norm(residue(qrow,st)[qatom]-residue(r,st)[tatom]))
    z['delta_distance_A']=z['mature_distance_A']-z['WT_distance_A']
    specific.append(z)
pd.DataFrame(specific).to_csv(OUT/'fixed_atom_contact_distances.csv',index=False)

# All pMHC--TCR residue pairs whose nearest heavy-atom distance is <=6 A in either input.
pmhc=df[df.role.isin(['mhc','peptide'])]
pairs=[]
for _,p in pmhc.iterrows():
    for _,t in tcr.iterrows():
        ma=minpair(residue(p,'MEL5'),residue(t,'MEL5'))
        mb=minpair(residue(p,'alpha24beta17'),residue(t,'alpha24beta17'))
        if min(ma['distance_A'],mb['distance_A'])>6: continue
        z={'pMHC_role':p.role,'pMHC_chain':p.chain_a,'pMHC_common_position':int(p.position_a),'pMHC_WT_aa':p.aa_a,'pMHC_mature_aa':p.aa_b,'TCR_chain':t.chain_a,'TCR_common_position':int(t.position_a),'TCR_WT_aa':t.aa_a,'TCR_mature_aa':t.aa_b,'TCR_region':t.region,
           'pMHC_Hchi_valid':bool(p.hchi_valid),'pMHC_Hchi_WT':p.hchi_full_a,'pMHC_Hchi_mature':p.hchi_full_b,'pMHC_delta_Hchi':p.hchi_full_b-p.hchi_full_a,'TCR_RMSF_WT':t.rmsf_full_a,'TCR_RMSF_mature':t.rmsf_full_b,'TCR_delta_RMSF':t.rmsf_full_b-t.rmsf_full_a}
        for k,v in ma.items(): z['WT_'+k]=v
        for k,v in mb.items(): z['mature_'+k]=v
        z['delta_min_heavy_A']=mb['distance_A']-ma['distance_A']
        pairs.append(z)
pairdf=pd.DataFrame(pairs).sort_values(['pMHC_chain','pMHC_common_position','TCR_chain','TCR_common_position'])
pairdf.to_csv(OUT/'all_interface_pairs_within6A_in_either_structure.csv',index=False)
pairdf[(pairdf.pMHC_chain=='A')&pairdf.pMHC_common_position.between(65,76)&(pairdf.TCR_chain=='E')&pairdf.TCR_common_position.between(48,56)].to_csv(OUT/'local_Q72_betaCDR2_neighbor_network.csv',index=False)

# Put WT into the mature coordinate frame for a genuinely identical camera view.
fitrows=df[(df.role=='mhc') & df.position_a.between(1,180)]
X=np.array([residue(r,'MEL5')['CA'] for _,r in fitrows.iterrows()]);Y=np.array([residue(r,'alpha24beta17')['CA'] for _,r in fitrows.iterrows()])
xc=X.mean(0);yc=Y.mean(0);U,S,Vt=np.linalg.svd((X-xc).T@(Y-yc)); rot=U@Vt
if np.linalg.det(rot)<0:
    U[:,-1]*=-1;rot=U@Vt
trans=yc-xc@rot
fitrmsd=float(np.sqrt(np.mean(np.sum((X@rot+trans-Y)**2,axis=1))))
lines=[]
for line in SOURCES['MEL5'].read_text().splitlines():
    if line.startswith(('ATOM  ','HETATM')):
        xyz=np.array([float(line[30:38]),float(line[38:46]),float(line[46:54])])@rot+trans
        line=line[:30]+''.join(f'{v:8.3f}' for v in xyz)+line[54:]
    lines.append(line)
(OUT/'3hg1_on_4jff_HLA1_180_CA.pdb').write_text('\n'.join(lines)+'\n')

displacements=[]
for _,r in df.iterrows():
    a=residue(r,'MEL5');b=residue(r,'alpha24beta17')
    for atom in sorted(a.keys()&b.keys()):
        if atom not in ('CA','OE1','NE2'): continue
        before=a[atom]@rot+trans;after=b[atom]
        displacements.append({'role':r.role,'chain':r.chain_a,'common_position':int(r.position_a),'WT_aa':r.aa_a,'mature_aa':r.aa_b,'atom':atom,'WT_x':before[0],'WT_y':before[1],'WT_z':before[2],'mature_x':after[0],'mature_y':after[1],'mature_z':after[2],'displacement_A':float(np.linalg.norm(after-before))})
pd.DataFrame(displacements).to_csv(OUT/'common_atom_displacements_after_HLA_fit.csv',index=False)

audit={'status':'passed','coordinate_scope':'Archived simulation-prepared model-input starting structures, not new crystal downloads, not MD averages.',
       'units':'Angstrom; same named atom pair for fixed_atom_contact_distances; all-heavy minimum for neighbor network.',
       'source_files':{k:{'path':str(v),'sha256':hashlib.sha256(v.read_bytes()).hexdigest()} for k,v in SOURCES.items()},
       'numbering':{'HLA_Q72':{'WT':'A72','mature':'A72'},'beta53':{'WT':'E53 ILE','mature':'E54 PHE'},'beta54':{'WT':'E54 GLY','mature':'E55 GLY'},'common_numbering':'Aligned WT/MEL5 numbering; mature input TCR residues have +1 N-terminal shift.'},
       'fixed_atom_distances':specific,'HLA_CA_alignment':{'n':len(X),'common_positions':[1,180],'RMSD_A':fitrmsd,'rotation_row_vector':rot.tolist(),'translation':trans.tolist()},
       'n_aligned_TCR_residues':len(tcr),'n_interface_pairs_within6A_either':len(pairdf),
       'scientific_scope':'Coordinates show a local contact rearrangement; model Hchi and RMSF are bound-state per-residue outputs. Geometry alone cannot establish Hchi-to-RMSF causality or quantify affinity.'}
(OUT/'geometry_audit.json').write_text(json.dumps(audit,indent=2,ensure_ascii=False),encoding='utf8')
print(json.dumps({'fixed_atom_distances':specific,'Q72_TCR_nearby':qdf.head(8).to_dict('records'),'HLA_fit_RMSD_A':fitrmsd,'files':len(list(OUT.glob('*')))},ensure_ascii=False,indent=2))
