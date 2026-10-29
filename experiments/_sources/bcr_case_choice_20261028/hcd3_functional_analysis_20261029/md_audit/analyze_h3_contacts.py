"""Exploratory functional audit of the complete H-CDR3, using existing MD only.

Before trajectory results: all 15 H3 positions, all protein residue partners,
heavy-atom distances, 4.0/4.5/5.0 A thresholds; no affinity/hotspot claims.
Frames are correlated samples of one trajectory, never independent replicates.
"""
from pathlib import Path
import csv, hashlib, json
import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array

ROOT = Path(r'E:\Project\FlexPPI\experiments\_sources\bcr_case_choice_20261028')
OUT = ROOT / 'hcd3_functional_analysis_20261029/md_audit'
PDB = Path(r'E:\Project\FlexPPI\experiments\core_corrected_pipeline\training_payload\dataset\structures\A00LZ.pdb')
XTC = Path(r'D:\data\dynarepo_non_nucleic_full\A00LZ\trajectory.xtc')
OUT.mkdir(parents=True, exist_ok=True)
THRESHOLDS = (4.0, 4.5, 5.0)
contract = dict(scope='all 15 H-CDR3 positions, authors 95..109, chain C',
                input=str(PDB), trajectory=str(XTC), thresholds_A=THRESHOLDS,
                pairs='all protein residue partners; same-chain |author difference|<=2 excluded from nonlocal summaries',
                contacts='minimum non-hydrogen atom distance with periodic box',
                fitting=['whole heavy chain CA (original label)', 'numbered VH CA', 'numbered VH framework CA', 'H-CDR3 CA itself'],
                inference='exploratory structural correspondence; no binding free energy, hotspot, causal or replicate-level inference')
(OUT/'analysis_scope.json').write_text(json.dumps(contract, ensure_ascii=False, indent=2), encoding='utf-8')
u = mda.Universe(str(PDB), str(XTC))
residues = list(u.residues)
keys = [(r.segid, int(r.resid), r.resname) for r in residues]
h3 = [i for i,k in enumerate(keys) if k[0]=='C' and 95<=k[1]<=109]
assert len(h3)==15 and len(residues)==519
heavy = u.atoms[np.char.upper(u.atoms.elements.astype(str))!='H']
starts = np.flatnonzero(np.r_[True, np.diff(heavy.resindices)!=0])
assert np.array_equal(heavy.resindices[starts], np.arange(519))
h3_atoms = [r.atoms[np.char.upper(r.atoms.elements.astype(str))!='H'] for r in (residues[i] for i in h3)]
distances = np.empty((len(u.trajectory),len(h3),len(residues)), dtype=np.float32)
ca = u.select_atoms('name CA')
assert len(ca)==519 and np.array_equal(ca.resindices,np.arange(519))
coords = np.empty((len(u.trajectory),519,3), dtype=np.float32)
for frame,ts in enumerate(u.trajectory):
    coords[frame] = ca.positions
    for j,atoms in enumerate(h3_atoms):
        d = distance_array(atoms.positions, heavy.positions, box=ts.dimensions).min(axis=0)
        distances[frame,j] = np.minimum.reduceat(d,starts)
    if frame%250==0:
        print('processed',frame,flush=True)
np.savez_compressed(OUT/'all_h3_partner_minimum_distances.npz', minimum_distances_A=distances,
                    h3_residue_indices=np.array(h3), chains=np.array([k[0] for k in keys]),
                    authors=np.array([k[1] for k in keys]), resnames=np.array([k[2] for k in keys]))
pred = list(csv.DictReader((ROOT/'6A77_all_positions_predictions_and_MD.csv').open(encoding='utf-8-sig')))
full = {int(r['author']):r for r in pred if r['variant']=='full' and r['chain']=='C'}
numbering = list(csv.DictReader((ROOT/'6A77_antibody_numbering.csv').open(encoding='utf-8-sig')))
vh_authors = {int(r['author']) for r in numbering if r['chain']=='C' and r['imgt_number'].strip()}
fw_authors = {int(r['author']) for r in numbering if r['chain']=='C' and r['imgt_number'].strip() and r['region']=='framework'}
fit_indices = {
    'whole_heavy_chain': [i for i,k in enumerate(keys) if k[0]=='C'],
    'VH': [i for i,k in enumerate(keys) if k[0]=='C' and k[1] in vh_authors],
    'VH_framework': [i for i,k in enumerate(keys) if k[0]=='C' and k[1] in fw_authors],
    'HCDR3_itself': h3,
}
fit_rmsf = {}
for fit_name,fit_idx in fit_indices.items():
    ref = coords[0,fit_idx].astype(float)
    ref_center = ref.mean(axis=0)
    ref_zero = ref-ref_center
    aligned = np.empty((len(coords),len(h3),3),dtype=np.float32)
    for frame in range(len(coords)):
        mobile = coords[frame,fit_idx].astype(float)
        center = mobile.mean(axis=0)
        left,_,right = np.linalg.svd((mobile-center).T @ ref_zero)
        rotation = left @ right
        if np.linalg.det(rotation)<0:
            left[:,-1]*=-1
            rotation=left@right
        aligned[frame] = (coords[frame,h3]-center)@rotation+ref_center
    fit_rmsf[fit_name] = np.sqrt(np.mean(np.sum((aligned-aligned.mean(axis=0))**2,axis=-1),axis=0))
rows=[]
pairs=[]
chain = np.array([k[0] for k in keys])
authors = np.array([k[1] for k in keys])
for j,idx in enumerate(h3):
    key=keys[idx]
    row={'chain':'C','author':key[1],'resname':key[2]}
    row.update({f'RMSF_fit_{name}':float(values[j]) for name,values in fit_rmsf.items()})
    for name,selector in [('antigen',chain=='A'),('light',chain=='B'),('heavy_nonlocal',(chain=='C')&(np.abs(authors-key[1])>2))]:
        closest = distances[:,j,selector].min(axis=1)
        row[f'{name}_mean_nearest_A']=float(closest.mean())
        for threshold in THRESHOLDS:
            row[f'{name}_contact_fraction_{threshold:g}A']=float(np.mean(closest<=threshold))
    for k in ['rmsf_pred','rmsf_MD','chi_entropy_pred','chi_entropy_MD','chi_entropy_valid']:
        row[k]=full[key[1]][k]
    rows.append(row)
    for partner,pkey in enumerate(keys):
        if pkey[0]=='C' and abs(pkey[1]-key[1])<=2:
            continue
        values=distances[:,j,partner]
        pair={'h3_author':key[1],'h3_resname':key[2],'partner_chain':pkey[0],
              'partner_author':pkey[1],'partner_resname':pkey[2],
              'minimum_A':float(values.min()),'mean_A':float(values.mean())}
        for threshold in THRESHOLDS:
            pair[f'fraction_{threshold:g}A']=float(np.mean(values<=threshold))
        pairs.append(pair)
for name,records in [('all_h3_dynamic_context.csv',rows),('all_h3_nonlocal_partner_occupancies.csv',pairs)]:
    with (OUT/name).open('w',encoding='utf-8-sig',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(records[0]));writer.writeheader();writer.writerows(records)
max_error=max(abs(row['RMSF_fit_whole_heavy_chain']-float(row['rmsf_MD'])) for row in rows)
audit={'frames':len(coords),'residues':519,'h3_residues':15,'atom_count':len(u.atoms),
       'heavy_atoms':len(heavy),'fit_position_counts':{k:len(v) for k,v in fit_indices.items()},
       'max_RMSF_difference_from_frozen_label_A':max_error,
       'sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [PDB,XTC,Path(__file__)]},
       'scope_and_limitations':contract,'runtime':{'MDAnalysis':mda.__version__,'numpy':np.__version__}}
(OUT/'audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
assert max_error<5e-5
print(json.dumps(audit,indent=2),flush=True)
for row in rows:
    print(row,flush=True)
