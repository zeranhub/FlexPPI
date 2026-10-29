"""Existing model-input geometry and original per-residue predictions, read-only."""
from pathlib import Path
import csv,json,hashlib
import numpy as np
ROOT=Path(__file__).resolve().parent
EXP=ROOT.parent.parent
SRC=EXP/'tcr_expanded_screen_v41_20261017/mel5_validation'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
with (SRC/'corrected_pair_all_aligned_positions.csv').open(encoding='utf8',newline='') as f:rows=list(csv.DictReader(f))
sources=[SRC/'corrected_pair_all_aligned_positions.csv']
def atoms(p):
    out=[]
    for line in p.read_text().splitlines():
        if not line.startswith('ATOM  ') or line[76:78].strip() in ['H','D']:continue
        out.append({'chain':line[21],'resid':line[22:27].strip(),'resname':line[17:20].strip(),'atom':line[12:16].strip(),
          'xyz':np.array([float(line[30:38]),float(line[38:46]),float(line[46:54])])})
    return out
def distance(aa,bb):
    d=np.linalg.norm(np.array([a['xyz'] for a in aa])[:,None,:]-np.array([b['xyz'] for b in bb])[None,:,:],axis=-1)
    i,j=np.unravel_index(np.argmin(d),d.shape)
    return {'distance_A':float(d[i,j]),'HLA_atom':aa[i]['atom'],'TCR_atom':bb[j]['atom'],
      'n_atom_pairs_le4A':int((d<=4).sum()),'n_atom_pairs_le4p5A':int((d<=4.5).sum())}
states={}
for side,pdb in [('a','3hg1'),('b','4jff')]:
    path=SRC/f'{pdb}_model_input.pdb';sources.append(path);states[side]=atoms(path)
records=[]
q=next(r for r in rows if r['role']=='mhc' and r['id_a']=='72')
for r in [r for r in rows if r['region']=='CDR2β']:
    record={'WT_input_id':r['id_a'],'derivative_input_id':r['id_b'],'WT_aa':r['aa_a'],'derivative_aa':r['aa_b'],
      'same_amino_acid':r['same_aa']=='True','WT_model_RMSF_A':float(r['rmsf_full_a']),'derivative_model_RMSF_A':float(r['rmsf_full_b']),
      'WT_existing_MD_RMSF_A':float(r['rmsf_md_a']),'derivative_existing_MD_RMSF_A':float(r['rmsf_md_b']),
      'WT_model_Q72_Hchi_nats':float(q['hchi_full_a']),'derivative_model_Q72_Hchi_nats':float(q['hchi_full_b'])}
    for side,label in [('a','WT'),('b','derivative')]:
        aa=[a for a in states[side] if a['chain']=='A' and a['resid']=='72']
        sc=[a for a in aa if a['atom'] not in ['N','CA','C','O','OXT']]
        bb=[a for a in states[side] if a['chain']==r['chain_'+side] and a['resid']==r['id_'+side]]
        backbone=[a for a in bb if a['atom'] in ['N','CA','C','O']]
        record[label+'_Q72_to_residue_all_heavy']=distance(aa,bb)
        record[label+'_Q72_sidechain_to_residue_backbone']=distance(sc,backbone)
    records.append(record)
report={'status':'passed','scope':'Complete five-position beta CDR2, correct global alignment; geometry from existing static model input only.',
 'sources':[{'path':str(p),'sha256':sha(p)} for p in sources],
 'model_observables':'Per-residue predictions of already-bound states, not time-resolved causal interventions.',
 'records':records}
(ROOT/'paired_Q72_beta_CDR2_static_and_model.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps(records,ensure_ascii=False,indent=2))
