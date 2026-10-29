import os
os.environ['OPENBLAS_NUM_THREADS']='1'
from pathlib import Path
from difflib import SequenceMatcher
import json
import numpy as np,pandas as pd,pyarrow.parquet as pq
from scipy.spatial.distance import cdist
from scipy.stats import spearmanr
R=Path('E:/Project/FlexPPI/paper/revision_v14_20261010/section33_expansion_20261025');B=R/'feature_baseline';O=Path(__file__).resolve().parent
sid='A0094';methods=['full','intra','aa_sasa_contacts','static47_sasa']
f=pq.read_table(B/'features_all.parquet',filters=[('sample_id','=',sid)]).to_pandas();p=pq.read_table(B/'oof_predictions_all.parquet',filters=[('sample_id','=',sid)]).to_pandas();pc=[c for c in p if c.startswith('pred_')];f=f.drop(columns=[c for c in pc if c in f]);f=f.merge(p[['sample_id','residue_index']+pc],on=['sample_id','residue_index'],validate='one_to_one')
aa3=dict(zip('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split(),'ARNDCQEGHILKMFPSTWYV'))
orig={}
for l in (O/'5ET1_original.pdb').read_text().splitlines():
 if l.startswith('ATOM') and l[16] in ' A' and l[12:16].strip()[0] not in 'HD':
  k=(l[21],l[22:27].strip());v=orig.setdefault(k,{'aa':aa3[l[17:20]],'coords':[]});v['coords'].append([float(l[30:38]),float(l[38:46]),float(l[46:54])])
mapping={};audit=[]
for modelch,originalchoices in [('A','AB'),('B','CD')]:
 part=f[f.chain==modelch].sort_values('residue_index');seq=''.join(part.amino_acid);candidates=[]
 for origch in originalchoices:
  keys=[k for k in orig if k[0]==origch];refseq=''.join(orig[k]['aa'] for k in keys);matcher=SequenceMatcher(None,seq,refseq,autojunk=False);m={}
  for a,b,n in matcher.get_matching_blocks():
   for i in range(n):m[int(part.iloc[a+i].residue_index)]=keys[b+i]
  candidates.append((len(m),origch,m))
 n,choice,m=max(candidates,key=lambda t:t[0]);mapping.update(m);audit.append(dict(model_chain=modelch,chosen_original_chain=choice,matched=n,model_n=len(part),candidate_coverage=[(v[0],v[1]) for v in candidates]))
f['original_chain']=[mapping.get(int(i),('',''))[0] for i in f.residue_index];f['original_residue']=[mapping.get(int(i),('',''))[1] for i in f.residue_index]
pep=f[f.chain=='B'];pdbch=pep.original_chain.iloc[0]
print('Myo3b mapped peptide:',','.join(pep.amino_acid+pep.original_residue),flush=True)
# Three binding surfaces independently defined by the deposited partner segments.
# Article's N-terminal RR end, middle KxL stretch, C-terminal YY helix.
segments={'RR_end':(1250,1253),'KxL_middle':(1254,1262),'YY_helix':(1263,1275)}
for name,(start,stop) in segments.items():
 coords=np.vstack([v['coords'] for (ch,r),v in orig.items() if ch==pdbch and start<=int(r)<=stop])
 vals=[]
 for i in f.residue_index:
  k=mapping.get(int(i));vals.append(float(cdist(orig[k]['coords'],coords).min()) if k and k[0]!=pdbch else np.nan)
 f[name+'_distance_A']=vals
f.to_csv(O/'espin_all_residues_mapped.csv',index=False)
q=f[(f.chain=='A')&f.valid_hchi&f.valid_rmsf].copy();rows=[]
masks={'whole_Espin':np.ones(len(q),bool),'all_partner_contacts_5A':q[[n+'_distance_A' for n in segments]].min(axis=1)<=5,'deep_buried_sasa5':q.sidechain_sasa_complex_A2<=5,'surface_sasa30':q.sidechain_sasa_complex_A2>=30}
for name in segments:masks[name+'_contact5']=q[name+'_distance_A']<=5;masks[name+'_shell8']=q[name+'_distance_A']<=8
for name,ix in masks.items():
 t=q[ix];row=dict(region=name,n=len(t),rmsf_mean=t.md_rmsf.mean(),rmsf_sd=t.md_rmsf.std(),rmsf_min=t.md_rmsf.min(),rmsf_max=t.md_rmsf.max(),hchi_mean=t.md_hchi.mean(),hchi_sd=t.md_hchi.std(),hchi_min=t.md_hchi.min(),hchi_max=t.md_hchi.max(),hchi_iqr=t.md_hchi.quantile(.75)-t.md_hchi.quantile(.25))
 for m in methods:
  for v in ['rmsf','hchi']:
   y=t['md_'+v];yp=t['pred_'+m+'_'+v];row[m+'_'+v+'_mae']=float(abs(y-yp).mean());row[m+'_'+v+'_rho']=float(spearmanr(y,yp).statistic) if len(y)>=4 else np.nan;row[m+'_'+v+'_mean']=yp.mean()
 rows.append(row)
pd.DataFrame(rows).to_csv(O/'espin_complete_functional_region_metrics.csv',index=False)
(O/'espin_structure_mapping_audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
print(pd.DataFrame(rows)[['region','n','rmsf_mean','hchi_mean','hchi_sd','full_hchi_mae','intra_hchi_mae','aa_sasa_contacts_hchi_mae','static47_sasa_hchi_mae','full_hchi_rho']].to_string(index=False),flush=True)
print(q[q.original_residue.isin(['69','103','110','205'])][['amino_acid','original_residue','md_rmsf','pred_full_rmsf','md_hchi','pred_full_hchi','pred_intra_hchi','pred_aa_sasa_contacts_hchi','pred_static47_sasa_hchi']].to_string(index=False),flush=True)
