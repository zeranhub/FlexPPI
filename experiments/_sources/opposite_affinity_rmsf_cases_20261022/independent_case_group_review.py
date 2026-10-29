"""Independent source-level verification. Does not import inventory calculation code."""
from pathlib import Path
import json, hashlib
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from Bio.PDB import PDBParser

ROOT=Path(__file__).resolve().parent
I=ROOT/'inventory'
B=ROOT.parent
P=['3hg1','4jff','4jfd','4jfe']
Z={p:np.load(B/f'tcr_case_v40_20261017/family_predictions/full/itrace_{p}.npz',allow_pickle=False) for p in P}
A=pd.read_csv(I/'MEL5_three_contrasts_all_aligned_residues.csv')
S=pd.read_csv(I/'MEL5_three_contrasts_complete_scope_summaries.csv')
F=pd.read_csv(I/'MEL5_four_states_all_residues.csv')
bools=['paired','same_aa','four_state_pMHC_interface_8A','rmsf_valid','hchi_valid']
assert all(A[k].dtype==bool for k in bools)
errors=[]; source_values_checked=0
for r in A.itertuples(index=False):
    for suffix,p in [('before',r.pdb_before),('after',r.pdb_after)]:
        idx=getattr(r,'index_'+suffix)
        if pd.isna(idx):continue
        idx=int(idx);z=Z[p]
        for key,nkey in [('chain','chain'),('id','author_residue_id'),('aa','amino_acid')]:
            v=str(getattr(r,key+'_'+suffix));v=v[:-2] if key=='id' and v.endswith('.0') else v
            assert v==str(z[nkey][idx]),(p,idx,key,v,z[nkey][idx])
        for m,j in [('rmsf',0),('hchi',1)]:
            for field,nkey in [('full','mean'),('md','target'),('var','label_variance')]:
                expected=float(z[nkey][idx,j]);v=float(getattr(r,m+'_'+field+'_'+suffix))
                if field=='md' and not (z['valid'][idx,j] and np.isfinite(expected)):
                    assert np.isnan(v);continue
                errors.append(abs(v-expected));source_values_checked+=1
assert max(errors)<2e-14,max(errors)

summaries=[];numeric_residuals=[]
for r in S.itertuples(index=False):
    q=A[(A.pdb_before==r.pdb_before)&(A.pdb_after==r.pdb_after)]
    if r.scope=='complete_peptide':q=q[q.role=='peptide']
    elif r.scope=='all_six_CDR':q=q[q.role.isin(['tcr_alpha','tcr_beta'])&q.region_before.str.startswith('CDR')]
    elif r.scope=='four_state_union_pMHC_interface':q=q[q.four_state_pMHC_interface_8A]
    else:q=q[q.region_before==r.scope]
    m=r.metric;q=q[q.paired&q[m+'_valid']]
    if m=='hchi':q=q[q.same_aa]
    x=np.array([Z[a]['mean'][int(i),0 if m=='rmsf' else 1] for a,i in zip(q.pdb_before,q.index_before)],float)
    y=np.array([Z[a]['mean'][int(i),0 if m=='rmsf' else 1] for a,i in zip(q.pdb_after,q.index_after)],float)
    t0=np.array([Z[a]['target'][int(i),0 if m=='rmsf' else 1] for a,i in zip(q.pdb_before,q.index_before)],float)
    t1=np.array([Z[a]['target'][int(i),0 if m=='rmsf' else 1] for a,i in zip(q.pdb_after,q.index_after)],float)
    d=y-x;t=t1-t0
    rho=float(spearmanr(d,t).statistic) if len(d)>=3 and np.std(d)>0 and np.std(t)>0 else float('nan')
    row={'before':r.pdb_before,'after':r.pdb_after,'scope':r.scope,'metric':m,'n':len(d),
         'predicted_n_decrease':int((d<0).sum()),'predicted_n_increase':int((d>0).sum()),
         'MD_n_decrease':int((t<0).sum()),'MD_n_increase':int((t>0).sum()),
         'direction_agreement_n':int((np.sign(d)==np.sign(t)).sum()),
         'Spearman':rho,'MAE':float(abs(d-t).mean()) if len(d) else float('nan'),
         'zero_change_MAE':float(abs(t).mean()) if len(d) else float('nan')}
    if len(d):
        row['MAE_improvement_fraction_vs_zero']=1-row['MAE']/row['zero_change_MAE']
        for k,original in [('MAE',r.change_MAE_to_MD),('zero_change_MAE',r.zero_change_MAE_to_MD),('Spearman',r.change_Spearman_to_MD)]:
            if np.isfinite(row[k]) and np.isfinite(original):numeric_residuals.append(abs(row[k]-original))
        assert len(d)==r.n
        assert row['predicted_n_decrease']==r.n_decrease and row['predicted_n_increase']==r.n_increase
    summaries.append(row)
assert max(numeric_residuals)<2e-14,max(numeric_residuals)

# Independent static heavy-atom contact computation, full TCR chains D/E and HLA A.
contacts=[];parser=PDBParser(QUIET=True)
for p in P:
    mdl=parser.get_structure(p,f'E:/data/ITrace/trajectories/{p}_run1_topology.pdb')[0]
    for pos in [4,10]:
        residue=mdl['C'][(' ',pos,' ')]
        coords=np.array([a.coord.astype(float) for a in residue if a.element not in ['H','D']],float)
        out={'pdb':p,'peptide_position':pos,'resname':residue.resname}
        for role,chains in [('TCR',['D','E']),('HLA',['A'])]:
            min_d=float('inf');partner_n=0;atom_n=0
            for ch in chains:
                for rr in mdl[ch]:
                    if rr.id[0]!=' ':continue
                    cc=np.array([a.coord.astype(float) for a in rr if a.element not in ['H','D']],float)
                    if not len(cc):continue
                    ds=np.linalg.norm(coords[:,None,:]-cc[None,:,:],axis=-1)
                    min_d=min(min_d,float(ds.min())); n=int((ds<=4.5).sum())
                    atom_n+=n;partner_n+=int(n>0)
            out[role+'_minimum_heavy_distance_A']=min_d;out[role+'_heavy4_5A_partner_residue_n']=partner_n;out[role+'_heavy4_5A_atom_pair_n']=atom_n
        contacts.append(out)
C=pd.read_csv(I/'MEL5_four_states_all10_peptide_contact_summary.csv')
for r in contacts:
    q=C[(C.pdb==r['pdb'])&(C.peptide_id==r['peptide_position'])].iloc[0]
    for k,v in r.items():
        if k.startswith(('TCR_','HLA_')):
            assert abs(float(q[k])-v)<1e-12,(r['pdb'],r['peptide_position'],k,q[k],v)
assert all(r['TCR_heavy4_5A_partner_residue_n']==0 for r in contacts if r['peptide_position']==10)
assert all(r['TCR_heavy4_5A_partner_residue_n']>0 for r in contacts if r['peptide_position']==4)

key=[]
for p in P:
    z=Z[p]
    for name,ch,author_id in [('Q72','A','72'),('common_beta_G54','E','54' if p=='3hg1' else '55'),('P4','C','4'),('P5','C','5'),('P10','C','10')]:
        idx=np.flatnonzero((z['chain']==ch)&(z['author_residue_id']==author_id))
        assert len(idx)==1;idx=int(idx[0])
        key.append({'pdb':p,'name':name,'chain':ch,'author_id':author_id,'aa':str(z['amino_acid'][idx]),
                    'model_RMSF_A':float(z['mean'][idx,0]),'MD_RMSF_A':float(z['target'][idx,0]),
                    'model_Hchi_nats':float(z['mean'][idx,1]) if z['valid'][idx,1] else None,
                    'MD_Hchi_nats':float(z['target'][idx,1]) if z['valid'][idx,1] else None})

sources=[]
for p in P:
    for f in [B/f'tcr_case_v40_20261017/family_predictions/full/itrace_{p}.npz',Path(f'E:/data/ITrace/trajectories/{p}_run1_topology.pdb')]:
        sources.append({'path':str(f),'sha256':hashlib.sha256(f.read_bytes()).hexdigest()})
for f in I.glob('MEL5*.csv'):
    sources.append({'path':str(f),'sha256':hashlib.sha256(f.read_bytes()).hexdigest()})

verdict={
 'status':'numeric_verification_passed_with_scientific_scope_restrictions',
 'scope':'Exploratory reanalysis of four observed complex structures. Independent of inventory calculation code; not independent prospective biological validation.',
 'source_values_checked':source_values_checked,'source_maximum_absolute_error':max(errors),
 'summary_numeric_checks':len(numeric_residuals),'summary_maximum_absolute_error':max(numeric_residuals),
 'complete_region_sizes':{'peptide_RMSF_each_pair':10,'CDR_RMSF_each_pair':45,'four_state_union_pMHC_RMSF_each_pair':30,'pMHC_same_AA_Hchi_by_pair':[23,23,22]},
 'key_source_values':key,'P4_and_P10_static_contacts':contacts,'independently_recomputed_summaries':summaries,
 'scientific_verdict':{
   'main_claim_supported':False,
   'main_claim':'Joint RMSF and Hchi patterns explain or predict the sign of opposite affinity changes across this case group.',
   'narrow_claim_supported':'Complete spatial output preserves local differences that a peptide maximum-drop statistic removes. G4A P4 RMSF and complete peptide Hchi have some independently supported dynamic information; I7A broad RMSF increase is partly supported.',
   'body_placement_recommendation':'Keep as independent exploratory analysis now. Do not add as a main validation of affinity-mechanism discrimination. A concise worked example can be considered only if the paper explicitly claims local dynamic resolution and the complete adverse regional results remain reviewable.',
   'G4A_P10':'Maximum predicted drop occurs at P10, which has no direct TCR heavy-atom contacts <=4.5 A in any original state. Its MD RMSF increases, so terminal anchoring/actual rigidification is not established; this is also a prediction error.',
   'G4A_complete_peptide_RMSF':'Negative change correlation and MAE worse than zero change. P4 agreement cannot be generalized to the full peptide.',
   'G4A_complete_peptide_Hchi':'All 7 chemically defined unchanged peptide sites retained; rho about .821, modest MAE improvement. This supports local sidechain information but does not establish the origin of 60000-fold affinity loss.',
   'old_Q72_G54_probe':'New weak-state Q72 Hchi predictions are opposite to MD changes and G54 RMSF barely responds. The old positive cross-interface mechanism does not transfer to explain the weak states.',
   'I7A':'Predicted peptide RMSF increases at 9/10 positions and CDR RMSF at 39/45, with modest independent support. It is not an example of a large peptide RMSF decrease.',
   'water_and_entropy':'Water-bridge and temperature-dependent SPR literature can explain a physical layer not represented by model outputs. Do not present this as detected solvent entropy or as model validation.',
   'preorganization':'Bound mutation contrast cannot establish free-to-bound preorganization or its entropy cost.'
 },'sources':sources}
(ROOT/'final_case_group_review.json').write_text(json.dumps(verdict,ensure_ascii=False,indent=2,allow_nan=True),encoding='utf8')
print(json.dumps({'status':verdict['status'],'values_checked':source_values_checked,'max_error':max(errors),'summary_max_error':max(numeric_residuals),'contacts':contacts},ensure_ascii=False,indent=2))
