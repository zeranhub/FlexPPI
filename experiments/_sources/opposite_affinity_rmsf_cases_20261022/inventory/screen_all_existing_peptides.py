"""Inventory, not validation: retain every catalogued contrast and all peptide positions."""
from pathlib import Path
import sys,json,hashlib
import numpy as np,pandas as pd
from Bio.Align import PairwiseAligner
OUT=Path(__file__).resolve().parent;BASE=OUT.parent.parent
ATLAS=BASE/'candidate_case_atlas_20261021'
sys.path.insert(0,str(ATLAS/'data_and_general_figures'))
import prepare_all_data as prior
CAT=json.loads((ATLAS/'catalog/candidate_inventory.json').read_text(encoding='utf8'))
INTERVENTIONS={
'C01':'TCR engineering and separate same-TCR peptide alanine variants',
'C02':'Different unrelated natural TCRs; not a mutation parent/derivative',
'C03':'Same engineered TCR peptide G12D, plus separate TCR parent/derivative',
'C04':'Same TCR peptide F5L single mutation',
'C05':'Different natural TCRs; not a mutation parent/derivative',
'C06':'TCR beta E31D direct derivative',
'C07':'HLA single mutations, fixed TCR and peptide',
'C08':'Same TCR peptide C9V single mutation',
'C09':'Same TCR neoantigen K/R substitution; neither is wildtype',
'C10':'HLA A3 versus A11 multiple polymorphisms',
'C11':'TCR beta engineered derivative, construct differences',
'C12':'TCR CDR1alpha engineered derivative',
'C13':'Peptide decamer versus nonamer; length and anchor changed',
'C14':'Within each TCR peptide length variant; parallel DMF5 and DMF4',
'C15':'Same TCR different multi-substitution peptide antigens',
'C16':'TCR CDR2alpha/beta engineered derivative',
'C17':'TCR S32M, force-dependent endpoint not equilibrium KD',
'C18':'Different TCRs; endpoint not yet individually verified',
'C19':'Different TCRs; endpoint not yet individually verified',
'C20':'Different TCRs; endpoint not yet individually verified',
'C21':'Same TCR peptide point variants or cross-antigens; distinguish each contrast',
'C22':'Same TCR HBV single peptide variants; endpoint not yet individually verified',
'C23':'Same TCR gp100 P2 anchor variants; receptor-dependent effect',
'C24':'Same TCR HIV ligands; binding and signalling dissociation',
'C25':'Different multi-substitution MAGE/Titin antigens, crossreactivity',
'C26':'Same TCR SL9 escape variants; pMHC stability/presentation endpoint',
'C27':'Same TCR different multi-substitution antigens',
'C28':'Different TCR plus oligomeric state difference',
'C29':'Same TCR viral and self multi-substitution peptides',
'C30':'Different unrelated TCRs; not a mutation parent/derivative'}
KD={('C01','3hg1','4jff'):(18e-6,.6e-9,'exact SPR'),('C01','4jff','4jfd'):(.6e-9,36e-6,'exact SPR'),('C01','4jff','4jfe'):(.6e-9,31e-6,'exact SPR'),
('C03','7ow5','7ow6'):(3e-6,743e-12,'approximately 3uM WT; SPR'),('C04','1oga','5hhm'):(1.79e-6,200e-6,'after KD >200uM; lower bound weakening'),
('C07','1lp9','2j8u'):(9.3e-6,31.8e-6,'SPR'),('C07','1lp9','2jcc'):(9.3e-6,15.4e-6,'SPR'),('C07','1lp9','2uwe'):(9.3e-6,4.7e-6,'SPR'),
('C08','2bnr','2bnq'):(13.3e-6,5e-6,'SPR 25C; do not mix37C fixed-backbone pilot'),('C09','8yiv','8yj2'):(1.2e-6,3.4e-6,'SPR'),
('C10','8rym','8ryn'):(21e-6,72e-6,'SPR'),('C12','8rym','8ryo'):(21e-6,76e-9,'SPR'),
('C14','3qdg','3qdj'):(6e-6,40e-6,'KD as catalogued primary study'),('C14','3qdm','3qeq'):(29e-6,170e-6,'KD as catalogued primary study')}
S={};allrows=[];pep=[];errors=[]
for c in CAT['cases']:
    for a,b in c['contrasts']:
        row={'case_id':c['case_id'],'case_name':c['name'],'pdb_before':a,'pdb_after':b,'intervention':INTERVENTIONS[c['case_id']],
             'existing_limitations':c['model_cannot_support_or_limitations'],'KD_before_M':None,'KD_after_M':None,'affinity_fold_increase':None,'KD_status':'Not individually verified in this inventory; do not infer KD from cell activation',
             'source_type':'Observed experimental complex geometries; same 2430-record joint fold0 unless explicitly noted'}
        if (c['case_id'],a,b) in KD:
            k1,k2,sta=KD[(c['case_id'],a,b)];row.update(KD_before_M=k1,KD_after_M=k2,affinity_fold_increase=k1/k2,KD_status=sta)
        try:
            fa,sa=prior.load(a);fb,sb=prior.load(b)
            S[a]=sa;S[b]=sb
            aa=fa[fa.role.eq('peptide')].reset_index(drop=True);bb=fb[fb.role.eq('peptide')].reset_index(drop=True)
            sx=''.join(aa.aa);sy=''.join(bb.aa)
            hit=prior.make_aligner(peptide=True).align(sx,sy)[0]
            pairs=[]
            for col,(i,j) in enumerate(zip(*hit.indices)):
                pr={'case_id':c['case_id'],'pdb_before':a,'pdb_after':b,'alignment_column':col+1,'paired':i>=0 and j>=0,
                    'sequence_before':sx,'sequence_after':sy}
                ar=aa.iloc[i] if i>=0 else None;br=bb.iloc[j] if j>=0 else None
                for suffix,r in [('before',ar),('after',br)]:
                    for key in ['chain','id','aa']:pr[key+'_'+suffix]=r[key] if r is not None else ''
                    for m in ['rmsf','hchi']:
                        for key in ['full','md','valid']:pr[m+'_'+key+'_'+suffix]=r[m+'_'+key] if r is not None else (False if key=='valid' else np.nan)
                pr['same_aa']=bool(i>=0 and j>=0 and ar.aa==br.aa)
                pr['rmsf_drop_before_minus_after']=pr['rmsf_full_before']-pr['rmsf_full_after']
                pr['hchi_drop_before_minus_after']=pr['hchi_full_before']-pr['hchi_full_after'] if pr['same_aa'] and pr['hchi_valid_before'] and pr['hchi_valid_after'] else np.nan
                pep.append(pr)
                if pr['paired']:pairs.append(pr)
            drops=np.array([p['rmsf_drop_before_minus_after'] for p in pairs],float);j=int(np.argmax(drops));best=pairs[j]
            row.update(sequence_before=sx,sequence_after=sy,peptide_length_before=len(sx),peptide_length_after=len(sy),paired_position_n=len(pairs),
                       same_peptide_sequence=sx==sy,peptide_substitution_n=sum(not p['same_aa'] for p in pairs),
                       max_peptide_RMSF_drop_A=drops[j],max_drop_id_before=best['id_before'],max_drop_aa_before=best['aa_before'],
                       mean_peptide_RMSF_drop_A=float(drops.mean()),n_peptide_RMSF_decrease=int((drops>0).sum()),n_peptide_RMSF_increase=int((drops<0).sum()),
                       peptide_length_caution='Different length: alignment-based max is not same-position mutation metric' if len(sx)!=len(sy) else '',status='existing NPZ read')
        except Exception as e:
            row['status']='not computable';row['error']=repr(e);errors.append(row.copy())
        allrows.append(row)
pd.DataFrame(allrows).to_csv(OUT/'all_existing_contrasts_peptide_RMSF_and_comparability.csv',index=False,encoding='utf-8-sig',float_format='%.17g')
pd.DataFrame(pep).to_csv(OUT/'all_existing_contrasts_complete_peptide_RMSF_Hchi.csv',index=False,encoding='utf-8-sig',float_format='%.17g')
sources=[{'path':str(ATLAS/'catalog/candidate_inventory.json'),'sha256':hashlib.sha256((ATLAS/'catalog/candidate_inventory.json').read_bytes()).hexdigest()}]
for p,s in S.items():
    sources.append({'pdb':p,'path':s['prediction'],'sha256':s['prediction_sha256'],'structure_path':s['structure'],'structure_sha256':s['structure_sha256']})
(OUT/'all_existing_contrast_inventory_audit.json').write_text(json.dumps({'status':'passed' if not errors else 'partial source coverage','case_series_n':30,'contrasts_n':len(allrows),'read_predictions_n':len(S),'errors':errors,'new_inference':0,'new_MD':0,'sources':sources},ensure_ascii=False,indent=2),encoding='utf8')
print(pd.DataFrame(allrows)[['case_id','pdb_before','pdb_after','affinity_fold_increase','max_peptide_RMSF_drop_A','max_drop_id_before','status']].to_string(index=False))
