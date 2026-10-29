from pathlib import Path
import json,hashlib,datetime
import pandas as pd,numpy as np
R=Path(__file__).resolve().parent
O=R/'analysis_v1_20261026'
S=O/'shared'
S.mkdir(parents=True,exist_ok=True)
old=R/'inventory_agent'
cand=pd.read_csv(old/'short_peptide_5to40_candidates.csv')
cand=cand[cand.two_chain_unambiguous].copy().sort_values(['pdb_id','sample_id'])
assert len(cand)==56 and cand.sample_id.is_unique and cand.pdb_id.nunique()==54
cand.to_csv(S/'candidates_frozen56.csv',index=False)
pd.DataFrame([dict(pdb_id=pdb,n_records=len(g),sample_ids='|'.join(g.sample_id),split_groups='|'.join(sorted(set(g.family_id)))) for pdb,g in cand.groupby('pdb_id')]).to_csv(S/'pdb_groups_frozen54.csv',index=False)
source=Path(r'C:\Users\DELL\.codex\worktrees\2a52\FlexPPI\output\non_tcr_case_screen_20261025\archive_extract\family_folds\evaluation\oof_per_residue.csv.gz')
expected='0852469b8fb6dd2977796d040ce8bc7fe04f3de29f217e7fc563213c064ad7f3'
h=hashlib.sha256(source.read_bytes()).hexdigest()
assert h==expected
wanted=set(cand.sample_id)
pieces=[]
for chunk in pd.read_csv(source,chunksize=200000,dtype={'chain':str,'author_residue_id':str}):
    sub=chunk[chunk.sample_id.isin(wanted)]
    if len(sub):pieces.append(sub)
df=pd.concat(pieces,ignore_index=True)
keys=['sample_id','residue_index','variant']
assert not df.duplicated(keys).any()
rolemap=pd.read_csv(R/'short_peptide_preliminary'/'all56_records_two_roles.csv')[['sample_id','chain','role','pdb_id','component_id']]
assert not rolemap.duplicated(['sample_id','chain']).any()
df=df.merge(rolemap,on=['sample_id','chain'],validate='many_to_one')
df.to_parquet(S/'oof_candidates56.parquet',index=False)
manifest=dict(date='2026-10-26',analysis_status='Exploratory class extension after prior single-case and preliminary-class inspection. Not preregistered.',source_oof=str(source),source_oof_sha256=h,candidate_member_sha256=hashlib.sha256((S/'candidates_frozen56.csv').read_bytes()).hexdigest(),n_candidate_records=56,n_candidate_PDBs=54,n_split_groups=int(cand.family_id.nunique()),n_oof_rows=len(df),roles=sorted(df.role.unique()),n_full_residues=int((df.variant=='full').sum()),sources_read_only=True,model_quality_not_membership_criterion=True,manuscript_modified=False)
(S/'source_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(manifest,ensure_ascii=False))
