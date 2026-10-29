from pathlib import Path
import json, pandas as pd, numpy as np
B=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010\section33_expansion_20261025\feature_baseline')
R=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010\section33_expansion_20261025\model_first_story_search_20261025\class_level_feasibility_20261026')
fam=pd.read_csv(r'E:\Project\FlexPPI\experiments\_sources\family_folds\record_family_folds.csv')
cand=pd.read_csv(R/'inventory_agent/short_peptide_5to40_candidates.csv'); cand=cand[cand.two_chain_unambiguous]
f=pd.read_parquet(B/'features_all.parquet')
print('shape',f.shape,'memory MB',f.memory_usage(deep=True).sum()/1e6)
print('candidate fold counts',fam[fam.sample_id.isin(cand.sample_id)].fold.value_counts().sort_index().to_dict())
print('feature minmax categorical',[(c, f[c].min(),f[c].max()) for c in f if c.startswith('aa_offset_')])
print('label invariance',f[['valid_rmsf','valid_hchi']].sum().to_dict())
print('candidate IDs',len(cand),cand.sample_id.nunique(),'overlap',len(set(cand.sample_id)&set(f.sample_id)))
print('membership foldindex',fam['fold'].value_counts().sort_index().to_dict())
