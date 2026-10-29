from pathlib import Path
import sys,json
ROOT=Path(__file__).resolve().parent
EXP=ROOT.parent
sys.path.insert(0,str(EXP/'label_sensitivity'/'runtime'))
import pandas as pd
INV=json.loads((EXP/'tcr_case_v40_20261017'/'local_pair_screen'/'sequence_inventory.json').read_text(encoding='utf8'))
rows=[]
for pdb,r in INV.items():
    meta=r['metadata'];citation=meta.get('citation',{})
    row={'pdb':pdb,'title':meta.get('title',''),'citation_title':citation.get('title',''),'DOI':str(citation.get('doi','')).lower()}
    for role in ['mhc','peptide','tcr_alpha','tcr_beta']:
        row[role+'_chain']=r['roles'].get(role,{}).get('chain','')
        row[role+'_sequence']=r['roles'].get(role,{}).get('sequence','')
    rows.append(row)
df=pd.DataFrame(rows).sort_values(['DOI','pdb'])
df.to_csv(ROOT/'all244_ITrace_metadata_and_sequences.csv',index=False,encoding='utf-8-sig')
for doi,g in df.groupby('DOI',dropna=False):
    if len(g)>=2:
        print(f"{doi} | {','.join(g.pdb)} | {g.citation_title.iloc[0][:115]} | "+','.join(g.peptide_sequence.unique()))
print('Total',len(df),'unique DOI',df.DOI.nunique())
