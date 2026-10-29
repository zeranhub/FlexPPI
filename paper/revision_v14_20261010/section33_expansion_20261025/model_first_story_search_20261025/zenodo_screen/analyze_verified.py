"""Use recomputed chemical-correct MD labels; no legacy Hchi substitution."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
OUT=Path(__file__).resolve().parent
DATA=Path(r'E:\Project\FlexPPI\artifacts\zenodo_6638504\data')
def score(y,z):
    strong=np.abs(y)>=.1
    return dict(n=len(y),md_mean=float(y.mean()),pred_mean=float(z.mean()),zero_mae=float(np.abs(y).mean()),mae=float(np.abs(y-z).mean()),rho=float(spearmanr(y,z).statistic),
                strong_n=int(strong.sum()),strong_direction=float((np.sign(y[strong])==np.sign(z[strong])).mean()) if strong.any() else None)
def main():
    pairs=[json.loads(x) for x in (DATA/'zenodo_6638504_pairs.jsonl').read_text().splitlines()]
    orig=pd.read_csv(OUT/'all_pair_residues_provisional.csv');rows=[];metrics=[];windows=[]
    for p in pairs:
        a=OUT/'raw_recomputed'/p['apo']['sample_id']/'labels_and_windows.npz';h=OUT/'raw_recomputed'/p['holo']['sample_id']/'labels_and_windows.npz'
        if not a.exists() or not h.exists():continue
        a=np.load(a);h=np.load(h);ai=np.array(p['apo_indices']);hi=np.array(p['holo_indices'])
        assert np.array_equal(a['aa'][ai],h['aa'][hi]);part=orig[orig.pair==p['pair_id']].reset_index(drop=True)
        assert np.array_equal(part.apo_index,ai) and np.array_equal(part.holo_index,hi)
        for j in range(len(ai)):
            row=dict(pair=p['pair_id'],chain=part.chain.iloc[j],author=part.author.iloc[j],aa=part.aa.iloc[j],interface_ca8=bool(part.interface_ca8.iloc[j]))
            for m,metric in enumerate(['rmsf','hchi']):
                row['md_apo_'+metric]=float(a['mean'][ai[j],m]);row['md_holo_'+metric]=float(h['mean'][hi[j],m]);row['md_delta_'+metric]=float(h['mean'][hi[j],m]-a['mean'][ai[j],m])
                row['valid_'+metric]=bool(a['mask'][ai[j],m] and h['mask'][hi[j],m])
                for v in ['full','intrachain']:row[v+'_delta_'+metric]=float(part[v+'_delta_'+metric].iloc[j])
            rows.append(row)
        for m,metric in enumerate(['rmsf','hchi']):
            valid=a['mask'][ai,m]&h['mask'][hi,m]
            if m==1:valid&=~np.isin(a['aa'][ai],list('AGP'))
            for region in ['all','interface','noninterface']:
                keep=valid.copy()
                if region!='all':keep&=part.interface_ca8.to_numpy(bool)==(region=='interface')
                for v in ['full','intrachain']:
                    y=h['mean'][hi[keep],m]-a['mean'][ai[keep],m];z=part[v+'_delta_'+metric].to_numpy()[keep]
                    metrics.append(dict(pair=p['pair_id'],metric=metric,region=region,model=v,**score(y,z)))
                    for b in range(min(len(a['blocks']),len(h['blocks']))):
                        yw=h['blocks'][b,hi[keep],m]-a['blocks'][b,ai[keep],m]
                        windows.append(dict(pair=p['pair_id'],metric=metric,region=region,model=v,block=b+1,**score(yw,z)))
    pd.DataFrame(rows).to_csv(OUT/'verified_all_residues.csv',index=False)
    pd.DataFrame(metrics).to_csv(OUT/'verified_metrics.csv',index=False)
    pd.DataFrame(windows).to_csv(OUT/'verified_block_metrics.csv',index=False)
    print(pd.DataFrame(metrics).to_string(index=False))
    for pair,g in pd.DataFrame(rows).groupby('pair'):
        print('\nINTERFACE',pair);print(g[g.interface_ca8][['author','aa','md_apo_rmsf','md_holo_rmsf','full_delta_rmsf','md_apo_hchi','md_holo_hchi','full_delta_hchi']].to_string(index=False))
if __name__=='__main__':main()
