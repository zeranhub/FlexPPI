"""Measure deposited local contacts, distinct from dynamics and mutant inference."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from analyze_predictions import read_atoms, BASE, OUT

def minimum(a,b,names_a=None,names_b=None):
    ka=[k for k in a if names_a is None or k in names_a]
    kb=[k for k in b if names_b is None or k in names_b]
    d=np.linalg.norm(np.array([a[k] for k in ka])[:,None,:]-np.array([b[k] for k in kb])[None,:,:],axis=-1)
    i,j=np.unravel_index(np.argmin(d),d.shape)
    return ka[i],kb[j],float(d[i,j])

def main():
    specs={s['id']:s for s in json.loads((BASE/'states.json').read_text())}
    out=[]
    for sid,cha,chb in [('b4_5iw1_ab_full','A','B'),('b4_5iw1_cd_full','C','D'),('b4_5iw1_ef_full','E','F'),('b4_5ivx_full_full','E','F')]:
        ats=read_atoms(specs[sid]['path'])
        for r1,r2,atoms1,atoms2,label in [('121','134',['OE1','OE2'],['NZ'],'E121alpha-K134beta charged atoms'),
            ('125','133',['OH'],['OD1','ND2'],'Y125alpha-N133beta polar atoms'),('125','134',['OH'],['NZ'],'Y125alpha-K134beta polar atoms')]:
            n1,n2,d=minimum(ats[(cha,r1)],ats[(chb,r2)],atoms1,atoms2)
            out.append(dict(case='B4.2.3',state=sid,pair=label,chain_1=cha,author_1=r1,atom_1=n1,chain_2=chb,author_2=r2,atom_2=n2,distance_A=d))
    for p in ['8tw4','8tw6','9bbc','9c3e']:
        sid=f'1g4_{p}_full_full';ats=read_atoms(specs[sid]['path'])
        for label,names in [('S104-V182 CA distance',['CA']),('S104-V182 CB distance',['CB'])]:
            n1,n2,d=minimum(ats[('A','104')],ats[('A','182')],names,names)
            out.append(dict(case='1G4-CD3',state=sid,pair=label,chain_1='A',author_1='104',atom_1=n1,chain_2='A',author_2='182',atom_2=n2,distance_A=d))
    pd.DataFrame(out).to_csv(OUT/'local_experimental_geometry.csv',index=False)
    print(pd.DataFrame(out).to_string(index=False))

if __name__=='__main__':
    main()
