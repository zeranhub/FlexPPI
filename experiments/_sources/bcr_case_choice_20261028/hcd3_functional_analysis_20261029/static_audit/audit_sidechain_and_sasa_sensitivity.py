from pathlib import Path
import copy, json
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley

OUT=Path(__file__).parent
PDB=Path(r'E:\Project\FlexPPI\experiments\core_corrected_pipeline\training_payload\dataset\structures\A00LZ.pdb')
BACKBONE={'N','CA','C','O','OXT'}
s=PDBParser(QUIET=True).get_structure('6A77',str(PDB))
for r in s.get_residues():
    for a in list(r):
        if a.element.upper() in {'H','D'}:
            r.detach_child(a.id)
allres=list(s.get_residues())
rows=[]
for r in s[0]['C']:
    if not 95<=r.id[1]<=109: continue
    rec={'chain':'C','author':r.id[1],'aa3':r.resname}
    for atomtype,names in [('sidechain',None),('backbone',BACKBONE)]:
        target=[a for a in r if (a.name not in BACKBONE if atomtype=='sidechain' else a.name in BACKBONE)]
        rc=np.asarray([a.coord for a in target])
        for category,partner_chains in [('antigen',{'A'}),('light',{'B'})]:
            nearest=(float('inf'),None,None,None)
            hits={4.0:[],4.5:[],5.0:[]}
            for rr in allres:
                if rr.parent.id not in partner_chains: continue
                ds=cdist(rc,np.asarray([a.coord for a in rr]))
                d=float(ds.min())
                q=np.unravel_index(np.argmin(ds),ds.shape)
                if d<nearest[0]: nearest=(d,target[q[0]].name,f'{rr.parent.id}{rr.id[1]}{rr.resname}',list(rr)[q[1]].name)
                for threshold in hits:
                    if d<=threshold: hits[threshold].append(f'{rr.parent.id}{rr.id[1]}{rr.resname}')
            rec[f'{atomtype}_to_{category}_min_A']=nearest[0]
            rec[f'{atomtype}_to_{category}_nearest']=f'{nearest[1]}--{nearest[2]}:{nearest[3]}'
            for threshold in hits:
                rec[f'{atomtype}_{category}_contacts_le_{threshold:g}']=';'.join(hits[threshold])
    rows.append(rec)
table=pd.DataFrame(rows)
for points in [240,960]:
    for context,chains in [('complex',{'A','B','C'}),('Fab',{'B','C'}),('heavy',{'C'})]:
        st=copy.deepcopy(s)
        for c in list(st[0]):
            if c.id not in chains: st[0].detach_child(c.id)
        ShrakeRupley(n_points=points,probe_radius=1.4).compute(st,level='R')
        for k in table.index:
            r=st[0]['C'][(' ',int(table.loc[k,'author']),' ')]
            table.loc[k,f'{context}_residue_SASA_{points}']=sum(a.sasa for a in r)
            table.loc[k,f'{context}_sidechain_SASA_{points}']=sum(a.sasa for a in r if a.name not in BACKBONE)
table.to_csv(OUT/'complete_HCDR3_sidechain_contacts_SASA_sensitivity.csv',index=False,encoding='utf-8-sig')
print(table.to_string(index=False))
