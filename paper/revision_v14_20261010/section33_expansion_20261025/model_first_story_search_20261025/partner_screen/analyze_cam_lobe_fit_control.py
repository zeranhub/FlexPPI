"""Alignment control: fit each CaM lobe separately in each 100 ns block."""
import os
for k in ['OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS']:os.environ[k]='2'
from pathlib import Path
import sys
import numpy as np
import pandas as pd
OUT=Path(__file__).resolve().parent
sys.path.insert(0,r'E:\Project\FlexPPI\experiments\core_corrected_pipeline\training_payload\src')
from flexppi.data.labels import _chainwise_rmsf

def main():
    mapping=pd.read_csv(OUT/'cam_md_mature_number_mapping.csv')
    rows=[]
    for sid in ['A004C','A007U']:
        data=np.load(OUT/'verified_case_trajectories'/sid/'ca_coordinates.npz')
        mm=mapping[mapping.sample_id.eq(sid)].set_index('mature_cam_number')
        for name,positions in [('N_lobe',range(4,76)),('C_lobe',range(82,149))]:
            sub=mm.loc[list(positions)];indices=sub.residue_index.to_numpy(int)
            for block in range(5):
                x=data['coords_A'][block*1000:(block+1)*1000,indices]
                rmsf=_chainwise_rmsf(x,np.zeros(len(indices),dtype=int))
                for pos,value in zip(positions,rmsf):rows.append(dict(sample_id=sid,region=name,mature_cam_number=pos,block=block+1,lobe_fitted_rmsf=float(value)))
    data=pd.DataFrame(rows);data.to_csv(OUT/'cam_md_per_lobe_fitted_rmsf_five_blocks.csv',index=False)
    wide=data.pivot(index=['region','mature_cam_number','block'],columns='sample_id',values='lobe_fitted_rmsf').reset_index()
    wide['delta']=wide.A007U-wide.A004C
    wide.to_csv(OUT/'cam_md_per_lobe_fitted_rmsf_deltas.csv',index=False)
    summary=wide.groupby(['region','block'])[['A004C','A007U','delta']].mean().reset_index()
    summary.to_csv(OUT/'cam_md_per_lobe_fitted_rmsf_control_summary.csv',index=False)
    print(summary.to_string(index=False))
    print(summary.groupby('region')[['A004C','A007U','delta']].mean().to_string())
if __name__=='__main__':main()
