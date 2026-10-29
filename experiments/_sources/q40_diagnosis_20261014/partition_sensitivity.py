from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parent;OUT=ROOT/'partition';OUT.mkdir(exist_ok=True)
CASE=ROOT.parent/'free_energy_curves_v27_20261014/case_curves'
SRC=CASE/'angles_5000frames_all6GLN.csv';REF=CASE/'joint_entropy_original_five_windows.csv'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
d=pd.read_csv(SRC);ref=pd.read_csv(REF);raw=[];summaries=[];boundary=[]
def entropy(codes,n=27):
    counts=np.bincount(codes,minlength=n);p=counts[counts>0]/len(codes)
    return float(-np.sum(p*np.log(p)))
for site,g in d.groupby('residue',sort=False):
    a=g.sort_values('frame_zero_based')[['chi1_degrees','chi2_degrees','chi3_degrees']].to_numpy(dtype=np.float32)
    assert a.shape==(5000,3) and np.isfinite(a).all()
    vals=[]
    for offset in range(0,120,5):
        if offset==0:b=np.clip(np.floor((a+np.float32(180))/np.float32(120)).astype(int),0,2)
        else:b=np.floor(((a.astype(float)+180-offset)%360)/120).astype(int)
        h=[]
        for w in range(5):
            v=entropy((b[w*1000:(w+1)*1000]*(3**np.arange(3))).sum(axis=1));h.append(v)
            raw.append(dict(residue=site,grid_offset_degrees=offset,window=w+1,Hchi_joint_nats=v))
            if offset==0:assert abs(v-ref.query('residue==@site and window==@w+1').Hchi_nats.iloc[0])<1e-12
        vals.append(np.mean(h))
    summaries.append(dict(residue=site,original_5window_mean=vals[0],offset24_min_Hchi=min(vals),offset24_max_Hchi=max(vals),offset24_range_Hchi=max(vals)-min(vals),offset24_SD_Hchi=np.std(vals,ddof=0)))
    # Proximity to the original coarse-state boundaries, including the circular seam.
    af=a.astype(float);dist=np.min(np.stack([abs(((af-b+180)%360)-180) for b in [-180,-60,60]],axis=-1),axis=-1)
    for q in range(3):
        boundary.append(dict(residue=site,chi=q+1,within5deg_of_original_boundary_fraction=float((dist[:,q]<=5).mean()),within10deg_of_original_boundary_fraction=float((dist[:,q]<=10).mean())))
pd.DataFrame(raw).to_csv(OUT/'24_offsets_5windows_joint_entropy.csv',index=False,float_format='%.17g')
pd.DataFrame(summaries).to_csv(OUT/'partition_summary.csv',index=False,float_format='%.17g')
pd.DataFrame(boundary).to_csv(OUT/'original_boundary_proximity.csv',index=False,float_format='%.17g')
audit={'status':'passed','source_sha256':{str(p):sha(p) for p in [SRC,REF]},'definition':'Diagnostic sensitivity to uniformly shifting all three 120-degree state grids over24 offsets. Only offset0 is the trained/evaluated target; shifted labels must not be compared with the frozen model as new accuracy results. No bin-width change, no smoothing. Five windows retained.','sources_unchanged':True,'summary':summaries,'boundary':boundary}
(OUT/'audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf8')
print(pd.DataFrame(summaries).to_string(index=False))
print(pd.DataFrame(boundary).to_string(index=False))
