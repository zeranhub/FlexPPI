from pathlib import Path
import json,sys,types
ROOT=Path(__file__).resolve().parent
EXP=ROOT.parents[1]
sys.path.insert(0,str(EXP/'label_sensitivity/runtime'))
SRC=Path('E:/Project/FlexPPI/experiments/core_corrected_pipeline/training_payload/src/flexppi')
for name,sub in [('flexppi',''),('flexppi.data','data')]:
    mod=types.ModuleType(name);mod.__path__=[str(SRC/sub)];sys.modules[name]=mod
from flexppi.data.labels import _chainwise_rmsf
import numpy as np,pandas as pd
df=pd.read_csv(ROOT/'all408_per_replica_labels.csv')
summary=pd.read_csv(ROOT/'all408_ten_replica_mean_SD_SEM.csv')
assert len(df)==8160 and len(summary)==816
assert df.groupby(['condition','replica']).size().eq(408).all()
assert df.frames.eq(103).all()
assert df[df.aa.isin(['A','G'])].Hchi_nats.isna().all()
assert df[df.aa=='P'].Hchi_nats.notna().all()
maxerr=0
for condition in ['WT','S32M']:
    for rep in range(1,11):
        x=np.load(ROOT/f'{condition}_rep{rep:02d}_all408.npz')
        _,chain_index=np.unique(x['chain'],return_inverse=True)
        exact=_chainwise_rmsf(x['ca_coordinates_A'],chain_index)
        err=float(np.max(abs(exact-x['rmsf_A'])))
        maxerr=max(maxerr,err)
        np.testing.assert_allclose(exact,x['rmsf_A'],atol=3e-6,rtol=0)
        np.testing.assert_array_equal(x['frame_indices_zero_based'],np.arange(7,110))
audit={'status':'passed','rows':len(df),'summary_rows':len(summary),'independent_trajectories':20,
       'protein_residues_per_trajectory':408,'frames_per_trajectory':103,
       'vectorized_Kabsch_vs_original_training_label_chainwise_max_error_A':maxerr,
       'tolerance_A':3e-6,'AlaGly_Hchi_missing':True,'Proline_two_chi_computed':True}
(ROOT/'independent_numerical_QA.json').write_text(json.dumps(audit,indent=2))
print(json.dumps(audit,indent=2))
