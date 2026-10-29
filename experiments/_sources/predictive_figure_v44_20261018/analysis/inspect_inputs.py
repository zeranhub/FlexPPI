from pathlib import Path
import sys, json, gzip
BASE=Path('E:/Project/FlexPPI/experiments/_sources')
sys.path.insert(0,str(BASE/'completed_results_20261012/plot_runtime'))
import numpy as np
import pandas as pd
CORE=Path('E:/Project/FlexPPI/experiments/core_corrected_pipeline')
for p in [next((CORE/'neural_exports/evaluation/full').glob('*.npz')),CORE/'static_corrected/static_sequence_complex_oof.npz',CORE/'static_corrected/aa_only_oof.npz']:
 with np.load(p,allow_pickle=True) as z:
  print(str(p));print({k:{'shape':z[k].shape,'dtype':str(z[k].dtype),'first':z[k].reshape(-1)[:8].tolist()} for k in z.files})
truth=pd.read_csv(BASE/'baselines/cohort/truth_residues.csv.gz',nrows=6)
print(truth.to_string(index=False));print(truth.dtypes.to_string())
for p in [CORE/'neural_exports/evaluation/prediction_manifest.jsonl',CORE/'neural_exports/evaluation/per_record_metrics.csv',CORE/'static_corrected/five_fold_mean_sample_sd.csv']:
 print(str(p));print(p.read_text(encoding='utf8').splitlines()[:4])
