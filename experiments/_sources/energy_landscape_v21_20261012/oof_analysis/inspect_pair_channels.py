from pathlib import Path
import numpy as np
import json
root = Path('/root/flexppi_revision_20261002/corrected_training_v1')
for path in [root/'evaluation/full/16pk_A.npz',root/'dataset/labels/16pk_A.npz']:
    with np.load(path,allow_pickle=False) as value:
        print(path)
        for name in value.files:
            array=value[name]
            print(name,array.shape,array.dtype,str(array.tolist())[:900] if array.size<30 else '')
for path in [root/'configs/full_fold0.json',root/'configs/full_fold_0.json',root/'configs/full_fold0.yaml',root/'configs/full_fold_0.yaml']:
    if path.exists():
        print(path.read_text()[:3000])
print('configs',list((root/'configs').iterdir())[:8])
