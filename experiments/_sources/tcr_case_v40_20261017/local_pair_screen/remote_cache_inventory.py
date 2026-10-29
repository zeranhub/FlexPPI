from pathlib import Path
import json, zipfile
import numpy as np
root=Path('/root/flexppi_completion_20261010/family_folds/evaluation')
z=np.load(root/'full/itrace_1ao7.npz',allow_pickle=False)
print('SCHEMA',[(k,z[k].shape,str(z[k].dtype)) for k in z.files])
for key in ['metric_names','label_names','chain','author_residue_id','amino_acid']:
    if key in z: print(key,z[key][:15])
out=Path('/tmp/flexppi_tcr_family_predictions_v40_20261017.zip')
with zipfile.ZipFile(out,'w',compression=zipfile.ZIP_DEFLATED) as archive:
    for variant in ['full','intrachain']:
        for path in sorted((root/variant).glob('itrace_*.npz')):
            archive.write(path,path.relative_to(root))
    archive.write(root/'prediction_manifest.jsonl','prediction_manifest.jsonl')
print('ARCHIVE',out,out.stat().st_size)
