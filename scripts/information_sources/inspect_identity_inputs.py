import json
import collections
from pathlib import Path
import numpy as np

rows = [json.loads(s) for s in Path('artifacts/evaluation_atlas_hchi_v1/oof/prediction_manifest.jsonl').read_text().splitlines()]
print(collections.Counter(x['dataset'] for x in rows))
for dataset in ['dynarepo', 'itrace']:
    print(json.dumps(next(x for x in rows if x['dataset'] == dataset)['record'], indent=2))
arrays = np.load(rows[0]['cache_path'])
print(arrays.files)
print(str(arrays['metadata_json'])[:1800])
roots = [('/ranze/repos/FlexPPI', 'E:/Project/FlexPPI'),
         ('/ranze/data/ATLAS_extracted', 'E:/data/protein_md_datasets/ATLAS_extracted'),
         ('/ranze/data/dynarepo/dynarepo_non_nucleic_full', 'D:/data/dynarepo_non_nucleic_full'),
         ('/ranze/data/ITrace', 'E:/data/ITrace')]
counts = collections.Counter()
missing = []
for row in rows:
    path = row['record']['structure_path']
    for old, new in roots:
        path = path.replace(old, new)
    exists = Path(path).exists()
    counts[(row['dataset'], exists)] += 1
    if not exists and len(missing) < 5:
        missing.append(path)
print(counts)
print(missing)
