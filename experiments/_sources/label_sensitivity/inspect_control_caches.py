from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent/'runtime'))
import numpy as np,json
base=Path(r'E:\Project\FlexPPI\experiments\core_corrected_pipeline')
for relative in ['static_corrected/training_features.npz','static_corrected/static_sequence_complex_oof.npz','static_corrected/aa_only_oof.npz','backflip2_corrected_scoring/training_aligned.npz']:
    p=base/relative
    print(relative,p.stat().st_size)
    try:
        with np.load(p,allow_pickle=False) as data:print([(k,data[k].shape,str(data[k].dtype)) for k in data])
    except Exception as err:print(type(err).__name__,str(err))
schema=json.loads((base/'static_corrected/training_schema.json').read_text(encoding='utf8'))
print(schema.keys());print(schema['records'][0])
