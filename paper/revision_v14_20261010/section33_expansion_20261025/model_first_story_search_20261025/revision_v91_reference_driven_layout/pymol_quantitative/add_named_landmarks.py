from pathlib import Path
import json
out=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010\section33_expansion_20261025\model_first_story_search_20261025\revision_v91_reference_driven_layout\pymol_quantitative')
m=json.loads((out/'render_manifest.json').read_text())
for metric in ['RMSF','Hchi']:
 e=m['assets'][f'1MXE_{metric}_axis090']
 named={e['landmarks'][n]['name']:r for n,r in e['focus_landmarks'].items()}
 e['focus_named_landmarks']=named;m['delivery_assets'][metric]['focus_named_landmarks']=named
(out/'render_manifest.json').write_text(json.dumps(m,ensure_ascii=False,indent=2),encoding='utf8')