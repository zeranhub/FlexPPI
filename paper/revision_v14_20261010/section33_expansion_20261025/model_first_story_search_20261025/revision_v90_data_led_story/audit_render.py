from pathlib import Path
import json
S=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010\section33_expansion_20261025\model_first_story_search_20261025');V=S/'revision_v90_data_led_story';O=S/'revision_v89_mechanism_story'
same=[];changed=[]
for n in range(1,30):
 a=V/'render/release_01'/f'page-{n}.png';b=O/'render/release_02'/f'page-{n}.png'
 (same if a.read_bytes()==b.read_bytes() else changed).append(n)
assert changed==[10,11,12],changed
(V/'manuscript/render_page_identity_audit.json').write_text(json.dumps({'same_as_previously_visually_verified_v89_pages':same,'changed_pages_visually_reviewed_in_v90':changed,'page_count':29},ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({'same':same,'changed':changed}))