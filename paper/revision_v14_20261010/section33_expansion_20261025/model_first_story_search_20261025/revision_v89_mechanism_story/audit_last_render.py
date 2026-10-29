from pathlib import Path
import json
V=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010\section33_expansion_20261025\model_first_story_search_20261025\revision_v89_mechanism_story')
changed=[];same=[]
for n in range(1,30):
 a=V/'render/release_01'/f'page-{n}.png';b=V/'render/release_02'/f'page-{n}.png'
 (same if a.read_bytes()==b.read_bytes() else changed).append(n)
assert changed==[24],changed
(V/'manuscript/last_correction_render_audit.json').write_text(json.dumps({'change':'Add Y215 to complete six-member F298 contact list in paragraph163; no figure changed','changed_pages':changed,'identical_to_independently_reviewed_release01_pages':same,'final_render':'render/release_02'},ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({'changed_pages':changed,'identical_pages':len(same)}))