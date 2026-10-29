from pathlib import Path
import json
V=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010\section33_expansion_20261025\model_first_story_search_20261025\revision_v89_mechanism_story')
p=V/'narrative/paragraph_replacements.json';data=json.loads(p.read_text(encoding='utf-8'));assert 'Y215' not in data['paragraph_replacements']['163'];data['paragraph_replacements']['163']=data['paragraph_replacements']['163'].replace('F104、Y184和P216','F104、Y184、P216及Y215');p.write_text(json.dumps(data,ensure_ascii=False,indent=2),encoding='utf-8')
for name in ['prepare_manuscript.py','README_交付与复绘.md','prepare_release_audit.py']:
 p=V/name;t=p.read_text(encoding='utf-8').replace('F104、Y184和P216','F104、Y184、P216及Y215').replace('F104/Y184/P216','F104/Y184/P216/Y215');p.write_text(t,encoding='utf-8')