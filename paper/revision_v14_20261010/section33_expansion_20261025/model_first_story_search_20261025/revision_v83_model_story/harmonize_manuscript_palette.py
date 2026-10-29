from pathlib import Path
import json
v=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010\section33_expansion_20261025\model_first_story_search_20261025\revision_v83_model_story');f=v/'figure3';s=v/'supplement_cam'
for name in ['model_story_body.py','compose_fig3.py']:
    p=f/name;t=p.read_text(encoding='utf-8').replace("STATE_A='#0072B2'; STATE_B='#D55E00'; KEY='#009E73'","STATE_A='#B87936'; STATE_B='#2166AC'; KEY='#3A8D8C'");p.write_text(t,encoding='utf-8')
p=s/'compose_figS3.py';t=p.read_text(encoding='utf-8').replace("('NHE1', 'a', '#0072B2', 'o'), ('IQCG', 'b', '#D55E00', 's')","('NHE1', 'a', '#B87936', 'o'), ('IQCG', 'b', '#2166AC', 's')");p.write_text(t,encoding='utf-8')
p=f/'caption_zh.md';t=p.read_text(encoding='utf-8').replace('蓝色圆点/曲线表示NHE1，朱橙方点/曲线表示IQCG','暖赭圆点/曲线表示NHE1，蓝色方点/曲线表示IQCG').replace('关键Met为绿色','关键Met为青绿色');p.write_text(t,encoding='utf-8')
p=s/'caption_zh.md';t=p.read_text(encoding='utf-8').replace('蓝色为NHE1，朱橙色为IQCG','暖赭色为NHE1，蓝色为IQCG');p.write_text(t,encoding='utf-8')
p=v/'manuscript/v83_text_plan.json';d=json.loads(p.read_text(encoding='utf-8'));d['notes'].append('按图1/2/4/5及publication_style.py协调：NHE1暖赭#B87936、IQCG蓝#2166AC，结构Key Met青绿#3A8D8C；图3与S3统一。');p.write_text(json.dumps(d,ensure_ascii=False,indent=2),encoding='utf-8')
audit={'palette_source':'E:/Project/FlexPPI/paper/revision_v14_20261010/figure2_dense_20261024/source_internal/publication_style.py','inspected_existing_figures':[1,2,4,5],'state_roles':{'NHE1':'#B87936','IQCG':'#2166AC'},'structure_roles':{'CaM_N':'#8A9299','CaM_C':'#BAC1C7','partner':'#8B72A5','key_Met':'#3A8D8C'},'role_separation':'Molecular components do not reuse the state palette; state labels are black over structures.','existing_figure4_order':'warm baseline, blue alternate; same ordering retained for NHE1 to IQCG','state_marker_redundancy':{'NHE1':'circle','IQCG':'square'}}
(v/'palette_audit/palette_decision.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
print('Harmonized state palette and captions against the actual manuscript figures.')