from pathlib import Path
import json
v=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010\section33_expansion_20261025\model_first_story_search_20261025\revision_v82_shared_view');m=v/'manuscript';p=m/'v82_text_plan.json';plan=json.loads(p.read_text(encoding='utf-8'));old=plan['replacements']['75']
new='CaM两态的叶瓣运动与侧链状态呈区域性再分配。叶瓣独立拟合时，N叶瓣平均波动未随整链RMSF同步大幅下降，说明整链读数含叶间运动贡献。完整模型恢复了主要侧链变化，但低估主链变化幅度，连接区部分方向也不准，不能据此声称其主链预测优于链内模型。伙伴、构型与原始钙条件同时不同，因此不作更换伙伴或去钙的单因素归因。轨迹复算、叶瓣拟合与强静态基线见补充材料D。'
plan['replacements']['75']=new;p.write_text(json.dumps(plan,ensure_ascii=False,indent=2),encoding='utf-8')
b=m/'build_v82_manuscript.py';s=b.read_text(encoding='utf-8');needle="        for text in plan.get('insert_before_158', []):";assert needle in s;s=s.replace(needle,"        # Keep the full Figure 3 and caption on their own page, then resume the case narrative.\n        paragraph_setting(paragraphs[73], 'pageBreakBefore', val=1)\n"+needle);b.write_text(s,encoding='utf-8')
print({'old_paragraph75_characters':len(old),'new_paragraph75_characters':len(new)})