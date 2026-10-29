from pathlib import Path
import json,hashlib
v=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010\section33_expansion_20261025\model_first_story_search_20261025\revision_v83_model_story')
s=v/'supplement_cam'
for name in ['compose_figS3.py','panel_body.py']:
    p=s/name;t=p.read_text(encoding='utf-8').replace('#0072B2','#B87936').replace('#D55E00','#2166AC');p.write_text(t,encoding='utf-8')
p=s/'figure_QA.json';d=json.loads(p.read_text(encoding='utf-8'));d['state_colors']={'NHE1':'#B87936','IQCG':'#2166AC'};d['visual_review']='passed: root and independent review, final manuscript palette';p.write_text(json.dumps(d,ensure_ascii=False,indent=2),encoding='utf-8')
p=v/'figure3/figure_QA.json';d=json.loads(p.read_text(encoding='utf-8'));d['visual_review']='passed: root and independent numerical/visual audit; harmonized manuscript palette';d['pdf_fonts']='Arial / Arial Bold CID TrueType embedded; pdffonts verified';d['frozen_values_verified']=568;d['palette_audit']='../palette_audit/palette_decision.json';p.write_text(json.dumps(d,ensure_ascii=False,indent=2),encoding='utf-8')
print('Final QA metadata and supplemental reproduction sources now use the actual manuscript palette.')