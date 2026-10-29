from pathlib import Path
import shutil,json
r=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010\section33_expansion_20261025\model_first_story_search_20261025');v=r/'revision_v84_effect_audit';v.mkdir(exist_ok=True);p=v/'source_data';p.mkdir(exist_ok=True)
shutil.copy2(r/'revision_v83_model_story/figure3/source_data/cam_md_state_contrast_all_common_residues.csv',p/'cam_md_state_contrast_all_common_residues.csv')
s=(r/'revision_v83_model_story/figure3/compose_fig3.py').read_text(encoding='utf-8');h=s[s.index('# Academic Figure Skill Typography'):s.index('# --- panel code ---')]
header='# Academic Figure Skill Asset Confirmation (existing scientific profile and paired-point logic)\n# Full sequence -> LineTrend inherited; region/site absolute pairs -> existing paired-point logic inherited.\n# All frozen-model residue values retained. Fixed 9-residue regions are descriptive model-selected examples.\n'+h
(v/'baseline_header.py').write_text(header,encoding='utf-8')
print(str(v))