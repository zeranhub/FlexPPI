from pathlib import Path
import shutil
r=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010\section33_expansion_20261025\model_first_story_search_20261025');v=r/'revision_v83_model_story';f=v/'figure3';f.mkdir(parents=True,exist_ok=True)
shutil.copytree(r/'revision_v82_shared_view/figure3/source_data',f/'source_data',dirs_exist_ok=True)
a=f/'structure_assets';a.mkdir(exist_ok=True)
for name in ['render_manifest.json','whole_NHE1_trim.png','whole_IQCG_trim.png','detail_focus_M72_NHE1_trim.png','detail_focus_M72_IQCG_trim.png','detail_focus_M109_NHE1_trim.png','detail_focus_M109_IQCG_trim.png']:
    shutil.copy2(v/'pymol'/name,a/name)
s=(r/'revision_v82_shared_view/figure3/compose_fig3.py').read_text(encoding='utf-8');h=s[:s.index('# --- panel code ---')]
h=h.replace('# (c) shared-sequence profiles + contact tracks -> LineTrend/plot_trend.py -> param inherit','# (c) all-position Full profiles + fixed key-site neighbourhood boxes -> LineTrend/plot_trend.py -> param inherit').replace('# (d) group means + all5 continuous windows -> existing paired-point logic -> param inherit','# (d) Full means of four complete contact groups -> existing paired-point logic -> param inherit').replace('# (e) all9 Met absolute states + all5 windows -> existing paired-line logic -> param inherit','# (e) Full values for all9 Met -> existing paired-line logic -> param inherit').replace('# (f) all244 absolute-state MD / Full points -> existing correlation style -> param inherit','# (f) two key sites, paired Full RMSF and Hchi -> existing paired-point logic -> param inherit')
(f/'baseline_header.py').write_text(h,encoding='utf-8')
print('Prepared sources and inherited figure baselines.')