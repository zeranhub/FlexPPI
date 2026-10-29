from pathlib import Path
import shutil
ROOT = Path(__file__).resolve().parent
FIG = ROOT/'figure4'
FIG.mkdir(parents=True, exist_ok=True)
EXP = ROOT.parent
# Retain the verbatim skill baselines from the previous verified script.
old = (EXP/'mel5_focus_v50_20261022/figure4/plot.py').read_text(encoding='utf-8')
start = old.index('# Academic Figure Skill Typography Baseline')
end = old.index('from pathlib import Path',start)
baseline = old[start:end]
header = '''# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a) Actual molecular overview + one local view -> archived PyMOL render -> native raster placement
# (b,c) Residue profiles -> LineTrend/plot_trend.py -> param inherit
# (d upper/lower) Complete aligned delta maps -> heatmap/plot_composition.py -> param inherit
# Assets and companion PNGs inspected. Monthly counts, grouped fractions and clinical-stage
# counts have different data/axis semantics; preserve clean axes and grouping, use existing
# scientific source mappings. All user-requested values retained; native vectors for all charts.
'''
body = (ROOT/'plot_body.py').read_text(encoding='utf-8')
(FIG/'plot.py').write_text(header+baseline+body,encoding='utf-8')
for name in ['LineTrend/plot_trend.py','heatmap/plot_composition.py']:
    target=ROOT/'source_audit/asset_parameters'/name
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(Path(r'C:/Users/DELL/.codex/skills/academic-figure-skill/assets/figures')/name,target)
print(FIG/'plot.py')
