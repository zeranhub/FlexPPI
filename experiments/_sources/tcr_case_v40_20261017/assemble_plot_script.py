from pathlib import Path
ROOT=Path(__file__).resolve().parent
source=ROOT.parent/'reader_revision_v39_20261016/english_figures/fig7/plot_english.py'
text=source.read_text(encoding='utf8')
baseline=text[:text.index('import hashlib,json,xml')]
baseline=baseline.replace("Path(__file__).resolve().parents[3]","Path(__file__).resolve().parent.parent")
(ROOT/'plot_case.py').write_text(baseline+(ROOT/'plot_case_body.py').read_text(encoding='utf8'),encoding='utf8')
