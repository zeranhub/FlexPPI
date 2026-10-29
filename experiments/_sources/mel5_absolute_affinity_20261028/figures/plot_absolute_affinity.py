# Asset Confirmation Table
# Peptide profiles: LineTrend production asset is cumulative-publication data,
# incompatible with residue differences; inherit clean axes and external legend.
# Complete paired scatter: basic-plots contains no production implementation.
# Regional replicate plots: basic-plots contains no production implementation.
# All panels retain all predefined positions; only chemically undefined Hchi is masked.

# Academic Figure Skill Typography Baseline — COPY VERBATIM, place at TOP of script
import matplotlib as mpl
mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans"],
    "font.size": 8,
    "axes.titlesize": 8,
    "axes.labelsize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 8,
    "figure.titlesize": 9,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.6,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "legend.frameon": False,
})

# Academic Figure Skill Nature/Cell/Science Color Palette -- COPY VERBATIM
CATEGORICAL = ["#2166AC", "#B2182B", "#1B7837", "#F1A340", "#762A83", "#666666"]
CATEGORICAL_EXTENDED = [
    "#2166AC", "#B2182B", "#1B7837", "#F1A340", "#762A83", "#666666",
    "#4393C3", "#D6604D", "#5AAE61", "#B35806", "#9970AB", "#999999",
]
DIVERGING   = ["#2166AC", "#F7F7F7", "#B2182B"]
SEQUENTIAL  = ["#F7FBFF", "#6BAED6", "#08306B"]
ACCENT_RED  = "#B2182B"
GREY        = "#999999"
BLACK       = "#222222"

# Academic Figure Skill Export Baseline — COPY VERBATIM
mpl.rcParams.update({
    "pdf.fonttype": 42,         # TrueType font embedding
    "svg.fonttype": "none",     # editable text in SVG
    "savefig.bbox": "tight",    # trim whitespace
    "savefig.dpi": 300,
})

def save_cns_figure(fig, filename):
    """Standard Academic Figure Skill export: vector PDF + 300dpi PNG preview."""
    fig.savefig(f"{filename}.pdf", bbox_inches="tight", dpi=300)
    fig.savefig(f"{filename}.png", bbox_inches="tight", dpi=300)

from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

ROOT = Path(__file__).resolve().parent
source = ROOT.parent / 'affinity' / 'absolute_affinity_results.csv'
data = pd.read_csv(source)
order = ['3hg1', '4jff', '4jfd', '4jfe']
assert set(data.pdb) == set(order)
assert np.isfinite(data.predicted_pKD).all()
assert np.allclose(data.predicted_pKD, -np.log10(data.predicted_KD_M), atol=1e-12)
measured = data.drop_duplicates('pdb').set_index('pdb').loc[order]
assert len(measured) == 4
methods = list(dict.fromkeys(data.model))
fig, ax = plt.subplots(figsize=(183 / 25.4, 75 / 25.4))
fig.subplots_adjust(left=0.10, right=0.985, bottom=0.29, top=0.82)
x = np.arange(4)
offsets = np.linspace(-0.19, 0.19, len(methods) + 1)
ax.scatter(x + offsets[0], measured.experimental_pKD, color=BLACK, marker='D', s=31,
           label='SPR experiment', zorder=4)
colors = [CATEGORICAL[0], CATEGORICAL[4], CATEGORICAL[2]]
markers = ['o', '^', 's']
for j, method in enumerate(methods):
    d = data[data.model == method].set_index('pdb').loc[order]
    assert len(d) == 4
    ax.scatter(x + offsets[j + 1], d.predicted_pKD, color=colors[j], marker=markers[j],
               s=32, label=method, zorder=3)
ax.set_xticks(x, ['MEL5 / ELA\n3HG1', 'Mature / ELA\n4JFF', 'Mature / G4A\n4JFD', 'Mature / I7A\n4JFE'])
ax.set_ylabel('pKD = −log10(KD / M)')
ax.set_xlim(-0.5, 3.5)
y = np.r_[data.predicted_pKD.to_numpy(), measured.experimental_pKD.to_numpy()]
ax.set_ylim(np.floor(y.min() - 0.3), np.ceil(y.max() + 0.3))
ax.legend(loc='lower center', bbox_to_anchor=(0.5, 1.05), ncol=min(3, len(methods)+1),
          frameon=False, fontsize=7.5, columnspacing=1.6)
fig.text(0.10, 0.04, 'Higher pKD indicates stronger binding. One static structure per state; no fitted calibration.', fontsize=7)
save_cns_figure(fig, ROOT / 'figS5_absolute_affinity')
fig.savefig(ROOT / 'figS5_absolute_affinity.svg', bbox_inches='tight', dpi=300)
plt.close(fig)
report = {'source': str(source), 'sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
          'structures': order, 'methods': methods, 'units': 'pKD = -log10(KD/M)',
          'independence': '4 related states from one TCR family, not 4 independent families',
          'statistical_tests': 'None; retrospective per-state comparison, no four-point performance p-value',
          'error_bars': 'None: deterministic structure predictions and published point KD estimates',
          'temperature': 'PRODIGY KD converted at 25 C, matching SPR; other model temperature per author training',
          'complete_results': data.to_dict('records')}
(ROOT / 'affinity_figure_statistics.json').write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')
print(json.dumps({'status':'rendered', 'models':methods, 'states':4}, ensure_ascii=False))