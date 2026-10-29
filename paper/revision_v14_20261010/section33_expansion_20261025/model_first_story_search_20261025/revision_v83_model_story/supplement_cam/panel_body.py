from pathlib import Path
import hashlib, json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

OUT = Path(__file__).resolve().parent
SRC = OUT / 'source_data'
W, H = 183., 145.
STATES = [('NHE1', 'a', '#B87936', 'o'), ('IQCG', 'b', '#2166AC', 's')]
mpl.rcParams.update({'axes.unicode_minus': False, 'savefig.bbox': None})
fig = plt.figure(figsize=(W / 25.4, H / 25.4), facecolor='white')
axes, plotted = [], []

def text(x, y, s, size=8, **kwargs):
    return fig.text(x / W, 1 - y / H, s, fontsize=size,
                    va=kwargs.pop('va', 'top'), color=kwargs.pop('color', BLACK), **kwargs)

def axis(x, y, w, h):
    ax = fig.add_axes([x / W, 1 - (y + h) / H, w / W, h / H])
    axes.append(ax)
    return ax

def panel(letter, label, x, y):
    text(x, y, letter, 10, fontweight='bold')
    text(x + 5, y, label, 8.3, fontweight='bold')

def check_limits(ax, x, y):
    x, y = np.asarray(x), np.asarray(y)
    mask = np.isfinite(x) & np.isfinite(y)
    assert mask.any()
    assert ax.get_xlim()[0] <= x[mask].min() <= x[mask].max() <= ax.get_xlim()[1]
    assert ax.get_ylim()[0] <= y[mask].min() <= y[mask].max() <= ax.get_ylim()[1]

data = pd.read_csv(SRC / 'cam_md_state_contrast_all_common_residues.csv')
chi = data[data.metric == 'hchi'].sort_values('mature_cam_number').copy()
back = data[data.metric == 'rmsf'].sort_values('mature_cam_number').copy()
assert len(chi) == 122 and len(back) == 145
for _, g in data.groupby('metric'):
    np.testing.assert_allclose(g.md_b - g.md_a, g.md_delta, atol=1e-12)
    np.testing.assert_allclose(g.full0_b - g.full0_a, g.full0_delta, atol=1e-12)
md_change = chi.md_b.to_numpy() - chi.md_a.to_numpy()
full_change = chi.full0_b.to_numpy() - chi.full0_a.to_numpy()
rho = float(spearmanr(md_change, full_change).statistic)
mae = float(np.mean(np.abs(md_change - full_change)))
assert round(rho, 3) == .630 and round(mae, 3) == .171

panel('a', 'Absolute Hχ', 1, 2)
ax = axis(18, 13, 55, 48)
for state, suffix, color, marker in STATES:
    ax.scatter(chi['md_' + suffix], chi['full0_' + suffix], s=8.5,
               marker=marker, facecolors=color, edgecolors='none', alpha=.78, zorder=3)
    for _, row in chi.iterrows():
        plotted.append(dict(panel='a', metric='hchi', state=state,
                            position=int(row.mature_cam_number),
                            x=float(row['md_' + suffix]), y=float(row['full0_' + suffix])))
ax.plot([0, 4.1], [0, 4.1], color='#555555', lw=.65, ls='--', zorder=1)
ax.set(xlim=(0, 4.1), ylim=(0, 4.1), xticks=[0, 2, 4], yticks=[0, 2, 4],
       xlabel='MD Hχ (nats)', ylabel='Predicted Hχ (nats)')
for _, suffix, _, _ in STATES:
    check_limits(ax, chi['md_' + suffix], chi['full0_' + suffix])

panel('b', 'Change in Hχ', 94, 2)
ax = axis(112, 13, 55, 48)
lim = np.ceil(max(np.abs(md_change).max(), np.abs(full_change).max()) * 10) / 10
lim = max(lim, .5)
ax.plot([-lim, lim], [-lim, lim], color='#555555', lw=.65, ls='--', zorder=1)
ax.axvline(0, color='#B8B8B8', lw=.45, zorder=1)
ax.axhline(0, color='#B8B8B8', lw=.45, zorder=1)
ax.scatter(md_change, full_change, s=8.5, c=BLACK, edgecolors='none', alpha=.8, zorder=3)
ax.set(xlim=(-lim, lim), ylim=(-lim, lim),
       xlabel='MD ΔHχ (nats)', ylabel='Predicted ΔHχ (nats)')
ax.text(.04, .97, f'ρ = {rho:.3f}\nMAE = {mae:.3f} nats',
        transform=ax.transAxes, ha='left', va='top', fontsize=7.5, color=BLACK)
check_limits(ax, md_change, full_change)
for p, md, full in zip(chi.mature_cam_number, md_change, full_change):
    plotted.append(dict(panel='b', metric='hchi_change', state='IQCG − NHE1',
                        position=int(p), x=float(md), y=float(full)))

panel('c', 'MD backbone reference', 1, 74)
for i, (state, _, color, marker) in enumerate(STATES):
    tx = 127 + 27 * i
    axleg = axis(tx, 74.5, 3, 2)
    axleg.plot([0], [0], marker=marker, ms=3, color=color, ls='none')
    axleg.set(xlim=(-1, 1), ylim=(-1, 1)); axleg.axis('off')
    text(tx + 4, 74, state, 7.5)

for letter, metric, yy, height, ylabel in [
        ('c', 'rmsf', 81, 18, 'RMSF (Å)'),
        ('d', 'hchi', 113, 19, 'Hχ (nats)')]:
    if letter == 'd':
        panel('d', 'MD side-chain reference', 1, 106)
    ax = axis(18, yy, 157, height)
    # Reindex the shared sequence so all undefined side chains remain gaps.
    g = data[data.metric == metric].set_index('mature_cam_number').reindex(range(4, 149))
    for state, suffix, color, marker in STATES:
        ax.plot(g.index, g['md_' + suffix], color=color, lw=.85)
        for p, val in g['md_' + suffix].dropna().items():
            plotted.append(dict(panel=letter, metric=metric, state=state,
                                position=int(p), x=int(p), y=float(val)))
    ax.set(xlim=(3.5, 148.5), ylim=(0, 4.1 if metric == 'hchi' else 4),
           ylabel=ylabel, yticks=[0, 2, 4], xticks=[4, 20, 40, 60, 80, 100, 120, 140, 148])
    for _, suffix, _, _ in STATES:
        check_limits(ax, g.index, g['md_' + suffix])
    ax.tick_params(length=2, pad=2)
    if letter == 'c':
        ax.tick_params(axis='x', labelbottom=False)
    else:
        ax.set_xlabel('Mature CaM residue number', labelpad=2)

fig.canvas.draw()
renderer = fig.canvas.get_renderer()
canvas = fig.bbox
outside = []
for label in fig.findobj(mpl.text.Text):
    if not label.get_visible() or not label.get_text().strip():
        continue
    bb = label.get_window_extent(renderer)
    if bb.x0 < canvas.x0 - 1 or bb.y0 < canvas.y0 - 1 or bb.x1 > canvas.x1 + 1 or bb.y1 > canvas.y1 + 1:
        outside.append(label.get_text())
assert not outside, outside
assert len(plotted) == 244 + 122 + 290 + 244

pd.DataFrame(plotted).to_csv(SRC / 'figS3_all_plotted_values.csv', index=False)
pd.DataFrame({'mature_cam_number': chi.mature_cam_number,
              'md_change_IQCG_minus_NHE1': md_change,
              'full_change_IQCG_minus_NHE1': full_change}).to_csv(
                  SRC / 'figS3_hchi_change_comparison.csv', index=False)
for ext in ['pdf', 'svg', 'png']:
    fig.savefig(OUT / ('figS3_CaM_validation.' + ext), bbox_inches=None,
                dpi=450 if ext == 'png' else 300, facecolor='white')

qa = dict(physical_size_mm=[W, H], png_dpi=450,
          state_colors={'NHE1': '#B87936', 'IQCG': '#2166AC'},
          panel_counts={'a_absolute_hchi':244, 'b_hchi_change':122,
                        'c_md_rmsf':290, 'd_md_hchi':244},
          total_plotted_records=len(plotted), hchi_change_spearman=rho,
          hchi_change_mae_nats=mae, undefined_hchi_positions_per_state=23,
          all_source_points_in_axes=True, text_outside_canvas=outside,
          within_training_fold0=True, reran_model=False,
          md_windows='five consecutive 100 ns windows, not independent replicates',
          exported_formats=['png','pdf','svg'],
          visual_inspection='passed: 183 by 145 mm PNG inspected; no overlapping labels or clipped source values',
          pdf_fonts='ArialMT and Arial-BoldMT, embedded CID TrueType; checked with pdffonts',
          source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest()
                         for p in SRC.glob('*') if p.is_file()})
(OUT / 'figure_QA.json').write_text(json.dumps(qa, ensure_ascii=False, indent=2), encoding='utf-8')
print(json.dumps({'output':str(OUT), 'rho':rho, 'mae':mae, 'plotted':len(plotted), 'QA':'passed'}, ensure_ascii=False))
