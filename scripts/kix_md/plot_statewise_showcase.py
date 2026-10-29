# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a) Paired residue maps -> heatmap/plot_comparison.py -> param inherit
# (b) Paired continuous ranks + shared-threshold outlines only -> heatmap -> param inherit
# (c) Molecular context -> render_statewise_showcase.py authentic PyMOL assets -> image loader
# (d) Residue strip + mean -> cross-type inherit from LineTrend/plot_sweep.py -> param inherit
# RULE: "native run" = load pre-rendered PNG via Image.open().ax.imshow().
#       "param inherit" = drawing function below that copies Class A/B/C values.
#       If a panel says "native run" and you write a drawing function, you broke the contract.
# The supplied templates have nested clinical counts / method optimization data,
# not paired residue coordinates and fixed molecular states. Semantic mismatch
# prevents native execution. Inherit frame-free cells, explicit sequential scale,
# direct annotations, and restrained point hierarchy. No template data are used.
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
import argparse
import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize, to_rgba
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
import matplotlib.patheffects as path_effects
from PIL import Image

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / 'artifacts/kix_md_20260914/pilot_results_20260916'
OUT = ROOT / 'output/flexppi_kix_statewise_showcase_20260916'
STATES = ['apo', 'myb', 'mll', 'both']
LABELS = ['Apo', 'cMyb', 'MLL', 'Both']
RESIDUES = np.arange(586, 673)
L12 = np.arange(614, 622)
WIDTH_MM, HEIGHT_MM = 183, 190
MD, MODEL, FOCUS = BLACK, CATEGORICAL[0], CATEGORICAL[3]
SHARED, MD_ONLY, MODEL_ONLY, NEITHER = '#08306B', '#77818A', '#6BAED6', '#EEF1F4'
MISSING = '#D5D9DD'
SEQ = LinearSegmentedColormap.from_list('within_state_rank', SEQUENTIAL)
SEQ.set_bad(MISSING)
mpl.rcParams.update({'savefig.pad_inches': 0, 'mathtext.fontset': 'custom',
                     'mathtext.rm': 'Arial', 'mathtext.it': 'Arial:italic',
                     'mathtext.bf': 'Arial:bold'})


def select(frame, **kwargs):
    mask = np.ones(len(frame), dtype=bool)
    for key, value in kwargs.items():
        mask &= frame[key].eq(value)
    return frame.loc[mask].copy()


def load_data():
    tables = {name: pd.read_csv(BASE / (name + '.csv'), dtype={'model_fold': str})
              for name in ['relative_residues', 'relative_metrics',
                           'relative_regions', 'top_quintile_overlap']}
    r = tables['relative_residues']
    assert len(r) == 4 * 6 * (87 + 74) * 4
    assert np.isfinite(r[['md_value', 'model_value']]).all().all()
    rows, summaries, membership = [], [], []
    for state in STATES:
        for metric, count in [('rmsf', 87), ('hchi', 74)]:
            raw = select(r, state=state, metric=metric, model_fold='ensemble_mean',
                         normalization='raw').sort_values('author_residue_number')
            pct = select(r, state=state, metric=metric, model_fold='ensemble_mean',
                         normalization='percentile').sort_values('author_residue_number')
            assert len(raw) == len(pct) == count
            assert not raw.author_residue_number.duplicated().any()
            assert np.array_equal(raw.author_residue_number, pct.author_residue_number)
            if metric == 'rmsf':
                assert np.array_equal(raw.author_residue_number, RESIDUES)
            else:
                assert not raw.amino_acid.isin(list('AGP')).any()
            for source in ['md', 'model']:
                expected = 100 * (rankdata(raw[source + '_value'], method='average') - 1) / (count - 1)
                np.testing.assert_allclose(expected, pct[source + '_value'], atol=1e-11)
            rho = float(spearmanr(raw.md_value, raw.model_value).statistic)
            recorded = select(tables['relative_metrics'], state=state, metric=metric,
                              model_fold='ensemble_mean', normalization='raw').iloc[0]
            assert abs(rho - recorded.spearman) < 1e-12
            data = raw[['state', 'metric', 'author_residue_number', 'amino_acid']].copy()
            data['md_raw'] = raw.md_value.to_numpy()
            data['flexppi_raw'] = raw.model_value.to_numpy()
            data['md_percentile'] = pct.md_value.to_numpy()
            data['flexppi_percentile'] = pct.model_value.to_numpy()
            rows.append(data)
            summary = dict(state=state, metric=metric, n_residues=count, spearman=rho,
                           md_independent_trajectories=1, model_checkpoints=5)
            if metric == 'rmsf':
                md_set = set(pct.loc[pct.md_value >= 80, 'author_residue_number'])
                model_set = set(pct.loc[pct.model_value >= 80, 'author_residue_number'])
                recorded_top = select(tables['top_quintile_overlap'], state=state,
                                      metric=metric, model_fold='ensemble_mean').iloc[0]
                assert len(md_set) == len(model_set) == 18
                assert len(md_set & model_set) == recorded_top.overlap
                assert md_set == {int(v) for v in recorded_top.md_residues.split(';')}
                assert model_set == {int(v) for v in recorded_top.model_residues.split(';')}
                region = pct[pct.author_residue_number.isin(L12)]
                assert len(region) == 8
                region_old = select(tables['relative_regions'], state=state, metric=metric,
                                    normalization='percentile', model_fold='ensemble_mean').iloc[0]
                np.testing.assert_allclose([region.md_value.mean(), region.model_value.mean()],
                                           [region_old.md_mean, region_old.model_mean], atol=1e-12)
                summary.update(md_high=18, flexppi_high=18, shared_high=len(md_set & model_set),
                               shared_fraction=len(md_set & model_set)/18,
                               l12_n=8, l12_md_mean_percentile=float(region.md_value.mean()),
                               l12_flexppi_mean_percentile=float(region.model_value.mean()),
                               l12_md_high=int((region.md_value >= 80).sum()),
                               l12_flexppi_high=int((region.model_value >= 80).sum()))
                for rid in RESIDUES:
                    a, b = rid in md_set, rid in model_set
                    category = 'Shared high' if a and b else 'MD only' if a else 'FlexPPI only' if b else 'Neither'
                    membership.append(dict(state=state, author_residue_number=int(rid),
                                           md_high=a, flexppi_high=b, category=category))
            summaries.append(summary)
    source = pd.concat(rows, ignore_index=True)
    summary = pd.DataFrame(summaries)
    assert len(source) == 4*(87+74)
    for metric, expected in [('rmsf', [.879, .897, .878, .905]),
                             ('hchi', [.961, .935, .949, .956])]:
        vals = [select(summary, state=s, metric=metric).spearman.iloc[0] for s in STATES]
        np.testing.assert_array_equal(np.round(vals, 3), expected)
    return tables, source, summary, pd.DataFrame(membership)


def text_mm(fig, x, y, text, **kw):
    return fig.text(x / WIDTH_MM, y / HEIGHT_MM, text, **kw)


def axes_mm(fig, x, y, w, h):
    return fig.add_axes([x / WIDTH_MM, y / HEIGHT_MM, w / WIDTH_MM, h / HEIGHT_MM])


def panel(fig, x, y, letter, title):
    text_mm(fig, x, y, letter, fontsize=10.5, fontweight='bold', va='bottom')
    text_mm(fig, x + 6, y, title, fontsize=8.5, fontweight='bold', va='bottom')


def plot_maps(fig, source, summary):
    panel(fig, 7, 168, 'a', 'Residue-level agreement across four states')
    for col, (metric, title, x) in enumerate([('rmsf', 'Main-chain RMSF', 26),
                                            ('hchi', r'Side-chain $H_\chi$', 105)]):
        w = 61
        text_mm(fig, x, 160, title, fontsize=8.2, fontweight='bold')
        text_mm(fig, x + w + 5.5, 160, 'Spearman ρ', fontsize=6.2, ha='center')
        ax = axes_mm(fig, x, 120, w, 36)
        for i, state in enumerate(STATES):
            group = select(source, state=state, metric=metric).set_index('author_residue_number')
            values = group[['md_percentile', 'flexppi_percentile']].reindex(RESIDUES).to_numpy().T
            # Blank group spacing is not missing data; grey cells are excluded Hchi positions.
            y_top = 4*i
            ax.pcolormesh(np.arange(585.5, 673.5), [y_top, y_top+1.65, y_top+3.3],
                          values, cmap=SEQ, vmin=0, vmax=100, shading='flat',
                          linewidth=0, edgecolors='face', rasterized=False)
            ax.hlines(y_top+1.65, 585.5, 672.5, color='white', linewidth=.35)
            if col == 0:
                y = 120 + 36*(1-(y_top+1.65)/15.3)
                text_mm(fig, 7, y, LABELS[i], fontsize=7.2, va='center')
                for dy, label in [(.825, 'MD'), (2.475, 'FlexPPI')]:
                    yy = 120 + 36*(1-(y_top+dy)/15.3)
                    text_mm(fig, 24.8, yy, label, ha='right', va='center', fontsize=6.1)
            rho = select(summary, state=state, metric=metric).spearman.iloc[0]
            ax.text(1.10, 1-(y_top+1.65)/15.3, f'{rho:.3f}', transform=ax.transAxes,
                    va='center', ha='center', fontsize=8.3, color=SHARED, fontweight='bold')
        ax.add_patch(Rectangle((613.5, 0), 8, 15.3, fill=False, ec=FOCUS, lw=.9))
        ax.text(617.5, -.5, 'L12', ha='center', va='bottom', fontsize=6.4, color='#A86816')
        ax.set_ylim(15.3, 0)
        ax.set_xlim(585.5, 672.5)
        ax.set_xticks([586, 610, 630, 650, 672])
        ax.tick_params(axis='x', length=0, pad=3, labelsize=6.2)
        ax.set_yticks([])
        ax.set_frame_on(False)
        ax.set_xlabel('KIX residue number', fontsize=6.5, labelpad=3)
    cax = axes_mm(fig, 140, 108, 26, 1.8)
    cb = fig.colorbar(mpl.cm.ScalarMappable(norm=Normalize(0, 100), cmap=SEQ), cax=cax,
                      orientation='horizontal')
    cb.outline.set_visible(False)
    cb.set_ticks([0, 50, 100])
    cb.ax.tick_params(labelsize=6, length=1.8, pad=1.3)
    text_mm(fig, 107, 108.4, 'Within-state percentile', fontsize=6.5)
    fig.patches.append(Rectangle((26/WIDTH_MM, 108/HEIGHT_MM), 2.8/WIDTH_MM, 1.8/HEIGHT_MM,
                                 transform=fig.transFigure, fc=MISSING, ec='none'))
    text_mm(fig, 30.5, 108.3, r'Excluded A/G/P in $H_\chi$', fontsize=6.3, color='#59616A')


def contiguous_ranges(residues):
    """Combine adjacent selected positions; never bridge an unselected residue."""
    residues = sorted(set(int(v) for v in residues))
    if not residues:
        return []
    ranges = []
    start = end = residues[0]
    for residue in residues[1:]:
        if residue == end + 1:
            end = residue
        else:
            ranges.append((start, end))
            start = end = residue
    ranges.append((start, end))
    assert {v for a,b in ranges for v in range(a,b+1)} == set(residues)
    return ranges


def draw_rank_pair(ax, source, membership, state, y_top):
    data = select(source, state=state, metric='rmsf').set_index('author_residue_number').loc[RESIDUES]
    values = data[['md_percentile','flexppi_percentile']].to_numpy().T
    assert values.shape == (2,87) and np.isfinite(values).all()
    ax.pcolormesh(np.arange(585.5,673.5), [y_top,y_top+1.65,y_top+3.3],
                  values, cmap=SEQ, vmin=0, vmax=100, shading='flat',
                  linewidth=0, edgecolors='face', rasterized=False)
    ax.hlines(y_top+1.65,585.5,672.5,color='white',linewidth=.45,zorder=2)
    calls = select(membership,state=state).set_index('author_residue_number').loc[RESIDUES]
    np.testing.assert_array_equal(calls.md_high,values[0]>=80)
    np.testing.assert_array_equal(calls.flexppi_high,values[1]>=80)
    # Colour encodes continuous rank. Only the shared threshold set is outlined;
    # a source-only threshold crossing is not a flexible-versus-rigid conflict.
    for category,dy,height,color,style in [
            ('Shared high',0,3.3,ACCENT_RED,'-')]:
        positions = calls.index[calls.category.eq(category)]
        for first,last in contiguous_ranges(positions):
            box=Rectangle((first-.5+.065,y_top+dy+.055),last-first+1-.13,height-.11,
                          fill=False,edgecolor=color,linewidth=.85,linestyle=style,zorder=4)
            box.set_path_effects([path_effects.Stroke(linewidth=1.6,foreground='white'),
                                  path_effects.Normal()])
            ax.add_patch(box)


def plot_top_sets(fig, source, summary, membership):
    panel(fig, 7, 98, 'b', 'RMSF | Continuous flexibility ranks and shared top-ranked sites')
    for x,label,color,style in [(26,'Both ≥80th percentile',ACCENT_RED,'-')]:
        fig.patches.append(Rectangle((x/WIDTH_MM,92/HEIGHT_MM),5/WIDTH_MM,2.1/HEIGHT_MM,
                                     transform=fig.transFigure,fill=False,ec=color,lw=.85,ls=style))
        text_mm(fig,x+7,92.2,label,fontsize=6.6)
    text_mm(fig,113,92.2,'Colour scale as a',fontsize=6.4,color='#59616A')
    text_mm(fig, 164, 92.2, 'Overlap', fontsize=7, ha='center', fontweight='bold')
    ax = axes_mm(fig, 26, 67, 121, 21)
    for i, state in enumerate(STATES):
        y_top=4*i
        draw_rank_pair(ax,source,membership,state,y_top)
        s = select(summary, state=state, metric='rmsf').iloc[0]
        yy = 67+21*(1-(y_top+1.65)/15.3)
        text_mm(fig, 7, yy, LABELS[i], fontsize=7.2, va='center')
        for dy,label in [(.825,'MD'),(2.475,'FlexPPI')]:
            text_mm(fig,24.8,67+21*(1-(y_top+dy)/15.3),label,fontsize=6,
                    ha='right',va='center')
        text_mm(fig, 161, yy+.25, f'{int(s.shared_high)}/18', fontsize=9, color=SHARED,
                fontweight='bold', ha='center', va='center')
        text_mm(fig, 176, yy+.2, f'{s.shared_fraction:.1%}', fontsize=6.7,
                color='#59616A', ha='right', va='center')
    ax.set_ylim(15.3,0)
    ax.set_xlim(585.5,672.5)
    ax.set_xticks([586,600,620,640,660,672])
    ax.tick_params(axis='x',labelsize=6.2,length=0,pad=3)
    ax.set_yticks([])
    ax.set_frame_on(False)


def export_panel_b_detail(source, summary, membership):
    """Standalone readable residue index, using exactly the same values/outlines."""
    fig=plt.figure(figsize=(183/25.4,91/25.4),facecolor='white')
    fig.text(.04,.93,'KIX | RMSF percentile and high-flexibility agreement',fontsize=10,fontweight='bold')
    handles=[Rectangle((0,0),1,1,fill=False,ec=ACCENT_RED,lw=.9,
                       label='Both ≥80th percentile: solid box across MD and FlexPPI')]
    fig.legend(handles=handles,loc='upper left',bbox_to_anchor=(.032,.89),ncol=2,
               fontsize=6.8,handlelength=1.8,columnspacing=1.5)
    ax=fig.add_axes([.15,.235,.77,.52])
    for i,state in enumerate(STATES):
        draw_rank_pair(ax,source,membership,state,4*i)
        for dy,label in [(.825,'MD'),(2.475,'FlexPPI')]:
            ax.text(-.01,1-(4*i+dy)/15.3,label,transform=ax.transAxes,
                    ha='right',va='center',fontsize=6.5)
        ax.text(-.135,1-(4*i+1.65)/15.3,LABELS[i],transform=ax.transAxes,
                ha='left',va='center',fontsize=7.5)
        count=int(select(summary,state=state,metric='rmsf').shared_high.iloc[0])
        ax.text(1.018,1-(4*i+1.65)/15.3,f'{count}/18',transform=ax.transAxes,
                va='center',fontsize=8,fontweight='bold',color=SHARED)
    ax.set_xlim(585.5,672.5)
    ax.set_ylim(15.3,0)
    ax.set_yticks([])
    ax.set_xticks(RESIDUES)
    ax.set_xticklabels(RESIDUES,rotation=90,fontsize=6)
    ax.tick_params(axis='x',length=1.5,pad=2)
    ax.set_frame_on(False)
    ax.set_xlabel('KIX residue number',fontsize=7,labelpad=5)
    cax=fig.add_axes([.74,.095,.18,.025])
    cb=fig.colorbar(mpl.cm.ScalarMappable(norm=Normalize(0,100),cmap=SEQ),
                    cax=cax,orientation='horizontal')
    cb.outline.set_visible(False)
    cb.set_ticks([0,50,80,100]);cb.ax.tick_params(labelsize=6,length=2,pad=1.5)
    fig.text(.53,.098,'Within-state RMSF percentile',fontsize=6.8)
    fig.text(.04,.026,'Boxes mark shared threshold membership only; unboxed does not mean rigid. All 87 residues retained.',
             fontsize=6.4,color='#59616A')
    fig.add_artist(Rectangle((0,0),1,1,transform=fig.transFigure,fill=False,edgecolor='none'))
    target=OUT/'panel_b_continuous_heatmap_detail'
    save_cns_figure(fig,str(target))
    fig.savefig(str(target)+'.svg',bbox_inches='tight',dpi=300)
    plt.close(fig)


def white_canvas_bbox(path):
    image = np.array(Image.open(path).convert('RGB'))
    nonwhite = np.any(image < 248, axis=2)
    ys, xs = np.where(nonwhite)
    assert len(xs) > 0, str(path)
    pad = 18
    return [max(0,xs.min()-pad), max(0,ys.min()-pad),
            min(image.shape[1],xs.max()+pad+1), min(len(image),ys.max()+pad+1)]


def image_panel(fig, rect, path, crop=None):
    ax = axes_mm(fig, *rect)
    image = np.array(Image.open(path).convert('RGB'))
    crop = white_canvas_bbox(path) if crop is None else crop
    # Only outer white canvas is cropped. The structure pixels are not edited.
    image = image[crop[1]:crop[3],crop[0]:crop[2]]
    ax.imshow(image)
    ax.axis('off')
    ax.source_crop = crop
    return ax


def plot_structure(fig):
    panel(fig, 7, 56, 'c', 'L12 in the ternary complex (Both)')
    asset = OUT / 'structure_assets'
    whole = image_panel(fig, (6, 14, 44, 38), asset/'context_whole.png')
    provenance = json.loads((asset/'render_provenance.json').read_text(encoding='utf-8'))
    xy = provenance['L12_visible_screen_target']['normalized_top_left_origin']
    x0,y0,_,_ = whole.source_crop
    target = (xy[0]*1400-x0,xy[1]*1100-y0)
    whole.annotate('L12', xy=target, xycoords='data', xytext=(.03,.99),
                   textcoords='axes fraction', fontsize=7, color='#A86816',
                   fontweight='bold', ha='left', va='bottom',
                   arrowprops=dict(arrowstyle='-',color=FOCUS,lw=.8,
                                   connectionstyle='angle3,angleA=0,angleB=90'))
    # Identical crop and placement dimensions preserve MD/model scale exactly.
    boxes=np.array([white_canvas_bbox(asset/f'l12_{s}.png') for s in ['md','model']])
    common_crop=[boxes[:,0].min(),boxes[:,1].min(),boxes[:,2].max(),boxes[:,3].max()]
    image_panel(fig, (54, 34, 39, 20), asset/'l12_md.png',common_crop)
    image_panel(fig, (54, 12, 39, 20), asset/'l12_model.png',common_crop)
    text_mm(fig, 53.5, 52, 'MD', fontsize=7.4, fontweight='bold')
    text_mm(fig, 53.5, 30, 'FlexPPI', fontsize=7.4, fontweight='bold', color=MODEL)
    text_mm(fig, 7, 12, 'L12: 614–621', fontsize=7, color='#A86816', fontweight='bold')
    for x,label,color in [(7,'KIX','#CCD0D4'),(20,'cMyb','#AFA0BD'),(36,'MLL','#D0BB94')]:
        fig.patches.append(Rectangle((x/WIDTH_MM,8/HEIGHT_MM),1.8/WIDTH_MM,1.8/HEIGHT_MM,
                                     transform=fig.transFigure,fc=color,ec='none'))
        text_mm(fig,x+2.7,8,label,fontsize=6.1,color='#59616A')
    text_mm(fig, 53.5, 9.2, 'B607–628; sticks: L12', fontsize=6.1, color='#59616A')
    text_mm(fig, 53.5, 6.3, 'RMSF percentile scale as a', fontsize=6.1, color='#59616A')


def plot_l12(fig, source, summary):
    panel(fig, 103, 56, 'd', 'L12 is locally flexible in all four states')
    ax = axes_mm(fig, 115, 16, 58, 33)
    ax.axvspan(80,100,facecolor='#F2F6FA',zorder=0)
    ax.axvline(80,color='#8997A7',ls=(0,(3,3)),lw=.65,zorder=1)
    for i,state in enumerate(STATES):
        y=6-2*i
        g=select(source,state=state,metric='rmsf').set_index('author_residue_number').loc[L12]
        # Fixed residue-order offsets ensure all eight residues remain visible;
        # offsets convey identity/order only, not another measured variable.
        offsets=np.linspace(-.12,.12,8)
        for col,color,marker,dy in [('md_percentile',MD,'o',.28),
                                    ('flexppi_percentile',MODEL,'s',-.28)]:
            ax.scatter(g[col],y+dy+offsets,s=9,color=color,alpha=.35,
                       marker=marker,linewidths=0,zorder=2)
            mean=float(g[col].mean())
            ax.scatter(mean,y+dy,s=35,color=color,marker=marker,
                       edgecolor='white',linewidth=.55,zorder=4)
            ax.text(103,y+dy,f'{mean:.1f}',fontsize=7,ha='left',va='center',color=color,
                    clip_on=False)
    ax.set_xlim(68,105)
    ax.set_ylim(-.8,6.9)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_bounds(70,100)
    ax.set_yticks([6,4,2,0],LABELS)
    ax.tick_params(axis='y',length=0,pad=7,labelsize=7)
    ax.set_xticks([70,80,90,100])
    ax.tick_params(axis='x',length=2,labelsize=6.7)
    ax.set_xlabel('RMSF percentile within KIX',fontsize=7,labelpad=4)
    ax.legend(handles=[Line2D([],[],color=MD,ls='',marker='o',ms=4,label='MD'),
                       Line2D([],[],color=MODEL,ls='',marker='s',ms=4,label='FlexPPI')],
              loc='lower left',bbox_to_anchor=(-.16,1.01),ncol=2,fontsize=6.6,
              handletextpad=.4,columnspacing=1.4,borderaxespad=0)
    text_mm(fig,171.5,50,'Mean',fontsize=6.3,ha='center',color='#59616A')
    text_mm(fig,104,6.3,'Small points: 8 residues; large points: region mean.',fontsize=6.1,color='#59616A')


def compose_figure(source, summary, membership, archetype='asymmetric_mixed'):
    assert archetype == 'asymmetric_mixed'
    fig=plt.figure(figsize=(WIDTH_MM/25.4,HEIGHT_MM/25.4),facecolor='white')
    text_mm(fig,7,183,'FlexPPI recovers the local flexibility landscape of KIX',
            fontsize=10.8,fontweight='bold',color=BLACK)
    text_mm(fig,7,177,'Four molecular states  |  MD reference and a frozen five-model ensemble',
            fontsize=7.2,color='#59616A')
    plot_maps(fig,source,summary)
    plot_top_sets(fig,source,summary,membership)
    plot_structure(fig)
    plot_l12(fig,source,summary)
    # Single, unobtrusive run-design line; extended scope is in the companion caption.
    text_mm(fig,7,2,'One 100 ns MD trajectory/state. RMSF: 87 residues; Hχ: 74. Both = cMyb + MLL. See caption for validation scope.',
            fontsize=6,color='#69727C')
    fig.add_artist(Rectangle((0,0),1,1,transform=fig.transFigure,fill=False,edgecolor='none'))
    fig.canvas.draw()
    renderer=fig.canvas.get_renderer()
    # Detect accidental off-canvas text before export (excluding intentionally hidden artists).
    bounds=fig.bbox
    overflow=[]
    for obj in fig.findobj(mpl.text.Text):
        if not obj.get_visible() or not obj.get_text():
            continue
        bb=obj.get_window_extent(renderer)
        if bb.x0 < -1 or bb.y0 < -1 or bb.x1 > bounds.x1+1 or bb.y1 > bounds.y1+1:
            overflow.append(obj.get_text())
    assert not overflow, 'Text outside canvas: '+str(overflow)
    filename=OUT/'fig01_kix_statewise_validation'
    save_cns_figure(fig,str(filename))
    fig.savefig(str(filename)+'.svg',bbox_inches='tight',dpi=300)
    plt.close(fig)


def write_sources(tables, source, summary, membership):
    source_dir=OUT/'source_data'
    source_dir.mkdir(parents=True,exist_ok=True)
    source.to_csv(source_dir/'statewise_residues.csv',index=False)
    summary.to_csv(source_dir/'statewise_summary.csv',index=False)
    membership.to_csv(source_dir/'rmsf_high_flexibility_membership.csv',index=False)
    range_rows=[]
    for state in STATES:
        group=select(membership,state=state)
        for category in ['Shared high']:
            for first,last in contiguous_ranges(group.loc[group.category.eq(category),'author_residue_number']):
                range_rows.append(dict(state=state,category=category,start_residue=first,
                                       end_residue=last,n_residues=last-first+1))
    pd.DataFrame(range_rows).to_csv(source_dir/'panel_b_outline_ranges.csv',index=False)
    source[source.author_residue_number.isin(L12)&source.metric.eq('rmsf')].to_csv(source_dir/'l12_residues.csv',index=False)
    # All five checkpoints retained for reproducibility, not represented as MD replicates.
    r=tables['relative_residues']
    r[r.normalization.isin(['raw','percentile'])].to_csv(source_dir/'all_five_checkpoints_and_ensemble.csv',index=False)
    inputs={name+'.csv':hashlib.sha256((BASE/(name+'.csv')).read_bytes()).hexdigest() for name in tables}
    manifest=dict(scope='State-wise residue ranking, fixed-threshold detection and L12 localization only',
                  revision='Panel b: continuous MD/model percentile pairs; source-only outlines removed',
                  panel_b='All 87 residues × two sources × four states; common 0–100 scale; solid boxes only for both sources ≥80th percentile; unboxed does not mean rigid',
                  panel_design_confirmed_by_user=True,normalization='100*(average_rank-1)/(N-1), within state and source',
                  ensemble='Mean five frozen raw predictions first; calculate percentile afterwards',
                  md_independent_trajectories_per_state=1,md_production_ns_per_state=100,
                  rmsf_n=87,hchi_n=74,hchi_excluded='A/G/P: 13 positions retained as grey cells',
                  l12='B614–621: all eight RMSF positions',threshold=80,
                  threshold_selection='Fixed across all states; post-readout analysis, not preregistered',
                  input_sha256=inputs,source_directory=str(BASE),
                  script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  width_mm=WIDTH_MM,height_mm=HEIGHT_MM,formats=['PDF','SVG','300dpi PNG'],
                  uncertainty='No MD error bars, p values, or confidence intervals from one trajectory',
                  training_scope='Near-identical KIX sequence present in training; not an independent-family test')
    (OUT/'figure_manifest.json').write_text(json.dumps(manifest,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--data-only',action='store_true')
    args=parser.parse_args()
    OUT.mkdir(parents=True,exist_ok=True)
    tables,source,summary,membership=load_data()
    write_sources(tables,source,summary,membership)
    if not args.data_only:
        compose_figure(source,summary,membership)
        export_panel_b_detail(source,summary,membership)
    print(summary.to_string(index=False))
    print(OUT)
