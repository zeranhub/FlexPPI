# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a) Molecular interface plate -> existing v49 real-coordinate renderer + v51 same-view renderer -> param inherit
# (b) Complete CDR curves -> LineTrend/plot_trend.py and inspected trend_by_month.png -> param inherit
# Original asset is a cumulative-month count, incompatible with paired residue RMSF.
# Retained semantic line comparison, direct annotations, external legend, no top/right spines.
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
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
    "savefig.bbox": "tight",
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
from PIL import Image
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, ConnectionPatch
from matplotlib.text import Text

ROOT = Path(__file__).resolve().parent
BASE = ROOT.parent
OLD = BASE.parent / 'mel5_joint_case_v49_20261019'
META = BASE / 'structure/dual_state_molecular_view.json'
meta = json.loads(META.read_text(encoding='utf-8'))
geometry = pd.read_csv(BASE / 'analysis/fixed_atom_contact_distances.csv')
cdr = pd.read_csv(OLD / 'figure4/source_plotted_CDR_RMSF.csv', float_precision='round_trip')
interface = pd.read_csv(OLD / 'figure4/source_plotted_interface_Hchi.csv', float_precision='round_trip')
key = json.loads((OLD / 'case_sources/residue_story/trajectory_link/key_findings.json').read_text(encoding='utf-8'))
WT = '#CF8C27'
ENG = CATEGORICAL[4]
TEAL = '#297B82'
GRAY = '#6B6B6B'
assert len(cdr) == 45 and len(interface) == 23
assert np.allclose(geometry.WT_distance_A, [5.848769101272505, 10.584427240054127], atol=1e-8)
assert np.allclose(geometry.mature_distance_A, [2.942889736296623, 3.119551249779373], atol=1e-8)
assert meta['screen_projection']['holdout_error_pixels'] < 1

# One shared crop from the union of both unmodified camera renders. Neither
# state is independently zoomed, stretched, cropped around a selected atom,
# edited or recolored during composition.
ims = {state: Image.open(meta['structures'][state]['PNG']).convert('RGBA') for state in ['wt', 'eng']}
assert ims['wt'].size == ims['eng'].size == tuple(meta['image_size'])
bounds = []
for im in ims.values():
    arr = np.asarray(im)
    visible = (arr[:, :, 3] > 8) & np.any(arr[:, :, :3] < 242, axis=2)
    yy, xx = np.where(visible)
    bounds.append([xx.min(), yy.min(), xx.max(), yy.max()])
crop = [max(0, min(b[0] for b in bounds)-40), max(0, min(b[1] for b in bounds)-95),
        min(ims['wt'].width, max(b[2] for b in bounds)+40), min(ims['wt'].height, max(b[3] for b in bounds)+30)]
crop = [int(value) for value in crop]

def point(state, name):
    return np.asarray(meta['structures'][state]['atom_screen_pixels_top_left'][name], float)

def note(ax, text, target, xytext, color=BLACK, fontsize=7.1, ha='left', weight='normal'):
    return ax.annotate(text, xy=target, xycoords='data', xytext=xytext,
                       textcoords='axes fraction', ha=ha, va='center', fontsize=fontsize,
                       color=color, fontweight=weight,
                       arrowprops=dict(arrowstyle='-', lw=.55, color=color),
                       bbox=dict(facecolor='white', edgecolor='none', pad=1.1), zorder=10)

def local_view(fig, state, rect):
    ax = fig.add_axes(rect)
    ax.imshow(ims[state], origin='upper', interpolation='lanczos')
    ax.set_xlim(crop[0], crop[2])
    ax.set_ylim(crop[3], crop[1])
    ax.axis('off')
    mature = state == 'eng'
    color = ENG if mature else WT
    # Residue-labelled model outputs stay attached to the residues that carry
    # them. No connecting causal arrow is drawn between the two predictions.
    suffix = 'engineered' if mature else 'WT'
    h = float(key['model'][f'Q72_Hchi_nats_{suffix}'])
    r = float(key['model'][f'beta_Gly_RMSF_A_{suffix}'])
    note(ax, f'HLA Q72\nHχ {h:.3f} nats', point(state, 'Q72_CA'),
         (.01, .18), color=TEAL, fontsize=7.1, weight='bold')
    note(ax, f'βGly54\nRMSF {r:.3f} Å', point(state, 'beta54_CA'),
         (.98, .98), color=color, ha='right', fontsize=7.1, weight='bold')
    sidechain = 'beta53_sidechain_center' if 'beta53_sidechain_center' in meta['structures'][state]['atom_screen_pixels_top_left'] else 'beta53_CA'
    note(ax, 'βPhe53' if mature else 'βIle53', point(state, sidechain),
         (.17, .965), color=color, fontsize=7.5, weight='bold')
    anchor52 = 'beta52_sidechain_center' if 'beta52_sidechain_center' in meta['structures'][state]['atom_screen_pixels_top_left'] else 'beta52_CA'
    note(ax, 'βPro52' if mature else 'βGly52', point(state, anchor52),
         (.03, .57), color=color, fontsize=6.8)
    # Each distance label is anchored at its own named atom-pair midpoint.
    for i, (first, second, xoffset, yoffset) in enumerate([
        ('Q72_OE1', 'beta53_N', -110 if mature else -90, -20),
        ('Q72_NE2', 'beta54_O', 85 if mature else 110, 10),
    ]):
        midpoint = (point(state, first) + point(state, second)) / 2
        distance = float(geometry.iloc[i]['mature_distance_A' if mature else 'WT_distance_A'])
        ax.text(midpoint[0]+xoffset, midpoint[1]+yoffset, f'{distance:.2f} Å',
                fontsize=7.3, ha='center', va='center', color=BLACK,
                bbox=dict(facecolor='white', edgecolor='none', pad=1.0), zorder=11)
    for name, label, shift in [
        ('Q72_OE1', 'OE1', (-25, -28)),
        ('Q72_NE2', 'NE2', (30, 55) if mature else (-115, 12)),
        ('beta53_N', 'N', (-120, -32)),
        ('beta54_O', 'O', (38, 0)),
    ]:
        p = point(state, name)
        ax.text(p[0]+shift[0], p[1]+shift[1], label, fontsize=6.0, color=GRAY,
                bbox=dict(facecolor='white', edgecolor='none', pad=.3), zorder=11)
    return ax

def compose_figure():
    fig = plt.figure(figsize=(183/25.4, 140/25.4), facecolor='white')
    fig.text(.02, .975, 'a', fontsize=11, fontweight='bold')
    fig.text(.062, .975, 'Affinity maturation reorganizes the Q72–CDR2β contact network', fontsize=8.5)
    # Small biological locator. The actual before/after interface occupies the
    # majority of panel a and carries its main conclusion.
    ov = Image.open(OLD / 'figure4/4jff_overview.png').convert('RGBA')
    oa = fig.add_axes([.028, .615, .185, .285]); oa.imshow(ov); oa.axis('off')
    oa.text(.06, .96, 'TCR α', transform=oa.transAxes, color=GRAY, fontsize=6.5)
    oa.text(.62, .96, 'TCR β', transform=oa.transAxes, color=ENG, fontsize=6.5)
    oa.add_patch(Circle((.49, .56), .085, transform=oa.transAxes, fill=False, lw=.85, ec=TEAL))
    fig.text(.029, .59, 'HLA-A*02:01', fontsize=7.1, color=TEAL)
    fig.text(.029, .564, '+ ELA peptide', fontsize=6.9, color=GRAY)
    fig.text(.029, .515, 'SPR affinity [28]', fontsize=7.1, fontweight='bold')
    fig.text(.029, .484, '18 μM → 0.6 nM', fontsize=7.8)
    fig.text(.029, .456, 'MEL5 → α24β17', fontsize=6.8, color=GRAY)
    fig.text(.405, .924, 'MEL5', ha='center', fontsize=9, color=WT, fontweight='bold')
    fig.text(.79, .924, 'α24β17', ha='center', fontsize=9, color=ENG, fontweight='bold')
    fig.text(.405, .901, '3HG1-derived structure', ha='center', fontsize=6.6, color=GRAY)
    fig.text(.79, .901, '4JFF-derived structure', ha='center', fontsize=6.6, color=GRAY)
    w = local_view(fig, 'wt', [.223, .447, .361, .425])
    e = local_view(fig, 'eng', [.608, .447, .361, .425])
    fig.text(.598, .665, '→', ha='center', fontsize=16, color=GRAY)
    # Contextual labels summarize the common residue set and paired local
    # substitutions. They do not imply that these are the only mutations.
    fig.text(.595, .414, 'Local TCR substitutions: βGly52→Pro52; βIle53→Phe53',
             ha='center', fontsize=7.3)
    fig.text(.595, .387, 'Q72 and βGly54 are unchanged residues on opposite sides of the interface',
             ha='center', fontsize=6.6, color=GRAY)
    fig.text(.02, .335, 'b', fontsize=11, fontweight='bold')
    fig.text(.062, .335, 'Predicted backbone response across all six recognition loops', fontsize=8.5)
    ax = fig.add_axes([.095, .092, .866, .188])
    groups = ['CDR1α', 'CDR2α', 'CDR3α', 'CDR1β', 'CDR2β', 'CDR3β']
    ticks, ticklabels = [], []
    for region in groups:
        q = cdr[cdr.region == region].sort_values('plot_x'); x=q.plot_x.to_numpy()
        if region == 'CDR2β': ax.axvspan(x[0]-.6, x[-1]+.6, color=ENG, alpha=.055, lw=0)
        for field, color, marker, fill in [
            ('MEL5_model_RMSF_A', WT, 'o', 'white'),
            ('alpha24beta17_model_RMSF_A', ENG, 's', ENG),
        ]:
            ax.plot(x, q[field], color=color, marker=marker, ms=2.7, mfc=fill,
                    mec=color, lw=.95, zorder=3)
        ax.text(x.mean(), 1.035, region, transform=ax.get_xaxis_transform(), ha='center', fontsize=7.5)
        ax.plot([x[0]-.2, x[-1]+.2], [1.01, 1.01], transform=ax.get_xaxis_transform(), color='#B5B5B5', lw=.6, clip_on=False)
        choose = [0, len(q)-1] if region != 'CDR3β' else [0, 6, len(q)-1]
        for j in choose: ticks.append(x[j]); ticklabels.append(str(int(q.iloc[j].common_TCR_position)))
        if region == 'CDR2β':
            g = q[q.common_TCR_position == 54].iloc[0]
            ax.annotate('Gly54', xy=(g.plot_x, g.alpha24beta17_model_RMSF_A), xytext=(g.plot_x+2.3, 1.0),
                        fontsize=7.3, color=ENG, arrowprops=dict(arrowstyle='-', lw=.6, color=ENG))
    ax.set(xlim=(-1, cdr.plot_x.max()+1.5), ylim=(0, 1.42), ylabel='RMSF (Å)', xticks=ticks, yticks=[0, .4, .8, 1.2])
    ax.set_xticklabels(ticklabels); ax.tick_params(length=2.5, pad=3)
    ax.set_xlabel('Common TCR residue numbering', labelpad=3, color=GRAY, fontsize=7.5)
    handles = [Line2D([], [], color=WT, marker='o', mfc='white', ms=3.5, lw=1, label='MEL5'),
               Line2D([], [], color=ENG, marker='s', mfc=ENG, ms=3.3, lw=1, label='α24β17')]
    fig.legend(handles=handles, loc='upper right', bbox_to_anchor=(.969, .356), ncol=2,
               fontsize=7.6, columnspacing=1.2, handlelength=1.4, handletextpad=.4, borderaxespad=0)
    fig.text(.095, .017, 'Same view and scale in a; dashed lines measure atom-pair distances. Hχ and RMSF are FlexPPI predictions.',
             fontsize=6.2, color=GRAY)
    return fig

fig = compose_figure()
fig.canvas.draw(); renderer = fig.canvas.get_renderer(); texts=[]
for item in fig.findobj(match=Text):
    if not item.get_visible() or not item.get_text(): continue
    bb = item.get_window_extent(renderer).transformed(fig.transFigure.inverted())
    texts.append(dict(text=item.get_text(), font_pt=item.get_fontsize(), bounds=list(bb.extents),
                      inside_canvas=bool(bb.x0>=-.001 and bb.y0>=-.001 and bb.x1<=1.001 and bb.y1<=1.001)))
assert all(t['inside_canvas'] for t in texts), [t for t in texts if not t['inside_canvas']]
assert min(t['font_pt'] for t in texts) >= 6
assert not any('\u4e00' <= char <= '\u9fff' for t in texts for char in t['text'])
for ext in ('png', 'pdf', 'svg'):
    with mpl.rc_context({'savefig.bbox': None}):
        fig.savefig(ROOT / f'fig4_mel5_interface.{ext}', dpi=450, bbox_inches=None)
plt.close(fig)
cdr.to_csv(ROOT / 'source_CDR_RMSF_all45.csv', index=False, encoding='utf-8-sig')
interface.to_csv(ROOT / 'source_interface_Hchi_all23.csv', index=False, encoding='utf-8-sig')
geometry.to_csv(ROOT / 'source_atom_pair_distances.csv', index=False, encoding='utf-8-sig')
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
audit = dict(status='generated_pending_visual_QA', dimensions_mm=[183, 140], raster_dpi=450,
             all_CDR_positions=45, interface_source_positions=23, shared_crop_pixels=crop,
             camera_metadata=str(META), camera_metadata_sha256=sha(META),
             coordinate_scope='Existing 3HG1/4JFF-derived simulation-prepared starting structures',
             minimum_font_pt=min(t['font_pt'] for t in texts), visible_text=texts,
             removed_panels=['v50 c all23 ranking and TCR-removal context', 'v49 d independent MD occupancy'],
             new_training=0, new_inference=0, new_MD=0,
             no_Hchi_to_RMSF_causal_arrow=True,
             model_anchors=key['model'],
             model_anchor_source=str(OLD / 'case_sources/residue_story/trajectory_link/key_findings.json'),
             source_model_SHA256='ff2974406d5d1ee07e6dc7118e43e57c0acc9849cc8feba7f5290476fca15ec5')
(ROOT / 'figure_data_audit.json').write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding='utf-8')
print(json.dumps({k:audit[k] for k in ['status','dimensions_mm','minimum_font_pt','shared_crop_pixels']}, ensure_ascii=False))

