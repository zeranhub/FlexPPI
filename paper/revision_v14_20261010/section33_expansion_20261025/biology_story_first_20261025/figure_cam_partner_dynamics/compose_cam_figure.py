# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a) molecular overview -> cross-type inherit -> param inherit; real PyMOL PNGs
# (b) molecular close-ups -> cross-type inherit -> param inherit; real PyMOL PNGs
# (c,d) measured-order heatmaps -> cross-type inherit -> param inherit
#       Available heatmap/plot_comparison.py is a method-performance chart.
# (e,f) paired measurements -> LineTrend clean direct-label styling -> param inherit
#       Available LineTrend/plot_trend.py is a publication timeline, not S2 data.
# (g) signed per-probe changes -> cross-type inherit -> param inherit
# RULE: no panel is marked native run; no template example data are used.
"""Figure-3-style experimental CaM story plate. All common probes retained.

Visual adaptation: inherited clean lines, direct labels and restrained color
from academic-figure-skill LineTrend assets. Heatmaps show measured S2, not the
method-performance data in the available heatmap template; no semantic reuse.
"""
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
import csv
import json
import hashlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.colors import LinearSegmentedColormap, Normalize
from PIL import Image

ROOT = Path(__file__).resolve().parent
W, H = 183., 231.
BLUE, WARM = '#2166AC', '#B87936'
PARTNERS = ['CaMKKalpha', 'eNOS', 'nNOS', 'smMLCK']
COLORS = ['#762A83', '#328C8D', BLUE, WARM]
LABELS = ['CaMKKα', 'eNOS', 'nNOS', 'smMLCK']
PDBS = ['1CKK', '1NIW', '2O60', '1CDL']
KEYS = [('L18', 'CD1'), ('L39', 'CD2'), ('L105', 'CD1'), ('M109', 'CE'), ('M144', 'CE')]
READ = lambda name: pd.read_csv(ROOT/'source_data'/name)
TEXTS = []
AXES = []


def txt(fig, x, y, text, size=8, **kw):
    t = fig.text(x/W, 1-y/H, text, fontsize=size, color=kw.pop('color', BLACK), va=kw.pop('va', 'top'), **kw)
    TEXTS.append(t)
    return t


def ax_mm(fig, x, y, w, h):
    a = fig.add_axes([x/W, 1-(y+h)/H, w/W, h/H])
    AXES.append(a)
    return a


def panel(fig, letter, title, x, y):
    txt(fig, x, y, letter, 9, fontweight='bold')
    txt(fig, x+5, y, title, 8)


def image_panel(fig, info, x, y, w, h):
    a = ax_mm(fig, x, y, w, h)
    im = Image.open(ROOT/'structure_renders'/info['file'])
    if info.get('crop'):
        im = im.crop(info['crop'])
    a.imshow(im)
    a.axis('off')
    return a


def callout(a, info, key, xytext, color=BLACK, fs=6.5, mark=False):
    point = info['landmarks'][key]
    if info.get('crop'):
        point = [point[0]-info['crop'][0], point[1]-info['crop'][1]]
    a.annotate(key, xy=point, xytext=xytext, textcoords='axes fraction',
               ha='center', va='center', fontsize=fs, color=color,
               arrowprops=dict(arrowstyle='-', color='#757B80', lw=.5,
                               shrinkA=1.5, shrinkB=1.5, connectionstyle='arc3,rad=0'),
               bbox=dict(fc='white', ec='none', alpha=.88, pad=.4))
    if mark:
        a.plot(point[0], point[1], 'o', ms=2, mfc=ACCENT_RED, mec='white', mew=.3)


def structures(fig, manifest):
    panel(fig, 'a', 'Four partner-bound CaM structures', 0, 14)
    for i, (partner, lab, pid, col) in enumerate(zip(PARTNERS, LABELS, PDBS, COLORS)):
        x = i*46
        txt(fig, x+22, 20, lab, 8, color=col, ha='center', fontweight='bold')
        info = manifest['assets']['overview_'+partner]
        a = image_panel(fig, info, x+.3, 23.5, 44, 36)
        if partner in ['nNOS', 'smMLCK']:
            for key, loc in [('L18', (.13, .46)), ('L39', (.30, .88)),
                             ('L105', (.90, .64)), ('M109', (.91, .42)), ('M144', (.65, .10))]:
                callout(a, info, key, loc, fs=6.1)
        txt(fig, x+22, 59.8, pid+(' · SeMet' if pid == '1NIW' else ''), 6.7, ha='center', color='#676C70')
    txt(fig, 3, 64, 'CaM: grey    Peptide: partner colour    Highlighted side chains: red', 6.6, color='#61676B')
    panel(fig, 'b', 'Key-site environments in the matched nNOS / smMLCK comparison', 0, 70)
    for j, site in enumerate(['L39', 'L105', 'M144']):
        x = j*61
        txt(fig, x+30.5, 75.7, site, 8, ha='center', fontweight='bold')
        for k, (partner, col) in enumerate([('nNOS', BLUE), ('smMLCK', WARM)]):
            info = manifest['assets']['local_'+site+'_'+partner]
            a = image_panel(fig, info, x+1+k*29.5, 79, 28.5, 30)
            # One leader to the experimentally monitored methyl axis position.
            p = info['landmarks']['methyl']
            a.annotate(r'CH$_3$', xy=p, xytext=(.75, .93), textcoords='axes fraction',
                       ha='center', fontsize=6, color=ACCENT_RED,
                       arrowprops=dict(arrowstyle='-', lw=.45, color='#73797D'),
                       bbox=dict(fc='white', ec='none', alpha=.87, pad=.3))
            txt(fig, x+15+k*29.5, 109.4, partner, 6.8, color=col, ha='center')
    txt(fig, 182, 114, 'Deposited snapshots locate the probes; they do not depict motion.', 6.3, ha='right', color='#686D72')


def heatmaps(fig, backbone, methyl):
    panel(fig, 'c', 'Backbone N–H: all 87 common probes', 0, 120)
    panel(fig, 'd', 'Side-chain methyl: all 44 common probes', 93, 120)
    cmap = LinearSegmentedColormap.from_list('constraint', ['#FFFFFF', '#B7CCD8', '#294E68'])
    for start, df, count, ylabel in [(0, backbone, 87, 'N–H'), (93, methyl, 44, 'methyl')]:
        a = ax_mm(fig, start+16, 127, 68 if count == 87 else 64.5, 17)
        values = df[[p+'_s2' for p in PARTNERS]].to_numpy().T
        assert values.shape == (4, count)
        a.imshow(values, cmap=cmap, vmin=0, vmax=1, interpolation='nearest', aspect='auto')
        a.set_yticks(range(4), LABELS, fontsize=6.2)
        for label, color in zip(a.get_yticklabels(), COLORS):
            label.set_color(color)
        indexes = []
        for residue in [18, 39, 105, 109, 144]:
            sub = np.flatnonzero(df.seq_id.to_numpy() == residue)
            if len(sub):
                if not (count == 87 and residue == 109):
                    indexes.append((int(sub[0]), str(residue)))
                else:
                    a.annotate('109', xy=(int(sub[0]), -.5), xytext=(int(sub[0])+3, -1.05),
                               xycoords='data', textcoords='data', fontsize=5.8, ha='center',
                               va='bottom', annotation_clip=False,
                               arrowprops=dict(arrowstyle='-', lw=.4, color=WARM))
                for idx in sub:
                    a.add_patch(Rectangle((idx-.5, -.5), 1, 4, fill=False, ec=WARM, lw=.55))
        a.set_xticks([v[0] for v in indexes], [v[1] for v in indexes], fontsize=6)
        plt.setp(a.get_xticklabels(), rotation=45, ha='right', rotation_mode='anchor')
        a.tick_params(length=0, pad=2)
        for spine in a.spines.values():
            spine.set_visible(False)
        a.set_xlabel('CaM residue number; probes in sequence order', fontsize=6, labelpad=1)
        np.savetxt(ROOT/'source_data'/f'plotted_{ylabel.replace("–", "_")}_matrix.csv', values, delimiter=',')
    cax = ax_mm(fig, 175.5, 127, 2, 17)
    mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=Normalize(0,1), ticks=[0, .5, 1])
    cax.tick_params(labelsize=5.5, length=1, pad=1)
    cax.set_ylabel('S²: greater restriction →', fontsize=5.8, labelpad=1)


def local_data(fig, key):
    panel(fig, 'e', 'Same sites: backbone constraints', 0, 157)
    panel(fig, 'f', 'Same sites: partner-specific side-chain motion', 93, 157)
    axes = [ax_mm(fig, 15, 164, 69, 23), ax_mm(fig, 108, 164, 69, 23)]
    for a, metric, error in zip(axes, ['NH_s2', 'methyl_s2'], ['NH_fit_error', 'methyl_fit_error']):
        for i, (residue, atom) in enumerate(KEYS):
            site = key[(key.methyl_probe_id == residue+'-'+atom)]
            vv = []
            for k, (partner, col, marker) in enumerate([('nNOS', BLUE, 'o'), ('smMLCK', WARM, 's')]):
                row = site[site.partner == partner].iloc[0]
                vv.append(row[metric])
                a.errorbar(i+(k-.5)*.16, row[metric], yerr=row[error], fmt=marker,
                           color=col, ms=3.1, elinewidth=.65, capsize=1.5, capthick=.6, zorder=3)
                if metric == 'methyl_s2' and residue != 'M109':
                    greater = (partner == 'smMLCK') if residue in ['L18', 'L105'] else (partner == 'nNOS')
                    a.annotate(f'{row[metric]:.3f}', (i+(k-.5)*.16, row[metric]),
                               xytext=(0, 6 if greater else -7), textcoords='offset points',
                               ha='center', va='bottom' if greater else 'top', fontsize=5.8, color=col)
            a.plot([i-.08, i+.08], vv, color='#BCC1C4', lw=.65, zorder=2)
        a.set_ylim(0, 1.05);a.set_xlim(-.5, 4.5)
        a.set_yticks([0, .5, 1]);a.set_ylabel('N–H S²' if metric == 'NH_s2' else 'Methyl-axis S²', labelpad=1)
        a.set_xticks(range(5), [r for r, _ in KEYS], fontsize=6.7)
        a.tick_params(pad=2, length=2)
    handles = [Line2D([], [], marker='o', ls='', color=BLUE, ms=3, label='nNOS'),
               Line2D([], [], marker='s', ls='', color=WARM, ms=3, label='smMLCK')]
    leg = fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(.51, 1-152.8/H),
                     ncol=2, fontsize=6.5, handletextpad=.3, columnspacing=1., borderaxespad=0)
    TEXTS.extend(leg.get_texts())
    txt(fig, 182, 191.3, 'Points ± deposited fit error; no confidence intervals or significance tests.',
        6.1, ha='right', color='#676C70')


def all_changes(fig, delta):
    panel(fig, 'g', 'All 44 methyl probes: opposing local changes despite similar means', 0, 197)
    a = ax_mm(fig, 15, 204, 162, 15)
    values = delta.delta_smMLCK_minus_nNOS_s2.to_numpy()
    a.axhline(0, color='#A6ABB0', lw=.55)
    a.bar(np.arange(len(delta)), values, width=.72, color=[WARM if v>0 else BLUE for v in values], zorder=2)
    a.set_xlim(-1, 44);a.set_ylim(-.53, .53);a.set_yticks([-.4, 0, .4]);a.set_ylabel('ΔS²', labelpad=1)
    key_ticks = [int(np.flatnonzero(delta.probe_id.to_numpy() == res+'-'+atom)[0]) for res,atom in KEYS]
    ticks = sorted(set([0, len(delta)-1]+key_ticks))
    a.set_xticks(ticks, delta.iloc[ticks].probe_id.str.replace('-', ':', regex=False), fontsize=5.8, rotation=0)
    a.tick_params(length=2, pad=1.5)
    for res, atom in KEYS:
        idx = np.flatnonzero(delta.probe_id.to_numpy() == res+'-'+atom)[0]
        val = values[idx]
        a.plot(idx, val+(.035 if val>0 else -.035), marker='v' if val>0 else '^', ms=2.3, color=BLACK)
    txt(fig, 182, 224.5, 'smMLCK − nNOS: 21 increase / 23 decrease; mean S² 0.546 → 0.549',
        6.4, ha='right', color='#555D62')


def audit_layout(fig):
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    bbox = fig.bbox
    out = []
    for text in fig.findobj(mpl.text.Text):
        if not text.get_visible() or not text.get_text():
            continue
        box = text.get_window_extent(renderer)
        if box.width and box.height and (box.x0 < -2 or box.y0 < -2 or box.x1 > bbox.x1+2 or box.y1 > bbox.y1+2):
            out.append(dict(text=text.get_text(), bounds=list(box.bounds)))
    return out


def audit_axes_overlap():
    overlap = []
    for i,a in enumerate(AXES):
        for j,b in enumerate(AXES[i+1:], start=i+1):
            aa, bb = a.get_position(), b.get_position()
            xx, yy = min(aa.x1,bb.x1)-max(aa.x0,bb.x0), min(aa.y1,bb.y1)-max(aa.y0,bb.y0)
            if xx > .0001 and yy > .0001:
                overlap.append([i,j,xx*W,yy*H])
    return overlap


def main():
    manifest = json.loads((ROOT/'structure_renders/render_manifest.json').read_text())
    # Crop all four overviews by the same bounds, preserving their common scale.
    boxes = []
    for partner in PARTNERS:
        info = manifest['assets']['overview_'+partner]
        im = Image.open(ROOT/'structure_renders'/info['file']).convert('RGBA')
        mask = np.asarray(im)[:, :, 3] > 20
        yy, xx = np.where(mask)
        boxes.append((xx.min(), yy.min(), xx.max()+1, yy.max()+1))
    common_crop = [max(0, min(b[0] for b in boxes)-70), max(0,min(b[1] for b in boxes)-70),
                   min(1400,max(b[2] for b in boxes)+70), min(1300,max(b[3] for b in boxes)+70)]
    for partner in PARTNERS:
        manifest['assets']['overview_'+partner]['crop'] = [int(x) for x in common_crop]
    backbone, methyl = READ('backbone_87_NH_wide.csv'), READ('sidechain_44_methyl_wide.csv')
    key, delta = READ('key_sites_matched_NH_methyl_long.csv'), READ('sidechain_44_methyl_delta_smMLCK_minus_nNOS.csv')
    assert len(backbone) == 87 and len(methyl) == 44 and len(delta) == 44
    fig = plt.figure(figsize=(W/25.4, H/25.4), dpi=180, facecolor='white')
    txt(fig, 0, 1, 'Calmodulin | Partner-specific side-chain motion within a constrained backbone', 9.1, fontweight='bold')
    txt(fig, 0, 7, 'Experimental NMR S² and deposited structural contexts', 7.1, color='#62686D')
    structures(fig, manifest)
    heatmaps(fig, backbone, methyl)
    local_data(fig, key)
    all_changes(fig, delta)
    # Reserve a small print margin while keeping all fonts at their point sizes.
    for a in AXES:
        p = a.get_position()
        a.set_position([3/W+p.x0*(1-6/W), 2/H+p.y0*(1-4/H),
                        p.width*(1-6/W), p.height*(1-4/H)])
    for t in fig.texts:
        x, y = t.get_position()
        t.set_position((3/W+x*(1-6/W), 2/H+y*(1-4/H)))
    outside = audit_layout(fig)
    base = ROOT/'fig_cam_partner_specific_dynamics'
    # Exact physical dimensions; standard export block above is intentionally
    # retained, with no tight crop here because panels use absolute-mm layout.
    for extension, dpi in [('pdf', 600), ('svg', 600), ('png', 400)]:
        fig.savefig(str(base)+'.'+extension, dpi=dpi, bbox_inches=fig.bbox_inches, facecolor='white')
    report = dict(status='rendered', width_mm=W, height_mm=H,
                  all_common_backbone_probes=87, all_common_methyl_probes=44,
                  selected_residues=[r for r,a in KEYS],
                  panel_order='a four structures; b three local structure pairs; c 87 NH; d 44 methyl; e key NH; f key methyl; g all methyl deltas',
                  out_of_canvas_text=outside,
                  overlapping_axes=audit_axes_overlap(),
                  numeric_audit='numeric_audit.json; 524 raw values and associated fit errors independently matched',
                  templates='Cross-type visual inheritance; production performance heatmap and publication timeline are not semantically compatible. No template example numbers are used.',
                  limitations=['Structures are snapshots, not motion ensembles', 'NMR S2 is not RMSF or Hchi',
                               'Key-site display is explanatory; all common probes are shown in c,d,g',
                               'No current model predictions included', 'Error bars are deposited fit errors'])
    (ROOT/'figure_layout_audit.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    print(json.dumps(report, indent=2), flush=True)


if __name__ == '__main__':
    main()
