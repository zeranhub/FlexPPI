"""Layout-neutral, unchanged native v110 figure components.

All boxes use millimetres (left, top, width, height), origin at page top left.
No full-page layout is specified or rendered here. Use only after an actual
PaperBanana result is available; do not label an invented layout PaperBanana.
"""
from pathlib import Path
import json
import sys
import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.legend_handler import HandlerBase
from PIL import Image

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
from case_profiles import draw_case, draw_sasa, method_handles
from plot_overlay import draw_overlay as draw_density

P = ROOT / 'figure3/source_data'
R = ROOT / 'pymol'
from palette import BLUE,ORANGE
WARM = ORANGE
ACCENT = WARM
PURPLE = '#8B72A5'
GRAY = '#BAC1C7'
INK = '#262626'
LIGHT = '#CCD1D5'


def apply_style():
    """Match manuscript Fig2 native 9/10/11 pt Arial typography."""
    mpl.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Helvetica', 'Liberation Sans'],
        'font.size': 9, 'axes.labelsize': 9, 'axes.titlesize': 10, 'xtick.labelsize': 9,
        'ytick.labelsize': 9, 'legend.fontsize': 9, 'axes.linewidth': .6,
        'axes.spines.top': False, 'axes.spines.right': False, 'axes.unicode_minus': False,
        'xtick.major.width': .6, 'ytick.major.width': .6, 'xtick.direction': 'out', 'ytick.direction': 'out',
        'pdf.fonttype': 42, 'svg.fonttype': 'none', 'savefig.bbox': None, 'savefig.dpi': 450})


def text(fig, canvas_mm, x, y, value, size=9, **kwargs):
    width, height = canvas_mm
    return fig.text(x/width, 1-y/height, value, fontsize=size,
                    va=kwargs.pop('va', 'top'), color=kwargs.pop('color', INK), **kwargs)


def axis(fig, canvas_mm, box):
    width, height = canvas_mm
    left, top, w, h = box
    return fig.add_axes([left/width, 1-(top+h)/height, w/width, h/height])


def panel(fig, canvas_mm, letter, title, x, y):
    text(fig, canvas_mm, x, y, letter, 11, fontweight='bold')
    if title:
        text(fig, canvas_mm, x+6, y+.25, title, 10)


def style_axis(a):
    a.tick_params(length=2.5, pad=2, width=.6)
    a.xaxis.labelpad = 2.5
    a.yaxis.labelpad = 2.5


def molecule(fig, canvas_mm, key, box, labels=None):
    """Use genuine frozen PyMOL image, crop and landmark annotations only."""
    manifest = json.loads((R/'render_manifest.json').read_text(encoding='utf-8'))
    z = manifest['assets'][key]
    a = axis(fig, canvas_mm, box)
    a.imshow(Image.open(R/z['png']).convert('RGBA'))
    a.axis('off')
    crop = z.get('crop',{'cam_overview': (200,110,1430,1315), 'ub_overview': (90,205,1535,1195)}[key])
    a.set_xlim(crop[0], crop[2]); a.set_ylim(crop[3], crop[1]); a.set_aspect('equal')
    if labels is None:
        labels = {'cam_overview': {'W303': (.26,.40), 'M316': (.97,.41)},
                  'ub_overview': {'Q41': (.22,.22), 'Q62': (.98,.44)}}[key]
    for name, location in labels.items():
        position = z['landmarks'][name]
        a.annotate(name, xy=position, xytext=location, textcoords='axes fraction', color=ACCENT,
                   fontsize=9, ha='center', va='center', arrowprops=(dict(arrowstyle='-', color=ACCENT, lw=.7) if name=='W303' else None),
                   bbox=dict(facecolor='white', edgecolor='none', pad=.1, alpha=.92))
    return a


def molecule_pair(fig, canvas_mm, cam_box, ubiquitin_box, cam_labels=None, ubiquitin_labels=None):
    return (molecule(fig, canvas_mm, 'cam_overview', cam_box, cam_labels),
            molecule(fig, canvas_mm, 'ub_overview', ubiquitin_box, ubiquitin_labels))


def component_legend(fig, canvas_mm, entries, y):
    """Entries are (left_mm, label, colour); no page coordinates are assumed."""
    width, height = canvas_mm
    for left, label, colour in entries:
        fig.add_artist(Line2D([left/width, (left+4)/width], [1-y/height]*2,
                              transform=fig.transFigure, color=colour, lw=2.2))
        text(fig, canvas_mm, left+5.5, y-1.6, label)


class PairedMetricLines(HandlerBase):
    """Display both metric colours as two parallel lines for each method."""
    def create_artists(self, legend, orig_handle, xdescent, ydescent, width, height, fontsize, trans):
        artists=[]
        for i, original in enumerate(orig_handle):
            y=height*(.25+.50*i)-ydescent
            line=Line2D([xdescent,xdescent+width],[y,y])
            self.update_prop(line,original,legend)
            line.set_transform(trans)
            artists.append(line)
        return artists


def method_legend(fig, canvas_mm, center_x, center_y):
    width, height = canvas_mm
    handles=[tuple(Line2D([],[],color=colour,lw=lw,ls=style,alpha=alpha) for colour in [BLUE,WARM])
             for lw,style,alpha in [(1.2,'-',1.),(.95,(0,(3,2)),1.)]]
    return fig.legend(handles=handles,labels=['Full','MD'],handler_map={tuple:PairedMetricLines()},
                       loc='center',bbox_to_anchor=(center_x/width,1-center_y/height),
                       ncol=2,handlelength=1.9,handleheight=1.25,columnspacing=1.2,frameon=False)


def annotate_case_sites(a, case):
    """The same v110 residue-axis annotations; do not imply a measured y-value."""
    if case == 'CaMKI':
        for position, name in [(303,'W303'), (316,'M316')]:
            a.axvline(position, color=BLUE, lw=.6, ls=':', alpha=1., zorder=0)
            a.text(position, .98, name, transform=a.get_xaxis_transform(), ha='right' if name=='M316' else 'center', va='top', fontsize=9, fontweight='bold', color=ACCENT)
    else:
        for position in [41,62]:
            a.axvline(position, color=BLUE, lw=.6, ls=':', alpha=1., zorder=0)
            a.text(position, .99, 'Q'+str(position), transform=a.get_xaxis_transform(), ha='center', va='top', fontsize=9, fontweight='bold', color=ACCENT)


def chi_legend(fig, canvas_mm, center_x, center_y):
    width, height = canvas_mm
    handles = [Line2D([], [], ls='', marker=m, ms=3.6, color=c, label=name)
               for name,m,c in [('Other','o',BLUE), ('Q41','D',ACCENT), ('Q62','s',ACCENT)]]
    return fig.legend(handles=handles, loc='center', bbox_to_anchor=(center_x/width,1-center_y/height),
                       ncol=3, handletextpad=.25, columnspacing=.7, borderpad=0, frameon=False)


def draw_4chi_row(fig, canvas_mm, lefts, top, width, height, title_y):
    """Keep four angle plots horizontal. Return axes and all 302 native values."""
    assert len(lefts) == 4
    defined = pd.read_csv(P/'defined_angles_long.csv', float_precision='round_trip')
    stats = pd.read_csv(P/'per_chi_correlations.csv', float_precision='round_trip')
    stats = stats[stats.predictor.eq('predicted_joint_Hchi') & stats.metric.eq('fine_entropy')]
    rows = []; axes = []
    for chi, left in enumerate(lefts, start=1):
        a = axis(fig, canvas_mm, (left,top,width,height)); style_axis(a); axes.append(a)
        z = defined[defined.chi.eq('chi'+str(chi))]
        stat = stats[stats.chi.eq('chi'+str(chi))].iloc[0]
        assert len(z) == {1:66,2:52,3:23,4:10}[chi]
        text(fig, canvas_mm, left, title_y, 'χ'+str(chi)+'  n='+str(len(z))+'  ρ='+format(stat.rho_AA_centered,'.2f'))
        for name, marker, colour in [('Other','o',BLUE),('Q41','D',ACCENT),('Q62','s',ACCENT)]:
            sub = z[~z.author_residue_id.isin([41,62])] if name=='Other' else z[z.author_residue_id.eq(int(name[1:]))]
            a.scatter(sub.predicted_Hchi_nats_AA_centered, sub.fine_entropy_four_origins_mean_nats_AA_centered,
                      s=8.5 if name=='Other' else 19, marker=marker, color=colour,
                      alpha=1., linewidths=.25,
                      edgecolors='white' if name=='Other' else colour, zorder=3 if name=='Other' else 4)
        for row in z.itertuples():
            rows.extend([dict(panel='e',kind='predicted_joint_Hchi_AA_centered',identity=str(int(row.author_residue_id)),chi=chi,series='Full',value=float(row.predicted_Hchi_nats_AA_centered)),
                         dict(panel='e',kind='MD_single_chi_fine_entropy_AA_centered',identity=str(int(row.author_residue_id)),chi=chi,series='MD',value=float(row.fine_entropy_four_origins_mean_nats_AA_centered))])
        a.axhline(0,color=LIGHT,lw=.4,ls=':'); a.axvline(0,color=LIGHT,lw=.4,ls=':')
        a.set(xlim=(-.7,.5),ylim=(-1,.75)); a.set_xticks([-.5,0,.5]); a.set_yticks([-.5,0,.5])
        if chi != 1: a.tick_params(labelleft=False)
    assert len(rows) == 302
    return axes, pd.DataFrame(rows)


def chi_shared_labels(fig, canvas_mm, ylabel_x, ylabel_y, xlabel_x, xlabel_y):
    return (text(fig, canvas_mm, ylabel_x, ylabel_y, 'MD entropy deviation\n(nats)', rotation=90, ha='center', va='center',linespacing=1.0),
            text(fig, canvas_mm, xlabel_x, xlabel_y, 'Predicted joint Hχ deviation (nats)', ha='center'))


apply_style()
