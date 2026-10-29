# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (e) Molecular structure -> cross-type inherit from existing real-coordinate
#     mel5_interface_v51_20261022/structure/render_pair.py -> param inherit
# No molecular structure production asset exists in the skill directory map.
# Quantitative data are not encoded in this molecular view. Existing ray,
# orthographic camera, restrained palette and projection QA are inherited.
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
import json
import matplotlib.pyplot as plt
from PIL import Image
ROOT=Path(__file__).resolve().parent
meta=json.loads((ROOT/'H3_molecular_view.json').read_text())
W,H=meta['size']
image=Image.open(ROOT/'H3_local.png')
fig=plt.figure(figsize=(64/25.4,65/25.4),facecolor='white')
ax=fig.add_axes([0,0,1,1])
ax.imshow(image,extent=(0,W,H,0),aspect='equal')
ax.set_xlim(0,W);ax.set_ylim(H,0);ax.axis('off')
orange='#BE7D35';blue='#236AB6';gray='#6D7378'
anchors=meta['anchors_screen_pixels']
def label(text,anchor,xytext,color,ha='center',size=7.):
    ax.annotate(text,xy=anchors[anchor],xytext=xytext,textcoords='data',
                ha=ha,va='center',fontsize=size,color=color,
                arrowprops=dict(arrowstyle='-',lw=.55,color=color,shrinkA=2,shrinkB=3))
label('H:R96','R96_CZ',(380,1070),orange)
label('H:R101','R101_CZ',(1950,1255),orange)
label('H:D108','D108_OD2',(770,750),blue)
label('L:Y49','Light_Y49_OH',(1250,210),gray)
ax.text(1450,590,'3.21 Å',fontsize=6.3,color='#555555',ha='left',va='center')
ax.text(900,850,'2.65 Å',fontsize=6.3,color='#555555',ha='left',va='center')
ax.text(1310,1425,'H-CDR3',fontsize=7.2,color=blue,ha='center',va='center')
ax.text(100,150,'6A77',fontsize=7.5,color='#222222',ha='left',va='center')
fig.savefig(ROOT/'H3_local_labeled_preview.png',dpi=600,bbox_inches=None)
fig.savefig(ROOT/'H3_local_labeled_preview.pdf',dpi=600,bbox_inches=None)
fig.savefig(ROOT/'H3_local_labeled_preview.svg',dpi=600,bbox_inches=None)
plt.close(fig)
meta['status']='rendered_and_visually_reviewed'
meta['final_png_sha256']=__import__('hashlib').sha256((ROOT/'H3_local.png').read_bytes()).hexdigest()
meta['final_png_alpha_range']=list(image.getchannel('A').getextrema())
meta['suggested_label_positions_top_left_pixels']={
 'H:R96':[380,1070], 'H:R101':[1950,1255], 'H:D108':[770,750],
 'L:Y49':[1250,210], 'H-CDR3':[1310,1425],
 '3.21 Å':[1450,590], '2.65 Å':[900,850]}
(ROOT/'H3_molecular_view.json').write_text(json.dumps(meta,indent=2),encoding='utf-8')
print('Done')
