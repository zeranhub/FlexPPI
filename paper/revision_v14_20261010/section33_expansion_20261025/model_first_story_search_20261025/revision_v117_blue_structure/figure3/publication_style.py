"""Semantic style shared by the four v59 manuscript result figures.

Blue: FlexPPI / higher magnitude; warm ochre: explicit comparison or focal site.
Other groups use neutral shades and redundant markers/line styles.
AP and percentile rank share an exact 0–1 colour mapping, but retain distinct labels.
Signed differences retain physical units and a symmetric, zero-centred colour map.
"""
from matplotlib.colors import LinearSegmentedColormap, Normalize, TwoSlopeNorm

BLUE = "#2166AC"
WARM = "#B87936"
TEAL = "#3A8D8C"
INK = "#262626"
MD = "#62676B"
MID = "#92999F"
LIGHT = "#CCD1D5"
PALE = "#F2F3F4"
SEQUENTIAL = ["#F7FBFF", "#6BAED6", "#08306B"]
SIGNED = ["#2166AC", "#F7F7F7", "#B87936"]
UNIT_TICKS = [0.0, 0.25, 0.5, 0.75, 1.0]

def unit_cmap():
    return LinearSegmentedColormap.from_list("publication_unit", SEQUENTIAL, N=256)

def unit_norm():
    return Normalize(vmin=0.0, vmax=1.0)

def signed_cmap():
    return LinearSegmentedColormap.from_list("publication_signed", SIGNED, N=256)

def signed_norm(limit):
    return TwoSlopeNorm(vmin=-limit, vcenter=0.0, vmax=limit)
