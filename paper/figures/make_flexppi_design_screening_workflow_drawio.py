from __future__ import annotations

import base64
from pathlib import Path
import xml.etree.ElementTree as ET
from xml.dom import minidom


HERE = Path(__file__).resolve().parent
OUT = HERE / "fig_flexppi_design_screening_workflow.drawio"
ASSETS = HERE / "fig_flexppi_design_screening_workflow_assets"
DESIGN_CANDIDATES = ASSETS / "design_candidates_clean_2_target_tight.png"
CANDIDATE_STRUCTURE = ASSETS / "candidate_structure_2_target_tight.png"


def image_uri(path: Path) -> str:
    payload = base64.b64encode(path.read_bytes()).decode("ascii")
    # Draw.io cell styles use semicolons as separators, so the semicolon in
    # the data-URI header must be percent-encoded.
    return f"data:image/png%3Bbase64,{payload}"


mxfile = ET.Element(
    "mxfile",
    {
        "host": "Electron",
        "modified": "2026-08-30T00:00:00.000Z",
        "agent": "Codex",
        "version": "24.7.17",
    },
)
diagram = ET.SubElement(
    mxfile,
    "diagram",
    {"id": "flexppi-design-screening", "name": "FlexPPI design screening workflow"},
)
model = ET.SubElement(
    diagram,
    "mxGraphModel",
    {
        "dx": "828",
        "dy": "425",
        "grid": "1",
        "gridSize": "10",
        "guides": "1",
        "tooltips": "1",
        "connect": "1",
        "arrows": "1",
        "fold": "1",
        "page": "1",
        "pageScale": "1",
        "pageWidth": "828",
        "pageHeight": "425",
        "math": "0",
        "shadow": "0",
        "background": "#FFFFFF",
    },
)
root = ET.SubElement(model, "root")
ET.SubElement(root, "mxCell", {"id": "0"})
ET.SubElement(root, "mxCell", {"id": "1", "parent": "0"})


def vertex(cell_id: str, value: str, style: str, x: float, y: float, w: float, h: float):
    cell = ET.SubElement(
        root,
        "mxCell",
        {"id": cell_id, "value": value, "style": style, "vertex": "1", "parent": "1"},
    )
    ET.SubElement(
        cell,
        "mxGeometry",
        {"x": str(x), "y": str(y), "width": str(w), "height": str(h), "as": "geometry"},
    )
    return cell


def edge(
    cell_id: str,
    points: list[tuple[float, float]],
    *,
    color: str = "#2A2A2A",
    width: float = 1.35,
    arrow: bool = False,
    curved: bool = False,
):
    style_parts = [
        "edgeStyle=none",
        "html=1",
        f"strokeColor={color}",
        f"strokeWidth={width}",
        "rounded=0",
    ]
    if curved:
        style_parts.append("curved=1")
    if arrow:
        style_parts.extend(["endArrow=classic", "endFill=1", "endSize=8"])
    else:
        style_parts.append("endArrow=none")
    cell = ET.SubElement(
        root,
        "mxCell",
        {
            "id": cell_id,
            "value": "",
            "style": ";".join(style_parts) + ";",
            "edge": "1",
            "parent": "1",
        },
    )
    geom = ET.SubElement(cell, "mxGeometry", {"relative": "1", "as": "geometry"})
    ET.SubElement(geom, "mxPoint", {"x": str(points[0][0]), "y": str(points[0][1]), "as": "sourcePoint"})
    ET.SubElement(geom, "mxPoint", {"x": str(points[-1][0]), "y": str(points[-1][1]), "as": "targetPoint"})
    if len(points) > 2:
        arr = ET.SubElement(geom, "Array", {"as": "points"})
        for px, py in points[1:-1]:
            ET.SubElement(arr, "mxPoint", {"x": str(px), "y": str(py)})
    return cell


TEXT_BASE = "text;html=1;whiteSpace=wrap;strokeColor=none;fillColor=none;fontFamily=Arial;fontColor=#272727;"


def text_cell(
    cell_id: str,
    value: str,
    x: float,
    y: float,
    w: float,
    h: float,
    *,
    size: float,
    align: str = "center",
    valign: str = "middle",
    bold: bool = False,
    italic: bool = False,
    rotation: int | None = None,
    color: str = "#272727",
):
    font_style = (1 if bold else 0) + (2 if italic else 0)
    style = (
        TEXT_BASE
        + f"align={align};verticalAlign={valign};fontSize={size};fontStyle={font_style};fontColor={color};"
        + (f"rotation={rotation};" if rotation is not None else "")
    )
    return vertex(cell_id, value, style, x, y, w, h)


vertex(
    "pageBackground",
    "",
    "rounded=0;whiteSpace=wrap;html=1;fillColor=#FFFFFF;strokeColor=none;",
    0,
    0,
    828,
    425,
)

# Upper-left target-profile chart.
edge("target-y-axis", [(37, 29), (37, 160)], color="#A8ADB2", width=1.0)
edge("target-x-axis", [(37, 160), (212, 160)], color="#A8ADB2", width=1.0)
edge(
    "target-profile-curve",
    [
        (46, 35), (50, 48), (53, 71), (57, 101), (63, 128), (69, 144),
        (77, 149), (85, 149), (91, 141), (97, 139), (103, 148), (110, 150),
        (116, 137), (121, 143), (126, 125), (130, 99), (134, 85), (138, 88),
        (143, 115), (148, 138), (154, 144), (160, 136), (166, 120), (171, 113),
        (176, 119), (182, 133), (190, 143), (202, 150),
    ],
    color="#55C995",
    width=1.45,
    curved=True,
)
text_cell("target-y-label", "Flexibility", 0, 77, 48, 26, size=14, rotation=270)
text_cell("target-x-label", "Residue", 78, 158, 93, 25, size=14)

# Step 1 and the protein-design-model stack.
edge("arrow-1", [(216, 100), (260, 100)], arrow=True)
text_cell("label-1", "(1)", 214, 61, 41, 35, size=21)
vertex("design-sheet-back", "", "rounded=0;whiteSpace=wrap;html=1;fillColor=#EAF4EE;strokeColor=none;opacity=88;", 273, 14, 248, 132)
vertex("design-sheet-mid", "", "rounded=0;whiteSpace=wrap;html=1;fillColor=#D9ECDF;strokeColor=none;opacity=94;", 284, 25, 248, 132)
vertex("design-sheet-front", "", "rounded=0;whiteSpace=wrap;html=1;fillColor=#A8DEB8;strokeColor=none;", 295, 37, 247, 132)
text_cell("design-title", "Protein Structure<br>Design Models", 305, 43, 227, 54, size=18.5)
text_cell(
    "design-subtitle",
    "Diverse protein structure<br>generators",
    305,
    99,
    227,
    47,
    size=14.5,
    italic=True,
)

# Step 2 and cropped candidate structures.
edge("arrow-2", [(553, 100), (605, 100)], arrow=True)
text_cell("label-2", "(2)", 557, 61, 41, 35, size=21)
vertex(
    "design-candidates-image",
    "",
    f"shape=image;imageAspect=0;aspect=fixed;html=1;image={image_uri(DESIGN_CANDIDATES)};strokeColor=none;fillColor=none;",
    607,
    8,
    164,
    181,
)

# Step 3 and FlexPPI prediction stack.
edge("arrow-3", [(694, 196), (694, 250)], arrow=True)
text_cell("label-3", "(3)", 652, 203, 41, 35, size=21)
vertex("predict-sheet-back", "", "rounded=0;whiteSpace=wrap;html=1;fillColor=#E7EEF8;strokeColor=none;", 588, 263, 190, 99)
vertex("predict-sheet-mid", "", "rounded=0;whiteSpace=wrap;html=1;fillColor=#D7E2F2;strokeColor=none;", 600, 271, 190, 99)
vertex("predict-sheet-front", "", "rounded=0;whiteSpace=wrap;html=1;fillColor=#B1C5E6;strokeColor=none;", 611, 279, 189, 100)
text_cell("predict-title", "FlexPPI", 620, 284, 171, 32, size=21)
text_cell("predict-subtitle", "Flexibility prediction", 620, 322, 171, 32, size=16.5, italic=True)

# Step 4, candidate structure, and editable flexibility color scale.
edge("arrow-4", [(554, 320), (488, 320)], arrow=True)
text_cell("label-4", "(4)", 515, 284, 41, 35, size=21)
vertex(
    "candidate-structure-image",
    "",
    f"shape=image;imageAspect=0;aspect=fixed;html=1;image={image_uri(CANDIDATE_STRUCTURE)};strokeColor=none;fillColor=none;",
    326,
    252,
    96,
    140,
)
bar_colors = [
    "#E61E2D", "#F13D4E", "#F76A78", "#FA98A2", "#FBC9CE",
    "#F7F7F7", "#DDE2F0", "#B8C5DF", "#899CC8", "#5C72AE", "#3B4C91",
]
for idx, fill in enumerate(bar_colors):
    vertex(
        f"colorbar-{idx+1}",
        "",
        f"rounded=0;whiteSpace=wrap;html=1;fillColor={fill};strokeColor=none;",
        434,
        261 + idx * 11,
        14,
        11,
    )
text_cell("colorbar-max", "max", 451, 254, 35, 20, size=11, align="left")
text_cell("colorbar-flex", "flex", 451, 306, 35, 20, size=11, align="left")
text_cell("colorbar-min", "min", 451, 369, 35, 20, size=11, align="left")
text_cell("candidate-label", "Candidate", 333, 391, 111, 27, size=14)

# Step 5 and lower-left comparison chart.
edge("arrow-5", [(318, 320), (230, 320)], arrow=True)
text_cell("label-5", "(5)", 258, 284, 41, 35, size=21)
edge("compare-y-axis", [(37, 252), (37, 382)], color="#A8ADB2", width=1.0)
edge("compare-x-axis", [(37, 382), (212, 382)], color="#A8ADB2", width=1.0)
edge(
    "compare-target-curve",
    [
        (47, 258), (51, 269), (55, 290), (59, 318), (64, 346), (69, 364),
        (76, 372), (84, 369), (90, 358), (96, 369), (103, 374), (109, 362),
        (115, 367), (121, 351), (126, 325), (131, 305), (135, 304), (140, 331),
        (146, 362), (152, 369), (157, 350), (163, 329), (168, 310), (173, 329),
        (179, 359), (187, 369), (196, 367), (204, 351),
    ],
    color="#55C995",
    width=1.45,
    curved=True,
)
edge(
    "compare-flexppi-curve",
    [
        (47, 277), (51, 301), (55, 329), (60, 309), (64, 344), (69, 371),
        (76, 366), (84, 350), (90, 342), (96, 354), (103, 366), (109, 360),
        (116, 374), (121, 347), (126, 326), (132, 317), (138, 340), (144, 370),
        (150, 373), (157, 352), (164, 340), (170, 359), (178, 376), (185, 369),
        (191, 348), (197, 360), (204, 340),
    ],
    color="#51458B",
    width=1.35,
    curved=True,
)
edge("legend-target-line", [(77, 265), (98, 265)], color="#55C995", width=1.4)
text_cell("legend-target-label", "Target profile", 102, 255, 83, 18, size=9.5, align="left")
edge("legend-flexppi-line", [(77, 279), (98, 279)], color="#51458B", width=1.35)
text_cell("legend-flexppi-label", "FlexPPI prediction", 102, 269, 101, 18, size=9.5, align="left")
text_cell("compare-y-label", "Flexibility", 0, 300, 48, 26, size=14, rotation=270)
text_cell("compare-x-label", "Residue", 78, 380, 93, 25, size=14)


rough = ET.tostring(mxfile, encoding="utf-8")
pretty = minidom.parseString(rough).toprettyxml(indent="  ", encoding="utf-8")
OUT.write_bytes(pretty)
print(OUT)
