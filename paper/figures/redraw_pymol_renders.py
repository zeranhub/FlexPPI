"""PyMOL renderer for the FlexPPI case-study figures.

Replaces the 3Dmol.js backend inside make_case_study_figures.py by
monkey-patching its ``render()`` entry point.  Geometry, residue selections,
contact distances, labels and statistics stay identical; the visual output is
produced with PyMOL ray tracing (cartoon, sticks, dashed distances, labels
with background boxes, per-residue b-factor spectrum coloring).

Requires the conda env at D:\\Anaconda\\envs\\pvb_cpu (PyMOL 3.x).

Usage:  python redraw_pymol_renders.py
"""

from __future__ import annotations

import subprocess
import sys
import tempfile
import os
from pathlib import Path

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

import make_case_study_figures as m  # noqa: E402

PYMOL_BAT = r"D:\Anaconda\envs\pvb_cpu\Scripts\pymol.bat"

# Per-panel camera framing: (selection expression, distance factor).
# factor < 1 moves the camera closer; focus selections are framed first.
FOCUS = {
    "6uon_full_md": ("all", 0.52),
    "6uon_full_prediction": ("all", 0.52),
    "2ptt_full_md": ("all", 0.72),
    "2ptt_full_prediction": ("all", 0.72),
    "2ptt_apo_bound_overlay": (
        "(obj0 and chain B and resi 38-43) or (obj0 and chain A and resi 37-41) "
        "or (obj1 and chain A and resi 38-43) or (obj2 and chain A and resi 37-41)",
        1.18,
    ),
    "2ptt_interface_contacts": (
        "(obj0 and chain B and resi 35+38+39+40+41+42+43+91+93) "
        "or (obj0 and chain A and resi 37+38+39+40+41+87+93)",
        1.05,
    ),
    "6uon_peptide_overlay": (
        "(obj0 and chain E and resi 1-10) or (obj1 and chain C and resi 1-10)",
        1.12,
    ),
    "6uon_contacts": (
        "(obj0 and chain E and resi 1-10) or (obj0 and chain A and resi 97) "
        "or (obj0 and chain B and resi 95+97)",
        0.85,
    ),
}

# Per-panel base-cartoon transparency overrides (0 = opaque).
TRANSPARENCY_OVERRIDE = {
    "6uon_peptide_overlay": 0.45,
}


# ---------------------------------------------------------------------------
# RMSF colour helpers (same blue sequential scale as the manuscript)
# ---------------------------------------------------------------------------


def sample_rmsf_lut(n=48):
    """Sample the paper's RMSF_CMAP (blue sequential) into [r,g,b] stops."""
    stops = np.array(
        [[r / 255.0 for r in hex_to_rgb(c)] for c in m.SEQUENTIAL],
        dtype=float,
    )
    # LinearSegmentedColormap semantics: linear between the stops.
    xs = np.linspace(0.0, 1.0, len(stops))
    out = []
    for value in np.linspace(0.0, 1.0, n):
        right = np.searchsorted(xs, value)
        if right == 0:
            color = stops[0]
        elif right >= len(stops):
            color = stops[-1]
        else:
            left = right - 1
            t = (value - xs[left]) / (xs[right] - xs[left])
            color = stops[left] * (1 - t) + stops[right] * t
        out.append([int(round(c * 255)) for c in color])
    return out


def hex_to_rgb(color_hex):
    color_hex = color_hex.lstrip("#")
    return [int(color_hex[i : i + 2], 16) for i in (0, 2, 4)]


def rmsf_inverse_table():
    """Map hex colour -> RMSF value by inverting the manuscript colormap."""
    samples = []
    for value in np.linspace(m.RMSF_MIN, m.RMSF_MAX, 4096):
        samples.append((float(value), np.array(hex_to_rgb(m.rmsf_hex(value)), dtype=int)))
    table = {}
    for value, rgb in samples:
        table[tuple(rgb)] = value
    return samples, table


INV_SAMPLES, INV_TABLE = rmsf_inverse_table()


def rmsf_from_hex(color_hex):
    rgb = tuple(hex_to_rgb(color_hex))
    if rgb in INV_TABLE:
        return INV_TABLE[rgb]
    best = min(INV_SAMPLES, key=lambda item: int(np.sum((item[1] - np.array(rgb)) ** 2)))
    return best[0]


def pymol_color(hex_color):
    """PyMOL accepts hex colors as 0xRRGGBB."""
    hex_color = hex_color.lstrip("#").upper()
    return "0x" + hex_color


# ---------------------------------------------------------------------------
# PDB annotation: write RMSF into the B-factor column
# ---------------------------------------------------------------------------


def annotate_bfactor(pdb_path, residue_values, out_path):
    lines = []
    for atom in m.parse_pdb_atoms(pdb_path):
        line = atom["line"]
        value = residue_values.get((atom["chain"], atom["resi"]))
        if value is not None and atom["record"] == "ATOM":
            line = line[:60] + f"{float(value):6.2f}" + line[66:]
        lines.append(line)
    Path(out_path).write_text("\n".join(lines) + "\nEND\n", encoding="utf-8")
    return out_path


def spectrum_pml_lines(model, obj_name, bfactor_pdb, mapped_residues):
    """Colour a whole structure by b-factor (blue sequential scale)."""
    lut = sample_rmsf_lut()
    colors = " ".join(pymol_color("#%02X%02X%02X" % tuple(c)) for c in lut)
    lines = [
        f"load {bfactor_pdb}, {obj_name}",
        f"hide everything, {obj_name}",
        f"show cartoon, {obj_name}",
        f"spectrum b, {colors}, {obj_name}, minimum={m.RMSF_MIN:.2f}, maximum={m.RMSF_MAX:.2f}, byres=1",
    ]
    unmapped = []
    for chain in sorted({c for (c, _) in mapped_residues}):
        mapped = [str(r) for (c2, r) in mapped_residues if c2 == chain]
        if mapped:
            unmapped.append(
                f"color gray80, {obj_name} and chain {chain} and not resi {'+'.join(mapped)}"
            )
    lines.extend(unmapped)
    return lines


def semantic_pml_lines(model, obj_name):
    """Colour cartoon/sticks per residue with semantic colours."""
    base = model.get("baseStyle", {}).get("cartoon", {})
    base_color = pymol_color(base.get("color", "#D9D9D9"))
    base_transparency = int(round((1.0 - float(base.get("opacity", 1.0))) * 100))
    lines = [
        f"load {model['path']}, {obj_name}",
        f"hide everything, {obj_name}",
        f"show cartoon, {obj_name}",
        f"color {base_color}, {obj_name}",
    ]
    if base_transparency:
        lines.append(f"set cartoon_transparency, {base_transparency}, {obj_name}")

    for entry in model.get("residueStyles", []):
        selection = entry.get("selection", {})
        style = entry.get("style", {})
        chain = selection.get("chain")
        resis = selection.get("resi", [])
        if isinstance(resis, int):
            resis = [resis]
        resi_sel = "+".join(str(r) for r in resis)
        sel = f"{obj_name} and chain {chain} and resi {resi_sel}"
        cartoon = style.get("cartoon", {})
        stick = style.get("stick", {})
        if cartoon.get("color"):
            color = pymol_color(cartoon["color"])
            lines.append(f"color {color}, {sel}")
            lines.append(f"set cartoon_transparency, 0, {sel}")
        if stick.get("color"):
            stick_color = pymol_color(stick["color"])
            lines.append(f"show sticks, {sel}")
            lines.append(f"color {stick_color}, {sel}")
            lines.append(f"set stick_radius, {float(stick.get('radius', 0.18))}, {sel}")
            lines.append(f"set stick_ball, 1, {sel}")
            lines.append(f"set stick_ball_ratio, 1.8, {sel}")
            lines.append(f"set stick_ball_color, {stick_color}, {sel}")

    for entry in model.get("addStyles", []):
        selection = entry.get("selection", {})
        style = entry.get("style", {})
        chain = selection.get("chain")
        resis = selection.get("resi", [])
        if isinstance(resis, int):
            resis = [resis]
        stick = style.get("stick", {})
        if not stick.get("color"):
            continue
        stick_color = pymol_color(stick["color"])
        sel = f"{obj_name} and chain {chain} and resi {'+'.join(str(r) for r in resis)}"
        lines.append(f"show sticks, {sel}")
        lines.append(f"color {stick_color}, {sel}")
        lines.append(f"set stick_radius, {float(stick.get('radius', 0.15))}, {sel}")
        lines.append(f"set stick_ball, 1, {sel}")
        lines.append(f"set stick_ball_ratio, 1.8, {sel}")
        lines.append(f"set stick_ball_color, {stick_color}, {sel}")
    return lines


# ---------------------------------------------------------------------------
# PyMOL panel driver
# ---------------------------------------------------------------------------


def pymol_render(config_name, config):
    output = Path(config["output"])
    output.parent.mkdir(parents=True, exist_ok=True)
    width = int(config.get("width", 1200))
    height = int(config.get("height", 900))
    scale = int(config.get("deviceScaleFactor", 2))
    zoom = float(config.get("zoom", 1.0))

    with tempfile.TemporaryDirectory(prefix="pymol_render_") as tmp:
        tmp_dir = Path(tmp)
        lines = [
            "reinitialize",
            "bg_color white",
            "set ray_opaque_background, 1",
            "set ray_shadows, 1",
            "set ray_shadow_decay_factor, 0.08",
            "set ray_trace_mode, 1",
            "set antialias, 2",
        ]
        obj_names = []
        obj_index = 0
        for model in config["models"]:
            obj_name = f"obj{obj_index}"
            obj_names.append(obj_name)
            obj_index += 1
            residue_styles = model.get("residueStyles", [])
            if len(residue_styles) > 60:
                mapped = []
                residue_values = {}
                for entry in residue_styles:
                    selection = entry.get("selection", {})
                    chain = selection.get("chain")
                    resis = selection.get("resi", [])
                    if isinstance(resis, int):
                        resis = [resis]
                    style = entry.get("style", {}).get("cartoon", {})
                    for resi in resis:
                        mapped.append((chain, int(resi)))
                        residue_values[(chain, int(resi))] = rmsf_from_hex(style["color"])
                bfactor_pdb = tmp_dir / f"{obj_name}_bfactor.pdb"
                annotate_bfactor(model["path"], residue_values, bfactor_pdb)
                lines.extend(spectrum_pml_lines(model, obj_name, bfactor_pdb, mapped))
            else:
                lines.extend(semantic_pml_lines(model, obj_name))
            if config_name in TRANSPARENCY_OVERRIDE:
                lines.append(
                    f"set cartoon_transparency, {TRANSPARENCY_OVERRIDE[config_name]}, {obj_name}"
                )

        for line in config.get("lines", []):
            start = line["start"]
            end = line["end"]
            obj = f"d{len([l for l in lines if l.startswith('distance ')])}"
            lines.append(
                f"distance {obj}, ({start['x']:.3f},{start['y']:.3f},{start['z']:.3f}), "
                f"({end['x']:.3f},{end['y']:.3f},{end['z']:.3f})"
            )
            lines.append(f"set dash_gap, 0.45")
            lines.append(f"set dash_length, 0.45")
            lines.append(f"set dash_color, {pymol_color(line.get('color', '#555555'))}")
            lines.append(f"hide labels, {obj}")

        label_index = 0
        for label in config.get("labels", []):
            pos = label["position"]
            obj = f"lbl{label_index}"
            label_index += 1
            lines.append(
                f"pseudoatom {obj}, pos=({pos['x']:.3f},{pos['y']:.3f},{pos['z']:.3f})"
            )
            text = label["text"].replace("\n", " ")
            lines.append(f'label {obj}, "{text}"')
            lines.append(f"set label_size, 18")
            lines.append(f"set label_color, {pymol_color(label.get('fontColor', '#222222'))}, {obj}")
            lines.append("set label_bg_color, white")
            lines.append("set label_bg_transparency, 0.10")
            lines.append(f"set label_outline_color, white, {obj}")

        focus_sel, factor = FOCUS.get(config_name, ("all", 1.0))
        if focus_sel == "all":
            lines.append("orient")
        else:
            lines.append(f"orient ({focus_sel})")
        for rotation in config.get("rotations", []):
            axis = rotation.get("axis")
            angle = float(rotation.get("angle", 0.0))
            if axis == "y":
                lines.append(f"turn y, {angle}")
            elif axis == "x":
                lines.append(f"turn x, {angle}")
            elif axis == "z":
                lines.append(f"turn z, {angle}")
        factor = float(os.environ.get("PY_ZOOM_BOOST", "1.0")) * factor
        lines.append(
            "python\n"
            "dist = cmd.get_view()[11]\n"
            f"cmd.move('z', dist * ({factor:.3f} - 1.0))\n"
            "python end"
        )
        lines.append(f"ray {width * scale}, {height * scale}")
        lines.append(f"png {output}, dpi=300")
        lines.append("quit")

        script = tmp_dir / f"{config_name}.pml"
        script.write_text("\n".join(lines) + "\n", encoding="utf-8")
        result = subprocess.run(
            [PYMOL_BAT, "-cq", str(script)],
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=600,
        )
        if result.returncode != 0 or not output.exists():
            raise RuntimeError(
                f"PyMOL failed for {config_name}:\n{result.stdout[-3000:]}\n{result.stderr[-3000:]}"
            )


if __name__ == "__main__":
    m.render = pymol_render
    m.main()
    print("PyMOL renders complete: figures written to", m.OUTPUT_2PTT.parent)
