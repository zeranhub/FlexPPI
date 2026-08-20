"""Vedo (VTK) molecular renderer for the FlexPPI case-study figures.

Replaces the 3Dmol.js backend inside make_case_study_figures.py by
monkey-patching its ``render()`` entry point, so all geometry, residue
selections, contact distances, labels and statistics stay identical while the
visual output is produced with smooth, lit molecular tubes, sticks, dashed
contact lines and 3D labels.

Usage:  python redraw_vedo_renders.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

import make_case_study_figures as m  # noqa: E402

import vedo  # noqa: E402
from vedo import Line, Plotter, Sphere, Text3D, Tube  # noqa: E402


# ---------------------------------------------------------------------------
# PDB helpers
# ---------------------------------------------------------------------------


def build_chains(atoms):
    """Group parsed atoms into {chain: {resi: {resname, atoms}}}."""
    chains = {}
    for atom in atoms:
        if atom["record"] != "ATOM":
            continue
        chain = chains.setdefault(atom["chain"], {})
        residue = chain.setdefault(atom["resi"], {"resname": atom["resname"], "atoms": {}})
        residue["atoms"].setdefault(
            atom["atom"], np.array([atom["x"], atom["y"], atom["z"]], dtype=float)
        )
    return chains


def ca_segments(chains, chain, residues, gap=4.2):
    """Return [(residues, coords)] for a chain, splitting at C-alpha gaps."""
    ordered = sorted(residues)
    segments = []
    current_res, current_pts = [], []
    for resi in ordered:
        ca = chains.get(chain, {}).get(resi, {}).get("atoms", {}).get("CA")
        if ca is None:
            if current_pts:
                segments.append((current_res, np.array(current_pts)))
                current_res, current_pts = [], []
            continue
        if current_pts and np.linalg.norm(ca - current_pts[-1]) > gap:
            segments.append((current_res, np.array(current_pts)))
            current_res, current_pts = [], []
        current_res.append(resi)
        current_pts.append(ca)
    if current_pts:
        segments.append((current_res, np.array(current_pts)))
    return segments


def to_vec(value):
    """Accept either a {x,y,z} dict or a sequence and return [x, y, z]."""
    if isinstance(value, dict):
        return [value["x"], value["y"], value["z"]]
    return list(value)


def tube_actors(chains, chain, specs, default_color, default_alpha, r_default, r_highlight):
    """Build smooth lit tubes for one chain, grouped by opacity."""
    actors = []
    by_alpha = {}
    for resi in sorted(chains.get(chain, {})):
        spec = specs.get((chain, resi), {})
        color = spec.get("color", default_color)
        alpha = float(spec.get("alpha", default_alpha))
        by_alpha.setdefault(alpha, []).append((resi, color))

    for alpha, items in by_alpha.items():
        if not items:
            continue
        lut, slot = [], {}
        resi_slot = {}
        for resi, color in items:
            if color not in slot:
                slot[color] = len(lut)
                lut.append(color)
            resi_slot[resi] = slot[color]
        radius = r_highlight if alpha > 0.9 else r_default
        for seg_resis, coords in ca_segments(chains, chain, [r for r, _ in items]):
            if len(coords) < 2:
                continue
            line = Line(coords).smooth(20)
            tube = Tube(line, r=radius, res=12)
            n_sections = len(seg_resis)
            section = np.minimum(np.arange(tube.npoints) // 12, n_sections - 1)
            scalars = np.array([resi_slot[seg_resis[i]] for i in section], dtype=float)
            tube.pointdata["scalars"] = scalars
            tube.cmap(lut, on="points")
            tube.alpha(alpha)
            tube.lighting("plastic", ambient=0.30, diffuse=0.82, specular=0.28, specular_power=20)
            actors.append(tube)
    return actors


def stick_actors(chains, chain, residue_colors, radius=0.17):
    """Draw bonds of selected residues as thin colored tubes + end spheres."""
    actors = []
    for resi, color in residue_colors.items():
        residue = chains.get(chain, {}).get(resi)
        if residue is None:
            continue
        heavy = {
            name: xyz
            for name, xyz in residue["atoms"].items()
            if name and name[0] not in ("H", "h")
        }
        names = list(heavy)
        for i in range(len(names)):
            for j in range(i + 1, len(names)):
                dist = float(np.linalg.norm(heavy[names[i]] - heavy[names[j]]))
                if dist < 1.78:
                    actors.append(
                        Tube(Line(heavy[names[i]], heavy[names[j]], res=1), r=radius, c=color)
                    )
        for name in names:
            actors.append(Sphere(heavy[name], r=radius * 1.02, c=color))
    return actors


def dashed_line_actors(start, end, color, width=0.08, dashes=14):
    actors = []
    start = np.asarray(to_vec(start), dtype=float)
    end = np.asarray(to_vec(end), dtype=float)
    t = np.linspace(0.0, 1.0, dashes * 2 + 1)
    pts = start[None, :] + t[:, None] * (end - start)[None, :]
    for i in range(0, dashes * 2, 2):
        actors.append(
            Tube(Line(pts[i], pts[i + 1], res=1), r=width, c=color)
        )
    return actors


def label_actors(text, position, size, color):
    actors = []
    position = np.asarray(to_vec(position), dtype=float)
    for dx, dy in ((0.10, 0.0), (-0.10, 0.0), (0.0, 0.10), (0.0, -0.10)):
        actors.append(
            Text3D(
                text,
                pos=position + [dx, dy, 0.0],
                s=size,
                c="white",
                justify="center",
            )
        )
    actors.append(Text3D(text, pos=position, s=size, c=color, justify="center"))
    return actors


def model_specs(model):
    """Flatten residueStyles/addStyles into per-residue color/alpha/stick maps."""
    specs = {}
    for entry in model.get("residueStyles", []):
        selection = entry.get("selection", {})
        style = entry.get("style", {})
        chain = selection.get("chain")
        resis = selection.get("resi", [])
        if isinstance(resis, int):
            resis = [resis]
        cartoon = style.get("cartoon", {})
        stick = style.get("stick", {})
        for resi in resis:
            spec = specs.setdefault((chain, int(resi)), {})
            if cartoon.get("color"):
                spec["color"] = cartoon["color"]
            if cartoon.get("opacity") is not None:
                spec["alpha"] = float(cartoon["opacity"])
            if stick.get("color"):
                spec["stick"] = stick["color"]
    for entry in model.get("addStyles", []):
        selection = entry.get("selection", {})
        style = entry.get("style", {})
        chain = selection.get("chain")
        resis = selection.get("resi", [])
        if isinstance(resis, int):
            resis = [resis]
        stick = style.get("stick", {})
        for resi in resis:
            if stick.get("color"):
                specs.setdefault((chain, int(resi)), {})["stick"] = stick["color"]
    return specs


def all_coordinates(models):
    coords = []
    for model in models:
        for atom in m.parse_pdb_atoms(model["path"]):
            if atom["record"] == "ATOM":
                coords.append([atom["x"], atom["y"], atom["z"]])
    return np.array(coords) if coords else np.zeros((1, 3))


def vedo_render(config_name, config):
    """Backend replacement for make_case_study_figures.render()."""
    del config_name
    output = Path(config["output"])
    output.parent.mkdir(parents=True, exist_ok=True)
    width = int(config.get("width", 1200))
    height = int(config.get("height", 900))
    scale = int(config.get("deviceScaleFactor", 2))

    chains_list = []
    for model in config["models"]:
        chains_list.append(build_chains(m.parse_pdb_atoms(model["path"])))

    coords = all_coordinates(config["models"])
    lo = coords.min(axis=0)
    hi = coords.max(axis=0)
    bbox_diag = float(np.linalg.norm(hi - lo))
    text_size = max(0.7, bbox_diag * 0.022)

    actors = []
    for model, chains in zip(config["models"], chains_list):
        base = model.get("baseStyle", {}).get("cartoon", {})
        base_color = base.get("color", "#D9D9D9")
        base_alpha = float(base.get("opacity", 1.0))
        specs = model_specs(model)
        for chain in sorted(chains):
            actors.extend(
                tube_actors(chains, chain, specs, base_color, base_alpha, 0.34, 0.46)
            )
            stick_colors = {
                resi: spec["stick"]
                for (ch, resi), spec in specs.items()
                if ch == chain and spec.get("stick")
            }
            if stick_colors:
                actors.extend(stick_actors(chains, chain, stick_colors))

    for line in config.get("lines", []):
        actors.extend(
            dashed_line_actors(
                line["start"],
                line["end"],
                line.get("color", "#555555"),
                width=0.06 + 0.02 * line.get("linewidth", 2),
            )
        )

    for label in config.get("labels", []):
        actors.extend(
            label_actors(
                label["text"],
                label["position"],
                text_size,
                label.get("fontColor", "#222222"),
            )
        )

    azimuth = elevation = roll = 0.0
    for rotation in config.get("rotations", []):
        axis = rotation.get("axis")
        angle = float(rotation.get("angle", 0.0))
        if axis == "y":
            azimuth += angle
        elif axis == "x":
            elevation += angle
        elif axis == "z":
            roll += angle

    plt = Plotter(offscreen=True, bg="white", size=(width, height))
    plt.show(
        *actors,
        zoom=float(config.get("zoom", 1.0)),
        azimuth=azimuth,
        elevation=elevation,
        roll=roll,
    )
    plt.screenshot(str(output), scale=scale)
    plt.close()


if __name__ == "__main__":
    m.render = vedo_render
    m.main()
    print("Vedo renders complete: figures written to", m.OUTPUT_2PTT.parent)
