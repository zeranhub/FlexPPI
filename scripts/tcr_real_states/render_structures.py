# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a) experimental 3D overlay -> no matching asset; existing project PyMOL -> param inherit
# (b) experimental local structure -> no matching asset; existing project PyMOL -> param inherit
# (c) predicted delta on deposited structure -> existing project PyMOL -> param inherit
# RULE: "native run" = load pre-rendered PNG via Image.open().ax.imshow().
#       "param inherit" = drawing function below that copies Class A/B/C values.
#       If a panel says "native run" and you write a drawing function, you broke the contract.

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

"""Render deposited TCR coordinates, never invented conformers or mutations.

The only coordinate operation is a single rigid-body fit of the observed apo
TCR onto the observed holo TCR. Local views inherit that same transform. All
rendered atoms are deposited heavy atoms (highest-occupancy alternative only).
Missing coordinates are NOT completed. Structure differences are not RMSF.

Raster 3D layers are composed with vector annotations by the parent figure
script; their PSE scenes and projection metadata remain editable/reproducible.
Prior project styling: scripts/render_mechanistic_structures.py (ray lighting,
orthographic camera, smooth cartoon). Palette is the parent figure contract.
"""

import argparse
import csv
import hashlib
import json
import os
from collections import defaultdict
from pathlib import Path

import numpy as np
from Bio.PDB import MMCIFParser
from Bio.SeqUtils import seq1

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "artifacts/tcr_pmhc_reassessment_20260915/sources"
OUT = ROOT / "output/tcr_real_states_20260915/structure_assets"
ENV = Path("D:/Anaconda/envs/pvb_cpu")
os.environ["PATH"] = str(ENV / "Library/bin") + os.pathsep + os.environ.get("PATH", "")
if hasattr(os, "add_dll_directory"):
    DLL_HANDLES = [os.add_dll_directory(str(p)) for p in [ENV / "Library/bin", ENV / "DLLs"] if p.exists()]
import pymol2

PALETTE = {
    "apo": "#526581", "holo": "#69A69B", "pmhc": "#C3B69F",
    "cd3": "#D7DADD", "missing": "#B7BEC4", "distal": "#762A83",
    "delta_blue": DIVERGING[0], "delta_zero": DIVERGING[1], "delta_red": DIVERGING[2],
}
SPECS = {
    "b423": {
        "label": "B4.2.3", "apo": "5IW1", "holo": "5IVX",
        "apo_chains": ["A", "B"], "holo_chains": ["E", "F"],
        "holo_partners": ["A", "B", "P"], "holo_cd3": [],
        "common_states": [("5IW1", ["A", "B"]), ("5IW1", ["C", "D"]), ("5IW1", ["E", "F"]), ("5IVX", ["E", "F"])],
        "focus": "chain F and resi 125-140", "apo_focus": "chain B and resi 125-140",
        "landmarks": [("F", "133", "N133β", "ASN"), ("F", "134", "K134β", "LYS"), ("E", "121", "E121α", "GLU"), ("E", "125", "Y125α", "TYR")],
        "apo_landmarks": [("B", "133", "N133β", "ASN"), ("B", "134", "K134β", "LYS"), ("A", "121", "E121α", "GLU"), ("A", "125", "Y125α", "TYR")],
        "constant_start": [116, 119],
    },
    "1g4": {
        "label": "1G4-CD3", "apo": "8TW4", "holo": "9C3E",
        "apo_chains": ["A", "B"], "holo_chains": ["A", "B"],
        "holo_partners": ["H"], "holo_cd3": ["D", "E", "F", "G", "X", "Y"],
        "common_states": [(p, ["A", "B"]) for p in ["8TW4", "8TW6", "9C3E", "9BBC"]],
        "focus": "chain A and resi 98-110+176-188", "apo_focus": "chain A and resi 98-110+176-188",
        "landmarks": [("A", "104", "S104", "SER"), ("A", "182", "V182", "VAL")],
        "apo_landmarks": [("A", "104", "S104", "SER"), ("A", "182", "V182", "VAL")],
        "constant_start": [140, 150],
    },
}


def read_source(pdb):
    """Same residue-level alternate and complete-backbone policy as inference."""
    from prepare_inputs import read_cif
    structure, _ = read_cif(pdb)
    result = {}
    for chain, residues in structure.items():
        records = {}
        for residue in residues:
            author = residue['author_id']
            atoms = {}
            for atom in residue['atoms']:
                atoms[atom['label_atom_id']] = dict(xyz=atom['xyz'],element=atom['type_symbol'],occupancy=atom['occupancy_float'],
                    bfactor=atom['bfactor_float'],source_alt=atom['label_alt_id'])
            records[author] = {"resname": residue['residue_name'], "atoms": atoms, "number": residue['resseq'], "ins": residue['insertion_code']}
        if records:
            result[chain] = records
    return result


def write_atoms(path, data, chains, rotation=None, translation=None):
    """An exact atom subset in PDB format; no extra atoms or coordinate repair."""
    rows, serial = [], 1
    for chain in chains:
        for author, residue in data[chain].items():
            for atom_name, atom in residue["atoms"].items():
                xyz = np.asarray(atom["xyz"])
                if rotation is not None:
                    xyz = xyz @ rotation + translation
                field = (" " + atom_name).ljust(4) if len(atom_name) < 4 else atom_name
                rows.append(f"ATOM  {serial:5d} {field} {residue['resname']:>3s} {chain:1s}{residue['number']:4d}{residue['ins']:1s}   "
                            f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}{atom['occupancy']:6.2f}{atom['bfactor']:6.2f}          {atom['element']:>2s}  ")
                serial += 1
        rows.append("TER")
    rows.append("END")
    path.write_text("\n".join(rows) + "\n", encoding="ascii")


def rigid_fit(mobile, target):
    x, y = np.asarray(mobile), np.asarray(target)
    u, _, vh = np.linalg.svd((x - x.mean(0)).T @ (y - y.mean(0)))
    correction = np.eye(3)
    correction[2, 2] = np.linalg.det(u @ vh)
    rotation = u @ correction @ vh
    translation = y.mean(0) - x.mean(0) @ rotation
    distances = np.linalg.norm(x @ rotation + translation - y, axis=1)
    return rotation, translation, float(np.sqrt(np.mean(distances ** 2)))


def common_fit(spec, sources):
    """Use identical author-numbered residues with all four backbone atoms.

    These experimentally matched constructs have shared author numbering and
    identical amino acids at all selected positions. The full shared-residue
    mapping is exported so this claim does not rely on a heuristic align call.
    """
    mapping, mobile, target = [], [], []
    for chain_index, (apo_ch, holo_ch) in enumerate(zip(spec["apo_chains"], spec["holo_chains"])):
        sets = []
        for pdb, chains in spec["common_states"]:
            residues = sources[pdb][chains[chain_index]]
            sets.append({r for r, d in residues.items() if all(a in d["atoms"] for a in ["N", "CA", "C", "O"])})
        for author in sorted(set.intersection(*sets), key=lambda r: (int(''.join(c for c in r if c.isdigit() or c == '-')), r)):
            aas = {sources[pdb][chains[chain_index]][author]["resname"] for pdb, chains in spec["common_states"]}
            if len(aas) != 1:
                raise ValueError(f"Common author ID has different amino acids: {author}, {aas}")
            mapping.append({"apo_chain": apo_ch, "holo_chain": holo_ch, "author_id": author, "resname": aas.pop()})
            mobile.append(sources[spec["apo"]][apo_ch][author]["atoms"]["CA"]["xyz"])
            target.append(sources[spec["holo"]][holo_ch][author]["atoms"]["CA"]["xyz"])
    rot, trans, rmsd = rigid_fit(mobile, target)
    return mapping, rot, trans, rmsd


def chain_selection(chains):
    return "(" + " or ".join("chain " + ch for ch in chains) + ")"


def setup(cmd):
    cmd.reinitialize()
    cmd.bg_color("white")
    settings = dict(ray_opaque_background=1, ray_trace_mode=1, ray_trace_gain=.06,
                    ray_shadows=0, antialias=2, orthoscopic=1, depth_cue=0,
                    specular=.13, ambient=.55, direct=.45, cartoon_fancy_helices=1,
                    cartoon_fancy_sheets=1, cartoon_flat_sheets=1, cartoon_smooth_loops=1,
                    cartoon_sampling=14, max_threads=2, field_of_view=20,
                    cartoon_loop_radius=.24, stick_radius=.16)
    for key, value in settings.items():
        cmd.set(key, value)
    for name, color in PALETTE.items():
        cmd.set_color(name, mpl.colors.to_rgb(color))


def camera_basis(cmd, spec):
    """Orient sequence domains with the antigen-binding side above constant domains."""
    centers_v, centers_c = [], []
    for ch, start in zip(spec["holo_chains"], spec["constant_start"]):
        centers_v.append(cmd.get_coords(f"holo and chain {ch} and name CA and resi 1-{start - 1}").mean(0))
        centers_c.append(cmd.get_coords(f"holo and chain {ch} and name CA and resi {start}-999").mean(0))
    yy = np.mean(centers_v, 0) - np.mean(centers_c, 0)
    yy /= np.linalg.norm(yy)
    xx = centers_c[1] - centers_c[0]
    xx -= np.dot(xx, yy) * yy
    xx /= np.linalg.norm(xx)
    zz = np.cross(xx, yy)
    return np.array([xx, yy, zz])


def set_camera(cmd, coords, basis, width, height, pad=1.10):
    camera_coords = np.asarray(coords) @ basis.T
    center_camera = (camera_coords.max(0) + camera_coords.min(0)) / 2
    center_world = center_camera @ basis
    spans = np.ptp(camera_coords, axis=0)
    half_height = max(spans[1] / 2, spans[0] * height / width / 2) * pad + 1.2
    distance = half_height / np.tan(np.deg2rad(10))
    clip = max(spans[2] / 2 + 15, 20)
    # PyMOL stores this rotation in column-major order. Transpose the row-basis
    # before flattening so the actual raster and vector projection agree.
    view = [*basis.T.reshape(-1), 0., 0., -distance, *center_world, distance - clip, distance + clip, 20.]
    cmd.set_view(view)
    cmd.set("orthoscopic", 1)
    cmd.viewport(width, height)
    return {"view": list(cmd.get_view()), "basis_rows": basis.tolist(), "center_world": center_world.tolist(),
            "half_height_A": float(half_height), "width_px": width, "height_px": height,
            "projection": "Orthographic; camera_xyz = (world_xyz - center_world) @ basis_rows.T. x_px = width/2 + camera_x * height/(2*half_height_A); y_px = height/2 - camera_y * height/(2*half_height_A)."}


def project(xyz, camera):
    xyz_camera = (np.asarray(xyz) - camera["center_world"]) @ np.asarray(camera["basis_rows"]).T
    scale = camera["height_px"] / (2 * camera["half_height_A"])
    return [float(camera["width_px"] / 2 + xyz_camera[0] * scale), float(camera["height_px"] / 2 - xyz_camera[1] * scale)]


def landmark_data(cmd, spec, camera):
    rows = []
    for obj, key in [("holo", "landmarks"), ("apo", "apo_landmarks")]:
        for ch, author, label, aa in spec[key]:
            sel = f"{obj} and chain {ch} and resi {author} and name CA"
            if cmd.count_atoms(sel) != 1:
                raise ValueError(f"Cannot locate exact landmark {sel}")
            atom = cmd.get_model(sel).atom[0]
            if atom.resn != aa:
                raise ValueError(f"Wrong residue for {label}: {atom.resn}")
            rows.append({"object": obj, "chain": ch, "author_id": author, "label": label,
                         "resname": aa, "xyz_A": list(atom.coord), "pixel_xy_top_left": project(atom.coord, camera)})
    return rows


def export_scene(cmd, name, camera, spec, metadata):
    cmd.hide("labels", "all")
    path = OUT / (name + ".png")
    cmd.png(str(path), width=camera["width_px"], height=camera["height_px"], dpi=300, ray=1)
    cmd.scene(name, "store")
    metadata["assets"].append({"name": name, "path": str(path), "camera": camera,
                               "landmarks": landmark_data(cmd, spec, camera)})
    print("Rendered", path.name, flush=True)


def render_system(name, predictions=None, contrasts=None, only_delta=False):
    spec = SPECS[name]
    source_ids = sorted({p for p, _ in spec["common_states"]})
    sources = {p: read_source(p) for p in source_ids}
    mapping, rotation, translation, rmsd = common_fit(spec, sources)
    holo_file, apo_file = OUT / (name + "_holo_deposited.pdb"), OUT / (name + "_apo_rigid_aligned.pdb")
    write_atoms(holo_file, sources[spec["holo"]], spec["holo_chains"] + spec["holo_partners"] + spec["holo_cd3"])
    write_atoms(apo_file, sources[spec["apo"]], spec["apo_chains"], rotation, translation)
    metadata = {"system": name, "source_pdb": {"apo": spec["apo"], "holo": spec["holo"]},
                "source_sha256": {p: hashlib.sha256((SOURCE / (p + ".cif")).read_bytes()).hexdigest() for p in source_ids},
                "common_state_set": spec["common_states"], "matching": mapping, "n_CA_for_rigid_fit": len(mapping),
                "static_CA_RMSD_A": rmsd, "apo_rotation_row_vector": rotation.tolist(),
                "apo_translation_A": translation.tolist(), "local_refit": False,
                "geometry_policy": "All displayed atoms come from deposited model 1; same one-alternate-per-residue policy as inference; no missing-atom construction, mutation, invented bond, or morph interpolation.",
                "partner_context": {"pMHC": spec["holo_partners"], "CD3": spec["holo_cd3"], "apo_non_TCR_chains_shown": False},
                "palette": PALETTE, "assets": [], "delta_maps": []}
    with pymol2.PyMOL() as pm:
        cmd = pm.cmd
        setup(cmd)
        cmd.load(str(holo_file), "holo")
        cmd.load(str(apo_file), "apo")
        cmd.hide("everything")
        cmd.show("cartoon", "all")
        cmd.dss("all")
        cmd.color("apo", "apo")
        cmd.color("holo", "holo")
        cmd.color("pmhc", "holo and " + chain_selection(spec["holo_partners"]))
        if spec["holo_cd3"]:
            cmd.color("cd3", "holo and " + chain_selection(spec["holo_cd3"]))
        # Matched structure comparisons: grey unpaired holo residues, not zeros.
        matched_holo = set((row["holo_chain"], row["author_id"]) for row in mapping)
        for ch in spec["holo_chains"]:
            for author in sources[spec["holo"]][ch]:
                if (ch, author) not in matched_holo:
                    cmd.color("missing", f"holo and chain {ch} and resi {author}")
        marker_selection = " or ".join(f"(chain {c} and resi {r})" for c, r, _, _ in spec["landmarks"])
        cmd.create("literature_sites", f"holo and name CA and ({marker_selection})")
        cmd.hide("everything", "literature_sites")
        cmd.show("spheres", "literature_sites")
        cmd.set("sphere_scale", .60, "literature_sites")
        cmd.color("distal", "literature_sites")
        basis = camera_basis(cmd, spec)
        whole_camera = set_camera(cmd, cmd.get_coords("all and name CA"), basis, 1600, 1400)
        if not only_delta:
            export_scene(cmd, name + "_experimental_overlay", whole_camera, spec, metadata)
            cmd.select("focus", "holo and (" + spec["focus"] + ")")
            cmd.select("apo_focus", "apo and (" + spec["apo_focus"] + ")")
            cmd.select("local_context", "byres ((holo or apo) within 6 of (focus or apo_focus))")
            cmd.hide("cartoon", "all and not local_context")
            cmd.hide("spheres", "literature_sites")
            for obj, key, color in [("holo", "landmarks", "distal"), ("apo", "apo_landmarks", "apo")]:
                sel = obj + " and (" + " or ".join(f"(chain {ch} and resi {author})" for ch, author, _, _ in spec[key]) + ")"
                cmd.show("sticks", sel)
                cmd.color(color, sel)
                cmd.show("spheres", sel + " and name CA")
                cmd.set("sphere_scale", .28, sel)
            local_camera = set_camera(cmd, cmd.get_coords("local_context"), basis, 1600, 1200, pad=1.13)
            export_scene(cmd, name + "_distal_local", local_camera, spec, metadata)
            metadata["local_context_selection"] = "All deposited residues with any atom within 6 A of either state's predeclared focus, after whole-TCR alignment; no local refit."
        if predictions is not None:
            selected = contrasts[name]
            values = defaultdict(list)
            with Path(predictions).open(newline="", encoding="utf-8-sig") as stream:
                for row in csv.DictReader(stream):
                    if row.get("contrast") != selected or not row.get("delta_rmsf") or row.get('common_all','').lower() != 'true':
                        continue
                    values[(row["chain"], row["author_id"])].append(float(row["delta_rmsf"]))
            if not values:
                raise ValueError(f"No delta_rmsf records for contrast {selected}")
            cmd.hide("everything")
            cmd.show("cartoon", "holo")
            cmd.color("missing", "holo")
            cmd.color("cd3", "holo and " + chain_selection(spec["holo_partners"] + spec["holo_cd3"]))
            cmd.alter("holo", "b=-999")
            for (ch, author), vals in values.items():
                if ch not in spec["holo_chains"]:
                    continue
                cmd.alter(f"holo and chain {ch} and resi {author}", "b=" + str(float(np.mean(vals))))
            cmd.spectrum("b", "delta_blue delta_zero delta_red", "holo and b > -998", minimum=-1, maximum=1)
            delta_camera = set_camera(cmd, cmd.get_coords("holo and name CA"), basis, 1600, 1400)
            export_scene(cmd, name + "_delta_rmsf", delta_camera, spec, metadata)
            metadata["delta_maps"].append({"contrast": selected, "source_csv": str(predictions), "metric": "delta_rmsf",
                                           "aggregation": "arithmetic mean over five checkpoints; TCR positions shared by all experimental references only",
                                           "shared_range_A": [-1, 1], "n_residues": len(values),
                                           "n_clipped": int(sum(abs(np.mean(v)) > 1 for v in values.values())),
                                           "unmatched_or_missing": "grey; not assigned zero"})
        cmd.save(str(OUT / (name + "_structure_scenes.pse")))
    (OUT / (name + "_render_provenance.json")).write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    return metadata


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--system", choices=list(SPECS) + ["both"], default="both")
    parser.add_argument("--predictions", type=Path)
    parser.add_argument("--b423-contrast")
    parser.add_argument("--1g4-contrast", dest="ig4_contrast")
    parser.add_argument("--only-delta", action="store_true")
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    names = list(SPECS) if args.system == "both" else [args.system]
    records = [render_system(name, args.predictions, {"b423": args.b423_contrast, "1g4": args.ig4_contrast}, args.only_delta) for name in names]
    prior_path = OUT / "render_provenance.json"
    prior = json.loads(prior_path.read_text(encoding="utf-8")) if prior_path.exists() else []
    by_system = {row["system"]: row for row in prior}
    by_system.update({row["system"]: row for row in records})
    prior_path.write_text(json.dumps(list(by_system.values()), indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
