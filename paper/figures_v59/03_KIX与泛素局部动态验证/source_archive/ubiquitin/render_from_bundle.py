# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (e) Protein complex structure -> no matching molecular asset -> PyMOL native render
# (f) Local molecular environment / observed side chains -> no matching molecular asset -> PyMOL native render
# RULE: Structure panels are ray-traced in PyMOL; no generic drawing substitutes.
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

import argparse
import csv
import hashlib
import json
import math
import runpy
import shutil
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_BUNDLE = ROOT / "output/positive_case_validation_20261009/figure_bundle/ubiquitin"
ROTAMERS = ROOT / "output/positive_case_validation_20261009/rotamers"
ORIGINAL = Path("D:/data/dynarepo_non_nucleic_full/A003R")


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_rows(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_rows(path, rows):
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def extract(bundle, source, rotamers):
    """Extract real trajectory frames once. Rendering later needs bundled PDBs only."""
    import MDAnalysis as mda
    from MDAnalysis.lib.distances import distance_array

    bundle.mkdir(parents=True, exist_ok=True)
    inputs = bundle / "inputs"
    inputs.mkdir(exist_ok=True)
    frames = bundle / "representative_frames"
    frames.mkdir(exist_ok=True)
    source_pdb = source / "structure.pdb"
    source_xtc = source / "trajectory.xtc"
    shutil.copy2(source_pdb, inputs / "A003R_input_structure.pdb")
    for filename in [
        "A003R_Q41_Q62_joint_frequencies_5000_label_frames.csv",
        "A003R_Q41_Q62_validation_summary.csv",
        "A003R_Q41_Q62_provenance_and_QA.json",
        "A003R_Q41_Q62_angles_all_5001_frames.csv",
    ]:
        shutil.copy2(rotamers / filename, inputs / filename)
    # The analysis input copy already has its XTC index: avoids creating files beside originals.
    analysis_source = rotamers / "working_inputs"
    universe = mda.Universe(str(analysis_source / "structure.pdb"),
                             str(analysis_source / "trajectory.xtc"))
    reference = mda.Universe(str(source_pdb))
    heavy = universe.select_atoms("protein and not name H*")
    ref_heavy = reference.select_atoms("protein and not name H*")
    assert len(heavy) == len(ref_heavy)
    assert len(universe.trajectory) == 5001
    assert len(reference.select_atoms("segid A").residues) == 72
    mapping = []
    # Input chain B residue IDs are global sequence 73-119. Original author numbering
    # is not silently substituted: all PDBs in this bundle retain the input IDs.
    for residue in reference.residues:
        mapping.append(dict(chain=residue.segid, input_residue_id=residue.resid,
                            chain_local_index_1based=(residue.resid if residue.segid == "A" else residue.resid - 72),
                            zero_based_global_index=residue.ix, residue_type=residue.resname))
    write_rows(bundle / "residue_mapping.csv", mapping)

    freq_rows = read_rows(rotamers / "A003R_Q41_Q62_joint_frequencies_5000_label_frames.csv")
    angle_rows = read_rows(rotamers / "A003R_Q41_Q62_angles_all_5001_frames.csv")
    representatives = []
    contacts = []
    for site in [41, 62]:
        residue_label = f"Q{site}"
        site_ref = reference.select_atoms(f"segid A and resid {site} and not name H*")
        assert site_ref.residues[0].resname == "GLN"
        # Local backbone alignment suppresses rigid motion near the site, without
        # changing ANY internal side-chain coordinate or choosing an average conformation.
        local_resids = reference.select_atoms(
            f"segid A and byres (around 9.0 (segid A and resid {site}))").residues.resids.tolist()
        if site not in local_resids:
            local_resids.append(site)
        local_resids.sort()
        residue_string = " ".join(str(value) for value in local_resids)
        fit = universe.select_atoms(f"segid A and resid {residue_string} and name N CA C O")
        fit_ref = reference.atoms[fit.indices].positions.astype(np.float64)
        fixed_center = fit_ref.mean(axis=0)
        target_centered = fit_ref - fixed_center
        # Contact measurements are descriptive source-structure geometry, not a
        # causal conclusion about entropy or a claim of binding energy.
        other_heavy = reference.select_atoms(f"protein and not name H* and not (segid A and resid {site})")
        sidechain_ref = site_ref.select_atoms("not name N CA C O")
        dists = distance_array(sidechain_ref.positions, other_heavy.positions)
        nearest_per_atom = dists.min(axis=0)
        for residue in other_heavy.residues:
            indexes = [j for j, atom in enumerate(other_heavy) if atom.resindex == residue.ix]
            dist = float(nearest_per_atom[indexes].min())
            if dist <= 5.0:
                contacts.append(dict(site=residue_label, partner_chain=residue.segid,
                                     partner_input_residue_id=residue.resid,
                                     partner_residue_type=residue.resname,
                                     min_sidechain_to_partner_heavy_atom_A=round(dist, 3),
                                     definition="minimum source-structure heavy-atom distance; no hydrogens"))

        sorted_states = sorted([row for row in freq_rows if row["residue"] == residue_label],
                               key=lambda row: (-float(row["probability"]), int(row["joint_code"])))
        for rank, state in enumerate(sorted_states[:3], start=1):
            code = int(state["joint_code"])
            candidate_rows = [row for row in angle_rows
                              if row["residue"] == residue_label
                              and int(row["frame_zero_based"]) < 5000
                              and int(row["joint_code"]) == code]
            assert len(candidate_rows) == int(state["n_frames"])
            angles = np.array([[float(row[f"chi{i}_degrees"]) for i in [1, 2, 3]]
                               for row in candidate_rows])
            radians = np.deg2rad(angles)
            circular_mean = np.arctan2(np.sin(radians).mean(axis=0), np.cos(radians).mean(axis=0))
            # Actual observed frame nearest to its category circular mean; never mean coordinates.
            circular_dist = np.sum(1.0 - np.cos(radians - circular_mean), axis=1)
            representative_index = int(np.argmin(circular_dist))
            row = candidate_rows[representative_index]
            frame = int(row["frame_zero_based"])
            universe.trajectory[frame]
            moving_fit = fit.positions.astype(np.float64)
            moving_center = moving_fit.mean(axis=0)
            moving_centered = moving_fit - moving_center
            u, _, vt = np.linalg.svd(moving_centered.T @ target_centered)
            rotation = u @ vt
            if np.linalg.det(rotation) < 0:
                u[:, -1] *= -1
                rotation = u @ vt
            aligned = (heavy.positions.astype(np.float64) - moving_center) @ rotation + fixed_center
            backbone_rmsd = np.sqrt(np.mean(np.sum((moving_centered @ rotation - target_centered) ** 2, axis=1)))
            saved_positions = heavy.positions.copy()
            heavy.positions = aligned.astype(np.float32)
            output_pdb = frames / f"{residue_label}_rank{rank}_code{code:02d}_frame{frame:04d}.pdb"
            with mda.Writer(str(output_pdb), n_atoms=len(heavy)) as writer:
                writer.write(heavy)
            heavy.positions = saved_positions
            representatives.append(dict(
                residue=residue_label, chain="A", input_residue_id=site,
                rank=rank, joint_code=code,
                chi1_bin=state["chi1_bin"], chi2_bin=state["chi2_bin"], chi3_bin=state["chi3_bin"],
                n_frames=int(state["n_frames"]), probability=float(state["probability"]),
                percentage=round(100 * float(state["probability"]), 2),
                trajectory_frame_zero_based=frame, time_ns=float(row["time_ns"]),
                chi1_degrees=float(row["chi1_degrees"]), chi2_degrees=float(row["chi2_degrees"]),
                chi3_degrees=float(row["chi3_degrees"]),
                aligned_backbone_rmsd_A=float(backbone_rmsd),
                alignment_atom_count=len(fit), alignment_residue_ids="+".join(map(str, local_resids)),
                pdb_path=output_pdb.relative_to(bundle).as_posix(),
                selection_rule="actual frame nearest circular mean of its joint bin; top 3 bins by all 5000 label frames",
                coordinate_processing="one rigid rotation and translation using local N/CA/C/O; internal geometry unchanged",
            ))
    write_rows(bundle / "representative_frames.csv", representatives)
    write_rows(bundle / "reference_environment_contacts.csv", contacts)
    provenance = dict(
        system="2QHO ubiquitin-EDD1 UBA; DynaRepo A003R",
        source_structure=str(source_pdb), source_trajectory=str(source_xtc),
        source_structure_sha256=sha256(source_pdb), source_trajectory_sha256=sha256(source_xtc),
        n_frames_in_file=5001, n_frames_used_by_labels=5000,
        frame_spacing_ns=0.1, n_windows=5, n_frames_per_window=1000,
        modeled_chain="A ubiquitin, 72 residues in actual input",
        partner_chain="B EDD1 UBA, 47 residues, input global residue IDs 73-119",
        residue_mapping="Q41 = A:41 = global zero-based 40; Q62 = A:62 = global zero-based 61",
        representation="Only the three most occupied observed categories per site are shown as structure sticks. Full 27-category frequencies are retained.",
        claim_boundary="Side-chain structures shown here are MD observations, not model-predicted conformations. The model predicts RMSF and side-chain entropy.",
        category_boundary="Joint bins follow original label code; bins need not correspond to physically separated rotamers, especially across angle boundaries.",
        extraction_software=dict(MDAnalysis=mda.__version__, numpy=np.__version__),
    )
    (bundle / "provenance.json").write_text(json.dumps(provenance, indent=2), encoding="utf-8")
    print(json.dumps(dict(bundle=str(bundle), representatives=representatives, contacts=contacts), indent=2, default=int))


def save_ubiquitin_sessions(bundle, cmd, views, source_coordinates, sasa_rows):
    """Save full molecular representations and direct-open panel-f sessions."""
    scene_audit = {}
    for scene, expected in views.items():
        cmd.scene(scene, "recall", animate=0)
        camera_difference = float(np.max(np.abs(np.asarray(cmd.get_view()) - expected)))
        assert camera_difference < 1e-4, (scene, camera_difference)
        ca_by_chain = {chain: cmd.count_atoms(
            f"complex and chain {chain} and name CA and rep cartoon")
                       for chain in ["A", "B"]}
        assert ca_by_chain == {"A": 72, "B": 47}, (scene, ca_by_chain)
        current_coordinates = np.asarray(cmd.get_coords("complex"))
        assert np.array_equal(current_coordinates, source_coordinates)
        scene_audit[scene] = dict(
            enabled_objects=cmd.get_names("objects", enabled_only=1),
            view=list(cmd.get_view()), complete_complex_CA_cartoon=ca_by_chain,
            complete_complex_CA_cartoon_total=sum(ca_by_chain.values()),
            max_camera_difference=camera_difference,
            atomic_coordinates_unchanged=True)
    scene_audit["source_structure_sidechain_sasa_A2"] = {
        row["residue"].removeprefix("Q"): row["sidechain_SASA_A2"] for row in sasa_rows}
    scene_audit["secondary_structure_CA"] = {
        row["residue"].removeprefix("Q"): row["secondary_structure_CA"] for row in sasa_rows}
    scene_audit["solvent_accessibility_protocol"] = (
        "Full complex 935 heavy atoms; probe radius 1.4 A; "
        "PyMOL3.0.0 dot_solvent=1 dot_density=4")
    (bundle / "pse_scene_audit.json").write_text(
        json.dumps(scene_audit, indent=2) + "\n", encoding="utf-8")
    cmd.scene("01_complex", "recall", animate=0)
    cmd.viewport(1500, 1250)
    cmd.save(str(bundle / "ubiquitin_2QHO.pse"))
    direct_audit = {}
    for site in [41, 62]:
        scene = f"02_Q{site}_MD_observed"
        cmd.scene(scene, "recall", animate=0)
        cmd.viewport(1500, 1250)
        cmd.save(str(bundle / f"f_Q{site}_local.pse"))
        pml = [
            "# Direct panel-f entry with the complete complex represented.",
            "# File > Run Script; all paths are derived from this file.",
            "python", "from pathlib import Path", "from pymol import cmd",
            "_f_bundle = Path(__script__).resolve().parent",
            'cmd.load(str(_f_bundle / "ubiquitin_2QHO.pse"))',
            f'cmd.scene("{scene}", "recall", animate=0)', "python end"]
        (bundle / f"f_Q{site}_local.pml").write_text(
            "\n".join(pml) + "\n", encoding="utf-8")
        direct_audit[str(site)] = dict(scene=scene,
            complete_complex_CA_cartoon_total=cmd.count_atoms(
                "complex and name CA and rep cartoon"),
            enabled_objects=cmd.get_names("objects", enabled_only=1),
            max_camera_difference=float(np.max(np.abs(
                np.asarray(cmd.get_view())-views[scene]))),
            atomic_coordinates_unchanged=np.array_equal(
                np.asarray(cmd.get_coords("complex")), source_coordinates),
            saved_scene_list=cmd.get_scene_list())
    (bundle / "f_local_entry_QA.json").write_text(json.dumps(dict(
        source_session_sha256=sha256(bundle / "ubiquitin_2QHO.pse"),
        direct_open=direct_audit, complete_protein_representation=True), indent=2)+"\n",
        encoding="utf-8")
    cmd.scene("01_complex", "recall", animate=0)
    cmd.viewport(1500, 1250)


def render(bundle, previews=False, build_only=False, command=None):
    """Create portable PSE, relative-path PML, and matching PNGs from bundle only."""
    if command is None:
        import pymol
        from pymol import cmd
        pymol.finish_launching(["pymol", "-cq"])
    else:
        cmd = command
    images = bundle / "renders"
    images.mkdir(exist_ok=True)
    crop_api = runpy.run_path(str(bundle / "whole_protein_crop.py"))
    plan_crop = crop_api["plan_crop"]
    render_crop = crop_api["render_crop"]
    crop_manifest = dict(
        rule="Keep all 119 residues represented; crop only the complete rendered image",
        source_structure="inputs/A003R_input_structure.pdb", cases={})
    cached_crops = (json.loads((bundle / "crop_manifest.json").read_text(encoding="utf-8"))
                    if build_only else None)
    rows = read_rows(bundle / "representative_frames.csv")
    bundled_inputs = [bundle / "inputs/A003R_input_structure.pdb"] + [
        bundle / row["pdb_path"] for row in rows]
    input_hashes_before = {str(path.relative_to(bundle)): sha256(path)
                           for path in bundled_inputs}
    cmd.reinitialize()
    cmd.set("auto_zoom", 0)
    cmd.load(str(bundle / "inputs/A003R_input_structure.pdb"), "complex")
    cmd.remove("hydro")
    cmd.dss("complex")
    assert cmd.count_atoms("complex and chain A and name CA") == 72
    assert cmd.count_atoms("complex and chain B and name CA") == 47
    source_coordinates = np.asarray(cmd.get_coords("complex")).copy()
    # Structural accessibility is computed on the same complete input complex.
    cmd.set("dot_solvent", 1)
    cmd.set("solvent_radius", 1.4)
    cmd.set("dot_density", 4)
    cmd.get_area("complex", load_b=1)
    sasa_atoms = []
    cmd.iterate("complex and chain A and resi 41+62", "sasa_atoms.append((resi,name,b,ss))",
                space={"sasa_atoms": sasa_atoms})
    sasa_rows = []
    for site in [41, 62]:
        sasa_rows.append(dict(
            residue=f"Q{site}",
            secondary_structure_CA=next(atom[3] for atom in sasa_atoms if atom[0] == str(site) and atom[1] == "CA"),
            sidechain_SASA_A2=sum(atom[2] for atom in sasa_atoms if atom[0] == str(site) and atom[1] not in ["N", "CA", "C", "O"]),
            whole_residue_SASA_A2=sum(atom[2] for atom in sasa_atoms if atom[0] == str(site)),
            probe_radius_A=1.4, dot_density=4, all_complex_heavy_atoms=cmd.count_atoms("complex"),
            software="PyMOL3.0.0 get_area dot_solvent=1", structure="A003R input structure, not MD-averaged"))
    write_rows(bundle / "reference_solvent_accessibility.csv", sasa_rows)
    cmd.set_color("ubiquitin_grey", [0.64, 0.66, 0.68])
    cmd.set_color("uba_green", [0.1294118, 0.4, 0.6745098])
    cmd.set_color("q41_orange", [0.7215686, 0.4745098, 0.2117647])
    cmd.set_color("q62_purple", [0.2274510, 0.5529412, 0.5490196])
    cmd.set_color("context_grey", [0.78, 0.79, 0.80])
    cmd.bg_color("white")
    for key, value in dict(orthoscopic=1, ray_opaque_background=1, antialias=2,
                           ray_trace_mode=0, ray_trace_gain=0.08, ray_shadows=0,
                           ambient=0.58, direct=0.55, specular=0.18,
                           cartoon_fancy_helices=1, cartoon_flat_sheets=0,
                           cartoon_smooth_loops=0, cartoon_side_chain_helper=0,
                           cartoon_round_helices=1, cartoon_loop_radius=0.22,
                           ribbon_trace_atoms=1, ribbon_sampling=5,
                           ribbon_radius=0.07, ribbon_width=1.2,
                           stick_radius=0.14, sphere_scale=0.35,
                           label_color="black", label_size=50,
                           label_font_id=7, depth_cue=0, fog=0).items():
        cmd.set(key, value)
    cmd.viewport(1500, 1250)
    cmd.hide("everything", "all")
    cmd.show("cartoon", "complex")
    cmd.color("ubiquitin_grey", "complex and chain A")
    cmd.color("uba_green", "complex and chain B")
    cmd.show("sticks", "complex and chain A and resi 41+62")
    cmd.color("q41_orange", "complex and chain A and resi 41")
    cmd.color("q62_purple", "complex and chain A and resi 62")
    cmd.show("spheres", "complex and chain A and resi 41+62 and name CA")
    cmd.orient("complex and name CA")
    cmd.turn("x", 15)
    cmd.turn("y", -30)
    cmd.turn("z", -12)
    cmd.zoom("complex", 3.0)
    cmd.label("complex and chain A and resi 41 and name CA", '"Q41"')
    cmd.label("complex and chain A and resi 62 and name CA", '"Q62"')
    cmd.set("label_position", [2.5, -1.5, 1.0], "complex and chain A and resi 41")
    cmd.set("label_position", [1.5, 1.5, 1.0], "complex and chain A and resi 62")
    overview = list(cmd.get_view())
    # Enlarge 1.23 times within the unchanged orientation; retain both whole chains.
    shift_z = overview[11] * (1 / 1.23 - 1)
    overview[11] += shift_z
    overview[15] -= shift_z
    overview[16] -= shift_z
    cmd.set_view(overview)
    rotation = np.asarray(overview[:9]).reshape(3, 3)
    for chain, name, text in [("A", "ubiquitin_label", "Ub"), ("B", "uba_label", "UBA")]:
        center = np.asarray(cmd.get_coords(f"complex and chain {chain} and name CA")).mean(axis=0)
        # PyMOL/OpenGL stores rotation columns, so world-space screen-up is column 1.
        anchor = center + rotation[:, 1] * 14.0 + rotation[:, 2] * 3.0
        cmd.pseudoatom(name, pos=anchor.tolist(), label=text)
        cmd.hide("everything", name)
        cmd.show("labels", name)
        cmd.set("label_position", [0, 0, 0], name)
    cmd.scene("01_complex", "store")
    if not build_only:
        cmd.png(str(images / "ubiquitin_context.png"), width=1500, height=1250, dpi=600, ray=1)
    else:
        assert (images / "ubiquitin_context.png").exists()
    views = {"01_complex": list(overview)}
    for row in rows:
        name = f'{row["residue"]}_observed_rank{row["rank"]}'
        cmd.load(str(bundle / row["pdb_path"]), name)
        cmd.remove(f"{name} and hydro")
        cmd.hide("everything", name)
        cmd.disable(name)
    for site in [41, 62]:
        color = "q41_orange" if site == 41 else "q62_purple"
        cmd.hide("everything", "all")
        cmd.disable("all")
        cmd.enable("complex")
        cmd.color("context_grey", "complex")
        cmd.select(f"Q{site}_reference_environment",
                   f"complex and byres (complex within 9.0 of (complex and chain A and resi {site}))")
        # The complete protein remains represented. The small manuscript panel
        # is obtained by cropping a full-complex render, never by truncating
        # the molecular object or hiding the backbone outside this neighborhood.
        cmd.show("cartoon", "complex")
        assert cmd.count_atoms("complex and name CA and rep cartoon") == 119
        cmd.set("cartoon_transparency", 0.12, "complex")
        # Nearby side chains supply structural context; only the nearest contacts
        # are shown so MD conformers retain the visual emphasis at 55 mm print width.
        cmd.select(f"Q{site}_neighbors",
                   f"complex and byres (complex within 3.8 of (complex and chain A and resi {site} and not name N+CA+C+O)) and not (chain A and resi {site})")
        # Retain the true N-CA-C backbone and CA-CB roots. Omitting CA previously
        # left apparently floating side chains in the molecular representation.
        cmd.show("sticks", f"Q{site}_neighbors and not hydro")
        cmd.set("stick_radius", 0.10, f"Q{site}_neighbors")
        cmd.set("stick_transparency", 0.28, f"Q{site}_neighbors")
        cmd.color("ubiquitin_grey", f"Q{site}_neighbors")
        cmd.hide("everything", f"complex and chain A and resi {site} and not name N+CA+C+O")
        for row in rows:
            if row["residue"] != f"Q{site}":
                continue
            name = f'{row["residue"]}_observed_rank{row["rank"]}'
            cmd.enable(name)
            # Each observed side chain belongs to its own genuine MD backbone.
            # Full-chain CA traces and the residue's actual peptide bonds make
            # this connection explicit without grafting onto the static input.
            backbone = f"{name} and chain A and name CA"
            residue = f"{name} and chain A and resi {site}"
            cmd.show("ribbon", backbone)
            cmd.set("ribbon_color", "ubiquitin_grey", name)
            cmd.set("ribbon_transparency", [0.35, 0.55, 0.70][int(row["rank"])-1], name)
            cmd.show("sticks", residue)
            cmd.color("ubiquitin_grey", name)
            cmd.set("stick_radius", 0.07, name)
            cmd.color(color, residue + " and not name N+C+O")
            cmd.set("stick_radius", 0.19, residue + " and not name N+C+O")
            cmd.set("stick_transparency", [0.0, 0.25, 0.45][int(row["rank"]) - 1], name)
            cmd.show("spheres", residue + " and name CA")
            cmd.set("sphere_scale", 0.18, residue + " and name CA")
            assert cmd.count_atoms(backbone + " and rep ribbon") == 72
            assert cmd.count_atoms(residue + " and rep sticks") == 9
        cmd.orient(f"complex and chain A and resi {site} and name N+CA+C+O+CB")
        cmd.turn("y", 65 if site == 41 else 45)
        cmd.turn("x", 20 if site == 41 else -25)
        cmd.turn("z", -10 if site == 41 else 20)
        cmd.zoom(f"complex and chain A and resi {site}", 7.0)
        # Identical orthographic camera distance gives identical molecular scale.
        view = list(cmd.get_view())
        view[11] = -42.0
        view[15], view[16] = 25.0, 59.0
        # The preceding site's focus view changes field_of_view. Reset the
        # requested scale before planning each site independently.
        view[17] = 20.0
        cmd.set_view(view)
        # Site names are the composite panel's vector titles; keep sticks unobscured.
        plan = plan_crop(cmd, "complex", view,
                         full_size=(6000, 5000), crop_aspect=(6, 5), margin_A=4.0)
        assert plan["full_CA_count"] == plan["full_CA_cartoon_count"] == 119
        whole_scene = f"03_Q{site}_whole_for_crop"
        cmd.viewport(*plan["full_size"])
        cmd.set_view(plan["full_view"])
        cmd.scene(whole_scene, "store")
        views[whole_scene] = list(cmd.get_view())
        scene = f"02_Q{site}_MD_observed"
        cmd.viewport(*plan["crop_size"])
        cmd.set_view(plan["focus_view"])
        cmd.scene(scene, "store")
        views[scene] = list(cmd.get_view())
        if not build_only:
            render_crop(cmd, plan, images / f"q{site}_whole_source.png",
                        images / f"q{site}_zoom.png", dpi=600,
                        preview_path=bundle / f"q{site}_crop_box_preview.png")
        else:
            cached = cached_crops["cases"][f"Q{site}"]
            assert plan["crop_box"] == cached["crop_box"]
            assert plan["full_size"] == cached["full_size"]
            for camera in ["full_view", "focus_view"]:
                assert np.allclose(plan[camera], cached[camera], rtol=0, atol=1e-4)
            for key in ["full_image", "crop_image", "full_image_sha256",
                        "crop_image_sha256", "exact_pixel_crop_verified", "dpi"]:
                plan[key] = cached[key]
            assert sha256(images / plan["full_image"]) == plan["full_image_sha256"]
            assert sha256(images / plan["crop_image"]) == plan["crop_image_sha256"]
            plan["cached_identical_full_render_verified"] = True
        plan.update(whole_scene=whole_scene, focus_scene=scene,
                    complete_complex_CA_per_chain={"A": 72, "B": 47})
        crop_manifest["cases"][f"Q{site}"] = plan
    input_hashes_after = {str(path.relative_to(bundle)): sha256(path)
                          for path in bundled_inputs}
    assert input_hashes_before == input_hashes_after
    assert np.array_equal(source_coordinates, np.asarray(cmd.get_coords("complex")))
    crop_manifest.update(input_PDB_sha256_before=input_hashes_before,
                         input_PDB_sha256_after=input_hashes_after,
                         all_source_and_representative_PDBs_unchanged=True,
                         input_heavy_atom_coordinates_unchanged=True)
    (bundle / "crop_manifest.json").write_text(
        json.dumps(crop_manifest, indent=2) + "\n", encoding="utf-8")
    save_ubiquitin_sessions(bundle, cmd, views, source_coordinates, sasa_rows)
    (bundle / "views.json").write_text(json.dumps(views, indent=2), encoding="utf-8")
    portable = bundle / "render_from_bundle.py"
    if Path(__file__).resolve() != portable.resolve():
        shutil.copy2(Path(__file__), portable)
    # PML contains only relative file paths and a cwd derived from itself.
    pml = [
        "# Open with PyMOL: File > Run Script, or @absolute_path/rebuild_relative.pml.",
        "# Exact complete-protein scenes are rebuilt without repeating ray tracing.",
        "# To regenerate PNGs, run render_from_bundle.py without --build-only.",
        "python", "from pathlib import Path", "import runpy",
        '_bundle = Path(__script__).resolve().parent',
        '_render_ns = runpy.run_path(str(_bundle / "render_from_bundle.py"))',
        '_render_ns["render"](_bundle, build_only=True, command=cmd)', "python end",
    ]
    (bundle / "rebuild_relative.pml").write_text("\n".join(pml) + "\n", encoding="utf-8")
    metadata = dict(pymol_version=list(cmd.get_version())[:2],
                    context_png_size_px=[1500, 1250],
                    local_png_size_px={key: value["crop_size"]
                                       for key, value in crop_manifest["cases"].items()},
                    whole_for_crop_png_size_px=[6000, 5000],
                    png_dpi=600, scenes=list(views), parallel_camera=True,
                    local_requested_field_height_A=2*42*math.tan(math.radians(10)),
                    local_generation="Exact pixel crop from complete protein ray image",
                    full_complex_CA_cartoon=119, relative_paths=True,
                    sidechain_connectivity="Full heavy-atom residues show real N-CA-C and CA-CB bonds; each MD frame has its own whole-chain CA trace",
                    literature_style_reference="Local OPUS-Rota4 Figure4, IsAb2.0 Figure3, Sequeiros-Borja2020 Figure3")
    (bundle / "render_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    for filename in ["ubiquitin_context.png", "q41_zoom.png", "q62_zoom.png",
                     "q41_whole_source.png", "q62_whole_source.png"]:
        shutil.copy2(images / filename, bundle / filename)
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--extract", action="store_true")
    parser.add_argument("--build-only", action="store_true",
                        help="Rebuild complete-protein scenes using verified existing PNGs")
    parser.add_argument("--bundle", type=Path, default=None)
    parser.add_argument("--source", type=Path, default=ORIGINAL)
    parser.add_argument("--rotamers", type=Path, default=ROTAMERS)
    args = parser.parse_args()
    bundle = args.bundle or (Path(__file__).resolve().parent if Path(__file__).name == "render_from_bundle.py" else DEFAULT_BUNDLE)
    if args.extract:
        extract(bundle, args.source, args.rotamers)
    else:
        render(bundle, build_only=args.build_only)
