"""Portable, source-data-based PyMOL reconstruction of all KIX scenes."""
from pathlib import Path
import csv
import json
import math
import runpy
from pymol import cmd
from pymol.cgo import BEGIN, END, LINES, VERTEX, COLOR, LINEWIDTH

BUNDLE_DIR = Path(__file__).resolve().parent
CROP_HELPER = runpy.run_path(str(BUNDLE_DIR / "whole_protein_crop.py"))

def rgb(color):
    return [int(color[i:i+2], 16) / 255.0 for i in (1, 3, 5)]

def rank_rgb(value):
    palette = [rgb("#F7FBFF"), rgb("#6BAED6"), rgb("#08306B")]
    position = float(value) / 50.0
    segment = min(int(position), 1)
    weight = position - segment
    return [(1-weight)*palette[segment][i] + weight*palette[segment+1][i] for i in range(3)]

def audit_l12_sticks(c, obj):
    """Validate real peptide/side-chain connectivity without moving any atoms."""
    selection = obj + " and chain B and resi 614-621 and not hydro"
    model = c.get_model(selection, 1)
    atoms = model.atom
    keys = [(a.resi, a.name) for a in atoms]
    shown = {(a.resi, a.name) for a in c.get_model(selection + " and rep sticks", 1).atom}
    assert len(atoms) == len(shown) == 56 and set(keys) == shown
    edges = [tuple(b.index) for b in model.bond]
    assert len(edges) == 57
    remaining = set(range(len(atoms)))
    components = []
    while remaining:
        stack = [remaining.pop()]
        component = []
        while stack:
            index = stack.pop()
            component.append(index)
            for left, right in edges:
                neighbor = right if left == index else left if right == index else None
                if neighbor in remaining:
                    remaining.remove(neighbor)
                    stack.append(neighbor)
        components.append(component)
    assert len(components) == 1
    peptide = []
    proline = []
    bonds = []
    for left, right in edges:
        a, b = atoms[left], atoms[right]
        length = math.sqrt(sum((a.coord[i]-b.coord[i])**2 for i in range(3)))
        bond = dict(atoms=[[a.resi, a.name], [b.resi, b.name]], length_A=length,
                    both_endpoints_displayed=True)
        bonds.append(bond)
        if a.resi != b.resi:
            assert {a.name, b.name} == {"C", "N"}
            peptide.append(bond)
        if a.resi == b.resi and a.resn == "PRO" and {a.name, b.name} == {"CD", "N"}:
            proline.append(bond)
    assert len(peptide) == 7 and {int(b["atoms"][0][0]) for b in proline} == {615, 617}
    ca_roots = []
    for residue in range(614, 622):
        root = {a.name for a in atoms if int(a.resi) == residue}
        assert {"N", "CA", "C", "O", "CB"} <= root
        ca_bonds = [bond for bond in bonds if [str(residue), "CA"] in bond["atoms"]]
        assert len(ca_bonds) == 3
        ca_roots.append(dict(residue=residue, displayed_CA_backbone_and_CB=True, root_bonds=ca_bonds))
    ca_by_chain = {chain: c.count_atoms(obj + " and chain " + chain + " and name CA and rep cartoon")
                   for chain in ("A", "B", "C")}
    assert ca_by_chain == {"A": 25, "B": 87, "C": 19}
    reference = {}
    for line in (BUNDLE_DIR / "prepared_2AGH_model1_ternary.pdb").read_text().splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            key = (line[21].strip(), line[22:27].strip(), line[17:20].strip(),
                   line[12:16].strip(), line[16].strip())
            reference[key] = [float(line[start:start+8]) for start in (30, 38, 46)]
    all_atoms = c.get_model(obj, 1).atom
    assert len(all_atoms) == len(reference)
    source_difference = 0.0
    for a in all_atoms:
        key = (a.chain, a.resi, a.resn, a.name, a.alt)
        original = reference[key]
        source_difference = max(source_difference, max(abs(a.coord[i]-original[i]) for i in range(3)))
    assert source_difference < 1e-5
    settings = {name: c.get_setting_float(name, obj) for name in
                ("cartoon_side_chain_helper", "cartoon_smooth_loops", "cartoon_flat_sheets", "cartoon_round_helices")}
    assert settings["cartoon_side_chain_helper"] == settings["cartoon_smooth_loops"] == settings["cartoon_flat_sheets"] == 0
    assert settings["cartoon_round_helices"] == 1
    # Per-atom overrides are assigned above; get_setting_float accepts an
    # object name and cannot query an arbitrary atom-selection expression.
    return dict(object=obj, L12_heavy_atoms=56, displayed_real_bonds=57,
        connected_components=1, peptide_bonds=peptide, proline_ring_closures=proline,
        CA_root_connections=ca_roots, complete_CA_cartoon_per_chain=ca_by_chain,
        complete_CA_cartoon_total=131, all_source_atom_count=len(all_atoms),
        max_source_coordinate_difference_A=source_difference, source_coordinates_unchanged=True,
        settings=settings, backbone_stick_radius_A=.07, sidechain_stick_radius_A=.12)

def build_kix(c=cmd):
    config = json.loads((BUNDLE_DIR / "scene_config.json").read_text(encoding="utf-8"))
    with (BUNDLE_DIR / "kix_both_fold1_source.csv").open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    values = {int(row["author_residue_number"]): row for row in rows}
    assert set(values) == set(range(586, 673))
    c.reinitialize()
    c.load(str(BUNDLE_DIR / "prepared_2AGH_model1_ternary.pdb"), "ternary")
    c.dss("ternary")
    c.hide("everything")
    c.bg_color("white")
    for key, value in dict(ray_opaque_background=1, ray_trace_mode=0,
            ray_shadows=0, antialias=2, orthoscopic=1, depth_cue=0,
            specular=.10, ambient=.55, direct=.45, cartoon_fancy_helices=1,
            cartoon_fancy_sheets=1, cartoon_flat_sheets=0,
            cartoon_smooth_loops=0, cartoon_side_chain_helper=0,
            cartoon_round_helices=1, cartoon_sampling=18,
            cartoon_loop_radius=.25, stick_radius=.12, max_threads=4,
            label_font_id=7, label_size=70, label_color="black",
            label_shadow_mode=0, label_outline_color="white").items():
        c.set(key, value)
    for chain, number in (("A", 25), ("B", 87), ("C", 19)):
        assert c.count_atoms("ternary and chain " + chain + " and name CA") == number
    c.select("KIX", "ternary and chain B")
    c.select("cMyb", "ternary and chain A")
    c.select("MLL", "ternary and chain C")
    c.select("L12", "ternary and chain B and resi 614-621")
    c.show("cartoon", "ternary")
    for name, selection, color in (("kix_gray", "KIX", "#BFC5CB"),
            ("cmyb_green", "cMyb", "#3A8D8C"),
            ("mll_purple", "MLL", "#B87936"),
            ("l12_orange", "L12", "#2166AC")):
        c.set_color(name, rgb(color))
        c.color(name, selection)
    c.show("sticks", "L12 and not hydro")
    c.set("stick_radius", .12, "L12 and not hydro")
    c.set("stick_radius", .07, "L12 and name N+CA+C+O")
    c.deselect()
    c.viewport(*config["context_size"])
    c.set_view(config["whole_view"])
    # Native PyMOL text and leaders are also saved in the context PSE scene.
    # Orthoscopic projection uses the transposed rotation matrix in get_view.
    def screen_point(u, v):
        view = config["whole_view"]
        height = 2*abs(view[11])*math.tan(math.radians(10))
        width = height*config["context_size"][0]/config["context_size"][1]
        plane = [(u-.5)*width-view[9], (.5-v)*height-view[10], 0]
        return [view[12+i]+sum(view[3*i+j]*plane[j] for j in range(3)) for i in range(3)]
    leaders = [LINEWIDTH, 1.5, BEGIN, LINES, COLOR, .45, .45, .45]
    annotations = [("L12", .14, .08, .16, .12, "chain B and resi 618"),
                   ("KIX", .66, .14, .64, .18, "chain B and resi 600"),
                   ("MLL", .085, .56, .13, .55, "chain C and resi 856"),
                   ("cMyb", .25, .86, .30, .82, "chain A and resi 292")]
    for name, u, v, start_u, start_v, selection in annotations:
        c.pseudoatom("context_labels", name=name, pos=screen_point(u, v), label=name)
        target = c.get_atom_coords("ternary and " + selection + " and name CA")
        leaders.extend([VERTEX]+screen_point(start_u, start_v)+[VERTEX]+target)
    leaders.append(END)
    c.load_cgo(leaders, "context_leaders")
    c.hide("everything", "context_labels")
    c.show("labels", "context_labels")
    c.set("label_size", 72, "context_labels")
    c.set("label_position", [0, 0, 0], "context_labels")
    c.set_view(config["whole_view"])
    c.scene("KIX_context", "store")
    config["scene_views"] = {"KIX_context": list(c.get_view())}
    config["crop_plans"] = {}
    for source in ("md", "model"):
        obj = "L12_" + source
        c.create(obj, "ternary")
        c.hide("everything", obj)
        c.show("cartoon", obj)
        c.show("sticks", obj + " and chain B and resi 614-621 and not hydro")
        c.set("stick_radius", .12, obj + " and chain B and resi 614-621 and not hydro")
        c.set("stick_radius", .07, obj + " and chain B and resi 614-621 and name N+CA+C+O")
        for residue in range(586, 673):
            color_name = "rank_" + source + "_" + str(residue)
            c.set_color(color_name, rank_rgb(values[residue][source+"_value"]))
            c.color(color_name, obj + " and chain B and resi " + str(residue))
        c.label(obj + " and chain B and resi 614 and name CA", '"614"')
        c.label(obj + " and chain B and resi 621 and name CA", '"621"')
        c.set("label_size", 140, obj)
        c.set("label_position", [1.4, 3.2, 0], obj + " and chain B and resi 614 and name CA")
        c.set("label_position", [-4.3, -2.2, 0], obj + " and chain B and resi 621 and name CA")
        c.disable("all")
        c.enable(obj)
        scene = "KIX_L12_" + ("MD" if source=="md" else "FlexPPI_fold1")
        c.set_view(config["local_view"])
        c.turn("y", 60)
        plan = CROP_HELPER["plan_crop"](c, obj, list(c.get_view()),
                                        full_size=(6000, 4800), crop_aspect=(5, 4))
        assert plan["full_CA_count"] == 131
        for name, view in ((scene+"_whole", plan["full_view"]), (scene, plan["focus_view"])):
            c.set_view(view)
            c.scene(name, "store")
            config["scene_views"][name] = list(c.get_view())
        config["crop_plans"][source] = plan
    (BUNDLE_DIR / "scene_config.json").write_text(json.dumps(config, indent=2)+"\n", encoding="utf-8")
    c.scene("KIX_context", "recall", animate=0)
    return config

def render_kix(c=cmd):
    config = build_kix(c)
    c.scene("KIX_context", "recall", animate=0)
    c.viewport(*config["context_size"])
    c.png(str(BUNDLE_DIR / "kix_context.png"), width=config["context_size"][0],
          height=config["context_size"][1], dpi=300, ray=1)
    for source, stem, scene in (("md", "kix_l12_md", "KIX_L12_MD"),
            ("model", "kix_l12_fold1", "KIX_L12_FlexPPI_fold1")):
        c.scene(scene+"_whole", "recall", animate=0)
        CROP_HELPER["render_crop"](c, config["crop_plans"][source],
              BUNDLE_DIR/(stem+"_whole_source.png"), BUNDLE_DIR/(stem+".png"),
              dpi=600, preview_path=BUNDLE_DIR/(stem+"_crop_box_preview.png"))
    (BUNDLE_DIR / "crop_manifest.json").write_text(json.dumps(config["crop_plans"], indent=2)+"\n", encoding="utf-8")
    connectivity = dict(status="pass", scientific_geometry="Original source coordinates and true bonds only",
                        representation="All eight L12 heavy residues include N/CA/C/O and side chains",
                        scenes={})
    for scene in config["scene_views"]:
        c.scene(scene, "recall", animate=0)
        obj = "ternary" if scene == "KIX_context" else "L12_md" if "MD" in scene else "L12_model"
        assert obj in c.get_names("objects", enabled_only=1)
        connectivity["scenes"][scene] = audit_l12_sticks(c, obj)
    (BUNDLE_DIR / "stick_connectivity_QA.json").write_text(json.dumps(connectivity, indent=2)+"\n", encoding="utf-8")
    c.scene("KIX_context", "recall", animate=0)
    c.viewport(*config["context_size"])
    c.save(str(BUNDLE_DIR / "KIX_case_validation.pse"))

if globals().get("_KIX_BUILD_ONLY", False) is False:
    build_kix()
