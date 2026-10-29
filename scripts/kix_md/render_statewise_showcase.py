# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (d) Protein cartoon and continuous L12 detail -> cross-type inherit -> param inherit
# No compatible protein-rendering asset exists in the skill's directory map or other/.
# Geometry and camera handling inherit scripts/kix_md/render_pilot.py, not synthetic art.
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


from pathlib import Path
import csv
import hashlib
import json
import os
import numpy as np
from PIL import Image

ENV = Path("D:/Anaconda/envs/pvb_cpu")
os.environ["PATH"] = str(ENV / "Library/bin") + os.pathsep + os.environ.get("PATH", "")
DLL_HANDLES = [os.add_dll_directory(str(p)) for p in (ENV / "Library/bin", ENV / "DLLs") if p.exists()]
import pymol2

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "artifacts/kix_md_20260914/remote_audit/both_protein.pdb"
DATA = ROOT / "artifacts/kix_md_20260914/pilot_results_20260916/relative_residues.csv"
OUT = ROOT / "output/flexppi_kix_statewise_showcase_20260916/structure_assets"


def rgb(hex_color):
    return np.array([int(hex_color[i:i+2], 16) / 255 for i in (1, 3, 5)])


def rank_rgb(value):
    assert np.isfinite(value) and 0 <= value <= 100
    position = value / 50
    segment = min(int(position), 1)
    weight = position - segment
    return (rgb(SEQUENTIAL[segment]) * (1-weight) + rgb(SEQUENTIAL[segment+1]) * weight).tolist()


def style(c):
    c.bg_color("white")
    for key, value in dict(ray_opaque_background=1, ray_trace_mode=0, ray_shadows=0,
                           antialias=2, orthoscopic=1, depth_cue=0, specular=.10,
                           ambient=.55, direct=.45, cartoon_fancy_helices=1,
                           cartoon_fancy_sheets=1, cartoon_flat_sheets=1,
                           cartoon_smooth_loops=1, cartoon_sampling=18,
                           cartoon_loop_radius=.21, stick_radius=.105,
                           max_threads=4).items():
        c.set(key, value)


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def pixel_audit(path):
    array = np.asarray(Image.open(path).convert("RGB")) / 255
    ys, xs = np.where(np.any(array < .975, axis=2))
    h, w = array.shape[:2]
    return dict(width=w, height=h, nonwhite_pixel_fraction=float(len(xs)/(h*w)),
                nonwhite_bbox_xyxy=[int(xs.min()), int(ys.min()), int(xs.max()), int(ys.max())],
                touches_canvas_edge=bool(xs.min()==0 or xs.max()==w-1 or ys.min()==0 or ys.max()==h-1))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    with DATA.open(encoding="utf-8", newline="") as handle:
        all_rows = list(csv.DictReader(handle))
    data = [row for row in all_rows if row["state"] == "both" and row["metric"] == "rmsf"
            and row["normalization"] == "percentile" and row["model_fold"] == "ensemble_mean"]
    residue_rows = {int(row["author_residue_number"]): row for row in data}
    assert len(data) == len(residue_rows) == 87
    assert set(residue_rows) == set(range(586, 673))
    for row in data:
        for field in ("md_value", "model_value"):
            assert np.isfinite(float(row[field])) and 0 <= float(row[field]) <= 100

    manifest = dict(
        structure_file=str(SOURCE), structure_sha256=sha256(SOURCE),
        data_file=str(DATA), data_sha256=sha256(DATA),
        reference_geometry="Prepared 2AGH model 1 ternary structure; not an MD snapshot or separately predicted conformation",
        state="both (cMyb + MLL)", source_data_rows_loaded=len(all_rows), selected_all_state_residues=87,
        structure_counts=dict(KIX=87, cMyb=25, MLL=19),
        selection_rationale="The full ternary geometry provides structural context. Local detail is a predefined continuous 22-residue segment B607-628, not selected by agreement.",
        local_scope="Continuous chain B residues 607-628; full L12 is 614-621. Only L12 side chains are sticks.",
        color_scope="All 22 displayed local residues carry their actual whole-KIX RMSF percentile, normalized separately within this state. Color scale is shared 0-100.",
        palette=SEQUENTIAL, scale=[0,100], clipped_values=0,
        statistic="Rank percentile 100*(average rank-1)/(87-1); model is rank of the five-checkpoint mean prediction, not mean rank.",
        n="One MD trajectory; five model checkpoints averaged. Residues are not independent MD repeats.",
        context_colors=dict(KIX="#CCD0D4", L12="#F1A340", cMyb="#AFA0BD", MLL="#D0BB94"),
        no_coordinate_modification=True, no_image_postprocessing=True,
        font_export_note="These text-free raster structure assets are placed in a vector-labeled master by the main plotting script.",
        renders=[])

    with pymol2.PyMOL() as pm:
        c = pm.cmd
        c.load(str(SOURCE), "protein")
        c.hide("everything")
        c.dss("protein")
        style(c)
        assert c.count_atoms("chain B and name CA") == 87
        assert c.count_atoms("chain A and name CA") == 25
        assert c.count_atoms("chain C and name CA") == 19
        c.show("cartoon", "protein")
        for name, selection, color in (("KIX", "chain B", "#CCD0D4"),
                                       ("cMyb", "chain A", "#AFA0BD"),
                                       ("MLL", "chain C", "#D0BB94"),
                                       ("L12", "chain B and resi 614-621", "#F1A340")):
            c.set_color(name, rgb(color).tolist())
            c.color(name, selection)
        focus = "chain B and resi 614-621"
        c.show("sticks", focus + " and not hydro and not name N+C+O")
        c.viewport(1400, 1100)
        c.orient("chain B")
        c.zoom("protein", buffer=2.2, complete=1)
        whole_view = list(c.get_view())
        context_file = OUT / "context_whole.png"
        c.png(str(context_file), width=1400, height=1100, dpi=300, ray=1)
        c.save(str(OUT / "context_whole.pse"))
        manifest["whole_view"] = whole_view

        # Screen leader target comes from actually rendered orange L12 pixels,
        # not a guessed geometric projection. Other context colors cannot pass
        # the orange hue and saturation thresholds used here.
        array = np.asarray(Image.open(context_file).convert("RGB")) / 255
        rr, gg, bb = array[..., 0], array[..., 1], array[..., 2]
        orange = (rr > .35) & (rr-gg > .12) & (gg-bb > .12) & ((rr-bb)/np.maximum(rr, 1e-6) > .48)
        ys, xs = np.where(orange)
        assert len(xs) > 100
        manifest["L12_visible_screen_target"] = dict(
            method="Centroid of visible orange L12 pixels in context_whole.png; image has no post-render crop",
            pixel_count=int(len(xs)),
            normalized_top_left_origin=[float(xs.mean()/1400), float(ys.mean()/1100)],
            normalized_bottom_left_origin=[float(xs.mean()/1400), float(1-ys.mean()/1100)],
            bbox_normalized_top_left_xyxy=[float(xs.min()/1400),float(ys.min()/1100),float(xs.max()/1400),float(ys.max()/1100)])
        manifest["renders"].append(dict(name="context_whole", camera=whole_view, audit=pixel_audit(context_file)))

        c.create("detail", "protein and chain B and resi 607-628")
        c.hide("everything")
        c.show("cartoon", "detail")
        c.show("sticks", "detail and resi 614-621 and not hydro and not name N+C+O")
        assert c.count_atoms("detail and name CA") == 22
        c.viewport(1200, 900)
        c.orient("detail and name CA")
        c.turn("x", 12)
        c.turn("y", -10)
        c.zoom("detail", buffer=2.1, complete=1)
        local_view = list(c.get_view())
        manifest["local_view"] = local_view
        for source in ("md", "model"):
            for residue_id in range(607, 629):
                name = f"rank_{source}_{residue_id}"
                c.set_color(name, rank_rgb(float(residue_rows[residue_id][source+"_value"])))
                c.color(name, f"detail and resi {residue_id}")
            c.set_view(local_view)
            output_file = OUT / f"l12_{source}.png"
            c.png(str(output_file), width=1200, height=900, dpi=300, ray=1)
            c.save(str(OUT / f"l12_{source}.pse"))
            manifest["renders"].append(dict(name=f"l12_{source}", camera=list(c.get_view()),
                                            audit=pixel_audit(output_file),
                                            values=[dict(residue=i,percentile=float(residue_rows[i][source+"_value"])) for i in range(607,629)]))
            print(f"Rendered l12_{source}", flush=True)
    assert all(not r["audit"]["touches_canvas_edge"] for r in manifest["renders"])
    (OUT / "render_provenance.json").write_text(json.dumps(manifest, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(manifest["L12_visible_screen_target"], indent=2), flush=True)


if __name__ == "__main__":
    main()
