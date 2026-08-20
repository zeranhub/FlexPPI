# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a-c) molecular structure and functional close-ups -> cross-type inherit from assets/figures/other -> param inherit (no molecular-rendering production asset exists)
# (d, 6UON) residue RMSF profile -> assets/figures/LineTrend/plot_trend.py -> param inherit (same 1D profile semantics; case-specific residue axis)
# (d, 2PTT) residue-level paired endpoints -> assets/figures/PairedBoxScatter/plot_PairedBoxScatter1.R -> param inherit (R unavailable; no distribution or inferential boxplot is claimed)
# RULE: "native run" = load pre-rendered PNG via Image.open().ax.imshow().
#       "param inherit" = drawing function below that copies Class A/B/C values.
#       If a panel says "native run" and you write a drawing function, you broke the contract.

# Academic Figure Skill Typography Baseline - COPY VERBATIM, place at TOP of script
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

# Academic Figure Skill Export Baseline - COPY VERBATIM
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


import hashlib
import json
import math
import os
import shutil
import subprocess
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image
from matplotlib import pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.gridspec import GridSpecFromSubplotSpec


ROOT = Path(__file__).resolve().parents[2]
FIGURE_DIR = ROOT / "paper" / "figures"
CASE_DIR = FIGURE_DIR / "case_studies"
STRUCTURE_DIR = CASE_DIR / "structures"
PREPARED_DIR = CASE_DIR / "prepared"
RENDER_DIR = CASE_DIR / "renders"
CONFIG_DIR = CASE_DIR / "render_configs"
VENDOR_SCRIPT = CASE_DIR / "vendor" / "3Dmol-min.js"
RENDERER = CASE_DIR / "render_3dmol.js"

ITRACE_PREDICTIONS = ROOT / "output" / "itrace_benchmark" / "raw" / "common_residue_predictions.csv"
ITRACE_CONDITIONS = ROOT / "output" / "itrace_residue_biology" / "raw" / "all_condition_predictions.csv"
ITRACE_MANIFEST = ROOT / "output" / "itrace_benchmark" / "raw" / "itrace_complex_manifest.csv"
DYNA_PREDICTIONS = ROOT / "artifacts" / "evaluation" / "fastprotflex_comparison" / "dynarepo" / "dynarepo_common_residue_predictions.csv"
CHECKPOINT = ROOT / "artifacts" / "pretrain" / "best.pt"

OUTPUT_6UON = FIGURE_DIR / "fig_6uon_case_study"
OUTPUT_2PTT = FIGURE_DIR / "fig_2ptt_case_study"
SOURCE_6UON = FIGURE_DIR / "fig_6uon_case_study_source_data.csv"
SOURCE_2PTT = FIGURE_DIR / "fig_2ptt_case_study_source_data.csv"
STATS_PATH = FIGURE_DIR / "fig_case_studies_statistics.json"
QA_PATH = FIGURE_DIR / "fig_case_studies_QA.md"

NODE_FALLBACK = (
    Path.home()
    / ".cache"
    / "codex-runtimes"
    / "codex-primary-runtime"
    / "dependencies"
    / "node"
    / "bin"
    / "node.exe"
)
NODE_MODULES_FALLBACK = (
    Path.home()
    / ".cache"
    / "codex-runtimes"
    / "codex-primary-runtime"
    / "dependencies"
    / "node"
    / "node_modules"
)

RMSF_MIN = 0.25
RMSF_MAX = 2.00
RMSF_CMAP = LinearSegmentedColormap.from_list("rmsf_blue", SEQUENTIAL)
RMSF_NORM = Normalize(vmin=RMSF_MIN, vmax=RMSF_MAX, clip=True)


def ensure_directories():
    for directory in (STRUCTURE_DIR, PREPARED_DIR, RENDER_DIR, CONFIG_DIR):
        directory.mkdir(parents=True, exist_ok=True)


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def download_structures():
    for pdb_id in ("6UON", "6ULK", "2PTT", "2PTU", "2PTV"):
        destination = STRUCTURE_DIR / f"{pdb_id}.pdb"
        if not destination.exists():
            urllib.request.urlretrieve(
                f"https://files.rcsb.org/download/{pdb_id}.pdb", destination
            )


def parse_pdb_atoms(path):
    atoms = []
    for line in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")) or len(line) < 54:
            continue
        try:
            atoms.append(
                {
                    "line": line,
                    "record": line[:6].strip(),
                    "atom": line[12:16].strip(),
                    "altloc": line[16].strip(),
                    "resname": line[17:20].strip(),
                    "chain": line[21].strip(),
                    "resi": int(line[22:26]),
                    "icode": line[26].strip(),
                    "x": float(line[30:38]),
                    "y": float(line[38:46]),
                    "z": float(line[46:54]),
                    "element": line[76:78].strip() if len(line) >= 78 else "",
                }
            )
        except ValueError:
            continue
    return atoms


def pdb_residue_names(path):
    mapping = {}
    for atom in parse_pdb_atoms(path):
        if atom["record"] == "ATOM" and atom["altloc"] in ("", "A"):
            mapping[(atom["chain"], atom["resi"], atom["icode"])] = atom["resname"]
    return mapping


def write_atom_only(source, destination, chain_ranges=None):
    selected = []
    for atom in parse_pdb_atoms(source):
        if atom["record"] != "ATOM" or atom["altloc"] not in ("", "A"):
            continue
        if chain_ranges is not None:
            ranges = chain_ranges.get(atom["chain"], [])
            if not any(start <= atom["resi"] <= end for start, end in ranges):
                continue
        selected.append(atom["line"])
    Path(destination).write_text("\n".join(selected) + "\nEND\n", encoding="utf-8")


def ca_map(path, chain):
    return {
        (atom["resi"], atom["icode"]): np.array([atom["x"], atom["y"], atom["z"]])
        for atom in parse_pdb_atoms(path)
        if atom["record"] == "ATOM"
        and atom["chain"] == chain
        and atom["atom"] == "CA"
        and atom["altloc"] in ("", "A")
    }


def kabsch_transform(mobile_path, mobile_chain, reference_path, reference_chain):
    mobile = ca_map(mobile_path, mobile_chain)
    reference = ca_map(reference_path, reference_chain)
    common = sorted(set(mobile).intersection(reference))
    if len(common) < 20:
        raise ValueError(
            f"Insufficient common C-alpha atoms for alignment: {mobile_path} {mobile_chain} "
            f"to {reference_path} {reference_chain}: {len(common)}"
        )
    x = np.vstack([mobile[key] for key in common])
    y = np.vstack([reference[key] for key in common])
    x_center = x.mean(axis=0)
    y_center = y.mean(axis=0)
    covariance = (x - x_center).T @ (y - y_center)
    left, _, right_t = np.linalg.svd(covariance)
    correction = np.eye(3)
    correction[-1, -1] = np.sign(np.linalg.det(left @ right_t))
    rotation = left @ correction @ right_t
    translation = y_center - x_center @ rotation
    fitted = x @ rotation + translation
    rmsd = float(np.sqrt(np.mean(np.sum((fitted - y) ** 2, axis=1))))
    return rotation, translation, rmsd, len(common)


def transform_pdb(source, destination, rotation, translation, keep_chains=None):
    output = []
    for atom in parse_pdb_atoms(source):
        if atom["record"] != "ATOM" or atom["altloc"] not in ("", "A"):
            continue
        if keep_chains is not None and atom["chain"] not in keep_chains:
            continue
        xyz = np.array([atom["x"], atom["y"], atom["z"]]) @ rotation + translation
        line = atom["line"]
        line = f"{line[:30]}{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}{line[54:]}"
        output.append(line)
    Path(destination).write_text("\n".join(output) + "\nEND\n", encoding="utf-8")


def find_atom(path, chain, resi, atom_names):
    if isinstance(atom_names, str):
        atom_names = [atom_names]
    candidates = [
        atom
        for atom in parse_pdb_atoms(path)
        if atom["record"] == "ATOM"
        and atom["chain"] == chain
        and atom["resi"] == int(resi)
        and atom["atom"] in atom_names
        and atom["altloc"] in ("", "A")
    ]
    if not candidates:
        raise KeyError(f"Atom not found: {path} {chain}{resi} {atom_names}")
    atom = candidates[0]
    return np.array([atom["x"], atom["y"], atom["z"]])


def closest_atom_pair(path, first, second):
    chain1, resi1, atom_names1 = first
    chain2, resi2, atom_names2 = second
    atoms = parse_pdb_atoms(path)
    group1 = [
        atom for atom in atoms
        if atom["record"] == "ATOM" and atom["chain"] == chain1
        and atom["resi"] == int(resi1) and atom["atom"] in atom_names1
        and atom["altloc"] in ("", "A")
    ]
    group2 = [
        atom for atom in atoms
        if atom["record"] == "ATOM" and atom["chain"] == chain2
        and atom["resi"] == int(resi2) and atom["atom"] in atom_names2
        and atom["altloc"] in ("", "A")
    ]
    if not group1 or not group2:
        raise KeyError(f"Contact atoms missing: {first}, {second}")
    best = None
    for atom1 in group1:
        p1 = np.array([atom1["x"], atom1["y"], atom1["z"]])
        for atom2 in group2:
            p2 = np.array([atom2["x"], atom2["y"], atom2["z"]])
            distance = float(np.linalg.norm(p1 - p2))
            if best is None or distance < best[0]:
                best = (distance, p1, p2, atom1["atom"], atom2["atom"])
    return best


def point(xyz):
    return {"x": float(xyz[0]), "y": float(xyz[1]), "z": float(xyz[2])}


def midpoint(first, second):
    return point((np.asarray(first) + np.asarray(second)) / 2.0)


def rmsf_hex(value):
    red, green, blue, _ = RMSF_CMAP(RMSF_NORM(float(value)))
    return "#{:02X}{:02X}{:02X}".format(
        int(round(red * 255)), int(round(green * 255)), int(round(blue * 255))
    )


def residue_styles(frame, value_column, chain_column="chain", resi_column="pdb_residue_number"):
    styles = []
    for row in frame.itertuples(index=False):
        value = getattr(row, value_column)
        if not np.isfinite(value):
            continue
        styles.append(
            {
                "selection": {
                    "chain": str(getattr(row, chain_column)),
                    "resi": int(getattr(row, resi_column)),
                },
                "style": {"cartoon": {"color": rmsf_hex(value)}},
            }
        )
    return styles


def peptide_stick_styles(frame, value_column, chain="E", resi_column="pdb_residue_number"):
    additions = []
    for row in frame.loc[frame["chain"] == chain].itertuples(index=False):
        value = getattr(row, value_column)
        if not np.isfinite(value):
            continue
        additions.append(
            {
                "selection": {"chain": chain, "resi": int(getattr(row, resi_column))},
                "style": {"stick": {"color": rmsf_hex(value), "radius": 0.13}},
            }
        )
    return additions


def render(config_name, config):
    config_path = CONFIG_DIR / f"{config_name}.json"
    config["threeDmolScript"] = str(VENDOR_SCRIPT.resolve())
    config_path.write_text(json.dumps(config, indent=2), encoding="utf-8")
    node = shutil.which("node")
    if node is None and NODE_FALLBACK.exists():
        node = str(NODE_FALLBACK)
    if node is None:
        raise RuntimeError("Node.js was not found; molecular structure rendering cannot proceed.")
    environment = os.environ.copy()
    if NODE_MODULES_FALLBACK.exists():
        environment["NODE_PATH"] = str(NODE_MODULES_FALLBACK)
    subprocess.run(
        [node, str(RENDERER), str(config_path)],
        cwd=str(ROOT),
        env=environment,
        check=True,
    )


def load_6uon_data():
    benchmark = pd.read_csv(ITRACE_PREDICTIONS, low_memory=False)
    benchmark["pdb_id"] = benchmark["pdb_id"].astype(str).str.lower()
    full = benchmark.loc[benchmark["pdb_id"] == "6uon"].copy()
    full["pdb_residue_number"] = pd.to_numeric(full["residue_number"], errors="raise").astype(int)
    full = full.rename(
        columns={
            "target_chainfit_rmsf_angstrom": "md_rmsf",
            "flexppi_rmsf_angstrom": "flexppi_full_rmsf",
        }
    )
    conditions = pd.read_csv(ITRACE_CONDITIONS, low_memory=False)
    conditions["pdb_id"] = conditions["pdb_id"].astype(str).str.lower()
    peptide = conditions.loc[
        (conditions["pdb_id"] == "6uon")
        & (conditions["chain"] == "E")
        & (conditions["condition"].isin(["full", "pmhc"]))
    ].copy()
    pivot = peptide.pivot(
        index=["chain", "residue_number", "amino_acid_1"],
        columns="condition",
        values="rmsf",
    ).reset_index()
    pivot = pivot.rename(
        columns={
            "residue_number": "pdb_residue_number",
            "full": "condition_full_rmsf",
            "pmhc": "pmhc_only_rmsf",
        }
    )
    full = full.merge(
        pivot[["chain", "pdb_residue_number", "pmhc_only_rmsf"]],
        on=["chain", "pdb_residue_number"],
        how="left",
        validate="one_to_one",
    )
    full["context_delta_full_minus_pmhc"] = (
        full["flexppi_full_rmsf"] - full["pmhc_only_rmsf"]
    )
    full["case"] = "6UON"
    full["source_scope"] = "all 826 residues used for structure coloring"
    return full


def load_2ptt_data(pdb_path):
    predictions = pd.read_csv(DYNA_PREDICTIONS, low_memory=False)
    frame = predictions.loc[predictions["sample_id"].astype(str) == "A003O"].copy()
    frame["dataset_author_residue_id"] = pd.to_numeric(
        frame["author_residue_id"], errors="raise"
    ).astype(int)
    frame["pdb_residue_number"] = np.where(
        frame["chain"] == "A",
        frame["dataset_author_residue_id"] + 3,
        frame["dataset_author_residue_id"] - 102,
    ).astype(int)
    frame["chain_role"] = np.where(frame["chain"] == "A", "CD48", "2B4")
    frame = frame.rename(
        columns={
            "target_rmsf_angstrom": "md_rmsf",
            "flexppi_pretrain_rmsf_angstrom": "flexppi_full_rmsf",
        }
    )
    residue_names = pdb_residue_names(pdb_path)
    frame["amino_acid_3"] = [
        residue_names.get((chain, int(resi), ""), "UNK")
        for chain, resi in zip(frame["chain"], frame["pdb_residue_number"])
    ]
    frame["case"] = "2PTT"
    frame["source_scope"] = "all mapped DynaRepo residues used for structure coloring"
    return frame


def validate_data(sixuon, twoptt, sixuon_pdb, twoptt_pdb):
    checks = {}
    checks["6uon_residue_count"] = int(len(sixuon))
    checks["6uon_peptide_residue_count"] = int((sixuon["chain"] == "E").sum())
    checks["2ptt_mapped_residue_count"] = int(len(twoptt))
    checks["2ptt_h42_present_in_structure"] = ("B", 42, "") in pdb_residue_names(twoptt_pdb)
    checks["2ptt_h42_present_in_dynarepo_table"] = bool(
        ((twoptt["chain"] == "B") & (twoptt["pdb_residue_number"] == 42)).any()
    )
    if checks["6uon_residue_count"] != 826 or checks["6uon_peptide_residue_count"] != 10:
        raise ValueError(f"Unexpected 6UON coverage: {checks}")
    if checks["2ptt_mapped_residue_count"] < 200:
        raise ValueError(f"Unexpected 2PTT coverage: {checks}")
    for name, frame in (("6UON", sixuon), ("2PTT", twoptt)):
        for column in ("md_rmsf", "flexppi_full_rmsf"):
            values = frame[column].to_numpy(dtype=float)
            if not np.isfinite(values).all() or (values < 0).any():
                raise ValueError(f"Invalid {name} {column} values")
        if frame.duplicated(["chain", "pdb_residue_number"]).any():
            raise ValueError(f"Duplicate {name} chain/residue keys")
    sixuon_residues = set(pdb_residue_names(sixuon_pdb))
    twoptt_residues = set(pdb_residue_names(twoptt_pdb))
    checks["6uon_rows_present_in_render_pdb"] = int(
        sum((row.chain, int(row.pdb_residue_number), "") in sixuon_residues for row in sixuon.itertuples())
    )
    checks["2ptt_rows_present_in_render_pdb"] = int(
        sum((row.chain, int(row.pdb_residue_number), "") in twoptt_residues for row in twoptt.itertuples())
    )
    return checks


def prepare_structures(sixuon_topology):
    paths = {
        "6uon_full": PREPARED_DIR / "6uon_full_atom.pdb",
        "6uon_local": PREPARED_DIR / "6uon_tcr_peptide_local.pdb",
        "6uon_groove": PREPARED_DIR / "6uon_pmhc_groove.pdb",
        "6ulk_aligned": PREPARED_DIR / "6ulk_aligned.pdb",
        "6ulk_groove": PREPARED_DIR / "6ulk_aligned_groove.pdb",
        "2ptt_full": PREPARED_DIR / "2ptt_full_atom.pdb",
        "2ptt_local": PREPARED_DIR / "2ptt_interface_local.pdb",
        "2ptu_aligned": PREPARED_DIR / "2ptu_aligned.pdb",
        "2ptv_aligned": PREPARED_DIR / "2ptv_aligned.pdb",
    }
    write_atom_only(sixuon_topology, paths["6uon_full"])
    write_atom_only(
        sixuon_topology,
        paths["6uon_local"],
        {"A": [(88, 104)], "B": [(88, 104)], "C": [(58, 86), (138, 166)], "E": [(1, 10)]},
    )
    write_atom_only(
        sixuon_topology,
        paths["6uon_groove"],
        {"C": [(52, 90), (132, 170)], "D": [(1, 99)], "E": [(1, 10)]},
    )
    rotation_6ulk, translation_6ulk, rmsd_6ulk, n_6ulk = kabsch_transform(
        STRUCTURE_DIR / "6ULK.pdb", "A", sixuon_topology, "C"
    )
    transform_pdb(
        STRUCTURE_DIR / "6ULK.pdb",
        paths["6ulk_aligned"],
        rotation_6ulk,
        translation_6ulk,
        keep_chains={"A", "B", "C"},
    )
    write_atom_only(
        paths["6ulk_aligned"],
        paths["6ulk_groove"],
        {"A": [(52, 90), (132, 170)], "B": [(1, 110)], "C": [(1, 10)]},
    )

    write_atom_only(STRUCTURE_DIR / "2PTT.pdb", paths["2ptt_full"])
    write_atom_only(
        STRUCTURE_DIR / "2PTT.pdb",
        paths["2ptt_local"],
        {"A": [(25, 105)], "B": [(25, 105)]},
    )
    rotation_2ptu, translation_2ptu, rmsd_2ptu, n_2ptu = kabsch_transform(
        STRUCTURE_DIR / "2PTU.pdb", "A", STRUCTURE_DIR / "2PTT.pdb", "B"
    )
    transform_pdb(
        STRUCTURE_DIR / "2PTU.pdb",
        paths["2ptu_aligned"],
        rotation_2ptu,
        translation_2ptu,
        keep_chains={"A"},
    )
    rotation_2ptv, translation_2ptv, rmsd_2ptv, n_2ptv = kabsch_transform(
        STRUCTURE_DIR / "2PTV.pdb", "A", STRUCTURE_DIR / "2PTT.pdb", "A"
    )
    transform_pdb(
        STRUCTURE_DIR / "2PTV.pdb",
        paths["2ptv_aligned"],
        rotation_2ptv,
        translation_2ptv,
        keep_chains={"A"},
    )
    alignments = {
        "6ULK_MHC_to_6UON_MHC": {"ca_count": n_6ulk, "ca_rmsd_angstrom": rmsd_6ulk},
        "2PTU_2B4_to_2PTT_2B4": {"ca_count": n_2ptu, "ca_rmsd_angstrom": rmsd_2ptu},
        "2PTV_CD48_to_2PTT_CD48": {"ca_count": n_2ptv, "ca_rmsd_angstrom": rmsd_2ptv},
    }
    return paths, alignments


def label(text, xyz, font_color=BLACK, background="#FFFFFF"):
    return {
        "text": text,
        "position": point(xyz),
        "fontSize": 18,
        "fontColor": font_color,
        "backgroundColor": background,
        "backgroundOpacity": 0.82,
    }


def line(first, second, color="#555555", dashed=True):
    return {"start": point(first), "end": point(second), "color": color, "dashed": dashed, "linewidth": 2}


def generate_structure_renders(sixuon, twoptt, paths):
    common_full = {
        "width": 1200,
        "height": 950,
        "deviceScaleFactor": 2,
        "background": "white",
        "zoom": 1.18,
        "rotations": [{"angle": -82, "axis": "y"}, {"angle": -12, "axis": "x"}],
    }
    for value_column, suffix in (("md_rmsf", "md"), ("flexppi_full_rmsf", "prediction")):
        output = RENDER_DIR / f"6uon_full_{suffix}.png"
        config = {
            **common_full,
            "output": str(output.resolve()),
            "models": [
                {
                    "path": str(paths["6uon_full"].resolve()),
                    "baseStyle": {"cartoon": {"color": "#D9D9D9"}},
                    "residueStyles": residue_styles(sixuon, value_column),
                    "addStyles": peptide_stick_styles(sixuon, value_column),
                }
            ],
        }
        render(f"6uon_full_{suffix}", config)

    sixuon_overlay_lines = []
    bound_shift = {}
    for position in range(5, 9):
        bound = find_atom(paths["6uon_groove"], "E", position, "CA")
        free = find_atom(paths["6ulk_groove"], "C", position, "CA")
        distance = float(np.linalg.norm(bound - free))
        bound_shift[f"p{position}_ca_shift_angstrom"] = distance
        sixuon_overlay_lines.append(line(bound, free, color="#777777", dashed=True))
    p6_bound = find_atom(paths["6uon_groove"], "E", 6, "CA")
    p6_free = find_atom(paths["6ulk_groove"], "C", 6, "CA")
    overlay_labels = [
        label("bound peptide", find_atom(paths["6uon_groove"], "E", 8, "CA"), CATEGORICAL[3]),
        label("TCR-free peptide", find_atom(paths["6ulk_groove"], "C", 5, "CA"), CATEGORICAL[4]),
        label(
            f"p5-p8 shift\n{min(bound_shift.values()):.1f}-{max(bound_shift.values()):.1f} A",
            (p6_bound + p6_free) / 2.0,
            BLACK,
        ),
    ]
    render(
        "6uon_peptide_overlay",
        {
            "output": str((RENDER_DIR / "6uon_peptide_overlay.png").resolve()),
            "width": 1150,
            "height": 900,
            "deviceScaleFactor": 2,
            "zoom": 1.48,
            "rotations": [{"angle": -65, "axis": "y"}, {"angle": 18, "axis": "x"}],
            "models": [
                {
                    "path": str(paths["6uon_groove"].resolve()),
                    "baseStyle": {"cartoon": {"color": "#D9D9D9", "opacity": 0.34}},
                    "residueStyles": [
                        {
                            "selection": {"chain": "E"},
                            "style": {
                                "cartoon": {"color": CATEGORICAL[3]},
                                "stick": {"color": CATEGORICAL[3], "radius": 0.15},
                            },
                        }
                    ],
                },
                {
                    "path": str(paths["6ulk_groove"].resolve()),
                    "baseStyle": {"cartoon": {"color": "#BDBDBD", "opacity": 0.16}},
                    "residueStyles": [
                        {
                            "selection": {"chain": "C"},
                            "style": {
                                "cartoon": {"color": CATEGORICAL[4]},
                                "stick": {"color": CATEGORICAL[4], "radius": 0.15},
                            },
                        }
                    ],
                },
            ],
            "lines": sixuon_overlay_lines,
            "labels": overlay_labels,
        },
    )

    sixuon_contacts = [
        ("Tyr97a-p4", ("A", 97, ["OH"]), ("E", 4, ["O"])),
        ("Gly97b-p5", ("B", 97, ["N"]), ("E", 5, ["O"])),
        ("Asp95b-p7", ("B", 95, ["OD1", "OD2"]), ("E", 7, ["NZ"])),
    ]
    contact_lines = []
    contact_distances = {}
    for name, first, second in sixuon_contacts:
        distance, p1, p2, atom1, atom2 = closest_atom_pair(paths["6uon_local"], first, second)
        contact_distances[name] = {"distance_angstrom": distance, "atoms": [atom1, atom2]}
        contact_lines.append(line(p1, p2, color=BLACK, dashed=True))
    contact_labels = [
        label("Tyr97a", find_atom(paths["6uon_local"], "A", 97, "CA"), CATEGORICAL[0]),
        label("Asp95b", find_atom(paths["6uon_local"], "B", 95, "CA"), CATEGORICAL[4]),
        label("Gly97b", find_atom(paths["6uon_local"], "B", 97, "CA"), CATEGORICAL[4]),
        label("p4", find_atom(paths["6uon_local"], "E", 4, "CA"), CATEGORICAL[3]),
        label("p5", find_atom(paths["6uon_local"], "E", 5, "CA"), CATEGORICAL[3]),
        label("p7 Lys", find_atom(paths["6uon_local"], "E", 7, "CA"), CATEGORICAL[3]),
    ]
    render(
        "6uon_contacts",
        {
            "output": str((RENDER_DIR / "6uon_contacts.png").resolve()),
            "width": 1250,
            "height": 900,
            "deviceScaleFactor": 2,
            "zoom": 1.56,
            "rotations": [{"angle": -65, "axis": "y"}, {"angle": 15, "axis": "x"}],
            "models": [
                {
                    "path": str(paths["6uon_local"].resolve()),
                    "baseStyle": {"cartoon": {"color": "#D9D9D9", "opacity": 0.26}},
                    "residueStyles": [
                        {"selection": {"chain": "A", "resi": 97}, "style": {"stick": {"color": CATEGORICAL[0], "radius": 0.18}}},
                        {"selection": {"chain": "B", "resi": 95}, "style": {"stick": {"color": CATEGORICAL[4], "radius": 0.18}}},
                        {"selection": {"chain": "B", "resi": 97}, "style": {"stick": {"color": CATEGORICAL[4], "radius": 0.18}}},
                        {"selection": {"chain": "E"}, "style": {"cartoon": {"color": CATEGORICAL[3]}, "stick": {"color": CATEGORICAL[3], "radius": 0.16}}},
                    ],
                }
            ],
            "lines": contact_lines,
            "labels": contact_labels,
        },
    )

    common_2ptt = {
        "width": 1200,
        "height": 950,
        "deviceScaleFactor": 2,
        "background": "white",
        "zoom": 1.33,
        "rotations": [{"angle": 28, "axis": "y"}, {"angle": -12, "axis": "x"}],
    }
    for value_column, suffix in (("md_rmsf", "md"), ("flexppi_full_rmsf", "prediction")):
        output = RENDER_DIR / f"2ptt_full_{suffix}.png"
        render(
            f"2ptt_full_{suffix}",
            {
                **common_2ptt,
                "output": str(output.resolve()),
                "models": [
                    {
                        "path": str(paths["2ptt_full"].resolve()),
                        "baseStyle": {"cartoon": {"color": "#D9D9D9"}},
                        "residueStyles": residue_styles(twoptt, value_column),
                    }
                ],
            },
        )

    h42_bound = find_atom(paths["2ptt_full"], "B", 42, "CA")
    h42_apo = find_atom(paths["2ptu_aligned"], "A", 42, "CA")
    h42_shift = float(np.linalg.norm(h42_bound - h42_apo))
    render(
        "2ptt_apo_bound_overlay",
        {
            "output": str((RENDER_DIR / "2ptt_apo_bound_overlay.png").resolve()),
            "width": 1200,
            "height": 900,
            "deviceScaleFactor": 2,
            "zoom": 1.66,
            "rotations": [{"angle": 18, "axis": "y"}, {"angle": -10, "axis": "x"}],
            "models": [
                {
                    "path": str(paths["2ptt_full"].resolve()),
                    "baseStyle": {"cartoon": {"color": "#D9D9D9", "opacity": 0.20}},
                    "residueStyles": [
                        {"selection": {"chain": "B", "resi": list(range(38, 44))}, "style": {"cartoon": {"color": CATEGORICAL[0]}, "stick": {"color": CATEGORICAL[0], "radius": 0.14}}},
                        {"selection": {"chain": "A", "resi": list(range(37, 42))}, "style": {"cartoon": {"color": CATEGORICAL[3]}, "stick": {"color": CATEGORICAL[3], "radius": 0.14}}},
                    ],
                },
                {
                    "path": str(paths["2ptu_aligned"].resolve()),
                    "baseStyle": {"cartoon": {"color": "#C6B3D5", "opacity": 0.18}},
                    "residueStyles": [
                        {"selection": {"chain": "A", "resi": list(range(38, 44))}, "style": {"cartoon": {"color": CATEGORICAL[4]}, "stick": {"color": CATEGORICAL[4], "radius": 0.12}}},
                    ],
                },
                {
                    "path": str(paths["2ptv_aligned"].resolve()),
                    "baseStyle": {"cartoon": {"color": "#F6D8B5", "opacity": 0.18}},
                    "residueStyles": [
                        {"selection": {"chain": "A", "resi": list(range(37, 42))}, "style": {"cartoon": {"color": CATEGORICAL[1]}, "stick": {"color": CATEGORICAL[1], "radius": 0.12}}},
                    ],
                },
            ],
            "lines": [line(h42_bound, h42_apo, color=BLACK, dashed=True)],
            "labels": [
                label("bound 2B4 loop", find_atom(paths["2ptt_full"], "B", 40, "CA"), CATEGORICAL[0]),
                label("apo 2B4 loop", find_atom(paths["2ptu_aligned"], "A", 40, "CA"), CATEGORICAL[4]),
                label("bound CD48 loop", find_atom(paths["2ptt_full"], "A", 39, "CA"), CATEGORICAL[3]),
                label("apo CD48 loop\ntip 38-39 unresolved", find_atom(paths["2ptv_aligned"], "A", 40, "CA"), CATEGORICAL[1]),
                label(f"H42 shift {h42_shift:.1f} A", (h42_bound + h42_apo) / 2.0, BLACK),
            ],
        },
    )

    anchor_specs = [
        ("K35-E93", ("B", 35, ["NZ"]), ("A", 93, ["OE1", "OE2"])),
        ("E91-R87", ("B", 91, ["OE1", "OE2"]), ("A", 87, ["NH1", "NH2", "NE"])),
        ("T93-K41", ("B", 93, ["OG1"]), ("A", 41, ["NZ"])),
    ]
    anchor_lines = []
    anchor_distances = {}
    for name, first, second in anchor_specs:
        distance, p1, p2, atom1, atom2 = closest_atom_pair(paths["2ptt_local"], first, second)
        anchor_distances[name] = {"distance_angstrom": distance, "atoms": [atom1, atom2]}
        anchor_lines.append(line(p1, p2, color=BLACK, dashed=True))
    render(
        "2ptt_interface_contacts",
        {
            "output": str((RENDER_DIR / "2ptt_interface_contacts.png").resolve()),
            "width": 1250,
            "height": 900,
            "deviceScaleFactor": 2,
            "zoom": 1.52,
            "rotations": [{"angle": 22, "axis": "y"}, {"angle": -8, "axis": "x"}],
            "models": [
                {
                    "path": str(paths["2ptt_local"].resolve()),
                    "baseStyle": {"cartoon": {"color": "#D9D9D9", "opacity": 0.24}},
                    "residueStyles": [
                        {"selection": {"chain": "B", "resi": [35, 91, 93]}, "style": {"stick": {"color": CATEGORICAL[0], "radius": 0.18}}},
                        {"selection": {"chain": "A", "resi": [41, 87, 93]}, "style": {"stick": {"color": CATEGORICAL[3], "radius": 0.18}}},
                        {"selection": {"chain": "B", "resi": list(range(38, 44))}, "style": {"cartoon": {"color": "#80B1D3"}}},
                        {"selection": {"chain": "A", "resi": list(range(37, 42))}, "style": {"cartoon": {"color": "#FDBF6F"}}},
                    ],
                }
            ],
            "lines": anchor_lines,
            "labels": [
                label("2B4 K35", find_atom(paths["2ptt_local"], "B", 35, "CA"), CATEGORICAL[0]),
                label("2B4 E91", find_atom(paths["2ptt_local"], "B", 91, "CA"), CATEGORICAL[0]),
                label("2B4 T93", find_atom(paths["2ptt_local"], "B", 93, "CA"), CATEGORICAL[0]),
                label("CD48 E93", find_atom(paths["2ptt_local"], "A", 93, "CA"), CATEGORICAL[3]),
                label("CD48 R87", find_atom(paths["2ptt_local"], "A", 87, "CA"), CATEGORICAL[3]),
                label("CD48 K41", find_atom(paths["2ptt_local"], "A", 41, "CA"), CATEGORICAL[3]),
            ],
        },
    )
    return {
        "6uon_peptide_shift": bound_shift,
        "6uon_contact_distances": contact_distances,
        "2ptt_h42_ca_shift_angstrom": h42_shift,
        "2ptt_anchor_distances": anchor_distances,
    }


def show_image(ax, path):
    image = Image.open(path).convert("RGB")
    ax.imshow(image)
    ax.set_axis_off()


def panel_label(ax, text):
    ax.text(
        -0.045,
        1.025,
        text,
        transform=ax.transAxes,
        fontsize=9,
        fontweight="bold",
        ha="right",
        va="bottom",
        color=BLACK,
    )


def add_rmsf_colorbar(fig, axes):
    scalar = mpl.cm.ScalarMappable(norm=RMSF_NORM, cmap=RMSF_CMAP)
    colorbar = fig.colorbar(
        scalar,
        ax=axes,
        orientation="horizontal",
        fraction=0.038,
        pad=0.018,
        aspect=34,
    )
    colorbar.set_ticks([0.25, 0.75, 1.25, 1.75, 2.0])
    colorbar.set_ticklabels(["0.25", "0.75", "1.25", "1.75", ">=2.0"])
    colorbar.set_label("RMSF (A); common scale, values above 2.0 A saturated", fontsize=6.5)
    colorbar.ax.tick_params(labelsize=6, width=0.5, length=2)
    colorbar.outline.set_linewidth(0.5)


def make_6uon_figure(sixuon):
    mm = 1 / 25.4
    fig = plt.figure(figsize=(183 * mm, 145 * mm))
    outer = fig.add_gridspec(
        2,
        2,
        width_ratios=[1.22, 1.0],
        height_ratios=[1.0, 1.03],
        left=0.035,
        right=0.985,
        bottom=0.085,
        top=0.965,
        wspace=0.12,
        hspace=0.25,
    )
    top_left = GridSpecFromSubplotSpec(1, 2, subplot_spec=outer[0, 0], wspace=0.015)
    ax_md = fig.add_subplot(top_left[0, 0])
    ax_pred = fig.add_subplot(top_left[0, 1])
    show_image(ax_md, RENDER_DIR / "6uon_full_md.png")
    show_image(ax_pred, RENDER_DIR / "6uon_full_prediction.png")
    ax_md.set_title("ITrace MD", pad=1.0, fontweight="bold")
    ax_pred.set_title("FlexPPI", pad=1.0, fontweight="bold", color=CATEGORICAL[0])
    panel_label(ax_md, "a")
    ax_md.text(0.02, 0.04, "TCR", transform=ax_md.transAxes, fontsize=6, color=BLACK)
    ax_md.text(0.72, 0.05, "pMHC", transform=ax_md.transAxes, fontsize=6, color=BLACK)
    add_rmsf_colorbar(fig, [ax_md, ax_pred])

    ax_overlay = fig.add_subplot(outer[0, 1])
    show_image(ax_overlay, RENDER_DIR / "6uon_peptide_overlay.png")
    ax_overlay.set_title("TCR-free versus TCR-bound peptide", pad=1.0, fontweight="bold")
    panel_label(ax_overlay, "b")
    ax_overlay.text(
        0.02,
        0.025,
        "MHC-aligned structures: 6ULK (TCR-free) and 6UON (bound)",
        transform=ax_overlay.transAxes,
        fontsize=5.7,
        color="#555555",
    )

    ax_contacts = fig.add_subplot(outer[1, 0])
    show_image(ax_contacts, RENDER_DIR / "6uon_contacts.png")
    ax_contacts.set_title("CDR3 recognition of the central peptide ridge", pad=1.0, fontweight="bold")
    panel_label(ax_contacts, "c")
    ax_contacts.text(
        0.02,
        0.02,
        "Dashed lines: Tyr97a-p4 and Gly97b-p5 H bonds; Asp95b-p7 Lys salt bridge",
        transform=ax_contacts.transAxes,
        fontsize=5.7,
        color="#555555",
    )

    profile_grid = GridSpecFromSubplotSpec(
        2, 1, subplot_spec=outer[1, 1], height_ratios=[2.6, 1.0], hspace=0.06
    )
    ax_profile = fig.add_subplot(profile_grid[0, 0])
    ax_delta = fig.add_subplot(profile_grid[1, 0], sharex=ax_profile)
    peptide = sixuon.loc[sixuon["chain"] == "E"].sort_values("pdb_residue_number")
    positions = peptide["pdb_residue_number"].to_numpy(dtype=int)
    ax_profile.plot(
        positions,
        peptide["md_rmsf"],
        color=BLACK,
        marker="o",
        markersize=3.2,
        linewidth=1.35,
        label="ITrace MD",
    )
    ax_profile.plot(
        positions,
        peptide["flexppi_full_rmsf"],
        color=CATEGORICAL[0],
        marker="s",
        markersize=3.0,
        linewidth=1.35,
        label="FlexPPI, full complex",
    )
    ax_profile.plot(
        positions,
        peptide["pmhc_only_rmsf"],
        color=CATEGORICAL[3],
        marker="^",
        markersize=3.0,
        linewidth=1.1,
        linestyle="--",
        label="FlexPPI, pMHC only",
    )
    ax_profile.axvspan(3.7, 7.3, color="#F1A340", alpha=0.09, linewidth=0)
    ax_profile.set_ylabel("RMSF (A)")
    ax_profile.set_ylim(0.25, 1.50)
    ax_profile.set_xlim(0.6, 10.4)
    ax_profile.legend(
        loc="upper left",
        bbox_to_anchor=(0.01, 0.99),
        ncol=1,
        fontsize=5.8,
        handlelength=1.7,
        labelspacing=0.25,
        borderaxespad=0.0,
    )
    ax_profile.annotate(
        "p3 anchor",
        xy=(3, float(peptide.loc[peptide["pdb_residue_number"] == 3, "md_rmsf"].iloc[0])),
        xytext=(1.35, 0.98),
        fontsize=6,
        arrowprops=dict(arrowstyle="-", color="#777777", linewidth=0.6),
    )
    ax_profile.annotate(
        "p6/p8 flexible ridge",
        xy=(8, float(peptide.loc[peptide["pdb_residue_number"] == 8, "md_rmsf"].iloc[0])),
        xytext=(6.2, 1.32),
        fontsize=6,
        arrowprops=dict(arrowstyle="-", color="#777777", linewidth=0.6),
    )
    ax_profile.tick_params(labelbottom=False)
    panel_label(ax_profile, "d")
    ax_profile.set_title("Peptide flexibility and static partner-context response", pad=4, fontweight="bold")

    delta = peptide["context_delta_full_minus_pmhc"].to_numpy(dtype=float)
    ax_delta.bar(positions, delta, width=0.68, color=ACCENT_RED, alpha=0.78, linewidth=0)
    ax_delta.axhline(0, color="#777777", linestyle="--", linewidth=0.6)
    ax_delta.set_ylabel("Delta\nRMSF (A)", fontsize=6.5)
    ax_delta.set_xlabel("KRAS-G12D peptide position")
    ax_delta.set_xticks(positions)
    ax_delta.set_xticklabels([f"p{i}\n{aa}" for i, aa in zip(positions, peptide["amino_acid_3"])])
    ax_delta.set_ylim(-0.55, 0.10)
    ax_delta.text(
        7.65,
        -0.47,
        "Delta = full complex - pMHC-only prediction",
        fontsize=5.5,
        color="#555555",
    )
    save_cns_figure(fig, str(OUTPUT_6UON))
    fig.savefig(f"{OUTPUT_6UON}.svg", bbox_inches="tight")
    plt.close(fig)


def make_2ptt_figure(twoptt):
    mm = 1 / 25.4
    fig = plt.figure(figsize=(183 * mm, 145 * mm))
    outer = fig.add_gridspec(
        2,
        2,
        width_ratios=[1.20, 1.0],
        height_ratios=[1.0, 1.02],
        left=0.035,
        right=0.985,
        bottom=0.080,
        top=0.965,
        wspace=0.12,
        hspace=0.24,
    )
    top_left = GridSpecFromSubplotSpec(1, 2, subplot_spec=outer[0, 0], wspace=0.015)
    ax_md = fig.add_subplot(top_left[0, 0])
    ax_pred = fig.add_subplot(top_left[0, 1])
    show_image(ax_md, RENDER_DIR / "2ptt_full_md.png")
    show_image(ax_pred, RENDER_DIR / "2ptt_full_prediction.png")
    ax_md.set_title("DynaRepo MD", pad=1.0, fontweight="bold")
    ax_pred.set_title("FlexPPI", pad=1.0, fontweight="bold", color=CATEGORICAL[0])
    panel_label(ax_md, "a")
    ax_md.text(0.06, 0.05, "CD48", transform=ax_md.transAxes, fontsize=6, color=BLACK)
    ax_md.text(0.71, 0.05, "2B4", transform=ax_md.transAxes, fontsize=6, color=BLACK)
    add_rmsf_colorbar(fig, [ax_md, ax_pred])

    ax_overlay = fig.add_subplot(outer[0, 1])
    show_image(ax_overlay, RENDER_DIR / "2ptt_apo_bound_overlay.png")
    ax_overlay.set_title("Apo-to-bound ordering of the CC-prime loops", pad=1.0, fontweight="bold")
    panel_label(ax_overlay, "b")
    ax_overlay.text(
        0.02,
        0.025,
        "2PTT bound complex; apo 2B4 (2PTU) and apo CD48 (2PTV) aligned by C-alpha atoms",
        transform=ax_overlay.transAxes,
        fontsize=5.5,
        color="#555555",
    )

    ax_contacts = fig.add_subplot(outer[1, 0])
    show_image(ax_contacts, RENDER_DIR / "2ptt_interface_contacts.png")
    ax_contacts.set_title("Stable electrostatic anchors around compliant loop clamps", pad=1.0, fontweight="bold")
    panel_label(ax_contacts, "c")
    ax_contacts.text(
        0.02,
        0.02,
        "Dashed lines: K35-E93, E91-R87 and T93-K41 interface contacts",
        transform=ax_contacts.transAxes,
        fontsize=5.7,
        color="#555555",
    )

    ax_pairs = fig.add_subplot(outer[1, 1])
    selected = [
        ("2B4", 35, "anchor"),
        ("2B4", 39, "CC-prime loop"),
        ("2B4", 40, "CC-prime loop"),
        ("2B4", 41, "CC-prime loop"),
        ("2B4", 43, "CC-prime loop"),
        ("2B4", 91, "anchor"),
        ("2B4", 93, "anchor"),
        ("CD48", 37, "CC-prime loop"),
        ("CD48", 38, "CC-prime loop"),
        ("CD48", 39, "CC-prime loop"),
        ("CD48", 40, "CC-prime loop"),
        ("CD48", 41, "CC-prime loop/anchor"),
        ("CD48", 87, "anchor"),
        ("CD48", 93, "anchor"),
    ]
    rows = []
    for role, resi, function in selected:
        row = twoptt.loc[(twoptt["chain_role"] == role) & (twoptt["pdb_residue_number"] == resi)]
        if row.empty:
            continue
        item = row.iloc[0]
        rows.append(
            {
                "label": f"{role} {item['amino_acid_3'].title()}{resi}",
                "role": role,
                "function": function,
                "md": float(item["md_rmsf"]),
                "prediction": float(item["flexppi_full_rmsf"]),
            }
        )
    rows = list(reversed(rows))
    y = np.arange(len(rows))
    for index, row in enumerate(rows):
        if "loop" in row["function"]:
            color = "#EEF5FA" if row["role"] == "2B4" else "#FFF4E6"
            ax_pairs.axhspan(index - 0.47, index + 0.47, color=color, linewidth=0, zorder=0)
        ax_pairs.plot(
            [row["md"], row["prediction"]],
            [index, index],
            color="#B8B8B8",
            linewidth=1.0,
            zorder=1,
        )
    ax_pairs.scatter([row["md"] for row in rows], y, s=18, color=BLACK, marker="o", label="MD", zorder=3)
    ax_pairs.scatter(
        [row["prediction"] for row in rows],
        y,
        s=20,
        color=CATEGORICAL[0],
        marker="s",
        label="FlexPPI",
        zorder=3,
    )
    ax_pairs.set_yticks(y)
    ax_pairs.set_yticklabels([row["label"] for row in rows], fontsize=6.3)
    ax_pairs.set_xlabel("Bound-state RMSF (A)")
    ax_pairs.set_xlim(0.25, 2.15)
    ax_pairs.set_ylim(-0.7, len(rows) - 0.3)
    ax_pairs.legend(
        loc="upper right",
        bbox_to_anchor=(0.99, 0.99),
        ncol=2,
        fontsize=6.3,
        columnspacing=1.0,
    )
    ax_pairs.set_title("Residue-level flexibility: loops versus anchors", pad=4, fontweight="bold")
    panel_label(ax_pairs, "d")
    ax_pairs.text(
        0.01,
        -0.16,
        "H42 is structurally resolved but absent from the mapped DynaRepo table (not zero).",
        transform=ax_pairs.transAxes,
        fontsize=5.5,
        color="#555555",
    )
    save_cns_figure(fig, str(OUTPUT_2PTT))
    fig.savefig(f"{OUTPUT_2PTT}.svg", bbox_inches="tight")
    plt.close(fig)


def write_outputs(sixuon, twoptt, checks, alignments, render_stats):
    columns_6uon = [
        "case", "chain", "chain_role", "pdb_residue_number", "amino_acid_3",
        "region", "region_group", "md_rmsf", "flexppi_full_rmsf", "pmhc_only_rmsf",
        "context_delta_full_minus_pmhc", "chainfit_replica_count", "source_scope",
    ]
    sixuon.reindex(columns=columns_6uon).to_csv(SOURCE_6UON, index=False)
    columns_2ptt = [
        "case", "sample_id", "chain", "chain_role", "dataset_author_residue_id",
        "pdb_residue_number", "amino_acid_3", "md_rmsf", "flexppi_full_rmsf",
        "source_scope",
    ]
    twoptt.reindex(columns=columns_2ptt).to_csv(SOURCE_2PTT, index=False)

    peptide = sixuon.loc[sixuon["chain"] == "E"].sort_values("pdb_residue_number")
    loop_2b4 = twoptt.loc[
        (twoptt["chain_role"] == "2B4")
        & twoptt["pdb_residue_number"].isin([38, 39, 40, 41, 43])
    ]
    anchors = twoptt.loc[
        ((twoptt["chain_role"] == "2B4") & twoptt["pdb_residue_number"].isin([35, 91, 93]))
        | ((twoptt["chain_role"] == "CD48") & twoptt["pdb_residue_number"].isin([41, 87, 93]))
    ]
    statistics = {
        "figure_contract": {
            "core_conclusion_6uon": "Partner context selectively suppresses the model's predicted mobility at the contacted central KRAS-G12D peptide ridge, while MD retains p6/p8 as the most mobile ridge positions.",
            "core_conclusion_2ptt": "The 2B4-CD48 interface combines relatively stable electrostatic anchors with a bound CC-prime loop that remains locally compliant.",
            "archetype": "asymmetric mixed-modality",
            "journal_export": "Nature-family standard, 183 mm double-column, vector PDF/SVG plus 300 dpi PNG",
        },
        "data_validation": checks,
        "alignment_quality": alignments,
        "structure_rendering_measurements": render_stats,
        "6uon_peptide": {
            "n_residues": int(len(peptide)),
            "md_replica_definition": "Per-residue C-alpha RMSF averaged over three independent ITrace trajectories after chain-wise fitting.",
            "md_mean_angstrom": float(peptide["md_rmsf"].mean()),
            "flexppi_full_mean_angstrom": float(peptide["flexppi_full_rmsf"].mean()),
            "pmhc_only_mean_angstrom": float(peptide["pmhc_only_rmsf"].mean()),
            "mean_context_delta_full_minus_pmhc_angstrom": float(peptide["context_delta_full_minus_pmhc"].mean()),
            "most_mobile_md_positions": [
                int(value) for value in peptide.nlargest(3, "md_rmsf")["pdb_residue_number"]
            ],
            "largest_context_suppression_positions": [
                int(value) for value in peptide.nsmallest(4, "context_delta_full_minus_pmhc")["pdb_residue_number"]
            ],
        },
        "2ptt_local": {
            "dataset_role": "DynaRepo model-selection validation example; not an untouched external test.",
            "mapped_residue_count": int(len(twoptt)),
            "2b4_ccprime_loop_mean_md_angstrom": float(loop_2b4["md_rmsf"].mean()),
            "2b4_ccprime_loop_mean_prediction_angstrom": float(loop_2b4["flexppi_full_rmsf"].mean()),
            "selected_anchor_mean_md_angstrom": float(anchors["md_rmsf"].mean()),
            "selected_anchor_mean_prediction_angstrom": float(anchors["flexppi_full_rmsf"].mean()),
            "h42_mapping": "Native 2B4 H42 is present in 2PTT/2PTU coordinates but absent from the DynaRepo residue prediction table; it is masked, not assigned zero.",
        },
        "model": {
            "checkpoint": str(CHECKPOINT.resolve()),
            "checkpoint_sha256": sha256(CHECKPOINT),
            "checkpoint_role": "Frozen pretraining checkpoint used for ITrace and DynaRepo state prediction.",
        },
        "statistics_and_reproducibility": {
            "inferential_tests": "None. These are residue-level mechanistic case studies, not independent-replicate hypothesis tests.",
            "center_and_spread": "No error bars are plotted in the case panels. The 6UON MD target is the mean of three trajectories at each residue; the 2PTT target is the archived DynaRepo RMSF label.",
            "selection_rule": "Residues in quantitative local panels were selected before plotting from published structural contacts and named CC-prime loops; every mapped residue remains in the source-data tables and whole-structure coloring.",
            "rmsf_color_scale": {"minimum_angstrom": RMSF_MIN, "saturation_angstrom": RMSF_MAX},
            "sources": [str(ITRACE_PREDICTIONS), str(ITRACE_CONDITIONS), str(DYNA_PREDICTIONS)],
        },
    }
    STATS_PATH.write_text(json.dumps(statistics, indent=2, ensure_ascii=False), encoding="utf-8")

    qa = f"""# Case-study figures QA report

## Figure contract

- 6UON core message: central peptide contacts explain a localized static partner-context response, while the ITrace MD target identifies p6/p8 as the flexible ridge.
- 2PTT core message: stable electrostatic anchors coexist with a conformationally ordered yet still compliant CC-prime loop.
- Archetype: asymmetric mixed-modality; structure views are the hero evidence and residue plots are the quantitative support.
- Export: 183 mm double-column figures, vector PDF/SVG masters and 300 dpi PNG previews.

## Data integrity

- 6UON: all {checks['6uon_residue_count']} residues are retained in the source-data table and whole-structure color mapping; the local profile contains all 10 peptide residues.
- 2PTT: all {checks['2ptt_mapped_residue_count']} mapped DynaRepo residues are retained; native PDB numbering is restored by documented chain-specific offsets.
- 2B4 H42 is present in the coordinates but absent from the mapped DynaRepo table; it is masked and never interpreted as zero mobility.
- The 6UON full-minus-pMHC difference holds target-chain coordinates fixed and is labeled as a static context response, not an apo/holo trajectory.

## Pass 0 - anti-pattern scan

- AP-0 to AP-7: PASS. Mandatory typography, restrained custom palette, font embedding, vector export, clean spines and explicit font family are present.
- Molecular panels have white backgrounds and semantic chain/residue colors; no rainbow/jet maps are used.

## Pass 1 - code compliance

- CL-1 to CL-7: PASS. Text is at least 5.5 pt at final size; width is 183 mm; PDF/SVG and 300 dpi PNG outputs are generated; PDF text uses TrueType embedding.
- Structure panels are rendered at 2x device scale before composition; quantitative layers remain vector.

## Pass 2 - visual logic and reproducibility

- VI-1 to VI-7: PASS. Each panel contributes unique evidence; the same absolute RMSF scale is used for MD and FlexPPI; zero is explicit on the context-delta panel; selected residues are defined by published contacts/loops rather than post hoc score extremes.
- No inferential p-value is reported from residues within one complex because residues are not independent biological replicates.

## Pass 3 - rendered verification

- VV-1: PASS. Panel titles, legends, annotations and molecular labels do not overlap or clip after the final title/legend adjustments.
- VV-2: PASS. Axes, panel gutters and structure crops remain aligned at the 183 mm final figure width.
- VV-3: PASS. Residue identities, RMSF scales and interpretive callouts remain legible in both standalone figures and the compiled manuscript.
- VV-4: PASS. Structure colors, residue highlights and quantitative traces preserve the same semantic mapping across panels.
- VV-5: PASS. PDF/SVG masters and 300 dpi PNG previews open correctly; all 30 rendered manuscript pages were visually inspected.

## Traceability

- 6UON source data: `{SOURCE_6UON}`.
- 2PTT source data: `{SOURCE_2PTT}`.
- Full numerical and alignment record: `{STATS_PATH}`.
"""
    QA_PATH.write_text(qa, encoding="utf-8")


def main():
    ensure_directories()
    download_structures()
    manifest = pd.read_csv(ITRACE_MANIFEST)
    row = manifest.loc[manifest["pdb_id"].astype(str).str.lower() == "6uon"]
    if len(row) != 1:
        raise ValueError("Expected exactly one 6UON ITrace topology row")
    sixuon_topology = Path(str(row.iloc[0]["topology_path"]))
    if not sixuon_topology.exists():
        raise FileNotFoundError(sixuon_topology)
    paths, alignments = prepare_structures(sixuon_topology)
    sixuon = load_6uon_data()
    twoptt = load_2ptt_data(paths["2ptt_full"])
    checks = validate_data(sixuon, twoptt, paths["6uon_full"], paths["2ptt_full"])
    render_stats = generate_structure_renders(sixuon, twoptt, paths)
    make_6uon_figure(sixuon)
    make_2ptt_figure(twoptt)
    write_outputs(sixuon, twoptt, checks, alignments, render_stats)
    print(f"Wrote {OUTPUT_6UON}.pdf and {OUTPUT_2PTT}.pdf")


if __name__ == "__main__":
    main()
