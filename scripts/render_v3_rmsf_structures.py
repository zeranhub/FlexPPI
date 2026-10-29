"""Render authentic test structures colored by MD or FlexPPI RMSF.

Inputs are the selected_cases.json list and provenance-checked OOF NPZ files.
Residues are matched by chain AND author residue ID (including insertion code),
and the full cache order/sequence is checked against the actual input PDB.
Both panels use identical coordinates and camera; these are static structures,
not generated trajectories or an MD conformational ensemble.

Example::

    python scripts/render_v3_rmsf_structures.py --cases \
        artifacts/evaluation_v3/analysis/selected_cases.json --vmax 5

PyMOL runs in an isolated, headless child process. The default executable is
the existing local conda runtime; --pymol-python can override it.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

COLORS = ("#2D4EA1", "#F7F7F7", "#FF3026")
MISSING_COLOR = "#B8BDC5"


def file_hash(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_json(path: Path, payload: Any) -> None:
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2,
                               allow_nan=False) + "\n", encoding="utf-8")


def prepare_case(case: dict, out: Path, vmax: float, size: int) -> dict:
    from flexppi.data.structures import parse_pdb

    case_id = str(case["case_id"])
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", case_id):
        raise ValueError(f"case_id must be filename-safe: {case_id!r}")
    cache_path = Path(case["cache_path"]).resolve()
    with np.load(cache_path, allow_pickle=False) as archive:
        cache = {key: archive[key] for key in archive.files}
    metadata = json.loads(str(cache["metadata_json"].item()))
    if metadata["sample_id"] != case["sample_id"]:
        raise ValueError(f"Cache sample ID mismatch for {case_id}")
    pdb_path = Path(case.get("structure_path") or metadata["structure_path"]).resolve()
    structure = parse_pdb(str(pdb_path))
    cache_keys = list(zip(cache["chain"].astype(str),
                          cache["author_residue_id"].astype(str)))
    structure_keys = list(zip(structure.chain_ids, structure.author_residue_ids))
    if cache_keys != structure_keys:
        raise ValueError(f"Full PDB/cache residue order mismatch for {case_id}")
    if len(set(cache_keys)) != len(cache_keys):
        raise ValueError(f"Duplicate residue keys in {case_id}")
    if "amino_acid" in cache:
        sequence = "".join(cache["amino_acid"].astype(str).tolist())
        if sequence != structure.sequence:
            raise ValueError(f"PDB/cache amino-acid sequence mismatch for {case_id}")
    expected_hash = metadata.get("identity", {}).get("input_sha256")
    source_hash = file_hash(pdb_path)
    if expected_hash and expected_hash != source_hash:
        raise ValueError(f"Input structure changed since prediction: {case_id}")

    metric = cache["metric_names"].astype(str).tolist().index("rmsf")
    target = cache["target"][:, metric].astype(float)
    prediction = cache["mean"][:, metric].astype(float)
    valid = cache["valid"][:, metric].astype(bool)
    valid &= np.isfinite(target) & np.isfinite(prediction)
    chain = str(case["chain"])
    selected = cache["chain"].astype(str) == chain
    if not selected.any():
        raise ValueError(f"No chain {chain!r} in {case_id}")
    if np.any(target[selected & valid] < 0) or np.any(prediction[selected & valid] < 0):
        raise ValueError(f"Negative RMSF in {case_id}")
    keys = {cache_keys[index]: index for index in np.flatnonzero(selected)}
    selected_valid = selected & valid
    if not selected_valid.any():
        raise ValueError(f"No valid RMSF residues for {case_id}")

    # Keep actual selected-chain atoms from the first PDB model. Altloc handling
    # matches parse_pdb, and residue IDs are never renumbered or offset guessed.
    atom_lines = []
    seen_atoms = set()
    saw_model = False
    with pdb_path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            record = line[:6].strip()
            if record == "MODEL":
                if saw_model:
                    break
                saw_model = True
                continue
            if record == "ENDMDL":
                break
            if record not in {"ATOM", "HETATM"} or len(line) < 54:
                continue
            if line[16:17] not in {" ", "A", "1"}:
                continue
            key = (line[21:22].strip() or "_",
                   line[22:26].strip() + line[26:27].strip())
            if key in keys:
                atom_key = (*key, line[12:16].strip())
                if atom_key in seen_atoms:
                    continue
                if line[17:20].strip().upper() != structure.residue_names[keys[key]]:
                    continue
                seen_atoms.add(atom_key)
                # PyMOL must not choose a different alternate atom conformer.
                clean = line.rstrip("\r\n").ljust(80)
                atom_lines.append((clean[:16] + " " + clean[17:], keys[key]))
    if not atom_lines:
        raise ValueError(f"No selected-chain atoms in {case_id}")

    assets = {}
    for role, values in (("md", target), ("prediction", prediction)):
        annotated = out / f"{case_id}_{role}.pdb"
        with annotated.open("w", encoding="ascii", errors="replace", newline="\n") as handle:
            handle.write("REMARK 900 STATIC INPUT COORDINATES; NO CONFORMATIONAL ENSEMBLE\n")
            handle.write(f"REMARK 900 B-FACTOR = {role.upper()} RMSF IN ANGSTROM\n")
            handle.write("REMARK 900 NEGATIVE B-FACTOR (-1) = UNMAPPED/MASKED RESIDUE\n")
            for line, index in atom_lines:
                value = float(values[index]) if valid[index] else -1.0
                if value > 999.99:
                    raise ValueError("RMSF exceeds PDB B-factor field capacity")
                handle.write(line[:60] + f"{value:6.2f}" + line[66:] + "\n")
            handle.write("TER\nEND\n")
        assets[role] = {"pdb": str(annotated),
                        "png": str(out / f"{case_id}_{role}.png")}

    csv_path = out / f"{case_id}_mapped_residues.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["sample_id", "chain", "author_residue_id", "amino_acid",
                         "cache_residue_index", "valid", "md_rmsf_angstrom",
                         "flexppi_rmsf_angstrom"])
        for index in np.flatnonzero(selected):
            writer.writerow([case["sample_id"], *cache_keys[index],
                             structure.sequence[index], int(index), bool(valid[index]),
                             float(target[index]) if valid[index] else "",
                             float(prediction[index]) if valid[index] else ""])

    return {"case_id": case_id, "sample_id": case["sample_id"], "chain": chain,
            "selection_group": case.get("selection_group"), "fold": metadata.get("fold"),
            "cache_path": str(cache_path), "cache_sha256": file_hash(cache_path),
            "source_structure": str(pdb_path), "source_structure_sha256": source_hash,
            "mapping": "exact chain and author residue ID, including insertion code",
            "selected_residues": int(selected.sum()), "valid_residues": int(selected_valid.sum()),
            "masked_residues": int((selected & ~valid).sum()),
            "md_above_color_limit": int((selected_valid & (target > vmax)).sum()),
            "prediction_above_color_limit": int((selected_valid & (prediction > vmax)).sum()),
            "coordinate_provenance": "first model of actual inference input PDB, selected chain",
            "ensemble": False, "numeric_scale_angstrom": [0, vmax],
            "color_stops": list(COLORS), "missing_color": MISSING_COLOR,
            "width": size, "height": size, "assets": assets,
            "mapped_residue_csv": str(csv_path),
            "pymol_script": str(out / f"{case_id}.pml"),
            "pymol_session": str(out / f"{case_id}.pse")}


def worker(payload_path: Path) -> None:
    """Separate process so the caller need not import PyMOL or its DLLs."""
    import pymol2
    from PIL import Image

    payload = json.loads(payload_path.read_text(encoding="utf-8"))
    settings = [("ray_opaque_background", 1), ("ray_trace_mode", 1),
                ("ray_trace_gain", 0.08), ("ray_shadows", 0),
                ("antialias", 2), ("orthoscopic", 1), ("depth_cue", 0),
                ("specular", 0.15), ("ambient", 0.55), ("direct", 0.45),
                ("cartoon_fancy_helices", 1), ("cartoon_fancy_sheets", 1),
                ("cartoon_flat_sheets", 1), ("cartoon_smooth_loops", 1),
                ("cartoon_sampling", 14), ("max_threads", 4)]
    color_hex = [value.lstrip("#") for value in COLORS]
    color_rgb = [tuple(int(value[index:index+2], 16) / 255
                       for index in (0, 2, 4)) for value in color_hex]
    with pymol2.PyMOL() as instance:
        cmd = instance.cmd
        cmd.reinitialize()
        cmd.bg_color("white")
        for setting, value in settings:
            cmd.set(setting, value)
        for index, value in enumerate(color_rgb):
            cmd.set_color(f"rmsf{index}", value)
        cmd.set_color("missing_rmsf", tuple(int(MISSING_COLOR[i:i+2], 16) / 255
                                          for i in (1, 3, 5)))
        cmd.viewport(payload["width"], payload["height"])
        lines = ["reinitialize", "bg_color white"]
        lines += [f"set {key}, {value}" for key, value in settings]
        for index, rgb in enumerate(color_rgb):
            lines.append(f"set_color rmsf{index}, [{','.join(map(str, rgb))}]")
        lines.append(f"set_color missing_rmsf, [0.721569,0.741176,0.772549]")
        camera = None
        for role in ("md", "prediction"):
            assets = payload["assets"][role]
            name = "MD" if role == "md" else "FlexPPI"
            cmd.load(assets["pdb"], name)
            cmd.hide("everything", name)
            cmd.show("cartoon", name)
            cmd.dss(name)
            cmd.spectrum("b", "rmsf0 rmsf1 rmsf2", name,
                         minimum=0, maximum=payload["numeric_scale_angstrom"][1], byres=1)
            cmd.color("missing_rmsf", f"{name} and b < 0")
            cmd.disable("all")
            cmd.enable(name)
            if camera is None:
                cmd.orient(name)
                cmd.zoom(name, buffer=2.5)
                camera = list(cmd.get_view())
            else:
                cmd.set_view(camera)
            cmd.png(assets["png"], width=payload["width"], height=payload["height"],
                    dpi=300, ray=1, quiet=1)
            if role == 'md':
                # Cartoon widths can extend past the atom-based zoom bounds.
                # Widen the camera only when the actual rendered image touches
                # an edge, then reuse that exact camera for the prediction.
                for buffer in (5.0, 8.0, 12.0):
                    pixels = np.asarray(Image.open(assets['png']).convert('RGB'))
                    yy, xx = np.where(np.any(pixels < 230, axis=2))
                    if len(xx) and min(xx.min(), yy.min(), pixels.shape[1]-1-xx.max(),
                                       pixels.shape[0]-1-yy.max()) >= 8:
                        break
                    cmd.zoom(name, buffer=buffer)
                    camera = list(cmd.get_view())
                    cmd.png(assets['png'], width=payload['width'], height=payload['height'],
                            dpi=300, ray=1, quiet=1)
                    payload['expanded_camera_buffer_angstrom'] = buffer
            # This script reproduces the actual shared camera, including its
            # orthographic clipping plane, without reorienting either panel.
            lines += [f'load "{Path(assets["pdb"]).as_posix()}", {name}',
                      f"hide everything, {name}", f"show cartoon, {name}",
                      f"dss {name}",
                      f"spectrum b, rmsf0 rmsf1 rmsf2, {name}, minimum=0, "
                      f"maximum={payload['numeric_scale_angstrom'][1]}, byres=1",
                      f"color missing_rmsf, {name} and b < 0", "disable all",
                      f"enable {name}", f"set_view ({','.join(map(str, camera))})",
                      f"ray {payload['width']}, {payload['height']}",
                      f'png "{Path(assets["png"]).as_posix()}", dpi=300']
        cmd.save(payload["pymol_session"])
        lines.append(f'save "{Path(payload["pymol_session"]).as_posix()}"')
        Path(payload["pymol_script"]).write_text("\n".join(lines) + "\n", encoding="utf-8")
        payload["camera"] = camera
        payload["pymol_version"] = str(cmd.get_version()[0])
    for assets in payload["assets"].values():
        png_path = Path(assets["png"])
        if not png_path.is_file() or png_path.stat().st_size < 1000:
            raise RuntimeError(f"PyMOL render missing or empty: {png_path}")
    payload["status"] = "complete"
    write_json(payload_path, payload)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cases", type=Path,
                        default=ROOT / "artifacts/evaluation_v3/analysis/selected_cases.json")
    parser.add_argument("--output-dir", type=Path,
                        default=ROOT / "artifacts/evaluation_v3/structure_renders")
    parser.add_argument("--pymol-python", type=Path,
                        default=Path(r"D:\Anaconda\envs\pvb_cpu\python.exe"))
    parser.add_argument("--vmax", type=float, default=5.0)
    parser.add_argument("--size", type=int, default=1200)
    parser.add_argument("--case", action="append", help="Only render the given case_id (repeatable)")
    parser.add_argument("--prepare-only", action="store_true")
    parser.add_argument("--worker", type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.worker:
        worker(args.worker)
        return
    if args.vmax <= 0 or args.size < 200:
        raise ValueError("Require a positive RMSF scale and image size >=200 pixels")
    cases_payload = json.loads(args.cases.read_text(encoding="utf-8"))
    cases = cases_payload if isinstance(cases_payload, list) else cases_payload["cases"]
    if args.case:
        cases = [case for case in cases if case["case_id"] in set(args.case)]
        if set(args.case) != {case["case_id"] for case in cases}:
            raise ValueError("Requested case ID does not exist")
    if not cases:
        raise ValueError("No cases to render")
    if len({case["case_id"] for case in cases}) != len(cases):
        raise ValueError("Duplicate selected case IDs")
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    executable = args.pymol_python.resolve()
    env = os.environ.copy()
    env["PATH"] = os.pathsep.join([str(executable.parent / "Library/bin"),
                                  str(executable.parent / "DLLs"),
                                  str(executable.parent), env.get("PATH", "")])
    reports = []
    for case in cases:
        payload = prepare_case(case, out, args.vmax, args.size)
        payload["status"] = "prepared"
        payload_path = out / f"{case['case_id']}_render.json"
        write_json(payload_path, payload)
        print(f"Prepared {case['case_id']}: {payload['valid_residues']}/"
              f"{payload['selected_residues']} mapped residues", flush=True)
        if not args.prepare_only:
            subprocess.run([str(executable), "-u", str(Path(__file__).resolve()),
                            "--worker", str(payload_path)], check=True, env=env,
                           timeout=240, creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))
            payload = json.loads(payload_path.read_text(encoding="utf-8"))
            print(f"Rendered {case['case_id']}", flush=True)
        reports.append(payload)
    # Per-case reports remain available when this invocation renders a subset;
    # the manifest deliberately describes only this invocation's exact cases.
    manifest = {"description": "Static inference inputs colored by MD target and FlexPPI RMSF",
                "ensemble": False, "numeric_scale_angstrom": [0, args.vmax],
                "color_stops": list(COLORS), "cases": reports}
    write_json(out / "render_manifest.json", manifest)


if __name__ == "__main__":
    main()
