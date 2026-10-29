"""Read-only, all-frame quality summaries of the completed four-state KIX pilot.

Run next to analyze_job.py. Uses its exact alignment/entropy implementation.
No raw trajectory, original analysis, model prediction or MD schedule is changed.
Default output: <root>/pilot_quality_20260916. No independent-replicate statistics.
"""

from __future__ import annotations

import argparse
REPEAT_NAME = "rep01"
import csv
from datetime import datetime, timezone
import io
import json
import os
from pathlib import Path

# Prevent an analysis worker from competing for many BLAS cores with MD/training.
for _thread_variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_thread_variable] = "1"

import numpy as np

from analyze_job import (Moments, align_batch, angle_bins, dihedral_degrees,
                         entropy_from_counts, parse_topology, segment_paths,
                         sha256_file, chi_hash, EXPECTED_CHI_HASH)


STATES = ("apo", "myb", "mll", "both")
EXPECTED_FRAMES = 100000
WINDOW_FRAMES = 10000
HALF_FRAMES = 50000
TRACE_BIN_FRAMES = 100


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"Refusing an unlabelled empty CSV: {path}")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def histograms_for(topology: dict) -> list:
    return [np.zeros(3 ** r["n_chi_angles"], dtype=np.int64) if r["n_chi_angles"] else None
            for r in topology["residues"]]


def snapshot_rows(state: str, topology: dict, moments: Moments, histograms: list,
                  kind: str, start_ns: float, end_ns: float) -> list[dict]:
    rmsf = dict(zip(topology["kix_rows"].tolist(), moments.rmsf().tolist()))
    rows = []
    expected = round((end_ns - start_ns) * 1000)
    assert moments.n == expected
    for index, residue in enumerate(topology["residues"]):
        histogram = histograms[index]
        if histogram is not None:
            assert int(histogram.sum()) == expected
        rows.append(dict(state=state, replicate=REPEAT_NAME, window_kind=kind,
                         start_ns_exclusive=start_ns, end_ns_inclusive=end_ns,
                         n_frames=expected, chain=residue["chain"],
                         author_residue_number=residue["author_residue_number"],
                         amino_acid=residue["amino_acid"],
                         rmsf_angstrom=rmsf.get(index, ""),
                         chi_entropy_nats=entropy_from_counts(histogram) if histogram is not None else "",
                         n_chi_angles=residue["n_chi_angles"],
                         hchi_informative_evaluation=histogram is not None and residue["amino_acid"] not in "AGP"))
    return rows


def native_contacts(topology: dict, cutoff: float = 4.5) -> dict:
    """Fixed initially contacting heavy atom pairs, grouped by residue pair.

    A residue pair is retained only if >=1 of its *initially contacting* atom
    pairs remains within cutoff. This is not a general all-heavy-contact count.
    No distance sampling or Cartesian pair product is repeated during streaming.
    """
    contacts = {}
    positions = topology["reference_positions"]
    residues = topology["residues"]
    for partner in ("A", "C"):
        atom_pairs, starts, labels = [], [], []
        for left in residues:
            if left["chain"] != "B":
                continue
            left_atoms = sorted({index for name, index in left["atom_lookup"].items()
                                 if not name.lstrip("0123456789").startswith(("H", "D"))})
            for right in residues:
                if right["chain"] != partner:
                    continue
                right_atoms = sorted({index for name, index in right["atom_lookup"].items()
                                      if not name.lstrip("0123456789").startswith(("H", "D"))})
                distances = np.linalg.norm(positions[left_atoms, None] - positions[None, right_atoms], axis=-1)
                ix, iy = np.where(distances <= cutoff)
                if len(ix):
                    starts.append(len(atom_pairs))
                    atom_pairs.extend((left_atoms[i], right_atoms[j]) for i, j in zip(ix, iy))
                    labels.append(dict(kix_author_residue=left["author_residue_number"],
                                       kix_amino_acid=left["amino_acid"], partner_chain=partner,
                                       partner_author_residue=right["author_residue_number"],
                                       partner_amino_acid=right["amino_acid"], native_atom_pairs=len(ix)))
        if atom_pairs:
            contacts[partner] = dict(atom_pairs=np.asarray(atom_pairs, dtype=np.int64),
                                     starts=np.asarray(starts), labels=labels,
                                     occupancy_counts=np.zeros(len(starts), dtype=np.int64))
    return contacts


def contact_values(frames: np.ndarray, contact: dict, cutoff: float = 4.5) -> tuple:
    pairs = contact["atom_pairs"]
    square_distances = np.sum((frames[:, pairs[:, 0]] - frames[:, pairs[:, 1]]) ** 2, axis=-1)
    atom_retained = square_distances <= cutoff * cutoff
    residue_retained = np.logical_or.reduceat(atom_retained, contact["starts"], axis=1)
    return atom_retained.mean(axis=1), residue_retained.mean(axis=1), residue_retained.sum(axis=0)


def trace_rows(state: str, frame_values: dict) -> list[dict]:
    rows = []
    for metric, values in frame_values.items():
        assert values.shape == (EXPECTED_FRAMES,) and np.isfinite(values).all()
        binned = values.reshape(-1, TRACE_BIN_FRAMES)
        for index, block in enumerate(binned):
            rows.append(dict(state=state, replicate=REPEAT_NAME, metric=metric,
                             start_ns_exclusive=index / 10, end_ns_inclusive=(index + 1) / 10,
                             n_frames=TRACE_BIN_FRAMES, minimum=float(block.min()),
                             mean=float(block.mean()), maximum=float(block.max())))
    return rows


def thermodynamic_rows(state: str, job: Path) -> tuple:
    """Preserve all production reporter samples, also fixed 0.1 ns summaries."""
    rows, source_hashes = [], []
    for path in sorted((job / "production").glob("segment_*.csv")):
        contents = path.read_text(encoding="utf-8")
        if contents.startswith("#"):
            contents = contents[1:]
        reader = csv.DictReader(io.StringIO(contents))
        keys = reader.fieldnames or []
        def field(prefix):
            matches = [key for key in keys if key.startswith(prefix)]
            if len(matches) != 1:
                raise ValueError(f"Cannot resolve {prefix!r} in {path}: {keys}")
            return matches[0]
        names = {"step": field("Step"), "time_ps": field("Time (ps)"),
                 "temperature_kelvin": field("Temperature"),
                 "potential_energy_kj_mol": field("Potential Energy"),
                 "kinetic_energy_kj_mol": field("Kinetic Energy")}
        segment_rows = []
        for record in reader:
            row = dict(state=state, replicate=REPEAT_NAME, source_segment=path.name)
            row.update({name: float(record[key]) for name, key in names.items()})
            row["step"] = int(row["step"])
            # Step-based clock avoids tiny integrator floating point time drift.
            row["production_ns"] = (row["step"] * .002 - 10000) / 1000
            if not all(np.isfinite(row[key]) for key in names):
                raise ValueError(f"Nonfinite thermodynamic sample in {path}")
            if not np.isclose(row["time_ps"], row["step"] * .002, atol=.01, rtol=0):
                raise ValueError(f"Reporter step/time mismatch in {path}")
            segment_rows.append(row)
        if len(segment_rows) != 100:
            raise ValueError(f"Expected 100 thermodynamic samples in {path}, saw {len(segment_rows)}")
        rows.extend(segment_rows)
        source_hashes.append(dict(filename=path.name, sha256=sha256_file(path), n_samples=len(segment_rows)))
    expected_steps = np.arange(5005000, 55000001, 5000)
    np.testing.assert_array_equal([row["step"] for row in rows], expected_steps)
    summaries = []
    for metric in ("temperature_kelvin", "potential_energy_kj_mol", "kinetic_energy_kj_mol"):
        values = np.array([row[metric] for row in rows]).reshape(-1, 10)
        for index, block in enumerate(values):
            summaries.append(dict(state=state, replicate=REPEAT_NAME, metric=metric,
                                  start_ns_exclusive=index / 10, end_ns_inclusive=(index + 1) / 10,
                                  n_samples=10, minimum=float(block.min()),
                                  mean=float(block.mean()), maximum=float(block.max())))
    return rows, summaries, source_hashes


def analyze_state(root: Path, state: str, verify_sha: bool, contacts_enabled: bool,
                  batch_frames: int) -> tuple:
    from MDAnalysis.coordinates.DCD import DCDReader
    pdb = root / "prepared" / state / "protein.pdb"
    job = root / "runs" / state / REPEAT_NAME
    summary_path = job / "analysis" / "summary.json"
    previous = json.loads(summary_path.read_text())
    done = json.loads((job / "DONE.json").read_text())
    if done["repeat"] != int(REPEAT_NAME[3:]) or previous["replicate"] != REPEAT_NAME:
        raise ValueError("Wrong repeat identity in completed source")
    if done["production_ns"] != 100 or done["final_step"] != 55000000:
        raise ValueError(f"{state}: incomplete DONE marker")
    if previous["n_frames"] != EXPECTED_FRAMES or previous["frame_interval_ps"] != 1.0:
        raise ValueError(f"{state}: incompatible prior analysis")
    if previous["protein_pdb_sha256"] != sha256_file(pdb):
        raise ValueError(f"{state}: protein topology changed since prior analysis")
    original_csv = job / "analysis" / "residue_metrics.csv"
    if previous["output_csv_sha256"] != sha256_file(original_csv):
        raise ValueError(f"{state}: residue metrics hash mismatch")
    if not (.5 <= previous["checked_covalent_distance_min_angstrom"] <= previous["checked_covalent_distance_max_angstrom"] <= 3):
        raise ValueError(f"{state}: prior all-frame geometry validation failed")
    topology = parse_topology(pdb)
    paths = segment_paths(job / "production")
    if len(paths) != 100 or len(previous["segments"]) != len(paths):
        raise ValueError(f"{state}: expected exactly 100 committed segments")
    contacts = native_contacts(topology) if contacts_enabled else {}
    if contacts_enabled and len(contacts) != {"apo": 0, "myb": 1, "mll": 1, "both": 2}[state]:
        raise ValueError(f"{state}: unexpected number of native interfaces")
    histograms = histograms_for(topology)
    half_histograms = [histograms_for(topology), histograms_for(topology)]
    moments, half_moments = Moments(), [Moments(), Moments()]
    residue_rows, position, reference = [], 0, None
    frame_values = {"KIX_CA_RMSD_angstrom": np.empty(EXPECTED_FRAMES),
                    "L12_CA_RMSD_after_KIX_fit_angstrom": np.empty(EXPECTED_FRAMES)}
    for chain in contacts:
        frame_values[f"partner_{chain}_native_atom_contact_fraction"] = np.empty(EXPECTED_FRAMES)
        frame_values[f"partner_{chain}_native_residue_contact_fraction"] = np.empty(EXPECTED_FRAMES)
    l12 = np.array([614 <= topology["residues"][i]["author_residue_number"] <= 621 for i in topology["kix_rows"]])
    assert int(l12.sum()) == 8
    for path, previous_segment in zip(paths, previous["segments"]):
        before = path.stat()
        if path.name != previous_segment["filename"] or before.st_size != previous_segment["bytes"]:
            raise ValueError(f"{state}: trajectory segment identity changed: {path}")
        if verify_sha:
            if sha256_file(path) != previous_segment["sha256"]:
                raise ValueError(f"{state}: DCD SHA-256 mismatch: {path}")
        elif before.st_mtime_ns != previous_segment["mtime_ns"]:
            raise ValueError(f"{state}: unchecked DCD timestamp changed: {path}")
        with DCDReader(str(path)) as reader:
            if len(reader) != 1000 or reader.n_atoms != topology["atom_count"] or not np.isclose(reader.dt, 1.0, rtol=5e-5):
                raise ValueError(f"{state}: unexpected DCD header: {path}")
            for start in range(0, len(reader), batch_frames):
                frames = np.stack([ts.positions.copy() for ts in reader[start:start + batch_frames]]).astype(np.float64)
                if not np.isfinite(frames).all():
                    raise ValueError(f"{state}: nonfinite trajectory coordinates")
                ca = frames[:, topology["kix_ca"]]
                if reference is None:
                    reference = ca[0].copy()
                aligned = align_batch(ca, reference)
                bins = angle_bins(dihedral_degrees(frames[:, topology["chi_quads"]]))
                codes = {}
                for i, residue in enumerate(topology["residues"]):
                    columns = residue["chi_columns"]
                    if columns:
                        codes[i] = np.sum(bins[:, columns].astype(np.int32) * (3 ** np.arange(len(columns))), axis=1)
                squared = np.sum((aligned - reference) ** 2, axis=-1)
                stop_position = position + len(frames)
                frame_values["KIX_CA_RMSD_angstrom"][position:stop_position] = np.sqrt(squared.mean(axis=1))
                frame_values["L12_CA_RMSD_after_KIX_fit_angstrom"][position:stop_position] = np.sqrt(squared[:, l12].mean(axis=1))
                for chain, contact in contacts.items():
                    atom_q, residue_q, occupied = contact_values(frames, contact)
                    frame_values[f"partner_{chain}_native_atom_contact_fraction"][position:stop_position] = atom_q
                    frame_values[f"partner_{chain}_native_residue_contact_fraction"][position:stop_position] = residue_q
                    contact["occupancy_counts"] += occupied
                local_start = 0
                # Always cut at predeclared 10 ns boundaries, independent of batch size.
                while local_start < len(frames):
                    global_start = position + local_start
                    next_boundary = (global_start // WINDOW_FRAMES + 1) * WINDOW_FRAMES
                    local_stop = min(len(frames), next_boundary - position)
                    piece = slice(local_start, local_stop)
                    half = 0 if global_start < HALF_FRAMES else 1
                    moments.update(aligned[piece])
                    half_moments[half].update(aligned[piece])
                    for i, values in codes.items():
                        counts = np.bincount(values[piece], minlength=len(histograms[i]))
                        histograms[i] += counts
                        half_histograms[half][i] += counts
                    total = position + local_stop
                    if total % WINDOW_FRAMES == 0:
                        residue_rows.extend(snapshot_rows(state, topology, moments, histograms,
                                                          "cumulative", 0.0, total / 1000))
                    local_start = local_stop
                position = stop_position
        after = path.stat()
        if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
            raise ValueError(f"{state}: segment changed during analysis: {path}")
    assert position == EXPECTED_FRAMES
    for half in (0, 1):
        residue_rows.extend(snapshot_rows(state, topology, half_moments[half], half_histograms[half],
                                          "disjoint_half", half * 50.0, (half + 1) * 50.0))
    original = list(csv.DictReader(original_csv.open(encoding="utf-8", newline="")))
    original_by_key = {(r["chain"], int(r["author_residue_number"])): r for r in original}
    final = [r for r in residue_rows if r["window_kind"] == "cumulative" and r["end_ns_inclusive"] == 100]
    differences = {"rmsf_angstrom": [], "chi_entropy_nats": []}
    for row in final:
        reference_row = original_by_key[(row["chain"], row["author_residue_number"])]
        for metric in differences:
            if row[metric] != "":
                expected = float(reference_row[metric])
                np.testing.assert_allclose(row[metric], expected, rtol=1e-10, atol=1e-10)
                differences[metric].append(abs(row[metric] - expected))
    occupancy_rows = []
    contact_summary = {}
    for chain, contact in contacts.items():
        for label, count in zip(contact["labels"], contact["occupancy_counts"]):
            occupancy_rows.append(dict(state=state, replicate=REPEAT_NAME, **label,
                                       n_frames=EXPECTED_FRAMES, frames_contact_retained=int(count),
                                       fraction_frames_retained=float(count / EXPECTED_FRAMES)))
        contact_summary[chain] = dict(native_atom_pairs=len(contact["atom_pairs"]),
                                      native_residue_pairs=len(contact["labels"]),
                                      mean_atom_pair_retention=float(frame_values[f"partner_{chain}_native_atom_contact_fraction"].mean()),
                                      mean_residue_pair_retention=float(frame_values[f"partner_{chain}_native_residue_contact_fraction"].mean()))
    thermo, thermo_binned, thermo_hashes = thermodynamic_rows(state, job)
    details = dict(state=state, n_frames=position, n_protein_residues=len(final), KIX_residues=87,
                   prior_summary_sha256=sha256_file(summary_path), original_residue_metrics_sha256=sha256_file(original_csv),
                   protein_pdb_sha256=sha256_file(pdb), dcd_hashes_verified=verify_sha,
                   full_repeat_metric_max_abs_difference={key: max(values) for key, values in differences.items()},
                   prior_all_frame_covalent_distance_min_angstrom=previous["checked_covalent_distance_min_angstrom"],
                   prior_all_frame_covalent_distance_max_angstrom=previous["checked_covalent_distance_max_angstrom"],
                   contacts=contact_summary, thermodynamic_sample_count=len(thermo), thermodynamic_sources=thermo_hashes)
    return residue_rows, trace_rows(state, frame_values), occupancy_rows, thermo, thermo_binned, details


def self_test() -> dict:
    rng = np.random.default_rng(2938)
    coords = rng.normal(size=(100, 87, 3))
    reference = coords[0]
    fitted = align_batch(coords, reference)
    all_moments = Moments()
    for piece in np.array_split(fitted, 9):
        all_moments.update(piece)
    expected = np.sqrt(np.mean(np.sum((fitted - fitted.mean(axis=0)) ** 2, axis=-1), axis=0))
    np.testing.assert_allclose(all_moments.rmsf(), expected, atol=1e-12)
    np.testing.assert_array_equal(angle_bins(np.array([-180., -60., 60., 180.])), [0, 1, 2, 2])
    np.testing.assert_allclose(entropy_from_counts(np.array([25, 25, 50])), 1.0397207708399179)
    fake_contact = dict(atom_pairs=np.array([[0, 1], [0, 2], [3, 4]]), starts=np.array([0, 2]))
    frames = np.array([[[0, 0, 0], [1, 0, 0], [2, 0, 0], [0, 0, 0], [4, 0, 0]],
                       [[0, 0, 0], [8, 0, 0], [2, 0, 0], [0, 0, 0], [9, 0, 0]]], dtype=float)
    aq, rq, occupancy = contact_values(frames, fake_contact)
    np.testing.assert_allclose(aq, [1, 1 / 3])
    np.testing.assert_allclose(rq, [1, .5])
    np.testing.assert_array_equal(occupancy, [2, 1])
    return dict(self_test="passed", fixture_type="synthetic only, not MD evidence",
                streamed_RMSF="passed", entropy_and_bins="passed", grouped_native_contact_retention="passed")


def main() -> None:
    global REPEAT_NAME
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--states", nargs="+", choices=STATES, default=list(STATES))
    parser.add_argument("--batch-frames", type=int, default=256)
    parser.add_argument("--skip-dcd-sha256", action="store_true", help="Use prior bytes+mtime instead of rehashing; recorded in output")
    parser.add_argument("--skip-contacts", action="store_true")
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--repeat", choices=["rep01", "rep02", "rep03"], default="rep01")
    args = parser.parse_args()
    REPEAT_NAME = args.repeat
    if args.self_test:
        print(json.dumps(self_test(), indent=2))
        return
    if args.root is None:
        parser.error("--root is required")
    if not 1 <= args.batch_frames <= 4096:
        parser.error("--batch-frames must be 1..4096")
    if len(args.states) != len(set(args.states)):
        parser.error("Duplicate states are not allowed")
    if chi_hash() != EXPECTED_CHI_HASH:
        raise ValueError("Frozen chi atom definition changed")
    root = args.root.resolve(strict=True)
    output = (args.output_dir or root / "quality" / REPEAT_NAME).resolve()
    if output.exists() and any(output.iterdir()):
        raise FileExistsError(f"Output must be empty to protect existing results: {output}")
    output.mkdir(parents=True, exist_ok=True)
    collections = [[], [], [], [], []]
    state_summaries = []
    for state in args.states:
        print(f"Analyzing all 100000 production frames: {state}", flush=True)
        result = analyze_state(root, state, not args.skip_dcd_sha256, not args.skip_contacts, args.batch_frames)
        for destination, source in zip(collections, result[:5]):
            destination.extend(source)
        state_summaries.append(result[5])
        print(json.dumps({"state": state, "status": "passed", "frames": result[5]["n_frames"],
                          "full_repeat_metric_max_abs_difference": result[5]["full_repeat_metric_max_abs_difference"]}), flush=True)
    outputs = {}
    names = ("residue_window_metrics.csv", "trajectory_trace_binned.csv", "native_contact_occupancy.csv",
             "thermodynamics_all_samples.csv", "thermodynamics_binned.csv")
    for name, rows in zip(names, collections):
        if rows:
            write_csv(output / name, rows)
            outputs[name] = dict(rows=len(rows), sha256=sha256_file(output / name))
    summary = dict(protocol="kix_rebuilt_single_repeat_allframe_quality_v1", replicate=REPEAT_NAME, created_utc=datetime.now(timezone.utc).isoformat(),
                   script_sha256=sha256_file(Path(__file__)), analyzer_script_sha256=sha256_file(Path(__file__).with_name("analyze_job.py")),
                   input_root=str(root), states=state_summaries, outputs=outputs,
                   time_definition="Production frame i (zero based) = (i+1) ps. 0.1 ns trace bins include all 100 frames, not subsampling; logs retain all 10 ps reporter samples.",
                   windows="Predeclared cumulative 0-10,...,0-100 ns and disjoint 0-50/50-100 ns; no burn-in deletion or best-window selection.",
                   RMSF_definition="Identical analyze_job implementation: fit all KIX CA to the first full-production frame; each window fluctuates around its own aligned-coordinate mean.",
                   Hchi_definition="Identical analyze_job full joint 3x120-degree chi-state entropy (natural logarithm), separately counted within each window. AG have no chi; P retained with hchi_informative_evaluation=False.",
                   RMSD_definition="Root mean square displacement from first full-production frame after KIX-only CA fit. L12 uses residues 614-621 under the same whole-KIX fit.",
                   native_contacts_definition="Prepared/protein.pdb heavy-atom pairs between KIX B and each partner initially <=4.5 Angstrom. Atomic retention tests those same fixed pairs in every frame. A native residue-pair contact is retained when >=1 of its originally contacting atom pairs remains <=4.5 Angstrom; new atom pairs are not counted. No minimum-image correction, consistent with prior verified unwrapped covalent geometry.",
                   geometry_validation="Requires unchanged prior analysis CSV, protein PDB and all trajectory identities; by default re-verifies every DCD SHA-256. Prior summary verified selected backbone/peptide/chi covalent bonds in every frame; this pass does not silently repair geometry.",
                   limitations=["This file describes one selected independent repeat per state; combine all three new repeats only after preserving each result. No convergence certificate.",
                                "Window agreement is descriptive stability evidence only; cumulative windows overlap.",
                                "MD values are simulation references, not experimental measurements.",
                                "Native contact retention excludes new contacts and may decline despite a rearranged bound interface.",
                                "Fixed 0.1 ns min/mean/max bands describe within-bin values, not uncertainty intervals."])
    (output / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"status": "complete", "output_dir": str(output), "states": args.states, "outputs": outputs}, indent=2))


if __name__ == "__main__":
    main()
