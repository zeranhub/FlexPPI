"""Recover only omitted GROMOS ILE Hchi references for the CDR reassessment.

The deposited GROMOS53a6 structures call ILE's terminal CG1-connected carbon
``CD``; FlexPPI's PDB atom convention calls it ``CD1``. This explicit, local
alias repairs the reference lookup, without changing structures, predictions,
trained weights, the common label code, or any previous analysis output.
"""
from pathlib import Path
import hashlib
import json
import re
import sys
import warnings

import MDAnalysis as mda
from MDAnalysis.lib.distances import calc_dihedrals, minimize_vectors
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from flexppi.constants import CHI_ATOMS
from flexppi.data.labels import _entropy

OLD = ROOT / "output/tcr_prediction_reference_revision_20260915"
OUT = ROOT / "output/tcr_cdr_reassessment_20260917"
EVIDENCE = ROOT / "artifacts/tcr_cdr_reassessment_20260917/ile_atom_alias"
RTP = EVIDENCE / "gromacs_v2019.6_gromos53a6_aminoacids.rtp"
RTP_URL = ("https://raw.githubusercontent.com/gromacs/gromacs/v2019.6/"
           "share/top/gromos53a6.ff/aminoacids.rtp")
ILE_QUADS = [("N", "CA", "CB", "CG1"), ("CA", "CB", "CG1", "CD")]


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main():
    warnings.filterwarnings("ignore", module="MDAnalysis")
    OUT.mkdir(parents=True, exist_ok=True)
    assert RTP.is_file(), "Download the official tagged RTP source before running."
    rtp_text = RTP.read_text(encoding="utf-8")
    ile_section = re.search(r"(?ms)^\[ ILE \](.*?)(?=^\[ [A-Z0-9]+ \]|\Z)", rtp_text).group(1)
    assert re.search(r"(?m)^\s+CD\s+CH3\s", ile_section)
    assert re.search(r"(?m)^\s+CG1\s+CD\s+gb_27\s", ile_section)
    assert "CD1" not in ile_section
    assert CHI_ATOMS["ILE"] == [("N", "CA", "CB", "CG1"), ("CA", "CB", "CG1", "CD1")]

    protected_paths = [OLD / "1g4_run_paired_residues.csv", OLD / "1g4_validation_audit.json",
                       ROOT / "scripts/tcr_validation_revision/1g4_reference.py",
                       ROOT / "src/flexppi/constants.py", ROOT / "src/flexppi/data/labels.py"]
    protected_before = {str(p): sha(p) for p in protected_paths}
    original = pd.read_csv(protected_paths[0], dtype={"author_id": str})
    audit = json.loads(protected_paths[1].read_text(encoding="utf-8"))
    old_hchi = original[original.metric.eq("hchi")].copy()
    missing_non_agp = old_hchi[~old_hchi.evaluation_mask & ~old_hchi.amino_acid.isin(list("AGP"))]
    assert len(missing_non_agp) == 280 and set(missing_non_agp.amino_acid) == {"I"}
    assert missing_non_agp.reference.isna().all()
    assert np.isfinite(missing_non_agp.prediction).all()
    identity_columns = ["chain", "author_id", "amino_acid"]
    expected_keys = set(map(tuple, missing_non_agp[identity_columns].drop_duplicates().values))
    assert len(expected_keys) == 14
    reference_rows = []
    run_audits = []
    control_differences = []

    for source in audit["runs"]:
        state, run = source["state"], int(source["run"])
        pdb, xtc = Path(source["source_pdb"]), Path(source["source_xtc"])
        assert sha(pdb) == source["pdb_sha256"] and sha(xtc) == source["xtc_sha256"]
        universe = mda.Universe(str(pdb), str(xtc))
        chain_names = dict(zip(source["target_chains"], ["alpha", "beta"]))
        residues = [r for seg in universe.segments if seg.segid in chain_names for r in seg.residues]
        recovered = [r for r in residues if r.resname == "ILE"]
        actual_keys = {(chain_names[r.segid], str(r.resid), "I") for r in recovered}
        assert actual_keys == expected_keys

        # Controls span both TCR chains and one through four chi angles. Their
        # coordinates/names are unchanged, and must reproduce stored labels.
        control_residues = []
        represented = set()
        for residue in residues:
            if residue.resname in {"ILE", "ALA", "GLY", "PRO"}:
                continue
            definitions = CHI_ATOMS.get(residue.resname, [])
            names = set(residue.atoms.names)
            group = (chain_names[residue.segid], len(definitions))
            if definitions and group not in represented and all(set(q) <= names for q in definitions):
                represented.add(group)
                control_residues.append(residue)
        assert len(control_residues) == 8

        quad_indices, residue_quad_map, bond_indices = [], [], []
        selected = recovered + control_residues
        for residue in selected:
            atom_names = list(residue.atoms.names)
            definitions = ILE_QUADS if residue.resname == "ILE" else CHI_ATOMS[residue.resname]
            required = set(sum((list(q) for q in definitions), []))
            assert all(atom_names.count(n) == 1 for n in required)
            lookup = {str(atom.name): atom.index for atom in residue.atoms}
            if residue.resname == "ILE":
                assert "CD1" not in lookup and atom_names.count("CD") == 1
                bond_indices.append([lookup["CG1"], lookup["CD"]])
            column_indices = []
            for names in definitions:
                column_indices.append(len(quad_indices))
                quad_indices.append([lookup[n] for n in names])
            residue_quad_map.append(column_indices)

        quad_indices = np.asarray(quad_indices, dtype=int)
        bond_indices = np.asarray(bond_indices, dtype=int)
        codes, times, bond_lengths = [], [], []
        for frame in universe.trajectory:
            if frame.time < 20000.0 - 1e-3:
                continue
            assert frame.time <= 100000.0 + 1e-3
            q = universe.atoms.positions[quad_indices]
            angles = np.degrees(calc_dihedrals(q[:, 0], q[:, 1], q[:, 2], q[:, 3],
                                              box=universe.dimensions))
            assert np.isfinite(angles).all()
            codes.append(np.floor(((angles + 180.0) % 360.0) / 120.0).astype(np.int8))
            bond_xyz = universe.atoms.positions[bond_indices]
            bond_delta = bond_xyz[:, 1] - bond_xyz[:, 0]
            if universe.dimensions is not None:
                bond_delta = minimize_vectors(bond_delta, universe.dimensions)
            bond_lengths.append(np.linalg.norm(bond_delta, axis=1))
            times.append(float(frame.time))
        assert len(times) == 641 and times[0] == 20000.0 and times[-1] == 100000.0
        assert np.allclose(np.diff(times), 125.0, atol=1e-3, rtol=0)
        codes, bond_lengths = np.stack(codes), np.stack(bond_lengths)
        assert np.all((bond_lengths > 1.3) & (bond_lengths < 1.75)), "CG1-CD is not a reasonable C-C bond"
        entropies = []
        for columns in residue_quad_map:
            joint = sum(codes[:, j].astype(np.int32) * (3 ** p) for p, j in enumerate(columns))
            entropies.append(_entropy(joint))
        old_run = old_hchi[old_hchi.state.eq(state) & old_hchi.run.eq(run)]
        controls = []
        for i, residue in enumerate(selected):
            chain, author_id = chain_names[residue.segid], str(residue.resid)
            value = float(entropies[i])
            if i < 14:
                old = old_run[old_run.chain.eq(chain) & old_run.author_id.eq(author_id)]
                assert len(old) == 1 and old.iloc[0].amino_acid == "I" and pd.isna(old.iloc[0].reference)
                reference_rows.append(dict(state=state, run=run, chain=chain,
                                           author_id=author_id, amino_acid="I", reference=value))
            else:
                old = old_run[old_run.chain.eq(chain) & old_run.author_id.eq(author_id)]
                assert len(old) == 1 and old.iloc[0].evaluation_mask
                difference = abs(value - float(old.iloc[0].reference))
                assert difference < 1e-10, (state, run, chain, author_id, difference)
                control_differences.append(difference)
                controls.append(dict(chain=chain, author_id=author_id, resname=residue.resname,
                                     n_chi=len(residue_quad_map[i]), recomputed=value,
                                     original=float(old.iloc[0].reference), absolute_difference=difference))
        run_audits.append(dict(state=state, run=run, source_pdb=str(pdb), source_xtc=str(xtc),
            pdb_sha256=source["pdb_sha256"], xtc_sha256=source["xtc_sha256"],
            target_chains=source["target_chains"], n_recovered=14, n_frames=len(times),
            time_start_ps=times[0], time_end_ps=times[-1], frame_interval_ps=125.0,
            CG1_CD_bond_min_A=float(bond_lengths.min()), CG1_CD_bond_max_A=float(bond_lengths.max()),
            controls=controls))
        print(f"{state} run {run:02d}: recovered 14 ILE; 8 existing controls match; "
              f"CG1-CD {bond_lengths.min():.4f}–{bond_lengths.max():.4f} A", flush=True)

    recovered_table = pd.DataFrame(reference_rows)
    assert len(recovered_table) == 280
    assert recovered_table.groupby(["state", "run"]).size().eq(14).all()
    assert recovered_table.groupby(identity_columns).size().eq(20).all()
    assert not recovered_table.duplicated(["state", "run", "chain", "author_id"]).any()
    assert np.isfinite(recovered_table.reference).all()
    assert recovered_table.reference.between(-1e-12, np.log(9) + 1e-12).all()
    output_csv = OUT / "ile_reference_recovery.csv"
    recovered_table.to_csv(output_csv, index=False)
    protected_after = {str(p): sha(p) for p in protected_paths}
    assert protected_before == protected_after
    recovery_audit = dict(
        status="passed", scope="New supplemental reference only; no retraining, inference, source trajectory edit, or old-output edit",
        omitted_reason="Original CHI_ATOMS uses ILE CD1; GROMOS53a6 topology provides CD. Only chi1 was found, so the original all-chi completeness check omitted ILE Hchi.",
        alias="ILE only: standard PDB CD1 is looked up as GROMOS CD; same atom coordinates, no mutation or geometric change",
        official_evidence=dict(url=RTP_URL, release="GROMACS v2019.6", local_path=str(RTP), sha256=sha(RTP),
            verified="ILE CD has CH3 atom type and an explicit CG1-CD gb_27 bond in the official residue topology"),
        chi_atoms=ILE_QUADS, chi_bin_width_degrees=120.0, entropy="Joint rotamer Shannon entropy, natural logarithm; project _entropy",
        periodic_boundary="MDAnalysis calc_dihedrals with each frame's unit cell; minimize_vectors for CG1-CD bond checks",
        reference_window_ps=[20000.0, 100000.0], expected_frames_per_run=641,
        n_runs=20, n_unique_ile_positions=14, n_recovered_values=280,
        positions=[dict(zip(identity_columns, k)) for k in sorted(expected_keys)],
        n_unchanged_control_values=len(control_differences), max_control_absolute_difference=max(control_differences),
        protected_old_files_sha256_before=protected_before, protected_old_files_sha256_after=protected_after,
        runner_sha256=sha(__file__), recovered_csv_sha256=sha(output_csv),
        library_versions=dict(numpy=np.__version__, pandas=pd.__version__, MDAnalysis=mda.__version__),
        source_archive_doi=audit["archive_doi"], source_paper_doi=audit["source_doi"], runs=run_audits)
    (OUT / "ile_reference_recovery_audit.json").write_text(json.dumps(recovery_audit, indent=2), encoding="utf-8")
    print(json.dumps({"status": "passed", "values": 280, "controls": len(control_differences),
                      "max_control_difference": max(control_differences), "output": str(output_csv)}, indent=2))


if __name__ == "__main__":
    main()
