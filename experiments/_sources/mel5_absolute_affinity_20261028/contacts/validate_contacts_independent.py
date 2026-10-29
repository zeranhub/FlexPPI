"""Re-read selected residue pairs directly from raw trajectories.

This deliberately avoids the batch reduce-at calculation in the production
script and independently checks residue identifiers, periodic distances, and
all saved-frame occupancies for the narrative contact candidates.
"""
from pathlib import Path
import json
import sys
sys.path.insert(0, r"E:\Project\FlexPPI\experiments\_sources\label_sensitivity\runtime")
import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array

ROOT = Path(__file__).resolve().parent
RAW = Path(r"E:\data\ITrace\trajectories")
source = pd.read_csv(ROOT / "all_peptide_neighbor_contacts_per_replica.csv")
mapping = pd.read_csv(ROOT / "all_state_to_4jff_and_paper_mapping.csv", dtype={"state_author_id": str})
audit = json.loads((ROOT / "contact_analysis_audit.json").read_text(encoding="utf-8"))
assert len(source) == 7920
assert len(source.drop_duplicates(["state", "replica", "peptide_position", "neighbor_chain", "neighbor_reference_position"])) == 7920
assert source.groupby(["state", "replica"]).size().eq(880).all()
assert source.contact4_fraction.between(0, 1).all()
assert (source.contact4_fraction <= source.contact4p5_fraction).all()
neighbors = pd.read_csv(ROOT / "complete_7A_neighbor_union.csv")
assert len(neighbors) == 88 and neighbors.same_neighbor_aa_all_states.all()
assert all(not block["unmapped_state_positions"] for block in audit["sequence_alignment"]["chains"])

# P7-beta96, P8-beta97, P8-beta99, P1-alpha29 and P1-HLA66.
# Neighbour reference positions are 4JFF numbering, not paper TCR numbering.
pairs = [(7, "E", 97), (8, "E", 98), (8, "E", 100), (1, "D", 30), (1, "A", 66)]
validated = []
for state in ("4jff", "4jfd", "4jfe"):
    for rep in (1, 2, 3):
        u = mda.Universe(str(RAW / f"{state}_run{rep}_topology.pdb"), str(RAW / f"{state}_run{rep}_md_processed.xtc"))
        byid = {(str(r.atoms.chainIDs[0]), str(r.resid)): r for r in u.residues}
        groups = []
        for p, chain, refpos in pairs:
            row = mapping[(mapping.state == state) & (mapping.chain == chain) & (mapping.reference_4jff_position == refpos)].iloc[0]
            a = byid[("C", str(p))]
            b = byid[(chain, row.state_author_id)]
            ai = [x.index for x in a.atoms if x.element not in ("H", "D")]
            bi = [x.index for x in b.atoms if x.element not in ("H", "D")]
            groups.append((ai, bi))
        distances = [[] for _ in pairs]
        for ts in u.trajectory:
            for j, (a, b) in enumerate(groups):
                distances[j].append(float(distance_array(ts.positions[a], ts.positions[b], box=ts.dimensions).min()))
        u.trajectory.close()
        cache = np.load(ROOT / f"{state}_rep{rep}_all_peptide_neighbor_distances.npz")
        for j, (p, chain, refpos) in enumerate(pairs):
            ids = np.flatnonzero((cache["neighbor_chains"] == chain) & (cache["neighbor_reference_positions"] == refpos))
            assert len(ids) == 1
            stored = cache["distances_A"][:, p-1, ids[0]]
            direct = np.array(distances[j])
            np.testing.assert_allclose(stored, direct, atol=2e-6, rtol=0)
            row = source[(source.state == state) & (source.replica == rep) & (source.peptide_position == p)
                         & (source.neighbor_chain == chain) & (source.neighbor_reference_position == refpos)].iloc[0]
            for threshold, col in [(4, "contact4_fraction"), (4.5, "contact4p5_fraction")]:
                assert abs(float((direct <= threshold).mean()) - row[col]) < 1e-12
            validated.append({"state": state, "replica": rep, "peptide_position": p,
                              "neighbor_chain": chain, "neighbor_reference_position": refpos,
                              "frames": len(direct), "max_distance_float32_storage_error_A": float(np.max(abs(stored-direct)))})
        cache.close()

report = {"status": "passed", "independent_direct_pair_recalculations": len(validated),
          "method": "Raw topology and trajectory re-read; each selected residue-pair distance calculated independently, without batch atom-to-residue reduction. All frames and both thresholds verified.",
          "complete_output_shape_and_identity_checks": True,
          "all_88_candidate_neighbors_same_amino_acid_across_three_states": True,
          "all_state_residues_mapped_to_reference": True,
          "no_contact_occupancy_changed_by_float32_distance_storage_for_checked_pairs": True,
          "checks": validated}
(ROOT / "independent_contact_validation.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
print("PASS", len(validated), "pair-trajectory direct recalculations")
