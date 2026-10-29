"""Complete peptide-neighbour contact audit for three MEL5 pMHC states.

Only existing, processed ITrace protein trajectories are read.  No trajectory,
topology, prediction or manuscript source is changed.  Occupancies describe
sampled bound-state contacts, not interaction free energies or lifetimes.
"""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys
import time

for variable in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ[variable] = "1"

ROOT = Path(__file__).resolve().parent
EXP = Path(r"E:\Project\FlexPPI\experiments\_sources")
sys.path.insert(0, str(EXP / "label_sensitivity/runtime"))
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array
import numpy as np
import pandas as pd
from Bio.Align import PairwiseAligner

RAW = Path(r"E:\data\ITrace\trajectories")
STATES = ("4jff", "4jfd", "4jfe")
SEQUENCES = {"4jff": "ELAGIGILTV", "4jfd": "ELAAIGILTV", "4jfe": "ELAGIGALTV"}
ROLES = {"A": "HLA", "B": "beta2M", "C": "peptide", "D": "TCR_alpha", "E": "TCR_beta"}
MAP_FILE = EXP / "tcr_expanded_screen_v41_20261017/mel5_validation/corrected_pair_all_aligned_positions.csv"
THREE = dict(zip(
    "ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split(),
    "ARNDCQEGHILKMFPSTWYV",
))
THREE.update({name: "H" for name in ("HID", "HIE", "HIP", "HSD", "HSE", "HSP")})


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_csv(frame: pd.DataFrame, filename: str) -> None:
    frame.to_csv(ROOT / filename, index=False, encoding="utf-8-sig", float_format="%.17g")


def heavy_atoms(residue) -> np.ndarray:
    return np.array([a.index for a in residue.atoms if a.element not in ("H", "D")], dtype=int)


def protein_chains(universe) -> dict:
    chains = {}
    for residue in universe.residues:
        assert residue.resname in THREE, ("Non-protein residue in processed topology", residue.resname)
        chain = str(residue.atoms.chainIDs[0])
        assert chain in ROLES, chain
        chains.setdefault(chain, []).append(residue)
    assert list(chains) == ["A", "B", "C", "D", "E"]
    return chains


paper = pd.read_csv(MAP_FILE, dtype={"id_a": str, "id_b": str})
paper_lookup = {
    (str(row.chain_b), str(row.id_b)): (str(row.id_a), str(row.region))
    for row in paper.itertuples()
}
source_hashes = {str(MAP_FILE): digest(MAP_FILE)}
initial = {}
chains = {}
mapping = {}
mapping_rows = []
sequence_audit = []
aligner = PairwiseAligner(mode="global", match_score=2, mismatch_score=-1,
                         open_gap_score=-6, extend_gap_score=-0.5)

for state in STATES:
    top = RAW / f"{state}_run1_topology.pdb"
    initial[state] = mda.Universe(str(top))
    chains[state] = protein_chains(initial[state])
    assert "".join(THREE[r.resname] for r in chains[state]["C"]) == SEQUENCES[state]
    mapping[state] = {}
    for chain, residues in chains[state].items():
        reference = chains["4jff"][chain]
        reference_sequence = "".join(THREE[r.resname] for r in reference)
        current_sequence = "".join(THREE[r.resname] for r in residues)
        alignment = aligner.align(reference_sequence, current_sequence)[0]
        pairs = []
        for ref_block, cur_block in zip(*alignment.aligned):
            assert ref_block[1] - ref_block[0] == cur_block[1] - cur_block[0]
            pairs.extend(zip(range(*ref_block), range(*cur_block)))
        differences = []
        for ri, ci in pairs:
            key = (chain, int(ri + 1))
            mapping[state][key] = residues[ci]
            ref_res, cur_res = reference[ri], residues[ci]
            paper_id, region = paper_lookup.get((chain, str(ref_res.resid)), ("", "outside_common_mapping"))
            row = {"state": state, "chain": chain, "role": ROLES[chain],
                   "reference_4jff_position": ri + 1, "reference_4jff_author_id": str(ref_res.resid),
                   "state_position": ci + 1, "state_author_id": str(cur_res.resid),
                   "reference_aa": THREE[ref_res.resname], "state_aa": THREE[cur_res.resname],
                   "same_aa": THREE[ref_res.resname] == THREE[cur_res.resname],
                   "paper_common_author_id": paper_id, "paper_region": region}
            mapping_rows.append(row)
            if not row["same_aa"]:
                differences.append({"reference_position": ri + 1, "state_position": ci + 1,
                                    "reference_aa": row["reference_aa"], "state_aa": row["state_aa"]})
        sequence_audit.append({"state": state, "chain": chain, "role": ROLES[chain],
                               "reference_sequence": reference_sequence, "state_sequence": current_sequence,
                               "reference_length": len(reference), "state_length": len(residues),
                               "aligned_residues": len(pairs), "amino_acid_differences": differences,
                               "unmapped_reference_positions": sorted(set(range(1, len(reference) + 1)) - {a + 1 for a, b in pairs}),
                               "unmapped_state_positions": sorted(set(range(1, len(residues) + 1)) - {b + 1 for a, b in pairs})})

write_csv(pd.DataFrame(mapping_rows), "all_state_to_4jff_and_paper_mapping.csv")
meta_lookup = {(r["state"], r["chain"], r["reference_4jff_position"]): r for r in mapping_rows}
peptide_keys = [("C", p) for p in range(1, 11)]
candidate_keys = set()
candidate_reasons = []
static_distances = {}

# Candidate definition is a union over all three run-1 preparation structures.
# No MD contact occupancy or prediction value is used to select neighbours.
for state in STATES:
    universe = initial[state]
    peptide_atoms = [heavy_atoms(mapping[state][key]) for key in peptide_keys]
    pi = np.concatenate(peptide_atoms)
    pstart = np.r_[0, np.cumsum([len(x) for x in peptide_atoms])[:-1]]
    non_peptide = [(key, residue) for key, residue in mapping[state].items() if key[0] != "C"]
    for key, residue in non_peptide:
        ni = heavy_atoms(residue)
        matrix = distance_array(universe.atoms.positions[pi], universe.atoms.positions[ni])
        per_peptide = np.minimum.reduceat(matrix, pstart, axis=0).min(axis=1)
        for p, value in enumerate(per_peptide, 1):
            static_distances[(state, p, key)] = float(value)
        minimum = float(per_peptide.min())
        if minimum <= 7:
            candidate_keys.add(key)
            candidate_reasons.append({"state": state, "neighbor_chain": key[0],
                                      "neighbor_reference_position": key[1],
                                      "minimum_static_peptide_distance_A": minimum})

candidate_keys = sorted(candidate_keys)
assert all(key in mapping[state] for state in STATES for key in candidate_keys)
candidate_rows = []
for key in candidate_keys:
    row = {"neighbor_chain": key[0], "neighbor_reference_position": key[1], "neighbor_role": ROLES[key[0]]}
    for state in STATES:
        meta = meta_lookup[(state, *key)]
        row[f"{state}_author_id"] = meta["state_author_id"]
        row[f"{state}_aa"] = meta["state_aa"]
        row[f"{state}_min_start_distance_A"] = min(static_distances[(state, p, key)] for p in range(1, 11))
    row["paper_common_author_id"] = meta["paper_common_author_id"]
    row["paper_region"] = meta["paper_region"]
    row["same_neighbor_aa_all_states"] = len({row[f"{s}_aa"] for s in STATES}) == 1
    candidate_rows.append(row)
write_csv(pd.DataFrame(candidate_rows), "complete_7A_neighbor_union.csv")
write_csv(pd.DataFrame(candidate_reasons), "neighbor_union_inclusion_reasons.csv")
print(f"Candidate union: {len(candidate_keys)} neighbours; {10 * len(candidate_keys)} pairs per trajectory", flush=True)

contact_rows = []
replica_audit = []
for state in STATES:
    for replica in (1, 2, 3):
        started = time.time()
        top = RAW / f"{state}_run{replica}_topology.pdb"
        trajectory = RAW / f"{state}_run{replica}_md_processed.xtc"
        top_hash, trajectory_hash = digest(top), digest(trajectory)
        source_hashes[str(top)] = top_hash
        source_hashes[str(trajectory)] = trajectory_hash
        universe = mda.Universe(str(top), str(trajectory))
        current_chains = protein_chains(universe)
        by_id = {(chain, str(r.resid)): r for chain, residues in current_chains.items() for r in residues}
        # Every replica must preserve the complete chain sequences and residue identifiers.
        for chain, residues in current_chains.items():
            expected = chains[state][chain]
            assert [str(r.resid) for r in residues] == [str(r.resid) for r in expected]
            assert [THREE[r.resname] for r in residues] == [THREE[r.resname] for r in expected]
        peptide = [by_id[("C", str(mapping[state][key].resid))] for key in peptide_keys]
        neighbors = [by_id[(key[0], str(mapping[state][key].resid))] for key in candidate_keys]
        patoms = [heavy_atoms(r) for r in peptide]
        natoms = [heavy_atoms(r) for r in neighbors]
        pi, ni = np.concatenate(patoms), np.concatenate(natoms)
        pstart = np.r_[0, np.cumsum([len(x) for x in patoms])[:-1]]
        nstart = np.r_[0, np.cumsum([len(x) for x in natoms])[:-1]]
        distances = []
        times = []
        for ts in universe.trajectory:
            assert ts.dimensions is not None and np.all(ts.dimensions[:3] > 0)
            matrix = distance_array(ts.positions[pi], ts.positions[ni], box=ts.dimensions)
            nearest = np.minimum.reduceat(np.minimum.reduceat(matrix, pstart, axis=0), nstart, axis=1)
            distances.append(nearest.astype(np.float32))
            times.append(float(ts.time) / 1000)
        universe.trajectory.close()
        distances = np.stack(distances)
        times = np.asarray(times)
        assert len(times) >= 1001 and abs(times[0]) < 1e-8 and abs(times[-1] - 200) < 1e-6
        # Statistics are evaluated from the same stored float32 nearest-distance values.
        # Boundary precision is far smaller than the 0.5 A cutoff sensitivity interval.
        c4 = (distances <= 4).mean(axis=0)
        c45 = (distances <= 4.5).mean(axis=0)
        mean_distance = distances.mean(axis=0, dtype=np.float64)
        median_distance = np.median(distances, axis=0)
        for p in range(1, 11):
            for j, key in enumerate(candidate_keys):
                meta = meta_lookup[(state, *key)]
                contact_rows.append({"state": state, "replica": replica, "frames": len(times),
                    "peptide_position": p, "peptide_aa": SEQUENCES[state][p-1],
                    "neighbor_chain": key[0], "neighbor_reference_position": key[1],
                    "neighbor_state_author_id": meta["state_author_id"], "neighbor_role": ROLES[key[0]],
                    "neighbor_aa": meta["state_aa"], "paper_common_author_id": meta["paper_common_author_id"],
                    "paper_region": meta["paper_region"],
                    "contact4_fraction": float(c4[p-1, j]), "contact4p5_fraction": float(c45[p-1, j]),
                    "mean_minimum_distance_A": float(mean_distance[p-1, j]),
                    "median_minimum_distance_A": float(median_distance[p-1, j]),
                    "starting_minimum_distance_A": static_distances[(state, p, key)]})
        cache = ROOT / f"{state}_rep{replica}_all_peptide_neighbor_distances.npz"
        np.savez_compressed(cache, distances_A=distances, time_ns=times,
                            peptide_positions=np.arange(1, 11),
                            neighbor_chains=np.array([key[0] for key in candidate_keys]),
                            neighbor_reference_positions=np.array([key[1] for key in candidate_keys]))
        replica_audit.append({"state": state, "replica": replica, "topology": str(top),
            "topology_sha256": top_hash, "trajectory": str(trajectory), "trajectory_sha256": trajectory_hash,
            "frames": len(times), "time_ns_first": float(times[0]), "time_ns_last": float(times[-1]),
            "atoms": universe.atoms.n_atoms, "protein_residues": len(universe.residues),
            "water_residues": 0, "nonprotein_residues": 0,
            "nearest_distance_cache": str(cache), "cache_sha256": digest(cache)})
        assert digest(top) == top_hash and digest(trajectory) == trajectory_hash
        print(state, replica, "frames", len(times), "seconds", round(time.time() - started, 2), flush=True)

contacts = pd.DataFrame(contact_rows)
write_csv(contacts, "all_peptide_neighbor_contacts_per_replica.csv")
keys = ["state", "peptide_position", "peptide_aa", "neighbor_chain", "neighbor_reference_position",
        "neighbor_state_author_id", "neighbor_role", "neighbor_aa", "paper_common_author_id", "paper_region"]
summary = contacts.groupby(keys, dropna=False).agg(
    contact4_mean=("contact4_fraction", "mean"), contact4_SD=("contact4_fraction", "std"),
    contact4p5_mean=("contact4p5_fraction", "mean"), contact4p5_SD=("contact4p5_fraction", "std"),
    median_minimum_distance_replica_mean_A=("median_minimum_distance_A", "mean"),
    starting_minimum_distance_A=("starting_minimum_distance_A", "first"),
).reset_index()
write_csv(summary, "all_peptide_neighbor_contacts_three_replica_summary.csv")

pair_keys = ["peptide_position", "neighbor_chain", "neighbor_reference_position"]
contrast_rows = []
for mutant, mutation in (("4jfd", "G4A"), ("4jfe", "I7A")):
    for key, group in contacts.groupby(pair_keys, dropna=False):
        ref = group[group.state == "4jff"].sort_values("replica")
        cur = group[group.state == mutant].sort_values("replica")
        assert len(ref) == len(cur) == 3
        row = {"mutation": mutation, "reference_state": "4jff", "mutant_state": mutant,
               "peptide_position": key[0], "reference_peptide_aa": ref.peptide_aa.iloc[0],
               "mutant_peptide_aa": cur.peptide_aa.iloc[0], "peptide_same_aa": ref.peptide_aa.iloc[0] == cur.peptide_aa.iloc[0],
               "neighbor_chain": key[1], "neighbor_reference_position": key[2],
               "neighbor_reference_author_id": ref.neighbor_state_author_id.iloc[0],
               "neighbor_mutant_author_id": cur.neighbor_state_author_id.iloc[0],
               "neighbor_role": ref.neighbor_role.iloc[0], "neighbor_reference_aa": ref.neighbor_aa.iloc[0],
               "neighbor_mutant_aa": cur.neighbor_aa.iloc[0], "neighbor_same_aa": ref.neighbor_aa.iloc[0] == cur.neighbor_aa.iloc[0],
               "paper_common_author_id": ref.paper_common_author_id.iloc[0], "paper_region": ref.paper_region.iloc[0]}
        for prefix, column in (("contact4", "contact4_fraction"), ("contact4p5", "contact4p5_fraction")):
            a, b = ref[column].to_numpy(), cur[column].to_numpy()
            row.update({f"{prefix}_reference_mean": float(a.mean()), f"{prefix}_mutant_mean": float(b.mean()),
                        f"{prefix}_delta": float(b.mean() - a.mean()),
                        f"{prefix}_reference_SD": float(a.std(ddof=1)), f"{prefix}_mutant_SD": float(b.std(ddof=1)),
                        f"{prefix}_all3_mutant_above_all3_reference": bool(b.min() > a.max()),
                        f"{prefix}_all3_mutant_below_all3_reference": bool(b.max() < a.min())})
            for j in range(3):
                row[f"{prefix}_reference_rep{j+1}"] = float(a[j])
                row[f"{prefix}_mutant_rep{j+1}"] = float(b[j])
        contrast_rows.append(row)
contrasts = pd.DataFrame(contrast_rows)
write_csv(contrasts, "all_peptide_neighbor_contact_contrasts.csv")

# Contact counts are descriptive expected numbers of contacting residue pairs.
# Each trajectory, not each frame or residue pair, is the replicate.
role_counts = contacts.groupby(["state", "replica", "peptide_position", "neighbor_role"], dropna=False)[
    ["contact4_fraction", "contact4p5_fraction"]].sum().reset_index().rename(columns={
        "contact4_fraction": "mean_contacting_neighbors4", "contact4p5_fraction": "mean_contacting_neighbors4p5"})
write_csv(role_counts, "per_peptide_role_contact_counts_per_replica.csv")
whole_counts = role_counts.groupby(["state", "replica", "neighbor_role"], dropna=False)[
    ["mean_contacting_neighbors4", "mean_contacting_neighbors4p5"]].sum().reset_index()
write_csv(whole_counts, "whole_peptide_role_contact_counts_per_replica.csv")

prominent = contrasts[
    (contrasts.contact4_delta.abs() >= 0.20)
    & (contrasts.contact4_all3_mutant_above_all3_reference | contrasts.contact4_all3_mutant_below_all3_reference)
].copy()
prominent["absolute_delta"] = prominent.contact4_delta.abs()
prominent = prominent.sort_values(["mutation", "absolute_delta"], ascending=[True, False])
write_csv(prominent, "exploratory_changes_atleast20pp_with_nonoverlapping_replica_ranges.csv")
for mutation in ("G4A", "I7A"):
    block = prominent[prominent.mutation == mutation]
    print("PROMINENT", mutation, len(block), "pairs; top12 for inspection", flush=True)
    print(block[["peptide_position", "neighbor_role", "paper_common_author_id", "paper_region",
                 "contact4_reference_mean", "contact4_mutant_mean", "contact4_delta"]].head(12).to_string(index=False), flush=True)

audit = {"status": "passed", "new_MD": False, "new_model_inference": False,
    "states": list(STATES), "replicas_per_state": 3, "frame_window_ns": [0, 200], "all_saved_frames": True,
    "peptide_sequences": SEQUENCES, "neighbor_count": len(candidate_keys),
    "residue_pairs_per_state": len(candidate_keys) * 10,
    "candidate_definition": "Union of all mapped non-peptide protein residues with any heavy atom within 7 A of any of the complete 10 peptide residues in any of the three run1 preparation topologies. No MD occupancy or prediction value used for selection.",
    "distance_definition": "Minimum heavy-atom distance between each peptide residue and each candidate protein residue, minimum-image convention using each trajectory frame box.",
    "thresholds_A": [4, 4.5], "occupancy_definition": "Fraction of every saved frame with minimum heavy-atom distance <= threshold; not a lifetime or free energy.",
    "replicate_unit": "One independently run MD trajectory; three per state. Frames and residue pairs are not independent biological replicates.",
    "summary": "Equal-weight mean and sample SD of the three trajectory occupancies. A range separation means every mutant replica exceeds every reference replica, or vice versa; it is descriptive, not a significance test.",
    "exploratory_filter": "Absolute mean occupancy change >=0.20 and non-overlapping three-replica ranges. All candidates and pairs are separately retained, including all zero contacts and opposite changes. No p-values or multiple-test significance claims.",
    "water_limitation": "All processed topologies contain only standard protein residues. Water bridges and solvent entropy cannot be analysed from these trajectories.",
    "sequence_limitation": "The available complex inputs are retrospective state structures with additional TCR-alpha terminal modelling differences. They are not a controlled mutation-only experiment. Sequence alignment and every difference are recorded.",
    "sequence_alignment": {"reference": "4jff", "mode": "global", "match_score": 2, "mismatch_score": -1,
                           "open_gap_score": -6, "extend_gap_score": -0.5, "chains": sequence_audit},
    "software": {"python": sys.version, "MDAnalysis": mda.__version__, "numpy": np.__version__, "pandas": pd.__version__},
    "source_files": replica_audit, "source_file_hashes": source_hashes,
    "script_sha256": digest(Path(__file__))}
(ROOT / "contact_analysis_audit.json").write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8")
readme = """# 肽突变动态接触证据\n\n本目录对已有4JFF、4JFD和4JFE三态各三条0–200 ns轨迹的全部保存帧进行接触分析，没有运行新MD。完整10个肽位置与三态起始结构7 Å邻居并集的所有组合均保留，包括无接触、反向变化及模型未涉及的位置。接触定义为两个残基的最小重原子距离≤4 Å，并提供4.5 Å敏感性口径。\n\n每条轨迹是一个重复，汇总为三重复等权均值和样本标准差。全部三条突变轨迹的接触占据都高于或低于全部三条参考轨迹，称为重复范围分离，仅作为描述，不作为统计显著性。`exploratory_changes_atleast20pp_with_nonoverlapping_replica_ranges.csv`是查看完整结果后的探索性筛选，不是独立验证；其阈值为平均占据改变至少20个百分点且三重复范围分离，全部未筛选结果另行完整保留。\n\n输入轨迹只有蛋白，不能从中验证动态水桥或溶剂熵。TCRα的端部还存在结构准备差异，所有序列差异和编号对齐均见audit及全映射表。以4JFF序列位置统一候选，再经既有3HG1–4JFF映射转换论文共同编号；不能直接把原始MD resid写成论文TCR编号。\n\n主要文件：`all_peptide_neighbor_contacts_per_replica.csv`逐重复全结果；`all_peptide_neighbor_contact_contrasts.csv`两种肽替换全差值；`whole_peptide_role_contact_counts_per_replica.csv`各相互作用方的整肽接触数；`complete_7A_neighbor_union.csv`全部候选；`contact_analysis_audit.json`方法和原始文件SHA256；每条轨迹的npz保存全部逐帧最小距离以供核查。\n"""
(ROOT / "README_肽突变接触证据.md").write_text(readme, encoding="utf-8")
print("COMPLETE", len(contact_rows), "pair-replica observations", flush=True)
