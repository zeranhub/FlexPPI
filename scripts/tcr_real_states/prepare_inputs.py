"""Prepare experimental TCR states without generating an apo by partner deletion.

Coordinates and author residue IDs are copied, not relaxed, completed or mutated.
Full/common and TCR-only are explicitly separate context sensitivity inputs.
"""
from __future__ import annotations

import hashlib
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.SeqUtils import seq1
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "artifacts/tcr_pmhc_reassessment_20260915/sources"
OUT = ROOT / "artifacts/tcr_real_states_20260915"
BB = {"N", "CA", "C", "O"}
AA = set("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split())
NANODISC_CHAINS = list("ABDEFGXY")


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2), encoding="utf-8")


def read_cif(pdb_id):
    path = SOURCE / (pdb_id + ".cif")
    cif = MMCIF2Dict(str(path))
    polymer_entities = {
        ent for ent, kind in zip(cif["_entity_poly.entity_id"], cif["_entity_poly.type"])
        if kind.startswith("polypeptide")
    }
    names = dict(zip(cif["_entity.id"], cif["_entity.pdbx_description"]))
    columns = ["group_PDB", "type_symbol", "label_atom_id", "label_alt_id",
               "label_comp_id", "auth_asym_id", "auth_seq_id", "pdbx_PDB_ins_code",
               "label_entity_id", "label_seq_id", "Cartn_x", "Cartn_y", "Cartn_z",
               "occupancy", "B_iso_or_equiv", "pdbx_PDB_model_num"]
    grouped = {}
    excluded = defaultdict(int)
    for values in zip(*[cif["_atom_site." + c] for c in columns]):
        row = dict(zip(columns, values))
        if row["pdbx_PDB_model_num"] != "1":
            excluded["other_model_atoms"] += 1
            continue
        if row["label_entity_id"] not in polymer_entities:
            excluded["nonprotein_atoms"] += 1
            continue
        if row["type_symbol"].upper() in ("H", "D"):
            excluded["hydrogen_atoms"] += 1
            continue
        if row["label_comp_id"] not in AA:
            excluded["nonstandard_protein_atoms"] += 1
            continue
        ins = row["pdbx_PDB_ins_code"]
        ins = "" if ins in ("?", ".") else ins
        key = (row["auth_asym_id"], row["auth_seq_id"] + ins)
        row.update(ins=ins, xyz=[float(row[c]) for c in ("Cartn_x", "Cartn_y", "Cartn_z")],
                   occupancy_float=float(row["occupancy"]), bfactor_float=float(row["B_iso_or_equiv"]))
        grouped.setdefault(key, []).append(row)
    chains = defaultdict(list)
    rejected, altlocs = [], []
    original_indices = defaultdict(int)
    for (chain, author), rows in grouped.items():
        original_index = original_indices[chain]
        original_indices[chain] += 1
        alternate_ids = sorted({r["label_alt_id"] for r in rows} - {".", "?"})
        chosen = None
        if alternate_ids:
            # One alternate per residue: highest mean occupancy, lexical tie break.
            chosen = min(alternate_ids, key=lambda a: (
                -np.mean([r["occupancy_float"] for r in rows if r["label_alt_id"] == a]), a))
            altlocs.append(dict(chain=chain, author_id=author, alternatives=alternate_ids, selected=chosen))
        selected = [r for r in rows if r["label_alt_id"] in (".", "?", chosen)]
        atoms = {}
        for row in sorted(selected, key=lambda r: (-r["occupancy_float"], r["label_alt_id"], r["label_atom_id"])):
            atoms.setdefault(row["label_atom_id"], row)
        missing = sorted(BB - set(atoms))
        if missing:
            rejected.append(dict(chain=chain, author_id=author, reason="missing_backbone", missing_atoms=missing))
            continue
        if len({r["label_comp_id"] for r in atoms.values()}) != 1:
            raise ValueError(f"Microheterogeneous residue requires manual review: {pdb_id}:{chain}:{author}")
        first = next(iter(atoms.values()))
        chains[chain].append(dict(chain=chain, author_id=author, resseq=int(first["auth_seq_id"]),
            insertion_code=first["ins"], residue_name=first["label_comp_id"],
            amino_acid=seq1(first["label_comp_id"]), entity_id=first["label_entity_id"],
            label_seq_id=int(first["label_seq_id"]), source_observed_chain_index=original_index,
            full_valid_chain_index=len(chains[chain]), atoms=list(atoms.values())))
    audit = dict(pdb_id=pdb_id, source_path=str(path), source_url=f"https://files.rcsb.org/download/{pdb_id}.cif",
                 source_sha256=digest(path), source_model=0, cif_model_number=1,
                 protein_entities={e:names[e] for e in polymer_entities},
                 method=cif.get("_exptl.method"),
                 resolution=cif.get("_em_3d_reconstruction.resolution", cif.get("_refine.ls_d_res_high")),
                 excluded_atom_counts=dict(excluded), rejected_residues=rejected,
                 altloc_selection=altlocs, full_valid_chain_sizes={c:len(rs) for c,rs in chains.items()})
    return dict(chains), audit


def same_author_common(structures, specifications):
    """Common, sequence-identical author IDs across specified equivalent chains."""
    sets = []
    records = []
    for pdb, chain in specifications:
        lookup = {r["author_id"]:r for r in structures[pdb][chain]}
        sets.append(set(lookup))
        records.append(lookup)
    common = set.intersection(*sets)
    mismatched = [a for a in common if len({r[a]["amino_acid"] for r in records}) != 1]
    if mismatched:
        raise ValueError(f"Author-ID correspondence has sequence mismatches: {specifications} {mismatched}")
    return common


def residue_order_key(r):
    return r["resseq"], r["insertion_code"]


def write_state(state_id, case, pdb_id, chain_order, target_chains, context, coverage,
                structures, audit, masks=None, full_state_id=None):
    rows = []
    for chain in chain_order:
        chain_rows = structures[pdb_id][chain]
        if masks and chain in masks:
            chain_rows = [r for r in chain_rows if r["author_id"] in masks[chain]]
        rows.extend(chain_rows)
    if not rows or set(chain_order) != {r["chain"] for r in rows}:
        raise ValueError("Empty selected protein chain")
    path = OUT / "inputs" / (state_id + ".pdb")
    path.parent.mkdir(parents=True, exist_ok=True)
    header = [f"REMARK 900 SOURCE {pdb_id} EXPERIMENTAL MODEL 1; NO COORDINATE CHANGES",
              f"REMARK 900 STATE {state_id}",
              f"REMARK 900 CONTEXT {context}; COVERAGE {coverage}",
              "REMARK 900 PROTEIN ONLY; NONPROTEIN COMPONENTS NOT REPRESENTED BY FLEXPPI"]
    lines = list(header)
    serial = 1
    previous_chain = None
    for r in rows:
        if previous_chain is not None and previous_chain != r["chain"]:
            lines.append("TER")
        previous_chain = r["chain"]
        for atom in sorted(r["atoms"], key=lambda a: (a["label_atom_id"] not in BB,
                            ["N", "CA", "C", "O"].index(a["label_atom_id"]) if a["label_atom_id"] in BB else a["label_atom_id"])):
            name = atom["label_atom_id"]
            element = atom["type_symbol"].upper()
            atom_field = f" {name:<3}" if len(element) == 1 and len(name) < 4 else f"{name:<4}"
            x, y, z = atom["xyz"]
            lines.append(f"ATOM  {serial:5d} {atom_field:4s} {r['residue_name']:>3s} {r['chain']:1s}"
                         f"{r['resseq']:4d}{r['insertion_code'] or ' ':1s}   "
                         f"{x:8.3f}{y:8.3f}{z:8.3f}{atom['occupancy_float']:6.2f}{atom['bfactor_float']:6.2f}          "
                         f"{element:>2s}  ")
            serial += 1
    lines.extend(["TER", "END"])
    path.write_text("\n".join(lines) + "\n", encoding="ascii")
    parsed = PDBParser(QUIET=True).get_structure(state_id, str(path))[0]
    parsed_residues = [r for c in parsed for r in c]
    assert len(parsed_residues) == len(rows)
    for original, parsed_residue in zip(rows, parsed_residues):
        assert set(parsed_residue.child_dict) >= BB
        assert parsed_residue.id[1:] == (original["resseq"], original["insertion_code"] or " ")
        for atom in original["atoms"]:
            assert np.allclose(parsed_residue[atom["label_atom_id"]].coord, atom["xyz"], atol=0.00051)
    residue_info = {}
    mapping = []
    full_index = {(r["chain"], r["author_id"]):i for i, r in enumerate(
        r for c in chain_order for r in structures[pdb_id][c])}
    for chain in chain_order:
        chain_rows = [r for r in rows if r["chain"] == chain]
        residue_info[chain] = dict(sequence="".join(r["amino_acid"] for r in chain_rows),
            author_ids=[r["author_id"] for r in chain_rows],
            source_residue_indices=[r["source_observed_chain_index"] for r in chain_rows],
            full_valid_chain_indices=[r["full_valid_chain_index"] for r in chain_rows],
            label_seq_ids=[r["label_seq_id"] for r in chain_rows])
    for i, r in enumerate(rows):
        mapping.append(dict(input_index=i, chain=r["chain"], author_id=r["author_id"],
            amino_acid=r["amino_acid"], full_state_index=full_index[(r["chain"], r["author_id"])],
            full_valid_chain_index=r["full_valid_chain_index"], label_seq_id=r["label_seq_id"],
            source_observed_chain_index=r["source_observed_chain_index"], backbone_valid=True))
    mapping_path = OUT / "mappings" / (state_id + ".json")
    save(mapping_path, mapping)
    state = dict(id=state_id, case=case, pdb_id=pdb_id, path=str(path), chains=chain_order,
        target_chains=target_chains, context=context, coverage=coverage, source_model=0,
        residues=residue_info, n_residues=len(rows), mapping_path=str(mapping_path),
        full_state_id=full_state_id or state_id,
        provenance=dict(source_path=audit["source_path"], source_sha256=audit["source_sha256"],
            source_url=audit["source_url"], input_sha256=digest(path),
            coordinate_rounding_max_A=0.0005, coordinate_operations="None; PDB serialization only",
            sequence_mutations=0, missing_residue_completion=False,
            nonprotein_context_omitted=True, backbone_rejected=audit["rejected_residues"]),
        valid_backbone_mask=[True] * len(rows),
        warning=("Context sensitivity input, not an experimentally measured isolated equilibrium state."
                 if context == "tcr_only" else "Full observed protein context; membrane, glycans, solvent and unresolved protein omitted."))
    print(state_id, len(rows), flush=True)
    return state


def contrast(cid, case, reference, comparison, chain_map, lookup, kind, notes):
    a, b = lookup[reference], lookup[comparison]
    rows = []
    for ca, cb in chain_map.items():
        ra, rb = a["residues"][ca], b["residues"][cb]
        da = dict(zip(ra["author_ids"], ra["sequence"]))
        db = dict(zip(rb["author_ids"], rb["sequence"]))
        for author in ra["author_ids"]:
            if author not in db:
                continue
            if da[author] != db[author]:
                raise ValueError(f"Contrast mismatch {cid}:{ca}:{author}")
            rows.append(dict(reference_chain=ca, reference_author_id=author,
                comparison_chain=cb, comparison_author_id=author, amino_acid=da[author]))
    if not rows:
        raise ValueError("No matched TCR residues")
    return dict(id=cid, case=case, reference_state=reference, comparison_state=comparison,
        comparison_minus_reference=True, kind=kind, chain_map=chain_map, mapping=rows,
        n_matched_residues=len(rows), notes=notes)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    structures, audits = {}, {}
    for pdb in ["5IW1", "5IVX", "8TW4", "8TW6", "9C3E", "9BBC"]:
        structures[pdb], audits[pdb] = read_cif(pdb)
    b4_common = {
        "alpha":same_author_common(structures, [("5IW1", c) for c in "ACE"] + [("5IVX", "E")]),
        "beta":same_author_common(structures, [("5IW1", c) for c in "BDF"] + [("5IVX", "F")])}
    nd_common = {c:same_author_common(structures, [(p,c) for p in ["8TW4", "8TW6", "9C3E", "9BBC"]])
                 for c in NANODISC_CHAINS}
    states = []
    for pair in ["AB", "CD", "EF"]:
        full_id = "b4_5iw1_" + pair.lower() + "_full"
        for coverage in ["full", "common"]:
            masks = {pair[0]:b4_common["alpha"], pair[1]:b4_common["beta"]}
            states.append(write_state(full_id.replace("_full", "_" + coverage), "B4.2.3", "5IW1",
                list(pair), list(pair), "native_apo_tcr", coverage, structures, audits["5IW1"],
                masks if coverage == "common" else None, full_id))
    for context, chains in [("full", list("EFABP")), ("tcr", list("EF"))]:
        full_id = "b4_5ivx_" + context + "_full"
        for coverage in ["full", "common"]:
            states.append(write_state(full_id.rsplit("_", 1)[0] + "_" + coverage, "B4.2.3", "5IVX",
                chains, list("EF"), "full_holo" if context == "full" else "tcr_only", coverage,
                structures, audits["5IVX"], {"E":b4_common["alpha"], "F":b4_common["beta"]}
                if coverage == "common" else None, full_id))
    for pdb in ["8TW4", "8TW6", "9C3E", "9BBC"]:
        for context in ["full", "tcr"]:
            chains = NANODISC_CHAINS + (["H"] if pdb == "9C3E" else []) if context == "full" else list("AB")
            full_id = "1g4_" + pdb.lower() + "_" + context + "_full"
            for coverage in ["full", "common"]:
                states.append(write_state(full_id.rsplit("_",1)[0] + "_" + coverage, "1G4-CD3", pdb,
                    chains, list("AB"), ("full_holo" if pdb == "9C3E" else "native_apo_full")
                    if context == "full" else "tcr_only", coverage, structures, audits[pdb],
                    nd_common if coverage == "common" else None, full_id))
    lookup = {s["id"]:s for s in states}
    contrasts = []
    for coverage in ["full", "common"]:
        for pair in ["AB", "CD", "EF"]:
            contrasts.append(contrast(f"b4_{pair.lower()}_apo_holo_{coverage}", "B4.2.3",
                f"b4_5iw1_{pair.lower()}_{coverage}", f"b4_5ivx_full_{coverage}",
                {pair[0]:"E", pair[1]:"F"}, lookup, "real_apo_holo",
                "Independent experimental apo and holo; apo ASU copies are technical/crystal controls, not biological replicates."))
        for pair in ["CD", "EF"]:
            contrasts.append(contrast(f"b4_apo_ab_vs_{pair.lower()}_{coverage}", "B4.2.3",
                f"b4_5iw1_ab_{coverage}", f"b4_5iw1_{pair.lower()}_{coverage}",
                {"A":pair[0], "B":pair[1]}, lookup, "apo_crystal_copy_control",
                "Within one apo crystal; not an uncertainty estimate across independent experiments."))
        for pdb in ["8tw4", "8tw6", "9bbc"]:
            contrasts.append(contrast(f"1g4_{pdb}_holo_{coverage}", "1G4-CD3",
                f"1g4_{pdb}_full_{coverage}", f"1g4_9c3e_full_{coverage}", {"A":"A", "B":"B"}, lookup,
                "detergent_apo_holo_environment_confound" if pdb == "9bbc" else "real_apo_holo",
                "True independently determined coordinates; common inputs equalize observed TCR/CD3 coverage, not fusion constructs or missing-state physics."))
        contrasts.append(contrast(f"1g4_nd_i_vs_ii_{coverage}", "1G4-CD3", f"1g4_8tw4_full_{coverage}",
            f"1g4_8tw6_full_{coverage}", {"A":"A", "B":"B"}, lookup, "apo_state_control",
            "Two unbound nanodisc conformations; not biological replicates."))
        for pdb in ["8tw4", "8tw6"]:
            contrasts.append(contrast(f"1g4_{pdb}_vs_detergent_{coverage}", "1G4-CD3",
                f"1g4_{pdb}_full_{coverage}", f"1g4_9bbc_full_{coverage}", {"A":"A", "B":"B"}, lookup,
                "apo_environment_control", "No antigen in either state; detergent/nanodisc and construct confounded."))
    # Keep context/coverage comparisons separate from biological state contrasts.
    for s in states:
        chain_map = {c:c for c in s["target_chains"]}
        if s["coverage"] == "common":
            contrasts.append(contrast(s["id"] + "_coverage_sensitivity", s["case"], s["full_state_id"],
                s["id"], chain_map, lookup, "observed_coverage_sensitivity",
                "Same coordinates at shared residues; omitted unresolved/nonshared context is artificial. Not apo/holo."))
        if s["context"] == "tcr_only":
            full_context_id = s["id"].replace("_tcr_", "_full_")
            contrasts.append(contrast(s["id"] + "_context_sensitivity", s["case"], full_context_id,
                s["id"], chain_map, lookup, "partner_context_sensitivity",
                "TCR-only extraction preserves the state's coordinates. This is NOT a new measured apo state."))
    contacts = []
    for pdb, targets in [("5IVX", "EF"), ("9C3E", "AB")]:
        target_xyz = np.asarray([a["xyz"] for c in targets for r in structures[pdb][c] for a in r["atoms"]])
        tree = cKDTree(target_xyz)
        for chain, residues in structures[pdb].items():
            if chain in targets:
                continue
            xyz = [a["xyz"] for r in residues for a in r["atoms"]]
            contacts.append(dict(pdb_id=pdb, partner_chain=chain, n_residues=len(residues),
                minimum_heavy_atom_distance_A=float(tree.query(xyz)[0].min()),
                role=("pMHC" if pdb == "5IVX" or chain == "H" else "CD3")))
    save(OUT / "states.json", states)
    save(OUT / "contrasts.json", contrasts)
    save(OUT / "preparation_audit.json", dict(source_audits=audits, interface_chain_audit=contacts,
        common_coverage=dict(b4_alpha=len(b4_common["alpha"]), b4_beta=len(b4_common["beta"]),
            one_g4_by_chain={c:len(ids) for c,ids in nd_common.items()}),
        common_ids=dict(b4={k:sorted(v) for k,v in b4_common.items()}, one_g4={k:sorted(v) for k,v in nd_common.items()}),
        state_count=len(states), contrast_count=len(contrasts), script_sha256=digest(__file__),
        selection_policy="First deposited model, deposited protein assembly/copy only, standard-AA heavy atoms, all N/CA/C/O required, no mutation/relaxation/cropping.",
        numbering_warning="PDB preserves author IDs. FlexPPI parser compresses observed positions; common/TCR subsets must use mappings.full_valid_chain_index to retain the corresponding full-input within-chain index.",
        context_warning="Full means all observed protein components, NOT a membrane simulation. No lipids, glycans, waters, ions or unresolved residues are invented."))
    print(json.dumps(dict(states=len(states), contrasts=len(contrasts), contacts=contacts), indent=2))


if __name__ == "__main__":
    main()
