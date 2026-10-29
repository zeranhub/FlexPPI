"""Adapter for the strict apo/holo PPI trajectories in Zenodo record 6638504.

The archive stores GROMACS ``gro`` topologies and ``xtc`` trajectories.  FlexPPI
uses PDB structures because chain identifiers are part of its residue geometry
representation.  This module converts the first GRO snapshot without changing
atom order, assigns chains from the exact free-partner atom counts, normalises a
small set of force-field residue aliases, and writes state/pair manifests.

The atom-order guarantee matters: an XTC file contains coordinates only, so a
reordered PDB could still open while silently attaching coordinates to the wrong
atoms.  Every generated topology is therefore reopened together with its XTC
and checked for equal atom and protein-residue counts.
"""

from __future__ import annotations

from copy import deepcopy
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

from flexppi.constants import THREE_TO_ONE
from flexppi.data.pairing import match_apo_to_holo
from flexppi.data.records import PairRecord, StateRecord, write_jsonl
from flexppi.data.structures import parse_pdb


FORCE_FIELD_RESIDUE_ALIASES = {
    "ASH": "ASP",
    "CYM": "CYS",
    "CYX": "CYS",
    "GLH": "GLU",
    "HID": "HIS",
    "HIE": "HIS",
    "HIP": "HIS",
    "LYN": "LYS",
}

TERMINAL_CAPS = {"ACE", "NME"}


def normalise_residue_name(name: str) -> str:
    """Map common Amber/GROMACS protonation aliases to standard PDB names."""

    upper = str(name).strip().upper()
    return FORCE_FIELD_RESIDUE_ALIASES.get(upper, upper)


def protein_residue_names(names: Iterable[str]) -> List[str]:
    """Return standard protein residue names after alias normalisation."""

    normalised = [normalise_residue_name(name) for name in names]
    return [name for name in normalised if name in THREE_TO_ONE]


def chain_ids_from_atom_counts(
    atom_count: int,
    chain_atom_counts: Sequence[int],
    chain_names: Sequence[str],
) -> List[str]:
    """Expand exact per-chain atom counts into one PDB chain ID per atom."""

    if len(chain_atom_counts) != len(chain_names):
        raise ValueError("chain_atom_counts and chain_names must have equal length")
    if any(int(count) <= 0 for count in chain_atom_counts):
        raise ValueError("Every chain must contain at least one atom")
    if sum(int(count) for count in chain_atom_counts) != int(atom_count):
        raise ValueError(
            "Chain atom counts {} do not sum to topology atom count {}".format(
                list(chain_atom_counts), atom_count
            )
        )
    result: List[str] = []
    for count, name in zip(chain_atom_counts, chain_names):
        if len(str(name)) != 1:
            raise ValueError("PDB chain identifiers must be one character")
        result.extend([str(name)] * int(count))
    return result


def _load_mdanalysis() -> Any:
    try:
        import MDAnalysis as mda
    except ImportError as exc:
        raise RuntimeError(
            "Zenodo PPI preparation requires MDAnalysis: pip install -e .[md]"
        ) from exc
    return mda


def _write_pdb_from_gro(
    gro_path: Path,
    output_path: Path,
    chain_atom_counts: Sequence[int],
    chain_names: Sequence[str],
) -> Dict[str, Any]:
    mda = _load_mdanalysis()
    universe = mda.Universe(str(gro_path))
    chain_ids = chain_ids_from_atom_counts(
        universe.atoms.n_atoms, chain_atom_counts, chain_names
    )
    normalised_names = [normalise_residue_name(name) for name in universe.residues.resnames]
    universe.residues.resnames = normalised_names
    if not hasattr(universe.atoms, "chainIDs"):
        universe.add_TopologyAttr("chainIDs")
    universe.atoms.chainIDs = chain_ids
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with mda.Writer(str(output_path), universe.atoms.n_atoms, reindex=True) as writer:
        writer.write(universe.atoms)
    structure = parse_pdb(str(output_path))
    return {
        "atoms": int(universe.atoms.n_atoms),
        "all_residues": int(universe.residues.n_residues),
        "protein_residues": int(structure.num_residues),
        "chains": structure.unique_chains(),
        "sequence": structure.sequence,
    }


def _write_chain_subset(source_pdb: Path, output_pdb: Path, chain_id: str) -> None:
    """Write a partner-deleted PDB while preserving retained atom coordinates."""

    output_pdb.parent.mkdir(parents=True, exist_ok=True)
    kept: List[str] = []
    with source_pdb.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(("ATOM  ", "HETATM")) and line[21:22] == chain_id:
                kept.append(line)
    if not kept:
        raise ValueError("No atoms for chain {} in {}".format(chain_id, source_pdb))
    kept.append("END\n")
    output_pdb.write_text("".join(kept), encoding="utf-8")


def _trajectory_metadata(topology: Path, trajectory: Path) -> Dict[str, Any]:
    mda = _load_mdanalysis()
    universe = mda.Universe(str(topology), str(trajectory))
    frames = int(len(universe.trajectory))
    first_time = float(universe.trajectory[0].time)
    last_time = float(universe.trajectory[-1].time)
    return {
        "atoms": int(universe.atoms.n_atoms),
        "residues": int(universe.residues.n_residues),
        "frames": frames,
        "frame_interval_ps": float(universe.trajectory.dt),
        "production_ns": (last_time - first_time) / 1000.0,
    }


def _state(
    sample_id: str,
    structure_path: Path,
    trajectory_path: Path,
    target_chains: Sequence[str],
    source_group: str,
    metadata: Dict[str, Any],
    label_dir: Path,
) -> StateRecord:
    label_path = label_dir / (sample_id + ".npz")
    structure = parse_pdb(str(structure_path))
    return StateRecord(
        sample_id=sample_id,
        dataset="zenodo_6638504",
        structure_path=str(structure_path.resolve()),
        label_path=str(label_path.resolve()) if label_path.exists() else None,
        trajectory_paths=[str(trajectory_path.resolve())],
        temperature=300.0,
        window_ns=100.0,
        target_chains=list(target_chains),
        source_group=source_group,
        metadata={
            **metadata,
            "residue_count": int(structure.num_residues),
            "zenodo_record": "6638504",
        },
    )


def _copy_state(record: StateRecord) -> StateRecord:
    return StateRecord.from_dict(deepcopy(record.to_dict()))


def index_zenodo_ppi(
    root: str,
    output_dir: str,
    label_dir: Optional[str] = None,
    limit: Optional[int] = None,
) -> Dict[str, Any]:
    """Prepare structures and manifests for the eight strict apo/holo systems."""

    dataset_root = Path(root).resolve()
    extracted_root = dataset_root / "extracted" if (dataset_root / "extracted").is_dir() else dataset_root
    output_root = Path(output_dir).resolve()
    structure_root = output_root / "structures"
    labels = Path(label_dir).resolve() if label_dir else output_root / "labels"
    systems = sorted(path for path in extracted_root.iterdir() if path.is_dir())
    if limit is not None:
        systems = systems[:limit]

    states: List[StateRecord] = []
    pairs: List[PairRecord] = []
    audit: Dict[str, Any] = {}
    for system in systems:
        pdb_id = system.name.upper()
        apo_a_gro = system / "A" / "start.gro"
        apo_b_gro = system / "B" / "start.gro"
        apo_a_xtc = system / "A" / "md.xtc"
        apo_b_xtc = system / "B" / "md.xtc"
        holo_gro = system / "complex" / "without_water" / "start.gro"
        holo_xtc = system / "complex" / "without_water" / "md.xtc"
        required = [apo_a_gro, apo_b_gro, apo_a_xtc, apo_b_xtc, holo_gro, holo_xtc]
        missing = [str(path) for path in required if not path.is_file()]
        if missing:
            raise FileNotFoundError("{} is missing required files: {}".format(pdb_id, missing))

        mda = _load_mdanalysis()
        apo_a_universe = mda.Universe(str(apo_a_gro))
        apo_b_universe = mda.Universe(str(apo_b_gro))
        holo_universe = mda.Universe(str(holo_gro))
        atom_a = int(apo_a_universe.atoms.n_atoms)
        atom_b = int(apo_b_universe.atoms.n_atoms)
        if atom_a + atom_b != int(holo_universe.atoms.n_atoms):
            raise ValueError(
                "{} complex atom count {} differs from free A+B {}+{}".format(
                    pdb_id, holo_universe.atoms.n_atoms, atom_a, atom_b
                )
            )

        system_output = structure_root / pdb_id
        apo_a_pdb = system_output / "apo_A.pdb"
        apo_b_pdb = system_output / "apo_B.pdb"
        holo_pdb = system_output / "holo_complex.pdb"
        holo_a_pdb = system_output / "holo_A_only.pdb"
        holo_b_pdb = system_output / "holo_B_only.pdb"
        pdb_a = _write_pdb_from_gro(apo_a_gro, apo_a_pdb, [atom_a], ["A"])
        pdb_b = _write_pdb_from_gro(apo_b_gro, apo_b_pdb, [atom_b], ["B"])
        pdb_holo = _write_pdb_from_gro(holo_gro, holo_pdb, [atom_a, atom_b], ["A", "B"])
        _write_chain_subset(holo_pdb, holo_a_pdb, "A")
        _write_chain_subset(holo_pdb, holo_b_pdb, "B")

        if pdb_holo["sequence"] != pdb_a["sequence"] + pdb_b["sequence"]:
            raise ValueError("{} protein sequence/order differs between free and complex states".format(pdb_id))
        trajectory_a = _trajectory_metadata(apo_a_pdb, apo_a_xtc)
        trajectory_b = _trajectory_metadata(apo_b_pdb, apo_b_xtc)
        trajectory_holo = _trajectory_metadata(holo_pdb, holo_xtc)
        if trajectory_a["frames"] != trajectory_b["frames"] or trajectory_a["frames"] != trajectory_holo["frames"]:
            raise ValueError("{} free and complex trajectories have different frame counts".format(pdb_id))

        base_metadata = {"pdb_id": pdb_id, "cluster_id": pdb_id}
        state_a = _state(
            pdb_id + "_apo_A", apo_a_pdb, apo_a_xtc, ["A"], "apo",
            {**base_metadata, **trajectory_a, "side": "A"}, labels,
        )
        state_b = _state(
            pdb_id + "_apo_B", apo_b_pdb, apo_b_xtc, ["B"], "apo",
            {**base_metadata, **trajectory_b, "side": "B"}, labels,
        )
        state_holo = _state(
            pdb_id + "_holo_complex", holo_pdb, holo_xtc, ["A", "B"], "holo",
            {**base_metadata, **trajectory_holo, "side": "AB"}, labels,
        )
        states.extend([state_a, state_b, state_holo])

        holo_structure = parse_pdb(str(holo_pdb))
        for side, apo_state in (("A", state_a), ("B", state_b)):
            apo_structure = parse_pdb(apo_state.structure_path)
            apo_indices, holo_indices, holo_chains, identity, coverage = match_apo_to_holo(
                apo_structure, holo_structure
            )
            if holo_chains != [side] or identity < 0.999 or coverage < 0.999:
                raise ValueError(
                    "{} side {} mapping failed: chains={}, identity={}, coverage={}".format(
                        pdb_id, side, holo_chains, identity, coverage
                    )
                )
            apo_copy = _copy_state(apo_state)
            holo_copy = _copy_state(state_holo)
            holo_copy.target_chains = [side]
            pairs.append(
                PairRecord(
                    pair_id="{}__{}".format(pdb_id, side),
                    apo=apo_copy,
                    holo=holo_copy,
                    apo_indices=apo_indices,
                    holo_indices=holo_indices,
                    sequence_identity=identity,
                    coverage=coverage,
                    metadata={
                        **base_metadata,
                        "side": side,
                        "partner_chain": "B" if side == "A" else "A",
                        "holo_target_chain": side,
                        "holo_target_only_pdb": str(
                            (holo_a_pdb if side == "A" else holo_b_pdb).resolve()
                        ),
                    },
                )
            )

        audit[pdb_id] = {
            "apo_A": {**pdb_a, **trajectory_a},
            "apo_B": {**pdb_b, **trajectory_b},
            "holo": {**pdb_holo, **trajectory_holo},
            "atom_order_check": "complex atom count equals exact free A+B counts",
            "sequence_check": "complex protein sequence equals free A followed by free B",
        }

    state_manifest = output_root / "zenodo_6638504_states.jsonl"
    pair_manifest = output_root / "zenodo_6638504_pairs.jsonl"
    write_jsonl(str(state_manifest), (record.to_dict() for record in states))
    write_jsonl(str(pair_manifest), (record.to_dict() for record in pairs))
    return {
        "systems": len(systems),
        "states": len(states),
        "pairs": len(pairs),
        "labelled_states": sum(record.label_path is not None for record in states),
        "state_manifest": str(state_manifest),
        "pair_manifest": str(pair_manifest),
        "label_dir": str(labels),
        "audit": audit,
    }
