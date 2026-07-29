"""Create PyMOL-ready flexibility comparison bundles."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any, Dict, List, Mapping, Sequence, Tuple

import numpy as np
import torch

from flexppi.baselines.fastprotflex import parse_heavy_protein_atoms
from flexppi.comparison import _load_fast_outputs, _run_fastprotflex
from flexppi.constants import THREE_TO_ONE
from flexppi.data.dataset import collate_states
from flexppi.data.structures import parse_pdb
from flexppi.evaluation import _rankdata
from flexppi.inference import _structure_example, load_model
from flexppi.trainer import move_to_device


ResidueKey = Tuple[str, str, str]


def _chain_percentiles(values: np.ndarray, chains: Sequence[str]) -> np.ndarray:
    output = np.full(len(values), np.nan, dtype=np.float64)
    chain_array = np.asarray(chains, dtype=str)
    for chain in dict.fromkeys(chains):
        selected = (chain_array == chain) & np.isfinite(values)
        count = int(selected.sum())
        if count == 0:
            continue
        if count == 1:
            output[selected] = 50.0
        else:
            output[selected] = 100.0 * (_rankdata(values[selected]) - 1.0) / (count - 1)
    return output


def _write_bfactor_pdb(
    input_path: str,
    output_path: Path,
    values_by_key: Mapping[ResidueKey, float],
    remarks: Sequence[str],
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with Path(input_path).open("r", encoding="utf-8", errors="replace") as source:
        lines = source.readlines()
    with output_path.open("w", encoding="utf-8", newline="") as destination:
        for remark in remarks:
            destination.write("REMARK 900 " + remark[:68] + "\n")
        for line in lines:
            record = line[:6].strip()
            if record not in {"ATOM", "HETATM"} or len(line) < 54:
                destination.write(line)
                continue
            residue_name = line[17:20].strip().upper()
            key = (
                line[21:22].strip() or "_",
                line[22:26].strip(),
                line[26:27].strip(),
            )
            value = values_by_key.get(key)
            if residue_name not in THREE_TO_ONE or value is None or not np.isfinite(value):
                destination.write(line)
                continue
            stripped = line.rstrip("\r\n").ljust(66)
            destination.write(stripped[:60] + "{:6.2f}".format(float(value)) + stripped[66:] + "\n")


def _write_pymol_script(
    output_path: Path,
    flex_percentile_pdb: Path,
    fast_percentile_pdb: Path,
    session_path: Path,
    pdb_id: str,
) -> None:
    flex_path = flex_percentile_pdb.resolve().as_posix()
    fast_path = fast_percentile_pdb.resolve().as_posix()
    session = session_path.resolve().as_posix()
    script = """# FlexPPI versus FastProtFlex flexibility comparison for {pdb_id}
# B-factor contains within-chain flexibility percentile: 0=rigid, 100=flexible.
reinitialize
bg_color white
set ray_opaque_background, off
set antialias, 2
set cartoon_fancy_helices, 1
set cartoon_smooth_loops, 1
set stick_radius, 0.16

load {flex_path}, FlexPPI_full_complex
load {fast_path}, FastProtFlex_chainwise

hide everything, all
show cartoon, FlexPPI_full_complex
show cartoon, FastProtFlex_chainwise
show sticks, (FlexPPI_full_complex or FastProtFlex_chainwise) and (chain C or ((chain D+E) and resi 106+107+108+109+110+111+112+135+136+137+138+139))

spectrum b, blue_white_red, FlexPPI_full_complex, minimum=0, maximum=100
spectrum b, blue_white_red, FastProtFlex_chainwise, minimum=0, maximum=100

# Chain roles: A+B=MHC class II, C=peptide, D=TCR alpha, E=TCR beta.
select FlexPPI_MHC, FlexPPI_full_complex and chain A+B
select FlexPPI_peptide, FlexPPI_full_complex and chain C
select FlexPPI_TCR, FlexPPI_full_complex and chain D+E
select FlexPPI_CDR3_alpha, FlexPPI_full_complex and chain D and resi 106+107+108+109+110+111+112+135+136+137+138+139
select FlexPPI_CDR3_beta, FlexPPI_full_complex and chain E and resi 106+107+108+109+110+111+112+135+136+137+138+139
select FastProtFlex_MHC, FastProtFlex_chainwise and chain A+B
select FastProtFlex_peptide, FastProtFlex_chainwise and chain C
select FastProtFlex_TCR, FastProtFlex_chainwise and chain D+E
select FastProtFlex_CDR3_alpha, FastProtFlex_chainwise and chain D and resi 106+107+108+109+110+111+112+135+136+137+138+139
select FastProtFlex_CDR3_beta, FastProtFlex_chainwise and chain E and resi 106+107+108+109+110+111+112+135+136+137+138+139

set grid_mode, 1
set grid_slot, 1, FlexPPI_full_complex
set grid_slot, 2, FastProtFlex_chainwise
orient
zoom all, 5
set seq_view, 1

# Opening this .pml creates a reusable PyMOL session at the path below.
save {session}
""".format(
        pdb_id=pdb_id,
        flex_path=json.dumps(flex_path),
        fast_path=json.dumps(fast_path),
        session=json.dumps(session),
    )
    output_path.write_text(script, encoding="utf-8")


@torch.no_grad()
def create_pymol_comparison_bundle(
    pdb_path: str,
    pdb_id: str,
    flexppi_checkpoint: str,
    output_dir: str,
    rscript_path: str,
    fastprotflex_batch_script: str,
    chain_roles: Mapping[str, str],
    temperature: float = 300.0,
    window_ns: float = 100.0,
    dataset_code: float = 0.0,
    device_name: str = "auto",
) -> Dict[str, Any]:
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    structure = parse_pdb(pdb_path)
    chains = structure.unique_chains()
    unknown_roles = sorted(set(chains) - set(chain_roles))
    if unknown_roles:
        raise ValueError("Missing roles for chains: {}".format(", ".join(unknown_roles)))

    model, model_config, device = load_model(flexppi_checkpoint, device_name)
    rmsf_index = model_config.metrics.index("rmsf")
    batch = collate_states([
        _structure_example(
            pdb_path,
            chains,
            temperature,
            window_ns,
            dataset_code,
        )
    ])
    model_output = model.forward_state(move_to_device(batch, device))
    flex_rmsf = model_output["mean"][0, :, rmsf_index].cpu().numpy().astype(np.float64)
    flex_std = torch.exp(
        0.5 * model_output["log_variance"][0, :, rmsf_index]
    ).cpu().numpy().astype(np.float64)
    flex_valid = model_output["mask"][0].cpu().numpy()
    flex_rmsf[~flex_valid] = np.nan
    flex_std[~flex_valid] = np.nan

    manifest_path = output / (pdb_id + "_fastprotflex_chain_manifest.tsv")
    fast_prediction_path = output / (pdb_id + "_fastprotflex_chain_predictions.csv")
    fast_timing_path = output / (pdb_id + "_fastprotflex_chain_timings.csv")
    with manifest_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample_id", "pdb_path", "chain_filter"])
        for chain in chains:
            writer.writerow([pdb_id + "__chain__" + chain, pdb_path, chain])
    _run_fastprotflex(
        manifest_path,
        fast_prediction_path,
        fast_timing_path,
        rscript_path,
        fastprotflex_batch_script,
    )
    _, fast_keyed, fast_timings = _load_fast_outputs(
        fast_prediction_path, fast_timing_path
    )
    failures = [
        row for row in fast_timings.values() if row.get("status") != "ok"
    ]
    if failures:
        raise RuntimeError("FastProtFlex chain failures: {}".format(failures))

    heavy_atoms = parse_heavy_protein_atoms(pdb_path)
    residue_keys = list(dict.fromkeys(heavy_atoms.residue_keys))
    if len(residue_keys) != structure.num_residues:
        raise ValueError(
            "PDB standard residue keys {} != FlexPPI residues {}".format(
                len(residue_keys), structure.num_residues
            )
        )
    fast_score = np.full(structure.num_residues, np.nan, dtype=np.float64)
    for residue_index, key in enumerate(residue_keys):
        chain = structure.chain_ids[residue_index]
        fast_score[residue_index] = fast_keyed[pdb_id + "__chain__" + chain].get(
            key, np.nan
        )
    if not np.all(np.isfinite(fast_score)):
        missing = np.flatnonzero(~np.isfinite(fast_score)).tolist()
        raise ValueError("Missing FastProtFlex scores at residue indices {}".format(missing))

    flex_percentile = _chain_percentiles(flex_rmsf, structure.chain_ids)
    fast_percentile = _chain_percentiles(fast_score, structure.chain_ids)
    flex_by_key = dict(zip(residue_keys, flex_rmsf.tolist()))
    flex_std_by_key = dict(zip(residue_keys, flex_std.tolist()))
    fast_by_key = dict(zip(residue_keys, fast_score.tolist()))
    flex_percentile_by_key = dict(zip(residue_keys, flex_percentile.tolist()))
    fast_percentile_by_key = dict(zip(residue_keys, fast_percentile.tolist()))

    flex_percentile_pdb = output / (pdb_id + "_FlexPPI_percentile.pdb")
    fast_percentile_pdb = output / (pdb_id + "_FastProtFlex_percentile.pdb")
    flex_raw_pdb = output / (pdb_id + "_FlexPPI_raw_RMSF_A.pdb")
    fast_raw_pdb = output / (pdb_id + "_FastProtFlex_raw_score.pdb")
    _write_bfactor_pdb(
        pdb_path,
        flex_percentile_pdb,
        flex_percentile_by_key,
        (
            "B-factor = within-chain FlexPPI RMSF percentile (0 rigid, 100 flexible)",
            "FlexPPI sees the full multichain complex",
        ),
    )
    _write_bfactor_pdb(
        pdb_path,
        fast_percentile_pdb,
        fast_percentile_by_key,
        (
            "B-factor = within-chain FastProtFlex percentile (0 rigid, 100 flexible)",
            "FastProtFlex predicts each protein chain independently",
        ),
    )
    _write_bfactor_pdb(
        pdb_path,
        flex_raw_pdb,
        flex_by_key,
        ("B-factor = raw FlexPPI predicted RMSF in Angstrom",),
    )
    _write_bfactor_pdb(
        pdb_path,
        fast_raw_pdb,
        fast_by_key,
        ("B-factor = raw normalized FastProtFlex score; not Angstrom",),
    )

    residue_csv = output / (pdb_id + "_residue_predictions.csv")
    rows: List[Dict[str, Any]] = []
    for index, key in enumerate(residue_keys):
        chain, residue_number, insertion = key
        rows.append({
            "residue_index": index,
            "chain": chain,
            "chain_role": chain_roles[chain],
            "author_residue_id": residue_number + insertion,
            "amino_acid": structure.sequence[index],
            "flexppi_rmsf_angstrom": float(flex_rmsf[index]),
            "flexppi_rmsf_std_angstrom": float(flex_std[index]),
            "fastprotflex_normalized_score": float(fast_score[index]),
            "flexppi_within_chain_percentile": float(flex_percentile[index]),
            "fastprotflex_within_chain_percentile": float(fast_percentile[index]),
            "percentile_difference_flexppi_minus_fastprotflex": float(
                flex_percentile[index] - fast_percentile[index]
            ),
        })
    with residue_csv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    pymol_script = output / (pdb_id + "_compare_FlexPPI_vs_FastProtFlex.pml")
    session_path = output / (pdb_id + "_FlexPPI_vs_FastProtFlex.pse")
    _write_pymol_script(
        pymol_script,
        flex_percentile_pdb,
        fast_percentile_pdb,
        session_path,
        pdb_id,
    )
    metadata_path = output / (pdb_id + "_bundle_metadata.json")
    metadata: Dict[str, Any] = {
        "pdb_id": pdb_id,
        "input_pdb": str(Path(pdb_path).resolve()),
        "residue_count": structure.num_residues,
        "chains": [
            {
                "chain": chain,
                "role": chain_roles[chain],
                "residue_count": int(len(structure.indices_for_chain(chain))),
                "fastprotflex_end_to_end_seconds": float(
                    fast_timings[pdb_id + "__chain__" + chain]["end_to_end_seconds"]
                ),
            }
            for chain in chains
        ],
        "flexppi": {
            "checkpoint": flexppi_checkpoint,
            "device": str(device),
            "context": "full multichain complex",
            "temperature": temperature,
            "window_ns": window_ns,
            "dataset_code": dataset_code,
        },
        "fastprotflex": {"context": "each chain independently"},
        "visual_scale": {
            "field": "PDB B-factor",
            "value": "within-chain percentile",
            "minimum": 0,
            "maximum": 100,
            "colors": "blue-white-red",
            "warning": (
                "Percentiles compare relative flexibility within the same chain; they are not "
                "absolute Angstrom values and should not be compared across chains."
            ),
        },
        "artifacts": {
            "flexppi_percentile_pdb": str(flex_percentile_pdb.resolve()),
            "fastprotflex_percentile_pdb": str(fast_percentile_pdb.resolve()),
            "flexppi_raw_pdb": str(flex_raw_pdb.resolve()),
            "fastprotflex_raw_pdb": str(fast_raw_pdb.resolve()),
            "residue_predictions": str(residue_csv.resolve()),
            "pymol_script": str(pymol_script.resolve()),
            "pymol_session_after_running_script": str(session_path.resolve()),
            "metadata": str(metadata_path.resolve()),
        },
    }
    with metadata_path.open("w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2, ensure_ascii=False)
    return metadata
