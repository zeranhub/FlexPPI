"""Run one BackFlip checkpoint over a frozen baseline manifest.

This file intentionally depends only on the Python standard library plus the
BackFlip installation, because it is executed inside the official baseline
container rather than the FlexPPI environment.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
import tempfile
import time
from typing import Dict, List, Sequence, Tuple

import numpy as np

from backflip.deployment.inference_class import BackFlip


OUTPUT_FIELDS = (
    "cohort",
    "sample_id",
    "residue_index",
    "chain",
    "author_residue_id",
    "residue_number",
    "insertion_code",
    "amino_acid_3",
    "backflip_rmsf_angstrom",
    "checkpoint_tag",
    "prediction_context",
    "inference_seconds",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--tag", default="backflip-1.0-seq")
    parser.add_argument("--limit", type=int)
    parser.add_argument("--device", default="cpu")
    parser.add_argument(
        "--chain-wise",
        action="store_true",
        help="Predict each PDB chain independently (required for complexes).",
    )
    return parser.parse_args()


def read_manifest(path: Path, limit: int | None) -> List[Dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return rows if limit is None else rows[:limit]


def residue_rows(path: str) -> List[Tuple[str, str, str, str, str]]:
    """Return unique protein CA residues in file order.

    A blank chain is normalized to ``_`` exactly as in the frozen FlexPPI
    common-residue files.  Alternative conformer A is retained; other
    conformers are ignored.
    """

    rows: List[Tuple[str, str, str, str, str]] = []
    seen = set()
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith("ATOM  ") or len(line) < 27:
                continue
            if line[12:16].strip() != "CA" or line[16:17] not in (" ", "A"):
                continue
            chain = line[21:22].strip() or "_"
            number = line[22:26].strip()
            insertion = line[26:27].strip()
            residue = line[17:20].strip().upper()
            key = (chain, number, insertion)
            if key in seen:
                continue
            seen.add(key)
            rows.append((chain, number + insertion, number, insertion, residue))
    return rows


def completed_samples(path: Path, tag: str) -> set[str]:
    if not path.is_file() or path.stat().st_size == 0:
        return set()
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    counts: Dict[str, int] = {}
    for row in rows:
        if row.get("checkpoint_tag") == tag:
            counts[row["sample_id"]] = counts.get(row["sample_id"], 0) + 1
    return {sample_id for sample_id, count in counts.items() if count > 0}


def prediction_array(result: Dict[str, object]) -> np.ndarray:
    values = result["global_rmsf"]
    if hasattr(values, "detach"):
        values = values.detach().cpu().numpy()
    return np.asarray(values, dtype=np.float64).squeeze().reshape(-1)


def write_chain_pdb(source: str, destination: Path, chain: str) -> None:
    chain_character = " " if chain == "_" else chain
    with open(source, "r", encoding="utf-8", errors="replace") as source_handle:
        lines = [
            line
            for line in source_handle
            if line.startswith("ATOM  ")
            and len(line) >= 22
            and line[21:22] == chain_character
        ]
    if not lines:
        raise ValueError(f"No ATOM records found for chain {chain} in {source}")
    with destination.open("w", encoding="utf-8", newline="\n") as output_handle:
        output_handle.writelines(lines)
        output_handle.write("TER\nEND\n")


def predict(
    model: BackFlip,
    path: str,
    residues: Sequence[Tuple[str, str, str, str, str]],
    chain_wise: bool,
) -> np.ndarray:
    if not chain_wise:
        return prediction_array(model.predict_from_pdb(path))
    output = np.full(len(residues), np.nan, dtype=np.float64)
    chains = list(dict.fromkeys(identity[0] for identity in residues))
    with tempfile.TemporaryDirectory(prefix="backflip_chains_") as directory:
        for chain in chains:
            indices = [
                index for index, identity in enumerate(residues) if identity[0] == chain
            ]
            chain_path = Path(directory) / f"chain_{chain}.pdb"
            write_chain_pdb(path, chain_path, chain)
            values = prediction_array(model.predict_from_pdb(str(chain_path)))
            if len(values) != len(indices):
                raise ValueError(
                    f"chain {chain}: BackFlip returned {len(values)} values "
                    f"for {len(indices)} residues"
                )
            output[indices] = values
    return output


def main() -> None:
    args = parse_args()
    manifest = read_manifest(Path(args.manifest), args.limit)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    done = completed_samples(output, args.tag)
    append = output.is_file() and output.stat().st_size > 0
    model = BackFlip.from_tag(
        tag=args.tag, device=args.device, progress_bar=False
    )
    with output.open("a", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_FIELDS)
        if not append:
            writer.writeheader()
        for position, record in enumerate(manifest, start=1):
            sample_id = record["sample_id"]
            if sample_id in done:
                continue
            path = record["container_pdb_path"]
            residues = residue_rows(path)
            started = time.perf_counter()
            prediction = predict(model, path, residues, args.chain_wise)
            elapsed = time.perf_counter() - started
            if len(prediction) != len(residues):
                raise ValueError(
                    f"{sample_id}: BackFlip returned {len(prediction)} values "
                    f"but the PDB contains {len(residues)} unique CA residues"
                )
            if not np.all(np.isfinite(prediction)):
                raise ValueError(f"{sample_id}: BackFlip returned non-finite values")
            for residue_index, (identity, value) in enumerate(zip(residues, prediction)):
                chain, author_id, number, insertion, residue = identity
                writer.writerow(
                    {
                        "cohort": record["cohort"],
                        "sample_id": sample_id,
                        "residue_index": residue_index,
                        "chain": chain,
                        "author_residue_id": author_id,
                        "residue_number": number,
                        "insertion_code": insertion,
                        "amino_acid_3": residue,
                        "backflip_rmsf_angstrom": float(value),
                        "checkpoint_tag": args.tag,
                        "prediction_context": (
                            "independent_chain" if args.chain_wise else "whole_pdb"
                        ),
                        "inference_seconds": elapsed,
                    }
                )
            handle.flush()
            print(
                f"[{position}/{len(manifest)}] {sample_id}: "
                f"{len(residues)} residues, {elapsed:.3f} s",
                flush=True,
            )


if __name__ == "__main__":
    main()
