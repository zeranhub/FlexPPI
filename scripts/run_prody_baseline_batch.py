"""Run ProDy ANM and GNM over a frozen RMSF baseline manifest."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
import time
from typing import Dict, List

import numpy as np
from prody import ANM, GNM, calcSqFlucts, confProDy, parsePDB


OUTPUT_FIELDS = (
    "cohort",
    "sample_id",
    "residue_index",
    "chain",
    "author_residue_id",
    "residue_number",
    "insertion_code",
    "amino_acid_3",
    "anm_sq_fluct",
    "gnm_sq_fluct",
    "anm_seconds",
    "gnm_seconds",
    "mode_count",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--modes", type=int, default=20)
    parser.add_argument("--cutoff", type=float, default=15.0)
    parser.add_argument("--gamma", type=float, default=1.0)
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--shard-count", type=int, default=1)
    return parser.parse_args()


def read_manifest(path: Path, limit: int | None) -> List[Dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return rows if limit is None else rows[:limit]


def completed_samples(path: Path) -> set[str]:
    if not path.is_file() or path.stat().st_size == 0:
        return set()
    with path.open("r", encoding="utf-8", newline="") as handle:
        return {row["sample_id"] for row in csv.DictReader(handle)}


def main() -> None:
    args = parse_args()
    confProDy(verbosity="none")
    if args.shard_count < 1 or not 0 <= args.shard_index < args.shard_count:
        raise ValueError("Require 0 <= shard-index < shard-count")
    manifest = read_manifest(Path(args.manifest), args.limit)
    manifest = manifest[args.shard_index :: args.shard_count]
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    done = completed_samples(output)
    append = output.is_file() and output.stat().st_size > 0
    with output.open("a", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_FIELDS)
        if not append:
            writer.writeheader()
        for position, record in enumerate(manifest, start=1):
            sample_id = record["sample_id"]
            if sample_id in done:
                continue
            atoms = parsePDB(record["container_pdb_path"])
            calphas = None if atoms is None else atoms.select("protein and name CA")
            if calphas is None or calphas.numAtoms() < 4:
                raise ValueError(f"{sample_id}: fewer than four protein CA atoms")
            count = int(calphas.numAtoms())
            mode_count = max(1, min(args.modes, count - 1, 3 * count - 6))

            started = time.perf_counter()
            anm = ANM(sample_id + "_anm")
            anm.buildHessian(calphas, cutoff=args.cutoff, gamma=args.gamma)
            anm.calcModes(n_modes=mode_count)
            anm_values = np.asarray(calcSqFlucts(anm), dtype=np.float64)
            anm_seconds = time.perf_counter() - started

            started = time.perf_counter()
            gnm = GNM(sample_id + "_gnm")
            gnm.buildKirchhoff(calphas, cutoff=args.cutoff, gamma=args.gamma)
            gnm.calcModes(n_modes=mode_count)
            gnm_values = np.asarray(calcSqFlucts(gnm), dtype=np.float64)
            gnm_seconds = time.perf_counter() - started
            if not np.all(np.isfinite(anm_values)) or not np.all(np.isfinite(gnm_values)):
                raise ValueError(f"{sample_id}: ProDy returned non-finite values")

            for residue_index in range(count):
                chain = str(calphas.getChids()[residue_index]).strip() or "_"
                number = str(int(calphas.getResnums()[residue_index]))
                insertion = str(calphas.getIcodes()[residue_index]).strip()
                writer.writerow(
                    {
                        "cohort": record["cohort"],
                        "sample_id": sample_id,
                        "residue_index": residue_index,
                        "chain": chain,
                        "author_residue_id": number + insertion,
                        "residue_number": number,
                        "insertion_code": insertion,
                        "amino_acid_3": str(calphas.getResnames()[residue_index]).upper(),
                        "anm_sq_fluct": float(anm_values[residue_index]),
                        "gnm_sq_fluct": float(gnm_values[residue_index]),
                        "anm_seconds": anm_seconds,
                        "gnm_seconds": gnm_seconds,
                        "mode_count": mode_count,
                    }
                )
            handle.flush()
            print(
                f"[{position}/{len(manifest)}] {sample_id}: {count} residues, "
                f"ANM {anm_seconds:.3f} s, GNM {gnm_seconds:.3f} s",
                flush=True,
            )


if __name__ == "__main__":
    main()
