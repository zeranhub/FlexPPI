"""Merge residue prediction shards with strict duplicate/coverage checks."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputs", nargs="+", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--expected-samples", type=int, required=True)
    args = parser.parse_args()

    rows = []
    fieldnames = None
    for input_text in args.inputs:
        path = Path(input_text)
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            if fieldnames is None:
                fieldnames = reader.fieldnames
            elif reader.fieldnames != fieldnames:
                raise ValueError(f"Shard headers differ: {path}")
            rows.extend(reader)
    keys = [(row["sample_id"], int(row["residue_index"])) for row in rows]
    if len(keys) != len(set(keys)):
        raise ValueError("Shards contain duplicate sample/residue keys")
    sample_count = len({sample_id for sample_id, _ in keys})
    if sample_count != args.expected_samples:
        raise ValueError(
            f"Expected {args.expected_samples} samples, found {sample_count}"
        )
    rows.sort(key=lambda row: (row["cohort"], row["sample_id"], int(row["residue_index"])))
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Merged {len(args.inputs)} shards: samples={sample_count}, rows={len(rows)}")


if __name__ == "__main__":
    main()
