"""Make the five train/test manifest pairs from local structures and MD labels."""

import argparse
import csv
import hashlib
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--structures", required=True, type=Path, help="Directory of <sample_id>.pdb files")
    parser.add_argument("--labels", required=True, type=Path, help="Directory of <sample_id>.npz files")
    parser.add_argument("--scheme", choices=("system", "family"), default="system")
    parser.add_argument("--output", type=Path, default=Path("data/splits"), help="Parent of the scheme directory")
    parser.add_argument("--verify", action="store_true", help="Verify every PDB and NPZ against the study SHA-256 checksums")
    args = parser.parse_args()

    membership = Path(__file__).with_name("fold_membership.csv")
    checksums = {}
    if args.verify:
        with membership.with_name("input_checksums.csv").open(encoding="utf-8", newline="") as handle:
            checksums = {row["sample_id"]: row for row in csv.DictReader(handle)}
    records = []
    with membership.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            sample_id = row["sample_id"]
            structure = (args.structures / (sample_id + ".pdb")).resolve()
            label = (args.labels / (sample_id + ".npz")).resolve()
            for path, key in ((structure, "structure_sha256"), (label, "label_sha256")):
                if not path.is_file():
                    parser.error(f"Missing input for {sample_id}: {path}")
                if args.verify:
                    expected = checksums.get(sample_id, {}).get(key)
                    if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
                        parser.error(f"Study checksum mismatch for {sample_id}: {path}")
            record = {
                "sample_id": sample_id,
                "dataset": row["dataset"],
                "structure_path": str(structure),
                "label_path": str(label),
                "source_group": row["source_group"],
                "target_chains": row["target_chains"].split(",") if row["target_chains"] else [],
                "metadata": {"pdb_id": row["pdb_id"], "biological_system": "pdb:" + row["pdb_id"]},
            }
            fold = int(row[args.scheme + "_fold"])
            records.append((fold, json.dumps(record)))

    destination = args.output / args.scheme
    destination.mkdir(parents=True, exist_ok=True)
    for fold in range(5):
        partitions = {"train": [], "test": []}
        for held_out, line in records:
            split = "test" if held_out == fold else "train"
            partitions[split].append(line)
        for split, lines in partitions.items():
            path = destination / f"fold_{fold}_{split}.jsonl"
            with path.open("w", encoding="utf-8", newline="\n") as handle:
                handle.writelines(line + "\n" for line in lines)
            print(f"{path}: {len(lines)} records")


if __name__ == "__main__":
    main()
