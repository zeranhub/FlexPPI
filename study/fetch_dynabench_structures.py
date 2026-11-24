"""Restore DynaBench structure inputs from the Inria DynaRepo node."""

import argparse
import csv
import hashlib
import http.client
import os
import sys
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import urlopen


API = "https://dynarepo.inria.fr/api/rest/current/projects"
STUDY = Path(__file__).resolve().parent
HISTIDINES = {"HSD", "HSE", "HSP", "HID", "HIE", "HIP"}


def read_inputs():
    with (STUDY / "input_checksums.csv").open(encoding="utf-8", newline="") as handle:
        checksums = {row["sample_id"]: row["structure_sha256"] for row in csv.DictReader(handle)}
    with (STUDY / "fold_membership.csv").open(encoding="utf-8", newline="") as handle:
        return {row["sample_id"]: checksums[row["sample_id"]] for row in csv.DictReader(handle)
                if row["dataset"] == "dynarepo" and row["source_group"] == "dynabench_complex"}


def canonical_pdb(source):
    # The reference input files use CRLF on every platform. Preserve column spacing.
    text = source.decode("utf-8").replace("\r\n", "\n").replace("\r", "\n")
    lines = text.splitlines(keepends=True)
    for index, line in enumerate(lines):
        if len(line) < 54 or not line.startswith(("ATOM  ", "HETATM")):
            continue
        residue = line[17:20]
        if residue in HISTIDINES:
            lines[index] = line[:17] + "HIS" + line[20:]
        elif residue == "ILE" and line[12:16].strip() == "CD":
            lines[index] = line[:12] + " CD1" + line[16:]
    return "".join(lines).replace("\n", "\r\n").encode("utf-8")


def download_structure(sample_id):
    url = f"{API}/{sample_id}/files/structure.pdb"
    for attempt in range(3):
        try:
            with urlopen(url, timeout=45) as response:
                return response.read()
        except (URLError, TimeoutError, http.client.IncompleteRead) as error:
            if isinstance(error, HTTPError) and error.code < 500 and error.code != 429:
                raise
            if attempt == 2:
                raise
            time.sleep(2 ** attempt)


def restore(sample_id, expected, directory):
    target = directory / f"{sample_id}.pdb"
    if target.is_file() and hashlib.sha256(target.read_bytes()).hexdigest() == expected:
        return
    content = canonical_pdb(download_structure(sample_id))
    if hashlib.sha256(content).hexdigest() != expected:
        raise ValueError("upstream structure differs from the reference checksum")
    with tempfile.NamedTemporaryFile(dir=directory, delete=False) as handle:
        temporary = Path(handle.name)
    try:
        temporary.write_bytes(content)
        os.replace(temporary, target)
    finally:
        temporary.unlink(missing_ok=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--structures", type=Path, default=Path("data/structures"),
                        help="Directory for <sample_id>.pdb files (default: data/structures)")
    args = parser.parse_args()
    args.structures.mkdir(exist_ok=True, parents=True)
    inputs = read_inputs()
    failures = 0
    with ThreadPoolExecutor(max_workers=3) as pool:
        futures = {pool.submit(restore, sample, checksum, args.structures): sample
                   for sample, checksum in inputs.items()}
        for done, future in enumerate(as_completed(futures), 1):
            try:
                future.result()
            except (OSError, ValueError, http.client.HTTPException) as error:
                failures += 1
                print(f"{futures[future]}: {error}", file=sys.stderr, flush=True)
            if done % 25 == 0 or done == len(inputs):
                print(f"Verified {done - failures}/{len(inputs)} structures", flush=True)
    return int(bool(failures))


if __name__ == "__main__":
    sys.exit(main())
