"""Run the unchanged ATLAS extractor with per-trajectory local temporary copies.

On a rotating source disk, many concurrent small XTC reads cause seeking.
Sequentially copying one XTC to the Windows temporary SSD before MDAnalysis
opens it changes only I/O. At most one trajectory per worker is staged; the
copy is removed as soon as that trajectory closes. Labels and fingerprints
still come from the original scientifically tested extraction implementation.
"""
import hashlib
import os
from pathlib import Path
import shutil
import sys
import tempfile

sys.path.insert(0, str(Path(__file__).resolve().parent))
import backfill_atlas_hchi as core

_extract_replica = core.extract_replica


def buffered_extract(pdb_path, trajectory_path):
    original = Path(trajectory_path)
    with tempfile.TemporaryDirectory(prefix="flexppi_atlas_hchi_") as directory:
        copied = Path(directory) / original.name
        digest = hashlib.sha256()
        # Read once sequentially, instead of issuing interleaved frame reads.
        with original.open("rb") as source, copied.open("wb") as target:
            for chunk in iter(lambda: source.read(4 * 1024 * 1024), b""):
                digest.update(chunk)
                target.write(chunk)
        if copied.stat().st_size != original.stat().st_size:
            raise ValueError("Incomplete trajectory scratch copy")
        # Verify the copy from the fast scratch disk, not by rereading the HDD.
        if core.sha256(copied) != digest.hexdigest():
            raise ValueError("Trajectory scratch copy hash mismatch")
        rows, qa = _extract_replica(pdb_path, copied)
        qa["trajectory_sha256"] = digest.hexdigest()
        qa["io_mode"] = "verified_per_trajectory_temporary_copy"
        return rows, qa


core.extract_replica = buffered_extract

if __name__ == "__main__":
    raise SystemExit(core.main())
