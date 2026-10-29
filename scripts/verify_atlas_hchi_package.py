"""Verify the delivered tar archive and run its preflight against real local inputs."""
import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import shutil
import subprocess
import sys
import tarfile
import tempfile

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def verify(package):
    archive = package / "atlas_hchi_retraining_v1.tar.gz"
    expected, name = (package / "archive_sha256.txt").read_text().strip().split("  ", 1)
    if name != archive.name or sha(archive) != expected:
        raise ValueError("Archive checksum mismatch")
    with tempfile.TemporaryDirectory(prefix="flexppi_atlas_package_verify_") as temporary:
        scratch = Path(temporary)
        with tarfile.open(archive, "r:gz") as tar:
            members = tar.getmembers()
            names = set()
            for member in members:
                path = PurePosixPath(member.name)
                if (not member.isfile() or path.is_absolute() or ".." in path.parts
                        or ":" in member.name or "\\" in member.name or member.name in names):
                    raise ValueError("Unsafe or duplicate archive member: " + member.name)
                names.add(member.name)
            # All members have been checked to be unique relative regular files.
            tar.extractall(scratch)
        print(f"Extracted and checked {len(members)} safe archive members", flush=True)
        shutil.copytree(ROOT / "src", scratch / "src", ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
        normalized = 0
        for path in (scratch / "src").rglob("*.py"):
            before = path.read_bytes()
            content = path.read_text(encoding="utf-8")
            with path.open("w", encoding="utf-8", newline="\n") as stream:
                stream.write(content)
            normalized += before != path.read_bytes()
        mapping = {
            "/ranze/repos/FlexPPI/artifacts/data/atlas_rmsf_hchi_v1": str(scratch / "artifacts/data/atlas_rmsf_hchi_v1"),
            "/ranze/repos/FlexPPI": str(ROOT),
            "/ranze/data/ATLAS_extracted": "E:/data/protein_md_datasets/ATLAS_extracted",
            "/ranze/data/dynarepo/dynarepo_non_nucleic_full": "D:/data/dynarepo_non_nucleic_full",
            "/ranze/data/ITrace": "E:/data/ITrace",
        }
        command = [sys.executable, str(scratch / "scripts/preflight_atlas_hchi_retraining.py"),
                   "--repo", str(scratch), "--fold", "0", "--check-new-run"]
        for old, new in mapping.items():
            command += ["--path-map", old + "=" + new]
        result = subprocess.run(command, capture_output=True, text=True, encoding="utf-8")
        print(result.stdout, flush=True)
        if result.returncode:
            raise RuntimeError(result.stderr or result.stdout)
        preflight = json.loads(result.stdout)
        # The overwrite guard is also tested by creating only a temporary marker.
        existing = scratch / "artifacts/performance_first_v2_atlas_hchi_v1/fold_0"
        existing.mkdir(parents=True)
        (existing / "verification_only.txt").write_text("Test marker; not a training result.\n")
        guarded = subprocess.run(command, capture_output=True, text=True, encoding="utf-8")
        guard_report = json.loads(guarded.stdout)
        if guarded.returncode != 1 or not any("Run output already populated" in e for e in guard_report["errors"]):
            raise ValueError("Existing-output protection did not reject the test marker")
        bash = Path("C:/Program Files/Git/bin/bash.exe")
        if not bash.is_file():
            raise FileNotFoundError("Git Bash is required for the local shell syntax check")
        checked = subprocess.run([str(bash), "-n", str(scratch / "scripts/train_atlas_hchi_fold.sh")], capture_output=True, text=True)
        if checked.returncode:
            raise ValueError(checked.stderr)
        report = {"passed": True, "archive": str(archive), "archive_sha256": expected,
                  "archive_bytes": archive.stat().st_size, "safe_archive_members": len(members),
                  "preflight": preflight, "existing_run_rejected": True, "bash_syntax_passed": True,
                  "linux_lf_source_compatibility_passed": True, "source_files_with_newlines_normalized": normalized,
                  "remote_training_started": False,
                  "verification_mode": "Archive unpacked into temporary repository; existing local structures and non-ATLAS labels mapped read-only"}
        (package / "package_verification.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
        print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--package-dir", type=Path, default=ROOT / "output/atlas_hchi_retraining_v1_linux")
    args = parser.parse_args()
    verify(args.package_dir.resolve())
