"""Reproduce an accepted figure in a separate output directory.

The accepted PNG/PDF/SVG/Draw.io masters are never overwritten by this entry.
Python plotters run against temporary copies of their local source bundles.
"""
from pathlib import Path
import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parent


def check_masters():
    manifest = json.loads((ROOT / "manifest.json").read_text(encoding="utf-8"))
    for item in manifest.get("shared_files", []):
        path = ROOT / item["path"]
        if hashlib.sha256(path.read_bytes()).hexdigest() != item["sha256"]:
            raise ValueError(f"Shared style changed: {path}")
    for figure in manifest["figures"]:
        for item in figure["accepted_files"]:
            source = ROOT / figure["folder"] / item["path"]
            actual = hashlib.sha256(source.read_bytes()).hexdigest()
            if actual != item["sha256"]:
                raise ValueError(f"Master changed: {source}")
        for item in figure.get("source_files", []):
            source = ROOT / figure["folder"] / item["path"]
            if hashlib.sha256(source.read_bytes()).hexdigest() != item["sha256"]:
                raise ValueError(f"Source changed: {source}")
    print("PASS: all six accepted figure packages match the manifest.")


def reproduce(figure, destination):
    source = ROOT / figure["folder"]
    destination.mkdir(parents=True, exist_ok=True)
    log = destination / "reproduce.log"
    if figure["figure"] == "1":
        candidates = [
            shutil.which("draw.io"), shutil.which("drawio"),
            str(Path(os.environ.get("ProgramFiles", "C:/Program Files")) / "draw.io/draw.io.exe"),
            str(Path.home() / "bin/draw.io.exe"),
        ]
        executable = next((p for p in candidates if p and Path(p).is_file()), None)
        if not executable:
            raise RuntimeError("Install Draw.io desktop or open fig1_framework.drawio and export from its File menu.")
        output = []
        for extension in ["pdf", "svg", "png"]:
            command = [executable, "--export", "--format", extension, "--output",
                       str(destination / (figure["stem"] + "." + extension)),
                       str(source / "fig1_framework.drawio")]
            result = subprocess.run(command, capture_output=True, text=True,
                                    encoding="utf-8", errors="replace")
            output.append(result.stdout + result.stderr)
            if result.returncode:
                log.write_text("\n".join(output), encoding="utf-8")
                raise RuntimeError(f"Draw.io export failed; see {log}")
        log.write_text("\n".join(output), encoding="utf-8")
    else:
        with tempfile.TemporaryDirectory(prefix="flexppi_figure_") as temporary:
            copied = Path(temporary) / source.name
            shutil.copytree(source, copied, ignore=shutil.ignore_patterns("__pycache__", "reproduced"))
            # Every v59 panel imports the same semantic colour and scale mapping.
            shutil.copy2(ROOT / "publication_style.py", Path(temporary) / "publication_style.py")
            command = [sys.executable, "-X", "utf8", str(copied / figure["plot_entry"])]
            generated = copied
            if figure["figure"] == "3":
                generated = copied / "reproduced"
                command += ["--output-dir", str(generated)]
            result = subprocess.run(command, cwd=copied, capture_output=True,
                                    text=True, encoding="utf-8", errors="replace")
            log.write_text(result.stdout + result.stderr, encoding="utf-8")
            if result.returncode:
                raise RuntimeError(f"Figure {figure['figure']} failed; see {log}")
            for extension in ["pdf", "svg", "png"]:
                path = generated / (figure["stem"] + "." + extension)
                if not path.is_file():
                    raise FileNotFoundError(path)
                shutil.copy2(path, destination / path.name)
    print(f"Figure {figure['figure']}: {destination}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("figure", nargs="?", choices=["1", "2", "3", "4", "5", "S1", "all"])
    parser.add_argument("--output-dir", type=Path, default=ROOT / "reproduced")
    parser.add_argument("--check", "-check", action="store_true", help="Only verify accepted source/export hashes.")
    args = parser.parse_args()
    check_masters()
    if args.check:
        return
    if args.figure is None:
        parser.error("provide a figure number, all, or --check")
    manifest = json.loads((ROOT / "manifest.json").read_text(encoding="utf-8"))
    for figure in manifest["figures"]:
        if args.figure in ["all", figure["figure"]]:
            reproduce(figure, args.output_dir.resolve() / figure["folder"])
    check_masters()


if __name__ == "__main__":
    main()
