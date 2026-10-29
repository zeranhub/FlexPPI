"""Prepare auditable A6 TCR inputs and residue mappings, never simulated labels.

This script retains deposited backbones. A sequence-reverted model is explicitly
labelled as a fixed-backbone construction, not a relaxed wild-type structure.
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "artifacts/mechanistic_case_studies_20260913/tcr"
SOURCE = ROOT / "artifacts/apo_holo_mechanistic_literature_v2_20260913/sources"
sys.path.insert(0, str(ROOT / "src"))
from flexppi.constants import THREE_TO_ONE
from flexppi.data.structures import parse_pdb


def dump(name, value):
    path = OUT / name
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2), encoding="utf-8")


def fetch_sources():
    (OUT / "sources").mkdir(parents=True, exist_ok=True)
    manifest = []
    for entry in (27440, 27441):
        path = OUT / "sources" / f"bmr{entry}_3.str"
        url = f"https://bmrb.io/ftp/pub/bmrb/entry_directories/bmr{entry}/bmr{entry}_3.str"
        if not path.exists():
            response = requests.get(url, timeout=(10, 40))
            response.raise_for_status()
            if not response.text.startswith("data_"):
                raise ValueError("BMRB response is not an NMR-STAR entry")
            path.write_bytes(response.content)
        manifest.append(dict(path=str(path), url=url, sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    dump("source_manifest.json", manifest)


def inspect_sources():
    for pdb in ("3QH3", "1QRN", "4FTV"):
        s = parse_pdb(str(SOURCE / f"{pdb}.pdb"))
        print(pdb, s.num_residues)
        for chain in s.unique_chains():
            indices = s.indices_for_chain(chain)
            print(chain, len(indices), s.chain_sequence(chain))
            if (pdb == "3QH3" and chain in "BD") or (pdb != "3QH3" and chain == "E"):
                print("beta 210-233", [(s.author_residue_ids[i], s.sequence[i]) for i in indices if 210 <= int(s.author_residue_ids[i].rstrip("ABCDEFGHIJKLMNOPQRSTUVWXYZ")) <= 233])


if __name__ == "__main__":
    fetch_sources()
    inspect_sources()
