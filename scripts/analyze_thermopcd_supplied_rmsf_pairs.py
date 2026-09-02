"""Screen ThermoPCD supplied RMSF labels for binding-response candidates.

This is a label-level audit, not a FlexPPI training run.  ThermoPCD's published
``gmx rmsf -res`` outputs are paired only when PDB id, temperature and random
seed agree.  Four manually curated systems have explicit complex and isolated
component descriptions plus RMSF on both sides.  Supplied RMSF is in nm and is
converted to Angstrom.  Because the original least-squares fitting group may
differ between a whole complex and an isolated component, these results are a
screen that should later be confirmed by target-chain-fitted trajectory RMSF.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


CURATED = {
    "5MV4": "citrullinated collagen-II epitope versus ACC1 Fab complex",
    "5OCX": "CII-C-13-CIT epitope versus ACPA E4 complex",
    "6BFS": "GM-CSF versus antibody complex (metadata antibody-name caveat)",
    "6SF6": "COMP P6 epitope versus mAb 15A complex",
}
PDB_RE = re.compile(r"(?i)(?<![A-Z0-9])([1-9][A-Z0-9]{3})(?![A-Z0-9])")
TEMP_RE = re.compile(r"(?i)(298|309|310|311(?:\.5)?|312|313)\s*K\b")
SEED_RE = re.compile(r"(?i)(?:random\s*seed\s*#?|\brs)[-_ #]*(\d+)")
DURATION_RE = re.compile(r"(?i)(\d+)\s*ns\b")


@dataclass
class Label:
    pdb_id: str
    kind: str
    temperature_k: float
    seed: int
    duration_ns: int | None
    directory: Path
    rmsf_path: Path
    pdb_path: Path | None
    values: dict[int, float]


def classify(text: str) -> str:
    lower = text.lower().replace("_", " ")
    return "complex" if re.search(r"\bcomplex\b|\bbound\b", lower) else "isolated"


def first_number(pattern: re.Pattern[str], text: str, cast):
    match = pattern.search(text)
    return cast(match.group(1)) if match else None


def parse_xvg(path: Path) -> dict[int, float]:
    values: dict[int, float] = {}
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith(("#", "@", "&")):
                continue
            fields = stripped.split()
            if len(fields) < 2:
                continue
            try:
                residue = int(round(float(fields[0])))
                rmsf_angstrom = float(fields[1]) * 10.0
            except ValueError:
                continue
            if math.isfinite(rmsf_angstrom):
                values[residue] = rmsf_angstrom
    return values


def choose_rmsf(directory: Path, pdb_id: str, temperature: float, seed: int) -> Path | None:
    candidates = [
        path
        for path in directory.iterdir()
        if path.is_file() and "rmsf" in path.name.lower() and path.suffix.lower() == ".xvg"
    ]
    if not candidates:
        return None
    temp_token = f"{temperature:g}k".lower()
    seed_tokens = (f"rs{seed}", f"rs-{seed}")

    def score(path: Path) -> tuple[int, int]:
        lower = path.name.lower()
        value = 0
        value += 4 if pdb_id.lower() in lower else 0
        value += 3 if temp_token in lower else 0
        value += 2 if any(token in lower for token in seed_tokens) else 0
        value -= 2 if re.search(r"rmsf[-_]1\.xvg$", lower) else 0
        return value, -len(path.name)

    return max(candidates, key=score)


def choose_pdb(directory: Path, temperature: float, seed: int) -> Path | None:
    candidates = [path for path in directory.iterdir() if path.is_file() and path.suffix.lower() == ".pdb"]
    if not candidates:
        return None
    temp_token = f"{temperature:g}k".lower()
    seed_token = f"rs{seed}"

    def score(path: Path) -> tuple[int, int]:
        lower = path.name.lower()
        value = 0
        value += 3 if "0ns" in lower else 0
        value += 2 if temp_token in lower else 0
        value += 2 if seed_token in lower else 0
        value -= 2 if "50ns" in lower or "100ns" in lower else 0
        return value, -len(path.name)

    return max(candidates, key=score)


def discover_labels(root: Path) -> list[Label]:
    labels: list[Label] = []
    seen_hashes: set[tuple[str, str, float, int, str]] = set()
    for collection in (root / "data", root / "ThermoPCD_labels"):
        if not collection.exists():
            continue
        for directory in collection.iterdir():
            if not directory.is_dir():
                continue
            pdb_match = PDB_RE.search(directory.name)
            if not pdb_match:
                continue
            pdb_id = pdb_match.group(1).upper()
            if pdb_id not in CURATED:
                continue
            temperature = first_number(TEMP_RE, directory.name, float)
            seed = first_number(SEED_RE, directory.name, int)
            duration = first_number(DURATION_RE, directory.name, int)
            if temperature is None or seed is None:
                continue
            rmsf_path = choose_rmsf(directory, pdb_id, temperature, seed)
            if rmsf_path is None:
                continue
            digest = hashlib.sha256(rmsf_path.read_bytes()).hexdigest()
            key = (pdb_id, classify(directory.name), temperature, seed, digest)
            if key in seen_hashes:
                continue
            seen_hashes.add(key)
            values = parse_xvg(rmsf_path)
            if not values:
                continue
            labels.append(
                Label(
                    pdb_id=pdb_id,
                    kind=classify(directory.name),
                    temperature_k=temperature,
                    seed=seed,
                    duration_ns=duration,
                    directory=directory,
                    rmsf_path=rmsf_path,
                    pdb_path=choose_pdb(directory, temperature, seed),
                    values=values,
                )
            )
    return labels


def parse_pdb(path: Path) -> dict[int, dict[str, object]]:
    residues: dict[int, dict[str, object]] = {}
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith(("ATOM  ", "HETATM")) or len(line) < 54:
                continue
            atom = line[12:16].strip()
            if atom.startswith("H"):
                continue
            try:
                resid = int(line[22:26])
                coord = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            except ValueError:
                continue
            item = residues.setdefault(
                resid,
                {"chain": line[21].strip(), "resname": line[17:20].strip(), "atoms": [], "ca": None},
            )
            item["atoms"].append(coord)
            if atom == "CA":
                item["ca"] = coord
    return residues


def interface_distances(pdb_path: Path | None, target_residues: set[int]) -> dict[int, tuple[float, float]]:
    if pdb_path is None:
        return {}
    residues = parse_pdb(pdb_path)
    partner = [item for resid, item in residues.items() if resid not in target_residues]
    if not partner:
        return {}
    partner_atoms = np.vstack([np.vstack(item["atoms"]) for item in partner if item["atoms"]])
    partner_ca = np.vstack([item["ca"] for item in partner if item["ca"] is not None])
    result: dict[int, tuple[float, float]] = {}
    for resid in target_residues:
        item = residues.get(resid)
        if item is None or not item["atoms"]:
            continue
        atoms = np.vstack(item["atoms"])
        heavy = float(np.sqrt(((atoms[:, None, :] - partner_atoms[None, :, :]) ** 2).sum(axis=2)).min())
        if item["ca"] is None or partner_ca.size == 0:
            ca = math.nan
        else:
            ca = float(np.sqrt(((partner_ca - item["ca"]) ** 2).sum(axis=1)).min())
        result[resid] = heavy, ca
    return result


def match_labels(labels: Iterable[Label]) -> tuple[pd.DataFrame, pd.DataFrame]:
    labels = list(labels)
    rows: list[dict[str, object]] = []
    qc_rows: list[dict[str, object]] = []
    for pdb_id in CURATED:
        apo_labels = [label for label in labels if label.pdb_id == pdb_id and label.kind == "isolated"]
        holo_labels = [label for label in labels if label.pdb_id == pdb_id and label.kind == "complex"]
        for apo in apo_labels:
            candidates = [
                holo
                for holo in holo_labels
                if holo.temperature_k == apo.temperature_k and holo.seed == apo.seed
            ]
            if not candidates:
                continue

            def pair_score(holo: Label) -> tuple[int, int, int]:
                common = len(set(apo.values) & set(holo.values))
                duration_gap = abs((holo.duration_ns or 0) - (apo.duration_ns or 0))
                exact_duration = int(holo.duration_ns == apo.duration_ns)
                return exact_duration, common, -duration_gap

            holo = max(candidates, key=pair_score)
            common = sorted(set(apo.values) & set(holo.values))
            distances = interface_distances(holo.pdb_path, set(common))
            qc_rows.append(
                {
                    "pdb_id": pdb_id,
                    "temperature_k": apo.temperature_k,
                    "seed": apo.seed,
                    "apo_duration_ns": apo.duration_ns,
                    "holo_duration_ns": holo.duration_ns,
                    "duration_matched": apo.duration_ns == holo.duration_ns,
                    "apo_residues": len(apo.values),
                    "holo_residues": len(holo.values),
                    "common_residues": len(common),
                    "apo_coverage": len(common) / len(apo.values),
                    "structure_mapped_residues": len(distances),
                    "apo_rmsf_path": str(apo.rmsf_path),
                    "holo_rmsf_path": str(holo.rmsf_path),
                    "holo_pdb_path": str(holo.pdb_path) if holo.pdb_path else None,
                }
            )
            for resid in common:
                heavy, ca = distances.get(resid, (math.nan, math.nan))
                apo_value = apo.values[resid]
                holo_value = holo.values[resid]
                rows.append(
                    {
                        "pdb_id": pdb_id,
                        "temperature_k": apo.temperature_k,
                        "seed": apo.seed,
                        "residue_index": resid,
                        "apo_rmsf_angstrom": apo_value,
                        "holo_rmsf_angstrom": holo_value,
                        "delta_rmsf_angstrom": holo_value - apo_value,
                        "min_partner_heavy_distance_angstrom": heavy,
                        "min_partner_ca_distance_angstrom": ca,
                        "interface_heavy5": bool(heavy <= 5.0) if math.isfinite(heavy) else None,
                        "interface_calpha8": bool(ca <= 8.0) if math.isfinite(ca) else None,
                        "duration_matched": apo.duration_ns == holo.duration_ns,
                    }
                )
    return pd.DataFrame(rows), pd.DataFrame(qc_rows)


def summarize(residue: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, object]]:
    comparison_rows: list[dict[str, object]] = []
    if residue.empty:
        return pd.DataFrame(), pd.DataFrame(), {}
    for keys, group in residue.groupby(["pdb_id", "temperature_k", "seed", "duration_matched"]):
        for region, subset in (
            ("all", group),
            ("interface_calpha8", group.loc[group["interface_calpha8"] == True]),  # noqa: E712
            ("noninterface_calpha8", group.loc[group["interface_calpha8"] == False]),  # noqa: E712
        ):
            if subset.empty:
                continue
            delta = subset["delta_rmsf_angstrom"].to_numpy(float)
            comparison_rows.append(
                {
                    "pdb_id": keys[0],
                    "temperature_k": keys[1],
                    "seed": keys[2],
                    "duration_matched": keys[3],
                    "region": region,
                    "residue_count": len(delta),
                    "mean_delta_rmsf_angstrom": float(np.mean(delta)),
                    "median_delta_rmsf_angstrom": float(np.median(delta)),
                    "negative_fraction": float(np.mean(delta < 0)),
                    "positive_fraction": float(np.mean(delta > 0)),
                }
            )
    comparisons = pd.DataFrame(comparison_rows)
    pdb_rows: list[dict[str, object]] = []
    for (pdb_id, region), group in comparisons.groupby(["pdb_id", "region"]):
        pdb_rows.append(
            {
                "pdb_id": pdb_id,
                "region": region,
                "matched_state_comparisons": len(group),
                "mean_delta_rmsf_angstrom": float(group["mean_delta_rmsf_angstrom"].mean()),
                "mean_negative_fraction": float(group["negative_fraction"].mean()),
                "comparisons_with_negative_mean_fraction": float(
                    (group["mean_delta_rmsf_angstrom"] < 0).mean()
                ),
            }
        )
    per_pdb = pd.DataFrame(pdb_rows)
    interface = per_pdb.loc[per_pdb["region"] == "interface_calpha8"]
    overall = {
        "candidate_systems_with_supplied_pair_labels": int(residue["pdb_id"].nunique()),
        "matched_state_comparisons": int(
            residue[["pdb_id", "temperature_k", "seed"]].drop_duplicates().shape[0]
        ),
        "mapped_residue_observations": int(len(residue)),
        "interface_residue_observations_calpha8": int((residue["interface_calpha8"] == True).sum()),  # noqa: E712
        "pdb_equal_interface_mean_delta_rmsf_angstrom": (
            float(interface["mean_delta_rmsf_angstrom"].mean()) if not interface.empty else None
        ),
        "pdb_equal_interface_negative_fraction": (
            float(interface["mean_negative_fraction"].mean()) if not interface.empty else None
        ),
        "pdb_equal_fraction_systems_negative_interface_mean": (
            float((interface["mean_delta_rmsf_angstrom"] < 0).mean()) if not interface.empty else None
        ),
        "scope": "screening only: source-provided gmx rmsf labels, not consistently target-chain-refitted trajectories",
    }
    return comparisons, per_pdb, overall


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(r"D:\data\ThermoPCD_trajectories"))
    parser.add_argument("--output-dir", type=Path, default=Path("artifacts/thermopcd_audit"))
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    labels = discover_labels(args.root)
    residue, qc = match_labels(labels)
    comparisons, per_pdb, summary = summarize(residue)
    residue.to_csv(args.output_dir / "supplied_rmsf_pair_residues.csv", index=False)
    qc.to_csv(args.output_dir / "supplied_rmsf_pair_qc.csv", index=False)
    comparisons.to_csv(args.output_dir / "supplied_rmsf_pair_comparisons.csv", index=False)
    per_pdb.to_csv(args.output_dir / "supplied_rmsf_pair_pdb_summary.csv", index=False)
    (args.output_dir / "supplied_rmsf_pair_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
