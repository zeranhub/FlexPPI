"""Audit the local ThermoPCD download without modifying FlexPPI.

The script separates Dataverse catalogue entries from files that are really
present on disk.  It also distinguishes continuous production trajectories
from isolated XTC snapshots and identifies PDB identifiers for which both a
complex and an isolated component appear to exist.  The latter are candidates
for a binding-response analysis, not automatically strict apo/holo pairs.
"""

from __future__ import annotations

import argparse
import json
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable

import pandas as pd


PDB_EXPLICIT = re.compile(r"(?i)\bPDB\s*[-_:,]?\s*([0-9][A-Z0-9]{3})\b")
PDB_LEADING = re.compile(r"(?i)^(?:DB\s+)?([0-9][A-Z0-9]{3})(?:\b|[_ ,.-])")
PDB_FILENAME = re.compile(r"(?i)^([0-9][A-Z0-9]{3})(?:\b|[-_.])")
TEMP_RE = re.compile(r"(?i)(298|309|310|311(?:\.5)?|312|313)\s*K\b")
SEED_RE = re.compile(r"(?i)(?:random\s*seed\s*#?|\brs)\s*[-_#]?\s*(\d+)")
DURATION_RE = re.compile(r"(?i)(\d+)\s*ns\b")


def extract_pdb(text: str, filenames: Iterable[str] = ()) -> str | None:
    for pattern in (PDB_EXPLICIT, PDB_LEADING):
        match = pattern.search(text)
        if match:
            return match.group(1).upper()
    for filename in filenames:
        match = PDB_FILENAME.search(filename)
        if match:
            candidate = match.group(1).upper()
            if not candidate.endswith("NS"):
                return candidate
    return None


def one_or_join(values: Iterable[str]) -> str | None:
    unique = sorted({value for value in values if value})
    if not unique:
        return None
    return unique[0] if len(unique) == 1 else "|".join(unique)


def extract_metadata(text: str, filenames: Iterable[str]) -> dict[str, Any]:
    joined = " ".join([text, *filenames])
    temperatures = [match.group(1) for match in TEMP_RE.finditer(joined)]
    seeds = [match.group(1) for match in SEED_RE.finditer(joined)]
    durations = [match.group(1) for match in DURATION_RE.finditer(text)]
    if not durations:
        durations = [match.group(1) for match in DURATION_RE.finditer(joined)]
    return {
        "pdb_id": extract_pdb(text, filenames),
        "temperature_k": one_or_join(temperatures),
        "seed": one_or_join(seeds),
        "duration_ns": one_or_join(durations),
    }


def classify_title(title: str) -> str:
    lower = title.lower().replace("_", " ")
    if re.search(r"\bcomplex\b|\bbound\b|\bbinding\b", lower):
        return "complex"
    antibody_terms = (
        "antibody",
        " mab",
        "mab ",
        "fab",
        "immunoglobulin",
        "durvalumab",
        "pembrolizumab",
        "infliximab",
        "canakinumab",
        "tremelimumab",
        "ipilimumab",
    )
    if any(term in lower for term in antibody_terms):
        return "isolated_antibody"
    return "isolated_other"


def file_flags(filenames: Iterable[str]) -> dict[str, Any]:
    names = list(filenames)
    lower = [name.lower() for name in names]
    production = [
        name
        for name in names
        if name.lower() in {"dynamic-nopbc.xtc", "system_reduced.xtc"}
        or re.search(r"(?i)(?:production|dynamic|md)[-_].*\.xtc$", name)
    ]
    snapshot_xtc = [
        name
        for name in names
        if name.lower().endswith(".xtc") and name not in production
    ]
    rmsf = [
        name
        for name in names
        if "rmsf" in name.lower() and name.lower().endswith((".xvg", ".tab", ".csv"))
    ]
    return {
        "file_count": len(names),
        "has_production_xtc": bool(production),
        "production_xtc_count": len(production),
        "production_xtc_names": "|".join(sorted(production)),
        "snapshot_xtc_count": len(snapshot_xtc),
        "has_rmsf": bool(rmsf),
        "rmsf_file_count": len(rmsf),
        "rmsf_names": "|".join(sorted(rmsf)),
        "has_tpr": any(name.endswith(".tpr") for name in lower),
        "has_gro": any(name.endswith(".gro") for name in lower),
        "has_pdb": any(name.endswith(".pdb") for name in lower),
    }


def load_manifest(path: Path) -> pd.DataFrame:
    grouped: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                row = json.loads(line)
                grouped[(row["doi"], row["entry"])].append(row)
    records: list[dict[str, Any]] = []
    for (doi, entry), files in grouped.items():
        filenames = [item["filename"] for item in files]
        record = {
            "doi": doi,
            "entry": entry,
            "entry_kind": classify_title(entry),
            "declared_bytes": sum(int(item.get("size", 0)) for item in files),
        }
        record.update(extract_metadata(entry, filenames))
        record.update(file_flags(filenames))
        records.append(record)
    return pd.DataFrame(records)


def load_local(root: Path) -> pd.DataFrame:
    records: list[dict[str, Any]] = []
    for collection in ("data", "ThermoPCD_labels"):
        collection_path = root / collection
        if not collection_path.exists():
            continue
        for directory in sorted(path for path in collection_path.iterdir() if path.is_dir()):
            files = [path for path in directory.iterdir() if path.is_file()]
            filenames = [path.name for path in files]
            production_paths = [
                path
                for path in files
                if path.name.lower() in {"dynamic-nopbc.xtc", "system_reduced.xtc"}
                or re.search(r"(?i)(?:production|dynamic|md)[-_].*\.xtc$", path.name)
            ]
            record = {
                "collection": collection,
                "directory": directory.name,
                "path": str(directory),
                "entry_kind": classify_title(directory.name),
                "local_bytes": sum(path.stat().st_size for path in files),
                "production_xtc_bytes": sum(path.stat().st_size for path in production_paths),
            }
            record.update(extract_metadata(directory.name, filenames))
            record.update(file_flags(filenames))
            records.append(record)
    return pd.DataFrame(records)


def candidate_pairs(local: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    present = local.loc[local["pdb_id"].notna()].copy()
    for pdb_id, group in present.groupby("pdb_id", sort=True):
        complex_rows = group.loc[group["entry_kind"] == "complex"]
        isolated_rows = group.loc[group["entry_kind"] != "complex"]
        if complex_rows.empty or isolated_rows.empty:
            continue
        complex_temps = set(complex_rows["temperature_k"].dropna().astype(str))
        isolated_temps = set(isolated_rows["temperature_k"].dropna().astype(str))
        complex_seeds = set(complex_rows["seed"].dropna().astype(str))
        isolated_seeds = set(isolated_rows["seed"].dropna().astype(str))
        rows.append(
            {
                "pdb_id": pdb_id,
                "complex_directories": len(complex_rows),
                "isolated_directories": len(isolated_rows),
                "complex_production_xtc": int(complex_rows["has_production_xtc"].sum()),
                "isolated_production_xtc": int(isolated_rows["has_production_xtc"].sum()),
                "complex_rmsf": int(complex_rows["has_rmsf"].sum()),
                "isolated_rmsf": int(isolated_rows["has_rmsf"].sum()),
                "shared_temperatures": "|".join(sorted(complex_temps & isolated_temps)),
                "shared_seed_labels": "|".join(sorted(complex_seeds & isolated_seeds)),
                "isolated_kinds": "|".join(sorted(set(isolated_rows["entry_kind"]))),
                "complex_titles": " || ".join(sorted(set(complex_rows["directory"]))),
                "isolated_titles": " || ".join(sorted(set(isolated_rows["directory"]))),
            }
        )
    return pd.DataFrame(rows)


def pdb_summary(local: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    present = local.loc[local["pdb_id"].notna()].copy()
    for pdb_id, group in present.groupby("pdb_id", sort=True):
        rows.append(
            {
                "pdb_id": pdb_id,
                "local_directories": len(group),
                "complex_directories": int((group["entry_kind"] == "complex").sum()),
                "isolated_antibody_directories": int((group["entry_kind"] == "isolated_antibody").sum()),
                "isolated_other_directories": int((group["entry_kind"] == "isolated_other").sum()),
                "production_xtc_directories": int(group["has_production_xtc"].sum()),
                "rmsf_directories": int(group["has_rmsf"].sum()),
                "temperatures": "|".join(sorted(set(group["temperature_k"].dropna().astype(str)))),
                "seeds": "|".join(sorted(set(group["seed"].dropna().astype(str)))),
            }
        )
    return pd.DataFrame(rows)


def exact_training_overlap(pdb_ids: set[str], repo: Path) -> dict[str, Any]:
    dynarepo_path = repo / "artifacts" / "data" / "dynarepo_pairs.jsonl"
    dynarepo: set[str] = set()
    if dynarepo_path.exists():
        with dynarepo_path.open("r", encoding="utf-8") as handle:
            for line in handle:
                row = json.loads(line)
                pdb_id = row.get("metadata", {}).get("pdb_id")
                if pdb_id:
                    dynarepo.add(str(pdb_id).upper())
    zenodo = {"1AK4", "1AY7", "1BRS", "1EMV", "1GCQ", "1PVH", "2OOB", "3SGB"}
    return {
        "dynarepo_pdb_ids_detected": sorted(dynarepo),
        "thermopcd_dynarepo_exact_pdb_overlap": sorted(pdb_ids & dynarepo),
        "thermopcd_zenodo_exact_pdb_overlap": sorted(pdb_ids & zenodo),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(r"D:\data\ThermoPCD_trajectories"))
    parser.add_argument("--output-dir", type=Path, default=Path("artifacts/thermopcd_audit"))
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    manifest = load_manifest(args.root / "manifest.jsonl")
    local = load_local(args.root)
    candidates = candidate_pairs(local)
    by_pdb = pdb_summary(local)

    manifest.to_csv(args.output_dir / "manifest_entries.csv", index=False)
    local.to_csv(args.output_dir / "local_entries.csv", index=False)
    candidates.to_csv(args.output_dir / "candidate_component_pairs.csv", index=False)
    by_pdb.to_csv(args.output_dir / "pdb_inventory.csv", index=False)

    local_pdb_ids = set(local["pdb_id"].dropna().astype(str))
    verify_text = (args.root / "verify.log").read_text(encoding="utf-8", errors="replace")
    verify_match = re.search(r"VERIFY DONE ok=(\d+) missing=(\d+) sizeMismatch=(\d+)", verify_text)
    verify = (
        {
            "ok": int(verify_match.group(1)),
            "missing": int(verify_match.group(2)),
            "size_mismatch": int(verify_match.group(3)),
        }
        if verify_match
        else {}
    )
    binary_mismatches = [
        line
        for line in verify_text.splitlines()
        if line.startswith("SIZE ") and re.search(r"(?i)\.(xtc|trr|tpr|edr)\s", line)
    ]
    kind_counts = Counter(local["entry_kind"])
    ambiguous_seed = local["seed"].fillna("").astype(str).str.contains("\\|", regex=True)
    ambiguous_temperature = (
        local["temperature_k"].fillna("").astype(str).str.contains("\\|", regex=True)
    )
    production = local["has_production_xtc"].astype(bool)
    summary = {
        "root": str(args.root),
        "catalogue_file_rows": sum(1 for _ in (args.root / "manifest.jsonl").open("r", encoding="utf-8")),
        "catalogue_entries": int(len(manifest)),
        "catalogue_unique_pdb_ids": int(manifest["pdb_id"].nunique()),
        "local_directories": int(len(local)),
        "local_unique_pdb_ids": int(len(local_pdb_ids)),
        "local_entry_kind_counts": dict(kind_counts),
        "local_production_xtc_directories": int(local["has_production_xtc"].sum()),
        "local_production_xtc_files": int(local["production_xtc_count"].sum()),
        "local_production_xtc_bytes": int(local["production_xtc_bytes"].sum()),
        "local_rmsf_directories": int(local["has_rmsf"].sum()),
        "local_rmsf_files": int(local["rmsf_file_count"].sum()),
        "local_directories_with_merged_seed_evidence": int(ambiguous_seed.sum()),
        "production_xtc_directories_with_merged_seed_evidence": int((production & ambiguous_seed).sum()),
        "local_directories_with_merged_temperature_evidence": int(ambiguous_temperature.sum()),
        "production_xtc_directories_with_merged_temperature_evidence": int(
            (production & ambiguous_temperature).sum()
        ),
        "candidate_pdb_ids_with_complex_and_isolated_component": int(len(candidates)),
        "candidate_pdb_ids": candidates["pdb_id"].tolist() if not candidates.empty else [],
        "verify": verify,
        "binary_size_mismatch_count": len(binary_mismatches),
        "binary_size_mismatches": binary_mismatches,
        "exact_pdb_overlap": exact_training_overlap(local_pdb_ids, Path.cwd()),
        "interpretation_boundary": (
            "A shared PDB identifier and temperature make a candidate comparison, "
            "but do not prove a strict matched apo/holo construction."
        ),
    }
    (args.output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
