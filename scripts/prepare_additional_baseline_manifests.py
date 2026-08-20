"""Prepare immutable manifests for additional RMSF baselines.

The output keeps both the native Windows path (for auditing) and the path seen
inside the Linux baseline containers.  Cohort membership is derived from the
already frozen common-residue evaluations, so new baselines cannot silently
change the evaluated sample set.
"""

from __future__ import annotations

import csv
from pathlib import Path, PureWindowsPath
from typing import Dict, Iterable, List


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "artifacts" / "evaluation" / "additional_baselines" / "manifests"


def read_rows(path: Path, delimiter: str = ",") -> List[Dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def container_path(path_text: str) -> str:
    path = PureWindowsPath(path_text)
    drive = path.drive.upper()
    mount = {"E:": "/data_e", "D:": "/data_d"}.get(drive)
    if mount is None:
        raise ValueError(f"No container mount is defined for {path_text}")
    relative_parts = path.parts[1:]
    if not relative_parts or relative_parts[0].lower() != "data":
        raise ValueError(f"Expected a path below {drive}\\data, got {path_text}")
    return mount + "/" + "/".join(relative_parts[1:])


def assert_paths_exist(rows: Iterable[Dict[str, str]]) -> None:
    missing = [row["pdb_path"] for row in rows if not Path(row["pdb_path"]).is_file()]
    if missing:
        preview = "\n".join(missing[:10])
        raise FileNotFoundError(f"{len(missing)} manifest structures are missing:\n{preview}")


def write_rows(path: Path, rows: List[Dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def prepare_static() -> List[Dict[str, str]]:
    base = ROOT / "artifacts" / "evaluation" / "fastprotflex_comparison"
    atlas = read_rows(base / "atlas_common_manifest.tsv", delimiter="\t")
    dynarepo = read_rows(
        base / "dynarepo" / "dynarepo_common_manifest.tsv", delimiter="\t"
    )
    complex_ids = {
        row["sample_id"]
        for row in read_rows(base / "chain_context" / "complex_per_system_metrics.csv")
    }
    rows: List[Dict[str, str]] = []
    for row in atlas:
        rows.append(
            {
                "cohort": "atlas_single_chain",
                "sample_id": row["sample_id"],
                "pdb_path": row["pdb_path"],
                "container_pdb_path": container_path(row["pdb_path"]),
            }
        )
    for row in dynarepo:
        rows.append(
            {
                "cohort": (
                    "dynarepo_complex"
                    if row["sample_id"] in complex_ids
                    else "dynarepo_single_chain"
                ),
                "sample_id": row["sample_id"],
                "pdb_path": row["pdb_path"],
                "container_pdb_path": container_path(row["pdb_path"]),
            }
        )
    assert_paths_exist(rows)
    return rows


def prepare_itrace() -> List[Dict[str, str]]:
    source = ROOT / "output" / "itrace_benchmark" / "raw" / "itrace_complex_manifest.csv"
    rows: List[Dict[str, str]] = []
    for record in read_rows(source):
        pdb_id = record["pdb_id"].lower()
        if pdb_id == "3kxf":
            continue
        rows.append(
            {
                "cohort": "itrace_external_clean",
                "sample_id": pdb_id,
                "pdb_path": record["topology_path"],
                "container_pdb_path": container_path(record["topology_path"]),
            }
        )
    assert_paths_exist(rows)
    return rows


def main() -> None:
    static_rows = prepare_static()
    itrace_rows = prepare_itrace()
    write_rows(OUTPUT / "static.tsv", static_rows)
    write_rows(
        OUTPUT / "static_complex.tsv",
        [row for row in static_rows if row["cohort"] == "dynarepo_complex"],
    )
    write_rows(OUTPUT / "itrace.tsv", itrace_rows)
    print(
        "Prepared manifests: "
        f"static={len(static_rows)} "
        f"(ATLAS={sum(r['cohort'] == 'atlas_single_chain' for r in static_rows)}, "
        f"DynaRepo-single={sum(r['cohort'] == 'dynarepo_single_chain' for r in static_rows)}, "
        f"DynaRepo-complex={sum(r['cohort'] == 'dynarepo_complex' for r in static_rows)}); "
        f"ITrace-external-clean={len(itrace_rows)}"
    )


if __name__ == "__main__":
    main()
