"""Make rigid receptor PDBQT inputs, preserving all crystal heavy atoms."""
from collections import OrderedDict
from pathlib import Path
import hashlib
import json
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
OBABEL = Path(sys.executable).parent / "obabel.exe"


def atoms_by_key(lines):
    return {(line[21], line[22:27], line[12:16].strip()): line[30:54]
            for line in lines if line.startswith("ATOM  ") and line[77:79].strip() not in {"H", "HD", "HS"}
            and not line[12:16].strip().startswith("H")}


records = []
for pdb in ["3hg1", "4jff", "4jfd", "4jfe"]:
    source = ROOT / "structures" / f"{pdb}_crystal_protein.pdb"
    raw = ROOT / "structures" / f"{pdb}_openbabel_raw.pdbqt"
    out = ROOT / "structures" / f"{pdb}_crystal_protein.pdbqt"
    command = [str(OBABEL), str(source), "-O", str(raw), "-xr", "-xn", "-h", "--partialcharge", "gasteiger"]
    result = subprocess.run(command, capture_output=True, text=True, check=True)
    source_lines = source.read_text(encoding="ascii").splitlines()
    raw_lines = raw.read_text(encoding="ascii").splitlines()
    assert atoms_by_key(source_lines) == atoms_by_key(raw_lines), "Heavy atoms/coordinates changed"
    grouped = OrderedDict()
    for line in source_lines:
        if line.startswith("ATOM  "):
            grouped.setdefault((line[21], line[22:27]), [])
    for line in raw_lines:
        if line.startswith("ATOM  "):
            # AutoDockTools AutoDock4_AtomTyper assigns non-oxidized sulfur SA.
            # The only sulfur atoms in these standard proteins are CYS SG / MET SD.
            # OpenBabel writes these as S, outside the published GNN vocabulary.
            if line[77:79].strip() == "S":
                assert (line[17:20], line[12:16].strip()) in {("CYS", "SG"), ("MET", "SD")}
                line = line[:77] + "SA" + line[79:]
            grouped[(line[21], line[22:27])].append(line)
    # OpenBabel appends hydrogens after all heavy atoms. The official model
    # parser expects contiguous residue records and a TER at each chain end.
    # Restore this standard record order; no coordinates/types are changed.
    lines = []
    for index, ((chain, resid), atom_lines) in enumerate(grouped.items()):
        lines.extend(atom_lines)
        keys = list(grouped)
        if index + 1 == len(keys) or keys[index + 1][0] != chain:
            lines.append(f"TER   {index + 1:5d}      {atom_lines[-1][17:20]} {chain}{resid}")
    lines.append("END")
    out.write_text("\n".join(lines) + "\n", encoding="ascii")
    records.append({
        "pdb": pdb, "command": command, "stdout": result.stdout, "stderr": result.stderr,
        "output_path": str(out), "output_sha256": hashlib.sha256(out.read_bytes()).hexdigest(),
        "residues": len(grouped), "atoms": sum(len(v) for v in grouped.values()),
        "heavy_atom_identity_and_coordinate_check": "identical to crystal protein input",
        "atom_types": sorted({line[77:79].strip() for line in lines if line.startswith("ATOM  ")}),
    })
    print(pdb, len(grouped), records[-1]["atoms"], flush=True)
(ROOT / "pdbqt_audit.json").write_text(json.dumps({
    "preparation": "OpenBabel 3.1.1.23, rigid receptor PDBQT, preserve atom names, add hydrogens, Gasteiger charges; standard CYS SG/MET SD sulfur types assigned SA following AutoDockTools AutoDock4_AtomTyper (non-oxidized sulfur), not OpenBabel S; reorder to contiguous residues and append per-chain TER for the official parser; crystal heavy atoms unchanged.",
    "sulfur_typing_source": "https://raw.githubusercontent.com/Valdes-Tresanco-MS/AutoDockTools_py3/master/AutoDockTools/atomTypeTools.py, AutoDock4_AtomTyper.setAutoDockElements sulfur branch; non-oxidized standard-protein CYS and MET only.",
    "deviation_from_author_readme": "Author suggests ADFR prepare_receptor. OpenBabel was used as the available PDBQT writer. This preparation difference and simplified inference model must be disclosed.",
    "records": records,
}, indent=2), encoding="utf-8")
