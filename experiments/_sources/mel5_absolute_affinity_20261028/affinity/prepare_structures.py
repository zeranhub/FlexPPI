"""Prepare fixed public crystal inputs without changing heavy-atom coordinates."""
from collections import OrderedDict
from pathlib import Path
import hashlib
import json

ROOT = Path(__file__).resolve().parent
AA = dict(zip(
    "ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split(),
    "ARNDCQEGHILKMFPSTWYV",
))
EXPECTED_PEPTIDES = {
    "3hg1": "ELAGIGILTV", "4jff": "ELAGIGILTV",
    "4jfd": "ELAAIGILTV", "4jfe": "ELAGIGALTV",
}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prepare(pdb):
    raw = ROOT / "structures" / f"{pdb}_crystal_raw.pdb"
    out = ROOT / "structures" / f"{pdb}_crystal_protein.pdb"
    # Keep the first model, protein ATOM records in chains A-E, and one
    # highest-occupancy conformer for each atom (blank, then A, for ties).
    atoms = OrderedDict()
    removed = {"heteroatom_records": 0, "hydrogen_records": 0, "alternate_records": 0}
    for line in raw.read_text(encoding="ascii").splitlines():
        if line.startswith("ENDMDL"):
            break
        if line.startswith("HETATM"):
            removed["heteroatom_records"] += 1
        if not line.startswith("ATOM  ") or line[21] not in "ABCDE":
            continue
        if line[17:20] not in AA:
            raise ValueError(f"Unrecognized ATOM residue: {line}")
        if line[76:78].strip() in {"H", "D"} or line[12:16].strip().startswith("H"):
            removed["hydrogen_records"] += 1
            continue
        key = (line[21], line[22:27], line[12:16])
        priority = (float(line[54:60] or 0), line[16] == " ", line[16] == "A")
        if key in atoms:
            removed["alternate_records"] += 1
        if key not in atoms or priority > atoms[key][0]:
            atoms[key] = (priority, line[:16] + " " + line[17:])
    lines = []
    chains = {}
    serial = 1
    for chain in "ABCDE":
        chain_atoms = [line for _, line in atoms.values() if line[21] == chain]
        assert chain_atoms, (pdb, chain)
        residues = OrderedDict()
        ca_residues = set()
        for line in chain_atoms:
            residues[line[22:27]] = line[17:20]
            if line[12:16].strip() == "CA":
                ca_residues.add(line[22:27])
            lines.append(line[:6] + f"{serial:5d}" + line[11:])
            serial += 1
        assert set(residues) == ca_residues, (pdb, chain, set(residues) - ca_residues)
        last = chain_atoms[-1]
        lines.append(f"TER   {serial:5d}      {last[17:20]} {chain}{last[22:27]}")
        serial += 1
        chains[chain] = {
            "role": {"A": "HLA-A*02:01", "B": "beta2-microglobulin", "C": "antigen_peptide", "D": "TCR_alpha", "E": "TCR_beta"}[chain],
            "residues": len(residues), "heavy_atoms": len(chain_atoms),
            "sequence": "".join(AA[res] for res in residues.values()),
            "author_residue_ids": list(residues),
        }
    assert chains["C"]["sequence"] == EXPECTED_PEPTIDES[pdb]
    assert chains["A"]["sequence"].startswith("GSHSMRYFF")
    assert chains["B"]["sequence"].startswith("MIQRTPKI")
    assert "EVEQNSGPLSVPEGA" in chains["D"]["sequence"]
    assert "SQTIHQWPATLVQPVG" in chains["E"]["sequence"]
    out.write_text("\n".join(lines) + "\nEND\n", encoding="ascii")
    return {
        "pdb": pdb, "source_url": f"https://files.rcsb.org/download/{pdb.upper()}.pdb",
        "raw_path": str(raw), "raw_sha256": sha(raw), "input_path": str(out),
        "input_sha256": sha(out), "chains": chains, "removed": removed,
    }


if __name__ == "__main__":
    records = [prepare(pdb) for pdb in EXPECTED_PEPTIDES]
    audit = {
        "temperature_C": 25.0,
        "interaction": "TCR alpha/beta (D,E) versus complete pMHC (A,B,C)",
        "preprocessing": "First crystal model; ATOM-only standard proteins; waters/ions/glycans/hydrogens removed; highest-occupancy atom conformer; no coordinate refinement or mutation modelling.",
        "records": records,
    }
    (ROOT / "structure_audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(json.dumps([{r["pdb"]: {c: v["residues"] for c, v in r["chains"].items()}} for r in records]))
