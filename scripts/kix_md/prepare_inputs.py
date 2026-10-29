"""Prepare audited KIX model-1 inputs without building atoms or running MD.

Only Python's standard library is needed.  Source coordinates and author
residue identifiers are preserved; the MD builder must repair missing heavy
atoms and add hydrogens after reading ``manifest.json``.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_SOURCE = ROOT / "artifacts/apo_holo_mechanistic_literature_v2_20260913/sources/2AGH.pdb"
DEFAULT_OUTPUT = ROOT / "artifacts/kix_md_20260914/inputs"
RANGES = {"A": (291, 315), "B": (586, 672), "C": (839, 857)}
STATES = {"apo": "B", "myb": "AB", "mll": "BC", "both": "ABC"}
NAMES = {"A": "cMyb", "B": "KIX", "C": "MLL"}
AA = dict(zip(
    "ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split(),
    "ARNDCQEGHILKMFPSTWYV",
))
SIDECHAINS = {
    "ALA": "CB", "ARG": "CB CG CD NE CZ NH1 NH2", "ASN": "CB CG OD1 ND2",
    "ASP": "CB CG OD1 OD2", "CYS": "CB SG", "GLN": "CB CG CD OE1 NE2",
    "GLU": "CB CG CD OE1 OE2", "GLY": "", "HIS": "CB CG ND1 CD2 CE1 NE2",
    "ILE": "CB CG1 CG2 CD1", "LEU": "CB CG CD1 CD2", "LYS": "CB CG CD CE NZ",
    "MET": "CB CG SD CE", "PHE": "CB CG CD1 CD2 CE1 CE2 CZ",
    "PRO": "CB CG CD", "SER": "CB OG", "THR": "CB OG1 CG2",
    "TRP": "CB CG CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2", "TYR": "CB CG CD1 CD2 CE1 CE2 CZ OH",
    "VAL": "CB CG1 CG2",
}


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def parse_model_one(source: Path) -> tuple[list[dict], list[dict]]:
    atoms, hydrogen_atoms, seen = [], [], set()
    active = False
    models = 0
    for line_number, raw_line in enumerate(source.read_text(encoding="ascii").splitlines(), 1):
        line = raw_line.ljust(80)
        if line.startswith("MODEL "):
            model_number = int(line[10:14])
            active = model_number == 1
            if active:
                models += 1
            continue
        if line.startswith("ENDMDL") and active:
            break
        if not active or not line.startswith("ATOM  "):
            continue
        if line[16] not in (" ", "A"):
            raise ValueError(f"Unexpected alternative location at source line {line_number}")
        name, resname, chain = line[12:16].strip(), line[17:20], line[21]
        if resname not in AA or chain not in RANGES:
            raise ValueError(f"Unexpected protein record at source line {line_number}")
        number, insertion = int(line[22:26]), line[26].strip()
        key = (chain, number, insertion, name)
        if key in seen:
            raise ValueError(f"Duplicate source atom {key}")
        seen.add(key)
        xyz = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
        if not all(math.isfinite(v) for v in xyz):
            raise ValueError(f"Nonfinite coordinate {key}")
        element = line[76:78].strip() or next(char for char in name if char.isalpha())
        record = {
            "source_line_number": line_number, "source_atom_serial": int(line[6:11]),
            "chain": chain, "author_residue_number": number, "insertion_code": insertion,
            "residue_name": resname, "atom_name": name, "element": element,
            "xyz_angstrom": xyz, "source_atom_line_sha256": digest(raw_line.encode("ascii")),
            "line": raw_line,
        }
        if element in ("H", "D"):
            hydrogen_atoms.append(record)
        else:
            if element not in ("C", "N", "O", "S"):
                raise ValueError(f"Unexpected heavy element {element}: {key}")
            atoms.append(record)
    if models != 1 or not atoms:
        raise ValueError("Exactly one nonempty MODEL 1 is required")
    return atoms, hydrogen_atoms


def by_residue(atoms: list[dict]) -> dict:
    grouped = {}
    for atom in atoms:
        key = (atom["chain"], atom["author_residue_number"], atom["insertion_code"])
        grouped.setdefault(key, []).append(atom)
    return grouped


def public_atom(atom: dict, output_index: int) -> dict:
    return {**{key: value for key, value in atom.items() if key != "line"},
            "output_atom_index_zero_based": output_index}


def pdb_bytes(atoms: list[dict], state: str) -> bytes:
    lines = [
        "REMARK 900 KIX MD INPUT; HEAVY ATOMS FROM 2AGH MODEL 1 ONLY",
        "REMARK 900 COORDINATES UNCHANGED; MISSING HEAVY ATOMS NOT REPAIRED",
        f"REMARK 900 STATE {state}; MLL AUTHOR RESIDUES 839-857 WHEN PRESENT",
    ]
    previous_chain = None
    for atom in atoms:
        if previous_chain is not None and atom["chain"] != previous_chain:
            lines.append("TER")
        lines.append(atom["line"])
        previous_chain = atom["chain"]
    lines.extend(["TER", "END"])
    return ("\n".join(lines) + "\n").encode("ascii")


def prepare(source: Path, output: Path) -> dict:
    source = source.resolve(strict=True)
    source_bytes = source.read_bytes()
    all_atoms, hydrogens = parse_model_one(source)
    all_residues = by_residue(all_atoms)
    selected = [a for a in all_atoms if RANGES[a["chain"]][0] <= a["author_residue_number"] <= RANGES[a["chain"]][1]]
    selected_residues = by_residue(selected)
    chains = {}
    mapping_rows = []
    missing_heavy = []
    for chain, (low, high) in RANGES.items():
        keys = [key for key in selected_residues if key[0] == chain]
        if keys != [(chain, number, "") for number in range(low, high + 1)]:
            raise ValueError(f"Residue numbering or coverage mismatch in chain {chain}: {keys}")
        sequence = "".join(AA[selected_residues[key][0]["residue_name"]] for key in keys)
        chains[chain] = {
            "molecule": NAMES[chain], "first_author_residue": low, "last_author_residue": high,
            "num_residues": len(keys), "sequence": sequence,
            "num_heavy_atoms": sum(len(selected_residues[key]) for key in keys),
        }
        for ordinal, key in enumerate(keys, 1):
            records = selected_residues[key]
            resname = records[0]["residue_name"]
            names = {a["atom_name"] for a in records}
            expected = set(("N CA C O " + SIDECHAINS[resname]).split())
            if ordinal == len(keys):
                expected.add("OXT")
            missing, extra = sorted(expected - names), sorted(names - expected)
            if extra:
                raise ValueError(f"Unexpected heavy atoms in {key}: {extra}")
            if missing:
                missing_heavy.append({"chain": chain, "author_residue_number": key[1],
                                      "residue_name": resname, "missing_atoms": missing})
            mapping_rows.append({
                "chain": chain, "molecule": NAMES[chain], "chain_sequence_index_one_based": ordinal,
                "author_residue_number": key[1], "insertion_code": key[2],
                "residue_name": resname, "amino_acid": AA[resname],
                "observed_heavy_atom_count": len(records),
                "observed_heavy_atoms": ";".join(a["atom_name"] for a in records),
                "missing_heavy_atoms": ";".join(missing),
                "present_in_states": ";".join(state for state, included in STATES.items() if chain in included),
            })
    if chains["C"]["sequence"] != "SDDGNILPSDIMDFVLKNT":
        raise ValueError("MLL main-protocol sequence mismatch")
    if {(m["chain"], m["author_residue_number"], tuple(m["missing_atoms"])) for m in missing_heavy} != {
        ("C", 841, ("CG", "OD1", "OD2")), ("C", 857, ("OXT",)),
    }:
        raise ValueError(f"Unexpected heavy-atom completeness findings: {missing_heavy}")
    h_by_residue = by_residue(hydrogens)
    his = []
    for key, records in selected_residues.items():
        if records[0]["residue_name"] == "HIS":
            ring_h = [a["atom_name"] for a in h_by_residue.get(key, []) if a["atom_name"] in ("HD1", "HE2")]
            if ring_h != ["HD1"]:
                raise ValueError(f"Unanticipated source histidine hydrogens at {key}: {ring_h}")
            his.append({"chain": key[0], "author_residue_number": key[1],
                        "source_ring_hydrogens": ring_h, "suggested_explicit_variant": "HID",
                        "status": "NMR coordinate convention, not independently established MD protonation"})
    manifest = {
        "protocol": "kix_2agh_model1_s1_mll839_857_heavy_inputs_v1",
        "source_pdb": str(source), "source_pdb_sha256": digest(source_bytes), "source_model": 1,
        "source_model1_full_heavy_atoms": len(all_atoms),
        "source_model1_full_hydrogen_atoms": len(hydrogens),
        "source_model1_full_chain_residue_counts": {chain: sum(key[0] == chain for key in all_residues) for chain in RANGES},
        "source_model1_full_heavy_atom_mapping": [public_atom(a, i) for i, a in enumerate(all_atoms)],
        "chains": chains, "states": {}, "missing_heavy_atoms_to_build": missing_heavy,
        "histidine_variants_from_source": his,
        "pdb_experimental_pH": 5.5,
        "s1_MD_pH": "not stated; do not infer from NMR sample pH",
        "termini": {
            "source": "Source A291/B586/C839 contain H1,H2,H3; full-chain C termini have OXT.",
            "proposed_MD": "Uncapped charged N and C termini, declared assumption. Add C857 OXT after truncation.",
            "paper_s1": "No terminal patch/cap or protonation details reported.",
            "experimental_synthetic_peptides": "Do not silently transfer experimental acetyl/amidated caps to this MD reconstruction.",
        },
        "sequence_conflict": "2AGH SEQADV C841 is deposited ASP versus reference UniProt CYS; retain ASP, do not auto-revert.",
        "transformations": ["select MODEL 1", "select author residue ranges", "remove all H/D atoms", "add chain TER records"],
        "not_performed": ["heavy-atom completion", "protonation", "capping", "coordinate relaxation", "MD", "model prediction"],
        "qa": {"all_selected_xyz_unchanged": True, "all_selected_heavy_lines_preserved": True,
               "full_backbone_coverage": True, "no_missing_residues_in_selected_ranges": True,
               "no_sidechain_reversion": True, "state_names": list(STATES)},
    }
    pending = {}
    for state, included in STATES.items():
        atoms = [a for a in selected if a["chain"] in included]
        data = pdb_bytes(atoms, state)
        output_lines = [line for line in data.decode("ascii").splitlines() if line.startswith("ATOM  ")]
        if output_lines != [a["line"] for a in atoms]:
            raise AssertionError("Input changed a source atom line")
        filename = state + ".pdb"
        pending[filename] = data
        manifest["states"][state] = {
            "pdb_filename": filename, "pdb_sha256": digest(data), "chains": list(included),
            "num_residues": sum(chains[ch]["num_residues"] for ch in included),
            "num_heavy_atoms": len(atoms),
            "source_pdb_model": 1, "source_coordinates_unchanged": True,
            "missing_heavy_atoms_to_build": [m for m in missing_heavy if m["chain"] in included],
            "atom_source_mapping": [public_atom(a, i) for i, a in enumerate(atoms)],
        }
    manifest_bytes = (json.dumps(manifest, ensure_ascii=False, indent=2) + "\n").encode("utf-8")
    pending["manifest.json"] = manifest_bytes
    import io
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=list(mapping_rows[0]))
    writer.writeheader()
    writer.writerows(mapping_rows)
    pending["residue_mapping.csv"] = buffer.getvalue().encode("utf-8-sig")
    pending["input_differences.md"] = (
        "# KIX MD 输入与既有预测的差异\n\n"
        "主方案遵照 Peiffer 2022 S1 的 MD 构建方法：2AGH 模型1，KIX B586–672（87残基），"
        "cMyb A291–315（25残基），MLL C839–857（19残基）。四状态为 apo、myb、mll、both。"
        "旧案例脚本使用 MLL840–858，序列 DDGNILPSDIMDFVLKNTP；本轮为839–857，序列 SDDGNILPSDIMDFVLKNT。"
        "不是相同序列，应针对新MD输入重新预测，不能直接套用旧图中的数值。\n\n"
        "2AGH作者存档MLL D841缺少CG、OD1、OD2，末端裁至T857后还缺OXT；这四个原子必须在MD构建阶段补齐。"
        "D841与参考UniProt的Cys存在SEQADV冲突，本输入保留实际存档的Asp，不执行野生型回退。"
        "本脚本没有补原子或改变任何保留原子的坐标，manifest逐原子追溯源行号、原子序号、坐标和哈希。\n\n"
        "源PDB五个His（B592、594、602、605、651）均仅含HD1环氢，支持把HID作为明确的初始假设，"
        "但不能证明原作者MD就用此质子化状态。源结构N端含H1/H2/H3，完整C端含OXT。"
        "拟采用带电未封端链端；S1未明确端基，不能把正文合成实验肽的乙酰化/酰胺化直接当成MD事实。"
        "PDB的pH5.5是NMR样本条件，不是S1确认的MD pH。\n\n"
        "必须由后续构建流程记录实际HIS模板、总电荷、重建原子与端基，并检查蛋白主链完整、"
        "没有意外链间共价键，再进入能量最小化及短测试。该输入产物本身不能直接视为MD就绪。\n"
    ).encode("utf-8")
    output.mkdir(parents=True, exist_ok=True)
    for filename, data in pending.items():
        target = output / filename
        if target.exists() and target.read_bytes() != data:
            raise FileExistsError(f"Refusing to overwrite differing existing input: {target}")
    for filename, data in pending.items():
        (output / filename).write_bytes(data)
    return {"output_directory": str(output.resolve()), "source_sha256": manifest["source_pdb_sha256"],
            "states": {key: {name: value[name] for name in ("num_residues", "num_heavy_atoms", "pdb_sha256")} for key, value in manifest["states"].items()},
            "missing_heavy_atoms_to_build": missing_heavy, "histidine_variants_from_source": his,
            "qa": manifest["qa"]}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    print(json.dumps(prepare(args.source, args.output), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
