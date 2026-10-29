"""Build an auditable, outcome-independent exploratory apo/holo cohort.

Downloaded structures are public RCSB files. Assembly models are flattened only
for biological-assembly files, never for NMR ensemble files.
"""
from __future__ import annotations

import concurrent.futures
import hashlib
import json
import string
import sys
from datetime import datetime, timezone
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "artifacts/apo_holo_exploration_20260913"
sys.path.insert(0, str(ROOT / "src"))
from flexppi.data.structures import parse_pdb
from flexppi.data.pairing import global_sequence_alignment

GRUNBERG = [
    ("c01", "1BRA", "*", "1AAP", "A", "1BRC", "E", "I"),
    ("c02", "2CGA", "A", "1HPT", "*", "1CGI", "E", "I"),
    ("c03", "2PKA", "AB", "5PTI", "*", "2KAI", "AB", "I"),
    ("c04", "1SUP", "*", "3SSI", "*", "2SIC", "E", "I"),
    ("c05", "1FGN", "LH", "1BOY", "*", "1AHW", "AB", "C"),
    ("c06", "1BVL", "AB", "3LZT", "*", "1BVK", "AB", "C"),
    ("c08", "1DQQ", "AB", "3LZT", "*", "1DQJ", "AB", "C"),
    ("c11", "1A19", "A", "1A2P", "A", "1BGS", "A", "E"),
    ("c13", "2BNH", "*", "7RSA", "*", "1DFJ", "E", "I"),
    ("c14", "1VXR", "*", "1FSC", "A", "1FSS", "A", "B"),
    ("c15", "1AVV", "*", "1SHF", "A", "1AVZ", "B", "C"),
    ("c16", "1AKZ", "*", "1UGI", "A", "1UGH", "E", "I"),
    ("c17", "1WER", "*", "5P21", "*", "1WQ1", "R", "G"),
    ("c19", "1PIF", "*", "2AIT", "*", "1BVN", "P", "T"),
    ("c20", "1HCL", "*", "1VIN", "*", "1FIN", "A", "B"),
    ("c21", "1B39", "A", "1FPZ", "A", "1FQ1", "A", "B"),
    ("c22", "1TBG", "AE", "1TAG", "*", "1GOT", "A", "BG"),
]
COLE = "1CDT 1G6N 1IL8 1MSB 1PP2 1UTG 1YPI 2CCY 2CTS 2GN5 2RUS 2RVE 2SOD 2TS1 2TSC 2WRP 3AAT 3GRS 3SDH 3SDP 3SSI 4MDH 5ADH 5HVP 8TIM".split()
CAM = [("smMLCK", "1CDL", "A", "E"), ("eNOS", "1NIW", "A", "B"),
       ("nNOS", "2O60", "A", "B"), ("CaMKI", "1MXE", "A", "E"),
       ("CaMKK", "1CKK", "A", "B")]


def dump(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False), encoding="utf-8")


def fetch(task):
    name, url = task
    path = OUT / "raw" / name
    path.parent.mkdir(parents=True, exist_ok=True)
    result = dict(name=name, url=url, retrieved_utc=datetime.now(timezone.utc).isoformat())
    try:
        if not path.exists():
            response = requests.get(url, timeout=(10, 35))
            response.raise_for_status()
            if len(response.content) < 100 or response.content.lstrip().startswith(b"<!DOCTYPE"):
                raise ValueError("Not a structure/data file")
            path.write_bytes(response.content)
        result.update(status="downloaded", path=str(path), sha256=hashlib.sha256(path.read_bytes()).hexdigest())
    except Exception as exc:
        result.update(status="failed", error=str(exc))
    return result


def flatten_assembly(pdb):
    source = OUT / "raw" / f"{pdb}.pdb1"
    text = source.read_text()
    # RCSB legacy assembly format uses MODEL for symmetry-related copies.
    model, mapping, records = 1, {}, []
    for line in text.splitlines():
        if line.startswith("MODEL"):
            model = int(line[10:14])
        if line.startswith(("ATOM  ", "HETATM")):
            key = (model, line[21])
            if key not in mapping:
                if len(mapping) >= 62:
                    raise ValueError("Assembly exceeds PDB chain namespace")
                mapping[key] = (string.ascii_uppercase + string.ascii_lowercase + string.digits)[len(mapping)]
            records.append(line[:21] + mapping[key] + line[22:])
    destination = OUT / "pdb" / f"{pdb}_assembly1.pdb"
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text("\n".join(records) + "\nEND\n")
    dump(destination.with_suffix(".provenance.json"), dict(source=str(source),
         source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
         transformation="Flatten deposited RCSB biological assembly 1; coordinates unchanged",
         chain_mapping=[dict(model=k[0], author_chain=k[1], output_chain=v) for k, v in mapping.items()]))
    return destination


def spec(pdb, chains):
    return dict(pdb_id=pdb, path=str(OUT / "raw" / f"{pdb}.pdb"),
                chains=list(chains) if chains != "*" else None)


def main():
    pdb_ids = {x for row in GRUNBERG for x in (row[1], row[3], row[5])}
    pdb_ids |= {"3LNX", "3LNY", "1CLL"} | {x[1] for x in CAM}
    tasks = [(f"{p}.pdb", f"https://files.rcsb.org/download/{p}.pdb") for p in sorted(pdb_ids)]
    tasks += [(f"{p}.pdb1", f"https://files.rcsb.org/download/{p}.pdb1") for p in COLE]
    tasks += [(f"bmr{p}_3.str", f"https://bmrb.io/ftp/pub/bmrb/entry_directories/bmr{p}/bmr{p}_3.str")
              for p in (15183, 15184, 15185, 15186, 15187, 15188, 15191)]
    with concurrent.futures.ThreadPoolExecutor(max_workers=6) as pool:
        downloads = []
        for item in pool.map(fetch, tasks):
            downloads.append(item)
            print(item["name"], item["status"], flush=True)
    dump(OUT / "download_manifest.json", downloads)
    cases = []
    for line in (ROOT / "artifacts/zenodo_6638504/data/zenodo_6638504_pairs.jsonl").read_text().splitlines():
        old = json.loads(line)
        apo, holo = old["apo"], old["holo"]
        cases.append(dict(case_id="md_" + old["pair_id"], direction="md_delta", group=old["metadata"]["pdb_id"],
            apo=dict(path=apo["structure_path"], chains=apo["target_chains"], pdb_id=old["metadata"]["pdb_id"], label_path=apo["label_path"]),
            holo=dict(path=holo["structure_path"], chains=None, pdb_id=old["metadata"]["pdb_id"], label_path=holo["label_path"]),
            target_chains=holo["target_chains"], exact_indices=[old["apo_indices"], old["holo_indices"]],
            state_type="MD-derived starting structures; paired trajectory labels; label protocol provenance pending",
            source_url="https://zenodo.org/records/6638504"))
    for cid, a, ac, b, bc, h, hc, pc in GRUNBERG:
        for side, target in enumerate((hc, pc)):
            cases.append(dict(case_id=f"grunberg_{cid}_{side+1}", direction="redistribution", group=h,
                apo_candidates=[spec(a, ac), spec(b, bc)], holo=spec(h, hc + pc), target_chains=list(target),
                state_type="experimentally unbound versus bound; sequence-based component assignment",
                notes="Table spells 3LTZ and 1BSG. Corrected to 3LZT and 1BGS based on RCSB protein identities and A:E author chains (1BSG is unrelated beta-lactamase). These are inferred typo repairs, not author-confirmed. Components assigned by sequence, not table column.",
                source_url="https://biskit.pasteur.fr/use/publications/pdf-files/Structure_06_with_errata.pdf"))
    for label, h, hc, pc in CAM:
        cases.append(dict(case_id="cam_" + label, direction="backbone_sidechain", group="CaM_" + label,
            apo=spec("1CLL", "A"), holo=spec(h, hc + pc), target_chains=[hc],
            state_type="calcium-loaded unliganded CaM versus peptide-bound CaM; calcium not encoded by model",
            source_url="https://doi.org/10.1038/nature05959"))
    cases.append(dict(case_id="cam_PDE", direction="backbone_sidechain", group="CaM_PDE",
                      status="unavailable", error="Matching PDE peptide complex structure not yet verified; not replaced by another ligand",
                      source_url="https://doi.org/10.1038/nature05959"))
    cases.append(dict(case_id="pdz2_RAGEF2", direction="remote_response", group="PTP1E_PDZ2",
        apo=spec("3LNX", "A"), holo=spec("3LNY", "AB"), target_chains=["A"],
        source_url="https://pubmed.ncbi.nlm.nih.gov/14698303/",
        responsive_residues="L18 I20 V22 V26 V30 A39 V40 V61 V64 L66 A69 L78 T81 V85".split(),
        response_list_source="https://pmc.ncbi.nlm.nih.gov/articles/PMC3510070/",
        state_type="experimentally unbound versus bound; response-list negatives unknown"))
    for p in COLE:
        base = dict(case_id="cole_" + p, direction="preorganization", group=p,
                    source_url="https://pmc.ncbi.nlm.nih.gov/articles/PMC2373749/",
                    state_type="partner-stripped holo conformation, NOT independently observed apo")
        try:
            path = flatten_assembly(p)
            structure = parse_pdb(str(path))
            chains = structure.unique_chains()
            if len(chains) != 2:
                raise ValueError(f"Assembly 1 contains {len(chains)} protein chains, not an unambiguous dimer")
            seq1, seq2 = [structure.chain_sequence(c) for c in chains]
            aln = global_sequence_alignment(seq1, seq2)
            if aln.identity < .95 or min(len(aln.indices_a)/len(seq1), len(aln.indices_b)/len(seq2)) < .85:
                raise ValueError("Assembly 1 is not a sequence-matched homodimer")
            base.update(apo=dict(path=str(path), chains=[chains[0]], pdb_id=p),
                        holo=dict(path=str(path), chains=chains, pdb_id=p), target_chains=[chains[0]])
        except Exception as exc:
            base.update(status="unavailable", error=str(exc))
        cases.append(base)
    dump(OUT / "cohorts.json", cases)
    print("Cohort:", len(cases), "cases; explicitly unavailable:", sum(x.get("status") == "unavailable" for x in cases))


if __name__ == "__main__":
    main()
