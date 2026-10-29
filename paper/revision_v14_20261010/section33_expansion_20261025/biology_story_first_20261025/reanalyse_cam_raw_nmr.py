"""Describe four matched CaM NMR datasets directly from deposited NMR-STAR files.

No model predictions, fitting, statistical significance tests, or figures are used.
Run with ordinary Python 3; the script uses only the standard library.
"""
from pathlib import Path
import csv
import hashlib
import json
import math
import re
import shlex
import statistics

SOURCE = Path(r"E:\Project\FlexPPI\artifacts\apo_holo_exploration_20260913")
OUT = Path(__file__).resolve().parent
PARTNERS = {"CaMKK": 15183, "eNOS": 15185, "nNOS": 15191, "smMLCK": 15186}
AA = dict(zip(
    "ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split(),
    "ARNDCQEGHILKMFPSTWYV"))


def write_csv(name, rows):
    with (OUT / name).open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_order_loops(text):
    for match in re.finditer(r"(?m)^\s*loop_\s*$([\s\S]*?)^\s*stop_\s*$", text):
        lines = [line.strip() for line in match[1].splitlines() if line.strip()]
        tags = [line for line in lines if line.startswith("_")]
        if not tags or not tags[0].startswith("_Order_param."):
            continue
        tokens = shlex.split("\n".join(lines[len(tags):]))
        assert len(tokens) % len(tags) == 0, "Incomplete NMR-STAR loop"
        for start in range(0, len(tokens), len(tags)):
            yield dict(zip([tag.split(".", 1)[1] for tag in tags], tokens[start:start + len(tags)]))


def parse_raw(label, number):
    path = SOURCE / "raw" / f"bmr{number}_3.str"
    text = path.read_text(encoding="utf-8")
    cam_entities = []
    for block in re.findall(r"(?m)^save_\S+\s*$([\s\S]*?)^save_\s*$", text):
        name = re.search(r"(?m)^\s*_Entity.Name\s+(.+)$", block)
        entity_id = re.search(r"(?m)^\s*_Entity.ID\s+(\S+)", block)
        sequence = re.search(r"_Entity.Polymer_seq_one_letter_code\s*\n;\s*\n([\s\S]*?)\n;", block)
        if name and entity_id and sequence and "calmodulin" in name[1].lower():
            cam_entities.append((entity_id[1], re.sub(r"\s+", "", sequence[1])))
    assert len(cam_entities) == 1, f"Expected one CaM entity in {path}"
    entity_id, sequence = cam_entities[0]
    rows = []
    for row in read_order_loops(text):
        if row["Entity_ID"] != entity_id or row["Order_param_val"] in (".", "?"):
            continue
        position = int(row["Comp_index_ID"])
        amino_acid = AA[row["Comp_ID"]]
        assert sequence[position - 1] == amino_acid, "Probe does not match CaM entity sequence"
        atom = row["Atom_ID"]
        expected_atom_name = atom == "N" or (amino_acid, atom) in {
            ("A", "CB"), ("I", "CD1"), ("I", "CG2"), ("L", "CD1"), ("L", "CD2"),
            ("M", "CE"), ("T", "CG2"), ("V", "CG1"), ("V", "CG2")}
        # Preserve deposited atom names exactly, including an eNOS THR29 "CG5" typo.
        # That unmatched probe cannot enter the four-way exact-site intersection.
        assert row["Order_parameter_list_ID"] == ("1" if atom == "N" else "2")
        error = float(row["Order_param_val_fit_err"]) if row["Order_param_val_fit_err"] not in (".", "?") else None
        rows.append(dict(bmrb=number, label=label, source_row_id=row["ID"], seq_id=position, amino_acid=amino_acid,
                         atom=atom, observable="backbone_NH" if atom == "N" else "methyl_axis",
                         s2=float(row["Order_param_val"]), fit_error=error,
                         expected_atom_name=expected_atom_name))
    keys = [(r["observable"], r["seq_id"], r["amino_acid"], r["atom"]) for r in rows]
    duplicate_keys = [key for key in set(keys) if keys.count(key) > 1]
    return rows, dict(label=label, bmrb=number, path=str(path), sha256=sha(path),
                      cam_entity_id=entity_id, cam_entity_length=len(sequence),
                      cam_sequence=sequence, parsed_rows=len(rows),
                      unrecognised_atom_names=[r for r in rows if not r["expected_atom_name"]],
                      duplicate_probe_rows=[r for r in rows if (r["observable"], r["seq_id"], r["amino_acid"], r["atom"]) in duplicate_keys])


def describe(values):
    values = [v for v in values if v is not None and math.isfinite(v)]
    return dict(n=len(values), mean=statistics.mean(values), median=statistics.median(values),
                minimum=min(values), maximum=max(values),
                population_sd=statistics.pstdev(values), zero_count=sum(v == 0 for v in values))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    match_path = SOURCE / "analysis" / "calmodulin_structure_experiment_match.csv"
    with match_path.open(encoding="utf-8-sig", newline="") as handle:
        matches = {r["case_id"]: r for r in csv.DictReader(handle)}
    gates = []
    for label in PARTNERS:
        row = matches[f"cam_{label}"]
        identity = float(row["peptide_identity"])
        coverage = float(row["peptide_bidirectional_coverage"])
        passed = identity >= 0.95 and coverage >= 0.85
        assert passed, f"Peptide matching gate failed: {label}"
        gates.append(dict(label=label, bmrb=PARTNERS[label], peptide_match_pass=passed,
                          peptide_identity=identity, peptide_bidirectional_coverage=coverage,
                          peptide_sequence_pdb=row["peptide_sequence_pdb"],
                          peptide_sequence_bmrb=row["peptide_sequence_bmrb"]))
    all_rows, inputs = [], []
    for label, number in PARTNERS.items():
        rows, audit = parse_raw(label, number)
        all_rows.extend(rows)
        inputs.append(audit)
    assert len({a["cam_sequence"] for a in inputs}) == 1, "CaM entity sequences differ"
    write_csv("four_partners_peptide_match_gate.csv", gates)
    write_csv("raw_nmr_four_partners_long.csv", all_rows)
    indexed = {}
    for row in all_rows:
        key = (row["observable"], row["seq_id"], row["amino_acid"], row["atom"])
        indexed.setdefault(key, {}).setdefault(row["label"], []).append(row)
    wide, excluded = [], []
    for key, partners in indexed.items():
        if any(len(values) != 1 for values in partners.values()):
            for values in partners.values():
                excluded.extend(values)
    for key, partners in sorted(indexed.items()):
        if set(partners) != set(PARTNERS):
            continue
        observable, position, amino_acid, atom = key
        ambiguous = any(len(partners[label]) != 1 for label in PARTNERS)
        if ambiguous:
            continue
        partners = {label: values[0] for label, values in partners.items()}
        out = dict(observable=observable, seq_id=position, amino_acid=amino_acid, atom=atom)
        for label in PARTNERS:
            out[f"{label}_s2"] = partners[label]["s2"]
            out[f"{label}_fit_error"] = partners[label]["fit_error"]
        lowest = min(PARTNERS, key=lambda label: partners[label]["s2"])
        highest = max(PARTNERS, key=lambda label: partners[label]["s2"])
        value_range = partners[highest]["s2"] - partners[lowest]["s2"]
        errors = [partners[lowest]["fit_error"], partners[highest]["fit_error"]]
        error_sum = sum(errors) if all(e is not None for e in errors) else None
        out.update(minimum_partner=lowest, maximum_partner=highest, across_partner_range=value_range,
                   extrema_fit_error_sum=error_sum,
                   range_over_extrema_fit_error_sum=value_range / error_sum if error_sum and error_sum > 0 else None,
                   range_exceeds_extrema_fit_error_sum=value_range > error_sum if error_sum is not None else None)
        wide.append(out)
    nh = [r for r in wide if r["observable"] == "backbone_NH"]
    methyl = [r for r in wide if r["observable"] == "methyl_axis"]
    assert len(nh) == 87 and len(methyl) == 44, (len(nh), len(methyl))
    write_csv("common_87_NH_probes_wide.csv", nh)
    write_csv("common_44_unambiguous_methyl_probes_wide.csv", methyl)
    write_csv("excluded_ambiguous_probe_all_partner_rows.csv", excluded)
    write_csv("common_all_probes_wide_with_ranges.csv", wide)
    partner_summary, range_summary = [], []
    for observable, rows in [("backbone_NH", nh), ("methyl_axis", methyl)]:
        residue_count = len({(r["seq_id"], r["amino_acid"]) for r in rows})
        for label in PARTNERS:
            values = describe([r[f"{label}_s2"] for r in rows])
            errors = describe([r[f"{label}_fit_error"] for r in rows])
            partner_summary.append(dict(observable=observable, label=label,
                common_probe_count=len(rows), distinct_residue_count=residue_count,
                **{f"s2_{k}": v for k, v in values.items()},
                **{f"reported_fit_error_{k}": v for k, v in errors.items()}))
        ranges = describe([r["across_partner_range"] for r in rows])
        error_sums = describe([r["extrema_fit_error_sum"] for r in rows])
        ratios = describe([r["range_over_extrema_fit_error_sum"] for r in rows])
        range_summary.append(dict(observable=observable, common_probe_count=len(rows),
            distinct_residue_count=residue_count,
            **{f"range_{k}": v for k, v in ranges.items()},
            **{f"extrema_fit_error_sum_{k}": v for k, v in error_sums.items()},
            **{f"ratio_range_to_extrema_error_sum_{k}": v for k, v in ratios.items()},
            range_exceeds_extrema_error_sum_count=sum(r["range_exceeds_extrema_fit_error_sum"] is True for r in rows)))
    write_csv("common_probes_partner_summary.csv", partner_summary)
    write_csv("common_probes_range_and_fit_error_summary.csv", range_summary)
    audit = dict(analysis="Experimental NMR descriptive reanalysis; current FlexPPI predictions not tested here",
                 inputs=inputs, peptide_match_source=str(match_path), peptide_match_source_sha256=sha(match_path),
                 peptide_gate="identity >= 0.95 and bidirectional exact-residue coverage >= 0.85; all four pass",
                 key="observable + CaM entity sequence position + amino acid + atom name; exact intersection of all four partners; exclude any probe duplicated within a partner",
                 model_prediction_files_read=[], figures_created=[], significance_tests=[],
                 warnings=["nNOS BMRB15191 Ile100 CG2 is deposited twice, with S2 0.708 and 0.793. The probe is absent for smMLCK, so it never enters the four-way intersection. The primary 44-probe methyl analysis contains no ambiguous duplicates and does not average duplicate entries.",
                    "eNOS BMRB15185 Thr29 CG5 is a nonstandard deposited atom name. It is retained unchanged in the raw table and does not enter the exact-site intersection.",
                    "Fit errors are deposited model-fit errors, not full experimental uncertainties.",
                    "Reported zero fit error is retained and does not mean uncertainty is absent.",
                    "The maximum-minus-minimum range selects two extrema among four partners; error ratios are descriptive diagnostics, not statistical significance.",
                    "Methyl probes are not independent residues; one residue may contribute more than one methyl group.",
                    "NH and methyl probes measure different local observables; numerical range comparison is descriptive, not a calibrated motion-amplitude comparison.",
                    "S2 is not Hchi or RMSF, and is not converted to entropy, binding energy, or affinity.",
                    "Four peptide-bound states are compared here; no apo-minus-holo values are calculated.",
                    "Terminal peptide differences remain despite the matching gate; conditions must be checked before confirmatory model validation."])
    (OUT / "raw_nmr_input_audit.json").write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8")
    summary = dict(partners=list(PARTNERS), partner_summary=partner_summary, range_summary=range_summary)
    (OUT / "raw_nmr_descriptive_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
