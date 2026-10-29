"""Build source tables for a CaM experimental NMR figure; standard library only.

No model predictions, statistical tests, or conversions to Hchi/RMSF are used.
The independent QA parser reads deposited NMR-STAR order-parameter loops anew.
"""
from pathlib import Path
from collections import Counter, defaultdict
import csv
import hashlib
import json
import re
import shlex
import statistics

OUT = Path(__file__).resolve().parent
BASE = OUT.parent
RAW = Path(r"E:\Project\FlexPPI\artifacts\apo_holo_exploration_20260913\raw")
DATA = OUT / "source_data"
PARTNERS = ("CaMKKalpha", "eNOS", "nNOS", "smMLCK")
OLD = {"CaMKKalpha": "CaMKK", "eNOS": "eNOS", "nNOS": "nNOS", "smMLCK": "smMLCK"}
BMRB = {"CaMKKalpha": 15183, "eNOS": 15185, "nNOS": 15191, "smMLCK": 15186}
AA = dict(zip("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split(), "ARNDCQEGHILKMFPSTWYV"))
KEY_METHYL = ((18, "CD1"), (39, "CD2"), (105, "CD1"), (105, "CD2"), (109, "CE"), (144, "CE"))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_csv(filename, rows):
    rows = list(rows)
    assert rows, filename
    path = DATA / filename
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return path


def load_wide(filename, expected_observable, expected_n):
    path = BASE / filename
    with path.open(encoding="utf-8-sig", newline="") as handle:
        source = list(csv.DictReader(handle))
    assert len(source) == expected_n
    rows = []
    for row in source:
        assert row["observable"] == expected_observable
        new = {"observable": expected_observable, "seq_id": int(row["seq_id"]), "amino_acid": row["amino_acid"], "atom": row["atom"]}
        new["probe_id"] = f'{new["amino_acid"]}{new["seq_id"]}-{new["atom"]}'
        for partner in PARTNERS:
            new[f"{partner}_s2"] = float(row[f"{OLD[partner]}_s2"])
            new[f"{partner}_fit_error"] = float(row[f"{OLD[partner]}_fit_error"])
            assert 0 <= new[f"{partner}_s2"] <= 1
            assert new[f"{partner}_fit_error"] >= 0
        rows.append(new)
    assert len({r["probe_id"] for r in rows}) == expected_n
    return rows


def make_long(rows):
    for row in rows:
        for order, partner in enumerate(PARTNERS):
            yield dict(observable=row["observable"], seq_id=row["seq_id"], amino_acid=row["amino_acid"], atom=row["atom"],
                       probe_id=row["probe_id"], partner=partner, partner_order=order, bmrb_id=BMRB[partner],
                       s2=row[f"{partner}_s2"], fit_error=row[f"{partner}_fit_error"])


def make_deltas(rows):
    result = []
    for row in rows:
        delta = row["smMLCK_s2"] - row["nNOS_s2"]
        result.append(dict(observable=row["observable"], seq_id=row["seq_id"], amino_acid=row["amino_acid"], atom=row["atom"],
                           probe_id=row["probe_id"], nNOS_s2=row["nNOS_s2"], nNOS_fit_error=row["nNOS_fit_error"],
                           smMLCK_s2=row["smMLCK_s2"], smMLCK_fit_error=row["smMLCK_fit_error"],
                           delta_smMLCK_minus_nNOS_s2=delta, abs_delta_s2=abs(delta),
                           fit_error_sum=row["smMLCK_fit_error"] + row["nNOS_fit_error"],
                           direction="higher_S2_more_restricted" if delta > 0 else "lower_S2_less_restricted" if delta < 0 else "unchanged_reported_value"))
    return result


def raw_qa():
    """Independently read raw loops, index by exact observable/position/AA/atom.

    Split only order-parameter loop fields, rather than using the prior analysis
    module. Duplicates are kept so ambiguity cannot be concealed by dict overwrite.
    """
    indexed = defaultdict(lambda: defaultdict(list))
    provenance = []
    for partner, number in BMRB.items():
        path = RAW / f"bmr{number}_3.str"
        text = path.read_text(encoding="utf-8")
        entity_blocks = re.findall(r"(?m)^save_\S+\s*\n([\s\S]*?)^save_\s*$", text)
        ids = []
        for block in entity_blocks:
            entity_name = re.search(r"(?m)^\s*_Entity.Name\s+(.+)$", block)
            entity_id = re.search(r"(?m)^\s*_Entity.ID\s+(\S+)", block)
            if entity_name and entity_id and "calmodulin" in entity_name[1].lower():
                ids.append(entity_id[1])
        assert len(ids) == 1
        raw_count = 0
        for body in re.findall(r"(?m)^\s*loop_\s*\n([\s\S]*?)^\s*stop_\s*$", text):
            lines = [line.strip() for line in body.splitlines() if line.strip()]
            tags = [line for line in lines if line.startswith("_")]
            if not tags or not tags[0].startswith("_Order_param."):
                continue
            tokens = shlex.split("\n".join(lines[len(tags):]))
            assert len(tokens) % len(tags) == 0
            fields = [tag.split(".", 1)[1] for tag in tags]
            for start in range(0, len(tokens), len(tags)):
                raw = dict(zip(fields, tokens[start:start + len(tags)]))
                if raw["Entity_ID"] != ids[0] or raw["Order_param_val"] in {".", "?"}:
                    continue
                key = ("backbone_NH" if raw["Atom_ID"] == "N" else "methyl_axis", int(raw["Comp_index_ID"]), AA[raw["Comp_ID"]], raw["Atom_ID"])
                indexed[key][partner].append((float(raw["Order_param_val"]), float(raw["Order_param_val_fit_err"])))
                raw_count += 1
        provenance.append(dict(partner=partner, bmrb_id=number, path=str(path), sha256=sha(path), raw_observation_count=raw_count))
    common = {key: values for key, values in indexed.items() if set(values) == set(PARTNERS) and all(len(values[p]) == 1 for p in PARTNERS)}
    return common, provenance


def summary(name, rows):
    deltas = [r["smMLCK_s2"] - r["nNOS_s2"] for r in rows]
    return dict(scope=name, observable=rows[0]["observable"], probe_count=len(rows), distinct_residue_count=len({r["seq_id"] for r in rows}),
                **{f"{p}_mean_s2": statistics.mean(r[f"{p}_s2"] for r in rows) for p in PARTNERS},
                nNOS_to_smMLCK_positive_delta_count=sum(d > 0 for d in deltas),
                nNOS_to_smMLCK_negative_delta_count=sum(d < 0 for d in deltas),
                nNOS_to_smMLCK_zero_delta_count=sum(d == 0 for d in deltas),
                mean_delta_s2=statistics.mean(deltas), median_abs_delta_s2=statistics.median(abs(d) for d in deltas),
                across_four_partners_median_range=statistics.median(max(r[f"{p}_s2"] for p in PARTNERS) - min(r[f"{p}_s2"] for p in PARTNERS) for r in rows))


def main():
    DATA.mkdir(parents=True, exist_ok=True)
    nh = load_wide("common_87_NH_probes_wide.csv", "backbone_NH", 87)
    methyl = load_wide("common_44_unambiguous_methyl_probes_wide.csv", "methyl_axis", 44)
    common_raw, raw_provenance = raw_qa()
    for row in nh + methyl:
        key = (row["observable"], row["seq_id"], row["amino_acid"], row["atom"])
        for p in PARTNERS:
            assert common_raw[key][p][0] == (row[f"{p}_s2"], row[f"{p}_fit_error"]), (key, p)
    assert Counter(key[0] for key in common_raw) == {"backbone_NH": 87, "methyl_axis": 44}

    files = []
    for label, rows in (("backbone_87_NH", nh), ("sidechain_44_methyl", methyl)):
        files.append(write_csv(f"{label}_wide.csv", rows))
        files.append(write_csv(f"{label}_long.csv", make_long(rows)))
        files.append(write_csv(f"{label}_delta_smMLCK_minus_nNOS.csv", make_deltas(rows)))

    nh_by_site = {r["seq_id"]: r for r in nh}
    methyl_sites = {r["seq_id"] for r in methyl}
    nh_matched = [r for r in nh if r["seq_id"] in methyl_sites]
    methyl_nh_matched = [r for r in methyl if r["seq_id"] in nh_by_site]
    files.append(write_csv("NH_at_common_methyl_residues_wide.csv", nh_matched))
    files.append(write_csv("methyl_with_corresponding_NH_wide.csv", methyl_nh_matched))

    key_methyl = [next(r for r in methyl if (r["seq_id"], r["atom"]) == key) for key in KEY_METHYL]
    key_nh = [nh_by_site[pos] for pos in sorted({key[0] for key in KEY_METHYL})]
    files.append(write_csv("key_sites_NH_and_methyl_wide.csv", key_nh + key_methyl))
    files.append(write_csv("key_sites_NH_and_methyl_long.csv", make_long(key_nh + key_methyl)))
    combined = []
    for met in key_methyl:
        backbone = nh_by_site[met["seq_id"]]
        for order, p in enumerate(PARTNERS):
            combined.append(dict(seq_id=met["seq_id"], amino_acid=met["amino_acid"], methyl_atom=met["atom"], methyl_probe_id=met["probe_id"],
                                 partner=p, partner_order=order, bmrb_id=BMRB[p], NH_s2=backbone[f"{p}_s2"], NH_fit_error=backbone[f"{p}_fit_error"],
                                 methyl_s2=met[f"{p}_s2"], methyl_fit_error=met[f"{p}_fit_error"]))
    files.append(write_csv("key_sites_matched_NH_methyl_long.csv", combined))

    summaries = [summary("all_87_common_NH_probes", nh), summary("all_44_common_methyl_probes", methyl),
                 summary("NH_at_common_methyl_residues_one_NH_per_residue", nh_matched),
                 summary("methyl_probes_with_corresponding_common_NH", methyl_nh_matched)]
    files.append(write_csv("descriptive_summary.csv", summaries))
    conditions = []
    for order, p in enumerate(PARTNERS):
        conditions.append(dict(partner=p, partner_order=order, bmrb_id=BMRB[p], temperature_K=308, pH=6.5,
                               imidazole_mM=10 if p == "CaMKKalpha" else 20, KCl_mM=100, CaCl2_mM=6,
                               sodium_azide="0.04%" if p == "CaMKKalpha" else "3 mM",
                               NH_proton_frequencies_MHz="500;750" if p == "CaMKKalpha" else "500;600",
                               methyl_proton_frequencies_MHz="500;750" if p == "CaMKKalpha" else "500;600",
                               source_url=f"https://bmrb.io/data_library/summary/index.php?bmrbId={BMRB[p]}"))
    files.append(write_csv("experimental_conditions.csv", conditions))

    delta_methyl = make_deltas(methyl)
    assert sum(r["delta_smMLCK_minus_nNOS_s2"] > 0 for r in delta_methyl) == 21
    assert sum(r["delta_smMLCK_minus_nNOS_s2"] < 0 for r in delta_methyl) == 23
    assert round(statistics.mean(r["nNOS_s2"] for r in methyl), 4) == 0.5459
    assert round(statistics.mean(r["smMLCK_s2"] for r in methyl), 4) == 0.5486
    key_deltas = make_deltas(key_methyl)
    assert abs(next(r["delta_smMLCK_minus_nNOS_s2"] for r in key_deltas if r["probe_id"] == "L39-CD2") + 0.450) < 1e-12
    assert abs(next(r["delta_smMLCK_minus_nNOS_s2"] for r in key_deltas if r["probe_id"] == "L105-CD1") - 0.305) < 1e-12

    audit = dict(status="passed", figure_type="experimental_NMR_data_only", model_prediction_files_read=[], p_values=[],
                 partner_order=list(PARTNERS), delta_definition="smMLCK S2 minus nNOS S2; positive = more restricted fast local orientation",
                 input_tables=[dict(path=str(BASE / f), sha256=sha(BASE / f)) for f in ("common_87_NH_probes_wide.csv", "common_44_unambiguous_methyl_probes_wide.csv")],
                 independent_raw_NMR_STAR_verification=dict(all_524_partner_probe_values_and_fit_errors_match=True, comparison_count=(87+44)*4, raw_inputs=raw_provenance),
                 summaries=summaries, key_site_methyl_deltas=key_deltas,
                 selected_site_note="Six displayed methyl probes represent five residues. CD1/CD2 at L105 are two methyl axes of the same residue, not independent residues.",
                 sampling_unit_note="87 NH probes cover 87 residues; 44 methyl probes cover 34 residues. Probe counts are not biological replicate counts.",
                 errorbar_note="Deposited fit_error is a model-fit error, not a confidence interval or full experimental uncertainty. Zero fit error does not imply no uncertainty.",
                 error_difference_note="Delta tables retain both endpoint fit errors and their descriptive sum; no propagated confidence interval or significance test is calculated.",
                 condition_note="nNOS/smMLCK deposited conditions and frequencies match; CaMKKalpha differs in buffer concentration, sodium azide declaration and proton frequencies.",
                 outputs=[dict(path=str(path), sha256=sha(path)) for path in files])
    (OUT / "numeric_audit.json").write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({"status": "passed", "summaries": summaries, "source_table_count": len(files), "key_site_methyl_deltas": key_deltas}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
