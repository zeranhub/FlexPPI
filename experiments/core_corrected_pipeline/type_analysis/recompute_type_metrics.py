"""Reuse the original sidechain_rows function on the preserved final OOF arrays."""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True,
                        help="New result directory; existing results are never overwritten")
    args = parser.parse_args()
    if args.output.exists():
        parser.error("Choose a new --output directory")
    here = Path(__file__).resolve().parent
    core = here.parent
    export = core / "neural_exports"
    schema_path = core / "static_corrected" / "training_schema.json"
    schema = json.loads(schema_path.read_text(encoding="utf-8"))
    assert len(schema["records"]) == 3063
    expected = {}
    for row in schema["records"]:
        for variant in ("full", "intrachain"):
            cache = export / "evaluation" / variant / (row["sample_id"] + ".npz")
            expected[str(cache.resolve())] = sha256(cache)
    baseline = core / "backflip2_corrected_scoring" / "training_aligned.npz"
    baseline_hash = sha256(baseline)
    implementation = here / "prepare_and_plot_main_oof_typography_v2.py"
    spec = importlib.util.spec_from_file_location("preserved_type_production", implementation)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.PIPE = core
    context = {"schema": schema, "output_hashes": expected,
               "audit": {"input_and_consumed_output_sha256": {
                   str(baseline.resolve()): baseline_hash}}}
    evidence = {}
    # This is the preserved production function, including identity, common-mask,
    # coverage and paired-support assertions, with no new statistical formula.
    tables = module.sidechain_rows(context, export, evidence)
    names = ["all_record_amino_acid_metrics.csv.gz",
             "amino_acid_support_and_macro_metrics.csv",
             "all_record_amino_acid_paired_mae.csv.gz",
             "amino_acid_paired_mae_support_and_macros.csv"]
    args.output.mkdir(parents=True)
    for table, name in zip(tables, names):
        table.to_csv(args.output / name, index=False, na_rep="",
                     compression="gzip" if name.endswith(".gz") else None)
    audit = {"records": 3063, "function": "sidechain_rows",
             "production_code_sha256": sha256(implementation),
             "schema_sha256": sha256(schema_path),
             "input_sha256": evidence,
             "outputs": {name: {"rows": len(table), "sha256": sha256(args.output / name)}
                         for name, table in zip(names, tables)},
             "description": "Original final-cohort within-type MAE function; no training, inference or plotting"}
    (args.output / "recomputation_audit.json").write_text(
        json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8")
    print(args.output.resolve())


if __name__ == "__main__":
    main()
