"""Run single-chain and multichain-context FlexPPI/FastProtFlex experiments."""

from __future__ import annotations

import argparse
import json

from flexppi.context_comparison import compare_chain_context


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="configs/evaluate_state_local.yaml")
    parser.add_argument(
        "--atlas-report",
        default="artifacts/evaluation/fastprotflex_comparison/atlas_comparison_report.json",
    )
    parser.add_argument(
        "--dynarepo-flex-predictions",
        default=(
            "artifacts/evaluation/fastprotflex_comparison/dynarepo/"
            "dynarepo_common_residue_predictions.csv"
        ),
    )
    parser.add_argument(
        "--dynarepo-fast-complex-predictions",
        default=(
            "artifacts/evaluation/fastprotflex_comparison/dynarepo/"
            "dynarepo_fastprotflex_predictions.csv"
        ),
    )
    parser.add_argument(
        "--dynarepo-fast-complex-timings",
        default=(
            "artifacts/evaluation/fastprotflex_comparison/dynarepo/"
            "dynarepo_fastprotflex_timings.csv"
        ),
    )
    parser.add_argument(
        "--output-dir",
        default="artifacts/evaluation/fastprotflex_comparison/chain_context",
    )
    parser.add_argument(
        "--rscript", default=r"C:\Program Files\R\R-4.6.0\bin\Rscript.exe"
    )
    parser.add_argument(
        "--fastprotflex-batch-script", default="scripts/fastprotflex_batch.R"
    )
    parser.add_argument("--interface-cutoff", type=float, default=8.0)
    parser.add_argument("--bootstrap-repetitions", type=int, default=2000)
    parser.add_argument("--reuse-fastprotflex-chain", action="store_true")
    arguments = parser.parse_args()
    report = compare_chain_context(
        config_path=arguments.config,
        atlas_report_path=arguments.atlas_report,
        dynarepo_flex_predictions=arguments.dynarepo_flex_predictions,
        dynarepo_fast_complex_predictions=arguments.dynarepo_fast_complex_predictions,
        dynarepo_fast_complex_timings=arguments.dynarepo_fast_complex_timings,
        output_dir=arguments.output_dir,
        rscript_path=arguments.rscript,
        fastprotflex_batch_script=arguments.fastprotflex_batch_script,
        interface_cutoff=arguments.interface_cutoff,
        bootstrap_repetitions=arguments.bootstrap_repetitions,
        reuse_fastprotflex_chain=arguments.reuse_fastprotflex_chain,
    )
    print(json.dumps({
        "single_chain_atlas_systems": report["single_chain"]["atlas"]["system_count"],
        "single_chain_dynarepo_systems": report["single_chain"]
        ["dynarepo_single_protein_chain"]["system_count"],
        "protein_complex_systems": report["protein_complex"]["system_count"],
        "protein_complex_chains": report["protein_complex"]["chain_count"],
    }, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
