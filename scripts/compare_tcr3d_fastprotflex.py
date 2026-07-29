"""Run the shared TCR3D proxy comparison."""

from __future__ import annotations

import argparse
import json

from flexppi.comparison import compare_tcr3d_proxies


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--flexppi-predictions",
        default=(
            "artifacts/evaluation/tcr3d_full/tcr3d_residue_predictions.csv"
        ),
    )
    parser.add_argument("--tcr3d-root", default=r"E:\Project\Data\TCR3d_database")
    parser.add_argument(
        "--output-dir",
        default="artifacts/evaluation/fastprotflex_comparison/tcr3d",
    )
    parser.add_argument(
        "--rscript", default=r"C:\Program Files\R\R-4.6.0\bin\Rscript.exe"
    )
    parser.add_argument(
        "--fastprotflex-batch-script", default="scripts/fastprotflex_batch.R"
    )
    parser.add_argument("--limit", type=int)
    parser.add_argument("--bootstrap-repetitions", type=int, default=2000)
    parser.add_argument("--reuse-fastprotflex", action="store_true")
    arguments = parser.parse_args()
    report = compare_tcr3d_proxies(
        flexppi_residue_predictions=arguments.flexppi_predictions,
        tcr3d_root=arguments.tcr3d_root,
        output_dir=arguments.output_dir,
        rscript_path=arguments.rscript,
        fastprotflex_batch_script=arguments.fastprotflex_batch_script,
        limit=arguments.limit,
        bootstrap_repetitions=arguments.bootstrap_repetitions,
        reuse_fastprotflex=arguments.reuse_fastprotflex,
    )
    print(json.dumps({
        "dataset": report["dataset"],
        "artifacts": report["artifacts"],
    }, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
