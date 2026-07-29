"""Command-line runner for the paired FlexPPI/FastProtFlex experiments."""

from __future__ import annotations

import argparse
import json

from flexppi.comparison import compare_atlas_rmsf


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="configs/evaluate_state_local.yaml")
    parser.add_argument("--pretrain-checkpoint", default="artifacts/pretrain/best.pt")
    parser.add_argument("--finetune-checkpoint", default="artifacts/finetune/best.pt")
    parser.add_argument(
        "--output-dir", default="artifacts/evaluation/fastprotflex_comparison"
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
    parser.add_argument("--dataset", default="atlas", choices=("atlas", "dynarepo"))
    arguments = parser.parse_args()
    report = compare_atlas_rmsf(
        config_path=arguments.config,
        pretrain_checkpoint=arguments.pretrain_checkpoint,
        finetune_checkpoint=arguments.finetune_checkpoint,
        output_dir=arguments.output_dir,
        rscript_path=arguments.rscript,
        fastprotflex_batch_script=arguments.fastprotflex_batch_script,
        limit=arguments.limit,
        bootstrap_repetitions=arguments.bootstrap_repetitions,
        reuse_fastprotflex=arguments.reuse_fastprotflex,
        dataset_name=arguments.dataset,
    )
    print(json.dumps({
        "paired_evaluable_proteins": report["dataset"]["paired_evaluable_proteins"],
        "paired_evaluable_residues": report["dataset"]["paired_evaluable_residues"],
        "report": report["artifacts"],
    }, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
