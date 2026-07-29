"""Create a PyMOL comparison bundle for the 8PJG TCR-pMHC structure."""

from __future__ import annotations

import argparse
import json

from flexppi.pymol_visualization import create_pymol_comparison_bundle


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--pdb",
        default=(r"E:\Project\Data\TCR3d_database\TCR_complexes\TCR_complexes"
                 r"\8pjg.trunc.fit.pdb"),
    )
    parser.add_argument("--pdb-id", default="8PJG")
    parser.add_argument("--checkpoint", default="artifacts/pretrain/best.pt")
    parser.add_argument(
        "--output-dir",
        default="artifacts/evaluation/tcr_pmhc_pymol/8PJG",
    )
    parser.add_argument(
        "--rscript", default=r"C:\Program Files\R\R-4.6.0\bin\Rscript.exe"
    )
    parser.add_argument(
        "--fastprotflex-batch-script", default="scripts/fastprotflex_batch.R"
    )
    arguments = parser.parse_args()
    metadata = create_pymol_comparison_bundle(
        pdb_path=arguments.pdb,
        pdb_id=arguments.pdb_id,
        flexppi_checkpoint=arguments.checkpoint,
        output_dir=arguments.output_dir,
        rscript_path=arguments.rscript,
        fastprotflex_batch_script=arguments.fastprotflex_batch_script,
        chain_roles={
            "A": "MHC class II alpha",
            "B": "MHC class II beta",
            "C": "antigen peptide PKYVKQNTLKLAR",
            "D": "TCR alpha",
            "E": "TCR beta",
        },
    )
    print(json.dumps(metadata, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
