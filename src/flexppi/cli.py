"""Command-line interface for indexing, labels, training, evaluation, and inference."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Dict, Optional, Sequence

from flexppi.config import load_config
from flexppi.data.atlas import prepare_atlas
from flexppi.data.dynarepo import index_dynarepo
from flexppi.data.labels import compute_dynarepo_labels
from flexppi.data.records import load_state_records, read_jsonl
from flexppi.data.zenodo_ppi import index_zenodo_ppi
from flexppi.evaluation import evaluate
from flexppi.inference import predict
from flexppi.tcr3d import evaluate_tcr3d, write_tcr3d_outputs
from flexppi.trainer import train


DEFAULT_ATLAS = r"E:\data\protein_md_datasets\ATLAS_extracted"
DEFAULT_DYNAREPO = r"D:\data\dynarepo_non_nucleic_full"
DEFAULT_ZENODO_PPI = r"E:\data\protein_md_datasets\zenodo_6638504_ppi_md"
DEFAULT_TCR3D = r"E:\Project\Data\TCR3d_database"


def _print_json(value: Any) -> None:
    print(json.dumps(value, indent=2, ensure_ascii=False, default=str))


def _index_atlas(args: argparse.Namespace) -> None:
    _print_json(prepare_atlas(args.root, args.output, args.limit))


def _index_dynarepo(args: argparse.Namespace) -> None:
    _print_json(
        index_dynarepo(
            args.root,
            args.output,
            label_dir=args.label_dir,
            limit=args.limit,
            minimum_identity=args.minimum_identity,
            minimum_coverage=args.minimum_coverage,
        )
    )


def _index_zenodo_ppi(args: argparse.Namespace) -> None:
    _print_json(
        index_zenodo_ppi(
            args.root,
            args.output,
            label_dir=args.label_dir,
            limit=args.limit,
        )
    )


def _compute_labels(args: argparse.Namespace) -> None:
    records = load_state_records(args.manifest)
    if args.groups:
        wanted = set(args.groups)
        records = [record for record in records if record.source_group in wanted]
    if args.limit is not None:
        records = records[: args.limit]
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    results: Dict[str, Any] = {"processed": 0, "failed": 0, "systems": {}}
    for record in records:
        output_path = output_dir / (record.sample_id + ".npz")
        try:
            summary = compute_dynarepo_labels(
                record,
                str(output_path),
                window_frames=args.window_frames,
                frame_stride=args.stride,
                include_sasa=args.sasa,
                max_windows=args.max_windows,
                sasa_points=args.sasa_points,
            )
            results["processed"] += 1
            results["systems"][record.sample_id] = summary
        except Exception as exc:
            results["failed"] += 1
            results["systems"][record.sample_id] = {"error": str(exc)}
            if not args.continue_on_error:
                _print_json(results)
                raise
    results["next_step"] = (
        "Re-run the matching dataset index command with --label-dir {} so state and "
        "pair records receive label paths."
        .format(output_dir)
    )
    _print_json(results)


def _train(args: argparse.Namespace) -> None:
    _print_json(train(load_config(args.config)))


def _evaluate(args: argparse.Namespace) -> None:
    thresholds: Dict[str, float] = {}
    for value in args.delta_threshold or []:
        if "=" not in value:
            raise ValueError("--delta-threshold must use metric=value, for example rmsf=0.1")
        name, raw_threshold = value.split("=", 1)
        thresholds[name.strip()] = float(raw_threshold)
    result = evaluate(
        load_config(args.config),
        args.checkpoint,
        split=args.split,
        limit=args.limit,
        bootstrap_repetitions=args.bootstrap,
        delta_thresholds=thresholds,
        interface_cutoff=args.interface_cutoff,
    )
    if args.output:
        output = Path(args.output)
        output.parent.mkdir(parents=True, exist_ok=True)
        with output.open("w", encoding="utf-8") as handle:
            json.dump(result, handle, indent=2, ensure_ascii=False)
        _print_json({"output": str(output.resolve()), "task": result["task"]})
    else:
        _print_json(result)


def _predict(args: argparse.Namespace) -> None:
    result = predict(
        args.checkpoint,
        args.apo_pdb,
        holo_pdb=args.holo_pdb,
        target_chains=args.target_chains,
        temperature=args.temperature,
        window_ns=args.window_ns,
        device_name=args.device,
    )
    if args.output:
        output = Path(args.output)
        output.parent.mkdir(parents=True, exist_ok=True)
        with output.open("w", encoding="utf-8") as handle:
            json.dump(result, handle, indent=2, ensure_ascii=False)
        _print_json({"output": str(output.resolve()), "rows": len(result["rows"])})
    else:
        _print_json(result)


def _inspect_manifest(args: argparse.Namespace) -> None:
    rows = list(read_jsonl(args.manifest))
    groups: Dict[str, int] = {}
    labelled = 0
    for row in rows:
        if "apo" in row:
            group = str(row.get("metadata", {}).get("pdb_id", "pair"))
            if row["apo"].get("label_path") and row["holo"].get("label_path"):
                labelled += 1
        else:
            group = str(row.get("source_group", "unknown"))
            if row.get("label_path"):
                labelled += 1
        groups[group] = groups.get(group, 0) + 1
    _print_json({"records": len(rows), "labelled": labelled, "groups": groups})


def _evaluate_tcr3d(args: argparse.Namespace) -> None:
    report, rows = evaluate_tcr3d(
        args.checkpoint,
        args.root,
        device_name=args.device,
        exclude_manifests=args.exclude_manifest,
        limit=args.limit,
        dataset_code=args.dataset_code,
        alternate_dataset_code=(
            None if args.no_condition_sensitivity else args.alternate_dataset_code
        ),
        temperature=args.temperature,
        window_ns=args.window_ns,
        bootstrap_repetitions=args.bootstrap,
    )
    outputs = write_tcr3d_outputs(report, rows, args.output_dir)
    _print_json({**outputs, "evaluated_chains": report["dataset"]["evaluated_chains"]})


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="flexppi")
    subparsers = parser.add_subparsers(dest="command", required=True)

    atlas = subparsers.add_parser("index-atlas", help="prepare ATLAS RMSF labels")
    atlas.add_argument("--root", default=DEFAULT_ATLAS)
    atlas.add_argument("--output", default="artifacts/data")
    atlas.add_argument("--limit", type=int)
    atlas.set_defaults(function=_index_atlas)

    dynarepo = subparsers.add_parser("index-dynarepo", help="index DynaRepo and build test53 pairs")
    dynarepo.add_argument("--root", default=DEFAULT_DYNAREPO)
    dynarepo.add_argument("--output", default="artifacts/data")
    dynarepo.add_argument("--label-dir")
    dynarepo.add_argument("--limit", type=int)
    dynarepo.add_argument("--minimum-identity", type=float, default=0.90)
    dynarepo.add_argument("--minimum-coverage", type=float, default=0.80)
    dynarepo.set_defaults(function=_index_dynarepo)

    zenodo_ppi = subparsers.add_parser(
        "index-zenodo-ppi",
        help="prepare Zenodo 6638504 strict apo/holo PPI trajectories",
    )
    zenodo_ppi.add_argument("--root", default=DEFAULT_ZENODO_PPI)
    zenodo_ppi.add_argument("--output", default="artifacts/zenodo_6638504/data")
    zenodo_ppi.add_argument("--label-dir")
    zenodo_ppi.add_argument("--limit", type=int)
    zenodo_ppi.set_defaults(function=_index_zenodo_ppi)

    labels = subparsers.add_parser("compute-dynarepo-labels", help="compute windowed MD labels")
    labels.add_argument("--manifest", default="artifacts/data/dynarepo_states.jsonl")
    labels.add_argument("--output-dir", default="artifacts/data/dynarepo_labels")
    labels.add_argument("--groups", nargs="+")
    labels.add_argument("--limit", type=int)
    labels.add_argument("--window-frames", type=int, default=1000)
    labels.add_argument("--stride", type=int, default=1)
    labels.add_argument("--max-windows", type=int)
    labels.add_argument("--sasa", action="store_true")
    labels.add_argument("--sasa-points", type=int, default=96)
    labels.add_argument("--continue-on-error", action="store_true")
    labels.set_defaults(function=_compute_labels)

    generic_labels = subparsers.add_parser(
        "compute-md-labels",
        help="compute windowed labels for any FlexPPI state manifest",
    )
    generic_labels.add_argument("--manifest", required=True)
    generic_labels.add_argument("--output-dir", required=True)
    generic_labels.add_argument("--groups", nargs="+")
    generic_labels.add_argument("--limit", type=int)
    generic_labels.add_argument("--window-frames", type=int, default=1000)
    generic_labels.add_argument("--stride", type=int, default=1)
    generic_labels.add_argument("--max-windows", type=int)
    generic_labels.add_argument("--sasa", action="store_true")
    generic_labels.add_argument("--sasa-points", type=int, default=96)
    generic_labels.add_argument("--continue-on-error", action="store_true")
    generic_labels.set_defaults(function=_compute_labels)

    train_parser = subparsers.add_parser("train")
    train_parser.add_argument("--config", required=True)
    train_parser.set_defaults(function=_train)

    evaluation = subparsers.add_parser("evaluate")
    evaluation.add_argument("--config", required=True)
    evaluation.add_argument("--checkpoint", required=True)
    evaluation.add_argument(
        "--split", choices=["all", "train", "val", "test"], default="all"
    )
    evaluation.add_argument("--limit", type=int)
    evaluation.add_argument("--bootstrap", type=int, default=200)
    evaluation.add_argument("--delta-threshold", action="append")
    evaluation.add_argument("--interface-cutoff", type=float, default=8.0)
    evaluation.add_argument("--output")
    evaluation.set_defaults(function=_evaluate)

    prediction = subparsers.add_parser("predict")
    prediction.add_argument("--checkpoint", required=True)
    prediction.add_argument("--apo-pdb", required=True)
    prediction.add_argument("--holo-pdb")
    prediction.add_argument("--target-chains", nargs="+")
    prediction.add_argument("--temperature", type=float, default=300.0)
    prediction.add_argument("--window-ns", type=float, default=100.0)
    prediction.add_argument("--device", default="auto")
    prediction.add_argument("--output")
    prediction.set_defaults(function=_predict)

    inspection = subparsers.add_parser("inspect-manifest")
    inspection.add_argument("--manifest", required=True)
    inspection.set_defaults(function=_inspect_manifest)

    tcr3d = subparsers.add_parser(
        "evaluate-tcr3d",
        help="external TCR3D B-factor/CDR proxy evaluation (not MD ground truth)",
    )
    tcr3d.add_argument("--checkpoint", required=True)
    tcr3d.add_argument("--root", default=DEFAULT_TCR3D)
    tcr3d.add_argument("--exclude-manifest", action="append", default=[])
    tcr3d.add_argument("--limit", type=int)
    tcr3d.add_argument("--device", default="auto")
    tcr3d.add_argument("--dataset-code", type=float, default=0.0)
    tcr3d.add_argument("--alternate-dataset-code", type=float, default=1.0)
    tcr3d.add_argument("--no-condition-sensitivity", action="store_true")
    tcr3d.add_argument("--temperature", type=float, default=300.0)
    tcr3d.add_argument("--window-ns", type=float, default=100.0)
    tcr3d.add_argument("--bootstrap", type=int, default=200)
    tcr3d.add_argument("--output-dir", default="artifacts/evaluation/tcr3d")
    tcr3d.set_defaults(function=_evaluate_tcr3d)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> None:
    args = build_parser().parse_args(argv)
    args.function(args)


if __name__ == "__main__":
    main()
