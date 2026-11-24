"""Command-line tools for structure prediction, training, and evaluation."""

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Any, Dict

from flexppi.config import load_config
from flexppi.evaluation import evaluate
from flexppi.inference import predict, predict_structure
from flexppi.trainer import train


def _print_json(result: Any) -> None:
    print(json.dumps(result, indent=2, ensure_ascii=False))


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
    if args.pdb:
        if args.holo_pdb:
            raise ValueError("--holo-pdb is used with --apo-pdb")
        rows = predict_structure(args.pdb, args.checkpoint, args.target_chains, args.device)
        if args.output:
            output = Path(args.output)
            output.parent.mkdir(parents=True, exist_ok=True)
            with output.open("w", encoding="utf-8", newline="") as handle:
                _write_table(rows, handle)
        else:
            _write_table(rows, sys.stdout)
        return
    if not args.checkpoint:
        raise ValueError("--apo-pdb requires a checkpoint; use --pdb for the bundled ensemble")
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


def _write_table(rows, handle) -> None:
    writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(prog="flexppi")
    subparsers = parser.add_subparsers(dest="command", required=True)
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
    prediction.add_argument("--checkpoint", help="Use a single inference checkpoint")
    structure = prediction.add_mutually_exclusive_group(required=True)
    structure.add_argument("--pdb", help="Input PDB file for residue-level CSV predictions")
    structure.add_argument("--apo-pdb", help="Input PDB for the single-checkpoint JSON interface")
    prediction.add_argument("--holo-pdb")
    prediction.add_argument("--target-chains", nargs="+")
    prediction.add_argument("--temperature", type=float, default=300.0)
    prediction.add_argument("--window-ns", type=float, default=100.0)
    prediction.add_argument("--device", default="auto")
    prediction.add_argument("--output")
    prediction.set_defaults(function=_predict)
    args = parser.parse_args()
    args.function(args)
