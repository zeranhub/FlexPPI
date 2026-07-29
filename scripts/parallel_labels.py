"""Parallel DynaRepo label computation using multiprocessing."""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from flexppi.data.labels import compute_dynarepo_labels
from flexppi.data.records import load_state_records


def process_one(record, output_dir, window_frames, frame_stride, sasa, max_windows, sasa_points, overwrite):
    output_path = Path(output_dir) / (record.sample_id + ".npz")
    if output_path.exists() and not overwrite:
        return {"sample_id": record.sample_id, "skipped": True}
    try:
        summary = compute_dynarepo_labels(
            record,
            str(output_path),
            window_frames=window_frames,
            frame_stride=frame_stride,
            include_sasa=sasa,
            max_windows=max_windows,
            sasa_points=sasa_points,
        )
        return {"sample_id": record.sample_id, "summary": summary, "error": None}
    except Exception as exc:
        return {"sample_id": record.sample_id, "summary": None, "error": str(exc)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", default="artifacts/data/dynarepo_states.jsonl")
    parser.add_argument("--output-dir", default="artifacts/data/dynarepo_labels")
    parser.add_argument("--groups", nargs="+", default=None)
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--window-frames", type=int, default=1000)
    parser.add_argument("--stride", type=int, default=1)
    parser.add_argument("--max-windows", type=int, default=None)
    parser.add_argument("--sasa", action="store_true")
    parser.add_argument("--sasa-points", type=int, default=96)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    records = load_state_records(args.manifest)
    if args.groups:
        wanted = set(args.groups)
        records = [r for r in records if r.source_group in wanted]
    if args.limit:
        records = records[:args.limit]

    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    completed = sum(1 for r in records if (output_dir / (r.sample_id + ".npz")).exists())
    total = len(records)
    print(f"Total: {total}, already computed: {completed}, remaining: {total - completed}")
    if completed == total:
        print("All labels already computed.")
        return

    print(f"Starting parallel label computation with {args.workers} workers...")
    t0 = time.time()
    successful = 0
    failed = 0
    skipped = 0

    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures = {}
        for record in records:
            if (output_dir / (record.sample_id + ".npz")).exists() and not args.overwrite:
                skipped += 1
                continue
            future = executor.submit(
                process_one,
                record,
                str(output_dir),
                args.window_frames,
                args.stride,
                args.sasa,
                args.max_windows,
                args.sasa_points,
                args.overwrite,
            )
            futures[future] = record.sample_id

        for future in as_completed(futures):
            result = future.result()
            sid = result["sample_id"]
            if result.get("skipped"):
                skipped += 1
            elif result["error"]:
                failed += 1
                print(f"  FAILED {sid}: {result['error']}")
            else:
                successful += 1
                s = result["summary"]
                elapsed = time.time() - t0
                rate = successful / max(elapsed, 1)
                remaining = total - completed - successful - failed
                eta = remaining / max(rate, 0.001) if rate > 0 else 0
                print(f"  [{successful}] {sid}: {s} | rate={rate:.2f}/s | ETA={eta/60:.1f}min")

    elapsed = time.time() - t0
    print(f"\nDone. Processed {successful}, failed {failed}, skipped {skipped} in {elapsed/60:.1f} min")
    print(f"Next: flexppi index-dynarepo --root /ranze/data/dynarepo/dynarepo_non_nucleic_full --output artifacts/data --label-dir {output_dir}")


if __name__ == "__main__":
    main()
