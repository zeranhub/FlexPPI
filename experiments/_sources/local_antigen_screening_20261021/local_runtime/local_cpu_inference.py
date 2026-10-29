"""CPU inference for the audited no-ITrace historical strict_full checkpoint.

The JSON manifest holds a list named ``candidates``. Each candidate requires
``candidate_id`` and ``input_pdb``; optional ``chain_sequences`` makes a clearly
labelled fixed-backbone sequence substitution. Coordinates are never relaxed in
that mode. No downloads or remote execution occur in this program.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import sys
import time
import traceback

ROOT = Path(__file__).resolve().parent
EXPECTED_CHECKPOINT_SHA256 = "b3eed060919df15f075329bc892f959d9355232be5feb561417edd49ea82b39f"
MODEL_IDENTITY = "historical strict_full 2431 no-ITrace; not family-fold Full 2430"
DEFAULT_SOURCE = Path("E:/Project/FlexPPI/src")
DEFAULT_CHECKPOINT = ROOT.parent / "model" / "epoch_040.pt"
ONE_TO_THREE = dict(zip("ACDEFGHIKLMNPQRSTVWY", [
    "ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU",
    "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR"]))
_DLL_HANDLES = []


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def write_json(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def import_runtime(source_root=DEFAULT_SOURCE, threads=2):
    """Keep Windows DLL search active; import the explicitly selected source tree."""
    dll = Path(sys.prefix) / "Library" / "bin"
    if os.name == "nt" and dll.is_dir():
        _DLL_HANDLES.append(os.add_dll_directory(str(dll)))
    os.environ.setdefault("HF_HUB_OFFLINE", "1")
    os.environ.setdefault("TRANSFORMERS_OFFLINE", "1")
    sys.path.insert(0, str(Path(source_root).resolve()))
    import numpy as np
    import torch
    import flexppi
    from flexppi.config import ModelConfig
    from flexppi.constants import CHI_ATOMS, RESTYPES
    from flexppi.data.dataset import collate_states
    from flexppi.data.structures import parse_pdb
    from flexppi.inference import _structure_example
    from flexppi.model.information_ablation import build_information_ablation
    torch.set_num_threads(int(threads))
    actual = Path(flexppi.__file__).resolve()
    if not actual.is_relative_to(Path(source_root).resolve()):
        raise RuntimeError("Wrong flexppi source loaded: " + str(actual))
    return locals()


def prepare_fixed_backbone(template_pdb, chain_sequences, output_pdb, runtime):
    """Synchronize residue identities, preserving every input backbone value.

    Mutated residues contain backbone atoms only: old side chains are deleted,
    rather than being relabelled as physically rebuilt mutant side chains.
    This is sufficient for this checkpoint, whose coordinate input is N/CA/C/O.
    """
    np, parse_pdb = runtime["np"], runtime["parse_pdb"]
    template_pdb, output_pdb = Path(template_pdb), Path(output_pdb)
    before = parse_pdb(str(template_pdb))
    changes = {}
    expected_sequence = list(before.sequence)
    for chain, sequence in chain_sequences.items():
        sequence = str(sequence).upper().strip()
        if not sequence or any(aa not in ONE_TO_THREE for aa in sequence):
            raise ValueError("Sequence must contain only the 20 standard amino acids: " + chain)
        indices = before.indices_for_chain(chain).tolist()
        if not indices:
            raise ValueError("Unknown chain: " + chain)
        if len(sequence) != len(indices):
            raise ValueError("Only same-length substitutions are supported: %s, %d != %d" % (chain, len(sequence), len(indices)))
        for i, new_aa in zip(indices, sequence):
            expected_sequence[i] = new_aa
            if before.sequence[i] != new_aa:
                key = (chain, before.author_residue_ids[i])
                changes[key] = {"chain": chain, "author_residue_id": key[1],
                    "residue_index": i, "from": before.sequence[i], "to": new_aa}
    lines = ["REMARK 900 FIXED-BACKBONE SEQUENCE SUBSTITUTION; NO MUTANT SIDE-CHAIN REBUILD\n",
             "REMARK 900 BACKBONE COORDINATES UNCHANGED; MUTATED SIDE CHAINS REMOVED\n"]
    for line in template_pdb.read_text(encoding="utf-8").splitlines(keepends=True):
        record = line[:6].strip()
        if record == "ENDMDL":
            break
        if record not in {"ATOM", "HETATM", "TER"}:
            continue
        key = (line[21:22].strip() or "_", line[22:26].strip() + line[26:27].strip())
        if key in changes:
            if record in {"ATOM", "HETATM"} and line[12:16].strip() not in {"N", "CA", "C", "O", "OXT"}:
                continue
            line = line[:17] + ONE_TO_THREE[changes[key]["to"]] + line[20:]
        lines.append(line if line.endswith("\n") else line + "\n")
    lines.append("END\n")
    output_pdb.write_text("".join(lines), encoding="utf-8")
    after = parse_pdb(str(output_pdb))
    assert after.sequence == "".join(expected_sequence)
    for name in ["backbone", "backbone_mask", "chain_index", "residue_index"]:
        assert np.array_equal(getattr(before, name), getattr(after, name)), name
    assert before.chain_ids == after.chain_ids
    assert before.author_residue_ids == after.author_residue_ids
    return {"input_mode": "fixed-backbone", "backbone_rebuilt": False,
        "sidechains_rebuilt": False, "old_mutated_sidechain_atoms_removed": True,
        "all_backbone_coordinates_bit_identical": True,
        "template_pdb": str(template_pdb.resolve()), "template_sha256": sha256(template_pdb),
        "prepared_pdb": str(output_pdb.resolve()), "prepared_pdb_sha256": sha256(output_pdb),
        "requested_chain_sequences": chain_sequences, "mutations": list(changes.values()),
        "sequence_after": {chain: after.chain_sequence(chain) for chain in after.unique_chains()},
        "esm_recomputed": "not applicable: learned sequence embedding, no ESM input"}


def load_audited_model(checkpoint, runtime):
    checkpoint = Path(checkpoint)
    actual_sha = sha256(checkpoint)
    if actual_sha != EXPECTED_CHECKPOINT_SHA256:
        raise ValueError("Refusing a different checkpoint: " + actual_sha)
    torch = runtime["torch"]
    data = torch.load(str(checkpoint), map_location="cpu", weights_only=False)
    config = data["config"]["model"]
    assert config["sequence_backend"] == "learned"
    assert config["use_condition_features"] is False
    assert config["enable_response_head"] is False
    assert config["information_ablation"] == "full"
    assert config["ipa_blocks"] == 4
    assert config["metrics"] == ["rmsf", "chi_entropy"]
    assert data["epoch"] == 40
    model = runtime["build_information_ablation"](runtime["ModelConfig"](**config), "full")
    model.load_state_dict(data["model"], strict=True)
    model.to("cpu").eval()
    return model, config


def run_manifest(manifest, output_dir, checkpoint=DEFAULT_CHECKPOINT, source_root=DEFAULT_SOURCE,
                 threads=2, limit=None, prepare_only=False):
    started = time.perf_counter()
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    runtime = import_runtime(source_root, threads)
    np, torch = runtime["np"], runtime["torch"]
    manifest_path = Path(manifest).resolve()
    payload = json.loads(manifest_path.read_text(encoding="utf-8-sig"))
    candidates = payload["candidates"] if isinstance(payload, dict) else payload
    if limit is not None:
        candidates = candidates[:int(limit)]
    ids = [row["candidate_id"] for row in candidates]
    if len(set(ids)) != len(ids):
        raise ValueError("candidate_id values must be unique")
    for candidate_id in ids:
        if not re.fullmatch(r"[A-Za-z0-9_.-]+", candidate_id):
            raise ValueError("Unsafe candidate_id: " + candidate_id)
    report = {"status": "running", "device": "cpu", "python": sys.executable,
        "torch": torch.__version__, "numpy": np.__version__, "threads": threads,
        "source_root": str(Path(source_root).resolve()),
        "source_sha256": {str(p.relative_to(source_root)): sha256(p)
            for p in sorted((Path(source_root)/"flexppi").rglob("*.py"))},
        "manifest": str(manifest_path), "manifest_sha256": sha256(manifest_path),
        "checkpoint": str(Path(checkpoint).resolve()), "expected_checkpoint_sha256": EXPECTED_CHECKPOINT_SHA256,
        "sequence_backend": "learned", "esm_used": False, "remote_compute": False,
        "model_identity": MODEL_IDENTITY,
        "prepare_only": prepare_only, "candidates": []}
    write_json(output_dir / "run_audit.json", report)
    model = config = None
    if not prepare_only:
        load_started = time.perf_counter()
        model, config = load_audited_model(checkpoint, runtime)
        report["model_load_seconds"] = time.perf_counter() - load_started
        report["checkpoint_sha256"] = sha256(checkpoint)
        report["model_config"] = config
    for row in candidates:
        case_started = time.perf_counter()
        candidate_id = row["candidate_id"]
        folder = output_dir / candidate_id
        folder.mkdir(exist_ok=True)
        path = Path(row["input_pdb"]).resolve()
        if row.get("chain_sequences"):
            prepared = folder / "model_input_fixed_backbone.pdb"
            preparation = prepare_fixed_backbone(path, row["chain_sequences"], prepared, runtime)
            path = prepared
        else:
            if row.get("input_mode", "observed_structure") == "fixed-backbone":
                raise ValueError("fixed-backbone requires explicit chain_sequences")
            preparation = {"input_mode": row.get("input_mode", "observed_structure"),
                "input_pdb": str(path), "input_sha256": sha256(path),
                "locally_prepared_by_this_program": False}
        write_json(folder / "input_preparation_audit.json", preparation)
        result = {"candidate_id": candidate_id, "input_pdb": str(path),
            "input_sha256": sha256(path), "preparation": preparation}
        if not prepare_only:
            st = runtime["parse_pdb"](str(path))
            ex = runtime["_structure_example"](str(path), row.get("target_chains"), 300.0, 200.0, 2.0)
            batch = runtime["collate_states"]([ex])
            infer_start = time.perf_counter()
            with torch.inference_mode():
                prediction = model.forward_state(batch)
            infer_seconds = time.perf_counter() - infer_start
            mean = prediction["mean"][0].numpy()
            log_variance = prediction["log_variance"][0].numpy()
            geometry_valid = prediction["mask"][0].numpy()
            assert mean.shape == log_variance.shape == (st.num_residues, 2)
            assert np.isfinite(mean).all() and np.isfinite(log_variance).all()
            hchi_defined = np.asarray([bool(runtime["CHI_ATOMS"].get(name)) for name in st.residue_names])
            target = ex["target_mask"].numpy()
            np.savez_compressed(folder / "prediction.npz", mean=mean, log_variance=log_variance,
                chain=np.asarray(st.chain_ids), author_residue_id=np.asarray(st.author_residue_ids),
                amino_acid=np.asarray(list(st.sequence)), geometry_valid=geometry_valid,
                target_mask=target, hchi_chemically_defined=hchi_defined,
                metrics=np.asarray(config["metrics"]), checkpoint_sha256=np.asarray(EXPECTED_CHECKPOINT_SHA256))
            with (folder / "prediction.csv").open("w", encoding="utf-8", newline="") as handle:
                fields = ["candidate_id", "residue_index", "chain", "author_residue_id", "amino_acid",
                    "geometry_valid", "target", "hchi_chemically_defined", "rmsf", "chi_entropy",
                    "rmsf_log_variance", "chi_entropy_log_variance"]
                writer = csv.DictWriter(handle, fieldnames=fields)
                writer.writeheader()
                for i in range(st.num_residues):
                    writer.writerow(dict(candidate_id=candidate_id, residue_index=i, chain=st.chain_ids[i],
                        author_residue_id=st.author_residue_ids[i], amino_acid=st.sequence[i],
                        geometry_valid=bool(geometry_valid[i]), target=bool(target[i]),
                        hchi_chemically_defined=bool(hchi_defined[i]), rmsf=float(mean[i, 0]),
                        chi_entropy=float(mean[i, 1]), rmsf_log_variance=float(log_variance[i, 0]),
                        chi_entropy_log_variance=float(log_variance[i, 1])))
            result.update(residues=st.num_residues, geometry_valid_count=int(geometry_valid.sum()),
                hchi_chemically_defined_count=int((geometry_valid & hchi_defined).sum()),
                inference_seconds=infer_seconds, all_predictions_finite=True,
                sequences_by_chain={ch: st.chain_sequence(ch) for ch in st.unique_chains()},
                mean_range_by_metric={name: [float(mean[:, j].min()), float(mean[:, j].max())]
                    for j, name in enumerate(config["metrics"])})
            if row.get("reference_prediction"):
                refpath = Path(row["reference_prediction"])
                with np.load(refpath, allow_pickle=False) as reference:
                    assert reference["mean"].shape == mean.shape
                    maxerr = float(np.max(np.abs(mean - reference["mean"])))
                    averageerr = float(np.mean(np.abs(mean - reference["mean"])))
                    for key, current in [("chain", st.chain_ids), ("author_residue_id", st.author_residue_ids),
                                         ("amino_acid", list(st.sequence))]:
                        if key in reference:
                            assert np.array_equal(reference[key].astype(str), np.asarray(current).astype(str)), key
                result["reference_reproduction"] = {"path": str(refpath), "sha256": sha256(refpath),
                    "max_absolute_error": maxerr, "mean_absolute_error": averageerr,
                    "tolerance": 2e-5, "passed": maxerr < 2e-5}
                if maxerr >= 2e-5:
                    raise AssertionError("Reference reproduction failed: " + str(maxerr))
        result["total_seconds"] = time.perf_counter() - case_started
        result["status"] = "passed"
        write_json(folder / "candidate_audit.json", result)
        report["candidates"].append(result)
        write_json(output_dir / "run_audit.json", report)
        print(json.dumps({"candidate_id": candidate_id, "status": "passed",
            "inference_seconds": result.get("inference_seconds"), "residues": result.get("residues")}), flush=True)
    report["elapsed_seconds"] = time.perf_counter() - started
    report["status"] = "passed"
    write_json(output_dir / "run_audit.json", report)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, default=DEFAULT_CHECKPOINT)
    parser.add_argument("--source-root", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--prepare-only", action="store_true")
    args = parser.parse_args()
    try:
        run_manifest(args.manifest, args.output_dir, args.checkpoint, args.source_root,
            args.threads, args.limit, args.prepare_only)
    except Exception as exc:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        write_json(args.output_dir / "failure.json", {"status": "failed", "exception": str(exc),
            "traceback": traceback.format_exc()})
        raise


if __name__ == "__main__":
    main()
