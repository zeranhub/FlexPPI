"""Released PLANET-MD full/seq/mini, one exact ESM3 embedding per unique sequence.

Native CPU float32 only. This script has no training, label access, cropping,
temperature fitting, output clipping, or mixed precision. --mode inspect is
tensor-free. Smoke verifies the RMSF-only execution path against the complete
official forward before a full benchmark is permitted.
"""
import argparse
import csv
import gc
import hashlib
import json
import os
from pathlib import Path
import sys
import threading
import time
import traceback
from datetime import datetime, timezone

from inspect_checkpoint import inspect

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
COHORT = ROOT / "cohort"
SOURCE = Path("E:/Project/protein_flexibility_models/PLANET-MD")
HF_HOME = Path("E:/Project/protein_flexibility_models/hf_cache")
OLD_WEIGHTS = HF_HOME / "hub/models--samsl--rocketshp/snapshots/7831d318182d4f4fc1acbdafdbf23862d8bc1657/checkpoints"
CHECKPOINTS = {
    "full": OLD_WEIGHTS / "rocketshp_pretrained_20250505_v1.ckpt",
    "seq": OLD_WEIGHTS / "rocketshp_pretrained_20250507_v1_seq.ckpt",
    "mini": ROOT / "weights/planet/checkpoints/rocketshp_pretrained_20250505_v1_mini.ckpt",
}
ESM_SNAPSHOT = HF_HOME / "hub/models--EvolutionaryScale--esm3-sm-open-v1/snapshots/47f0545b2b6daf26a93439a3cd610f4f7f3d5478"
ESM_WEIGHTS = ESM_SNAPSHOT / "data/weights/esm3_sm_open_v1.pth"
STRUCTURE_WEIGHTS = ESM_SNAPSHOT / "data/weights/esm3_structure_encoder_v0.pth"
ESM_PACKAGE = Path("E:/Project/protein_flexibility_models/envs/planet_md_py312/Lib/site-packages/esm")
AA = dict(zip("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split(), "ARNDCQEGHILKMFPSTWYV"))


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024**2), b""):
            digest.update(block)
    return digest.hexdigest()


def json_write(path, value):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(str(path) + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, ensure_ascii=False), encoding="utf-8")
    temporary.replace(path)


def csv_rows(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def prepare_metadata():
    checkpoints = {}
    for variant, path in CHECKPOINTS.items():
        info = inspect(path)
        info["sha256"] = sha256(path)
        assert info["hyper_parameters"]["embedding_dim"] == 1536
        assert info["hyper_parameters"]["seq_features"] is True
        assert info["hyper_parameters"]["struct_features"] == (variant == "full")
        assert info["default_temp_K"] == 300.0
        if info["hyper_parameters"]["struct_features"]:
            assert info["hyper_parameters"]["struct_stage"] in {"encoded", "pre-quantized", "quantized"}
        checkpoints[variant] = info
    original_freeze = ROOT.parent / "baselines/cohort"
    freeze = json.loads((original_freeze / "cohort_freeze.json").read_text(encoding="utf-8"))
    assert sha256(original_freeze / "chain_manifest.csv") == freeze["files"]["chain_manifest.csv"]
    chain_entries = csv_rows(COHORT / "chains.csv")
    records = {row["sample_id"]: row for row in csv_rows(COHORT / "records.csv")}
    original_entries = {row["input_id"]: row for row in csv_rows(original_freeze / "chain_manifest.csv")}
    for entry in chain_entries:
        entry["dataset"] = records[entry["sample_id"]]["dataset"]
        entry["record_structure_path"] = records[entry["sample_id"]]["structure_path"]
        entry["record_structure_sha256"] = records[entry["sample_id"]]["structure_sha256"]
        original = original_entries[entry["input_id"]]
        assert entry["sequence"] == original["sequence"]
        assert entry["structure_sha256"] == original["input_sha256"]
        assert entry["chain"] == original["chain"]
        assert len(entry["sequence"]) == int(entry["n_residues"])
    assert len(chain_entries) == len(original_entries) == 4973
    source_files = [
        "planet_md/modeling/architectures.py", "planet_md/features.py",
        "planet_md/esm3.py", "planet_md/structure/protein_chain.py",
        "planet_md/structure/normalize_coordinates.py", "LICENSE",
    ]
    metadata = {
        "model_family": "PLANET-MD_released_full_seq_mini", "cpu_threads": 2,
        "device": "cpu", "dtype": "float32", "temperature_K": 300.0,
        "checkpoint_metadata": checkpoints,
        "cohort_chains_sha256": sha256(COHORT / "chains.csv"),
        "cohort_records_sha256": sha256(COHORT / "records.csv"),
        "original_frozen_chain_manifest_sha256": sha256(original_freeze / "chain_manifest.csv"),
        "upstream_source_sha256": {name: sha256(SOURCE / name) for name in source_files},
        "esm_package_path": str(ESM_PACKAGE),
        "esm_implementation_sha256": {str(path.relative_to(ESM_PACKAGE)).replace("\\", "/"): sha256(path) for path in sorted(ESM_PACKAGE.rglob("*.py"))},
        "esm3_weight_path": str(ESM_WEIGHTS), "esm3_weight_bytes": ESM_WEIGHTS.stat().st_size,
        "structure_encoder_weight_path": str(STRUCTURE_WEIGHTS), "structure_encoder_weight_bytes": STRUCTURE_WEIGHTS.stat().st_size,
        "input": "Each frozen canonical chain independently; all original chain/author-residue identities and coordinates verified against complete structure.",
        "sequence_embedding_cache": "Exact float32 .npy per amino-acid-sequence SHA256; shared without modification between released heads. Cache metadata binds ESM weight/source hashes.",
        "checkpoint_loading": "Official architectures, mmap read-only released tensor storages, strict load_state_dict(assign=True) onto meta-created models; regenerate only untrained nonpersistent rotary constants by their official reset_parameters.",
        "RMSF_execution": "Official model._transform followed by the identical temperature append and model.rmsf_head. Skip independent SHP and pairwise heads; smoke requires bitwise equality to model(features)['rmsf'].",
        "RMSF_units": "Released native nm multiplied by exactly 10 to Angstrom; no fitting or clipping.",
        "Hchi_supported": False,
        "limitations": "Independent-chain input does not encode partners. Released training data may overlap benchmark sources. SHP output is not joint chi-state entropy and is not evaluated as Hchi.",
        "expected_chains": len(chain_entries), "unique_sequences": len({e["sequence"] for e in chain_entries}),
        "unique_sequence_residues": sum(len(s) for s in {e["sequence"] for e in chain_entries}),
    }
    return chain_entries, metadata


def selected_entries(entries, mode):
    if mode == "full":
        # Group equal sequences to reuse ESM without repeated disk reads. Every
        # chain still receives its own structure feature and identity mapping.
        return sorted(entries, key=lambda item: (len(item["sequence"]), hashlib.sha256(item["sequence"].encode()).hexdigest(), item["input_id"]))
    selected = []
    for source in ("atlas", "dynarepo", "itrace"):
        candidates = [e for e in entries if e["dataset"] == source]
        if source == "atlas":
            # A released-weight CPU result already exists for this exact sequence.
            # It independently checks the memory-efficient ESM loading path;
            # source PDB formatting is irrelevant to the sequence-only head.
            selected.append(next(e for e in candidates if e["sample_id"] == "1a62_A"))
            continue
        if source == "dynarepo":
            candidates = [e for e in candidates if e["main_group_candidate"] == "other_multichain"]
        ordered = sorted(int(e["n_residues"]) for e in candidates)
        median_size = ordered[len(ordered) // 2]
        selected.append(min(candidates, key=lambda item: (abs(int(item["n_residues"]) - median_size), item["sample_id"], item["input_id"])))
    return selected


def read_atoms(path):
    rows = []
    for line in Path(path).read_text(encoding="ascii").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        residue = line[17:20].strip()
        if residue not in AA or line[76:78].strip().upper() in {"H", "D"}:
            continue
        assert not line[16:17].strip(), (str(path), "alternate atom location")
        rows.append({"chain": line[21:22].strip() or "_", "author_id": str(int(line[22:26])) + line[26:27].strip(),
                     "atom": line[12:16].strip(), "aa": AA[residue],
                     "xyz": (float(line[30:38]), float(line[38:46]), float(line[46:54]))})
    return rows


def validate_input(entry, np):
    assert sha256(entry["structure_path"]) == entry["structure_sha256"]
    assert sha256(entry["record_structure_path"]) == entry["record_structure_sha256"]
    canonical_atoms = read_atoms(entry["structure_path"])
    original_atoms = read_atoms(entry["record_structure_path"])
    canonical_ca = [atom for atom in canonical_atoms if atom["atom"] == "CA"]
    original_all_ca = [atom for atom in original_atoms if atom["atom"] == "CA"]
    original_ca = [atom for atom in original_all_ca if atom["chain"] == entry["chain"]]
    sequence = entry["sequence"]
    assert len(canonical_ca) == len(original_ca) == len(sequence)
    assert "".join(atom["aa"] for atom in canonical_ca) == sequence
    assert "".join(atom["aa"] for atom in original_ca) == sequence
    assert {atom["chain"] for atom in canonical_atoms} == {"A"}
    assert [atom["author_id"] for atom in canonical_ca] == [str(i) for i in range(1, len(sequence) + 1)]
    original_ids = [atom["author_id"] for atom in original_ca]
    assert len(set(original_ids)) == len(original_ids)
    original_index = {author: index for index, author in enumerate(original_ids)}
    canon_index = {str(index + 1): index for index in range(len(sequence))}
    source_atoms = {(original_index[a["author_id"]], a["atom"]): a["xyz"] for a in original_atoms if a["chain"] == entry["chain"]}
    new_atoms = {(canon_index[a["author_id"]], a["atom"]): a["xyz"] for a in canonical_atoms}
    assert len(source_atoms) == sum(atom["chain"] == entry["chain"] for atom in original_atoms)
    assert len(new_atoms) == len(canonical_atoms)
    assert source_atoms.keys() == new_atoms.keys(), (entry["input_id"], "atom mapping differs")
    assert all(max(abs(x - y) for x, y in zip(source_atoms[key], new_atoms[key])) < 1e-9 for key in source_atoms)
    return {
        "input_index": np.arange(len(sequence), dtype=np.int64),
        "chain": np.full(len(sequence), entry["chain"], dtype=str),
        "author_residue_id": np.asarray(original_ids, dtype=str),
        "amino_acid": np.asarray(list(sequence), dtype=str),
        "record_residue_index": np.asarray([i for i, atom in enumerate(original_all_ca) if atom["chain"] == entry["chain"]], dtype=np.int64),
    }


def restore_rotary_constants(model, torch):
    from esm.layers.rotary import RotaryEmbedding
    for module in model.modules():
        if isinstance(module, RotaryEmbedding):
            module.device = torch.device("cpu")
            module.reset_parameters()
    assert all(parameter.device.type == "cpu" and parameter.dtype == torch.float32 for parameter in model.parameters())
    assert all(buffer.device.type == "cpu" for buffer in model.buffers())


def load_head(path, torch, PlanetMD):
    checkpoint = torch.load(path, map_location="cpu", weights_only=True, mmap=True)
    hp = checkpoint["hyper_parameters"]
    with torch.device("meta"):
        model = PlanetMD(hp["embedding_dim"], hp["output_dim"], hp["d_model"], hp["n_heads"], hp["n_layers"],
                         seq_only=not hp["struct_features"], struct_stage=hp.get("struct_stage", "quantized"), struct_dim=hp.get("struct_dim"))
    state = {key.replace("child_model.", ""): value for key, value in checkpoint["state_dict"].items()}
    model.load_state_dict(state, strict=True, assign=True)
    restore_rotary_constants(model, torch)
    assert model.default_temp.item() == 300.0
    model.eval()
    del checkpoint, state
    gc.collect()
    return model


def load_esm(torch, tokenizers):
    from esm.models.esm3 import ESM3
    from esm.pretrained import ESM3_function_decoder_v0, ESM3_structure_decoder_v0, ESM3_structure_encoder_v0
    with torch.device("meta"):
        model = ESM3(d_model=1536, n_heads=24, v_heads=256, n_layers=48,
                     structure_encoder_fn=ESM3_structure_encoder_v0, structure_decoder_fn=ESM3_structure_decoder_v0,
                     function_decoder_fn=ESM3_function_decoder_v0, tokenizers=tokenizers)
    state = torch.load(ESM_WEIGHTS, map_location="cpu", weights_only=True, mmap=True)
    model.load_state_dict(state, strict=True, assign=True)
    restore_rotary_constants(model, torch)
    model.eval()
    del state
    gc.collect()
    return model


def rmsf_only(model, features, torch):
    # These operations and their order are exactly the RMSF path in the official
    # PlanetMDModel.forward. The other outputs never feed back into this branch.
    temperature = torch.ones(features["seq_feats"].shape[1], device=features["seq_feats"].device).unsqueeze(0) * model.default_temp
    transformed = model._transform(features)
    return model.rmsf_head(torch.cat([transformed, temperature.unsqueeze(-1)], dim=-1))


def wait_for_ram(required_gib, maximum_wait_seconds, mode):
    import psutil
    started = time.perf_counter()
    started_utc = datetime.now(timezone.utc).isoformat()
    assert required_gib > 0 and maximum_wait_seconds >= 0
    while psutil.virtual_memory().available < required_gib * 1024**3:
        waited = time.perf_counter() - started
        status = {"status": "waiting_before_ESM_load", "mode": mode, "started_utc": started_utc,
                  "required_available_gib": required_gib, "available_gib": psutil.virtual_memory().available / 1024**3,
                  "wait_seconds": waited, "maximum_wait_seconds": maximum_wait_seconds}
        if waited >= maximum_wait_seconds:
            status["status"] = "blocked_insufficient_RAM_before_model_load"
            status["reason"] = "Available physical RAM did not reach the scheduling threshold before the finite waiting deadline; no models were loaded."
            json_write(HERE / "memory_status.json", status)
            print(json.dumps(status), flush=True)
            return status
        json_write(HERE / "memory_status.json", status)
        print(json.dumps(status), flush=True)
        time.sleep(min(30.0, maximum_wait_seconds - waited))
    status = {"status": "RAM_gate_ready", "mode": mode, "started_utc": started_utc,
              "required_available_gib": required_gib, "available_gib": psutil.virtual_memory().available / 1024**3,
              "wait_seconds": time.perf_counter() - started, "maximum_wait_seconds": maximum_wait_seconds}
    json_write(HERE / "memory_status.json", status)
    return status


class ResourceSampler:
    """Observe actual RSS/available RAM without changing inference settings."""
    def __init__(self, mode):
        import psutil
        self.psutil = psutil
        self.process = psutil.Process()
        self.path = HERE / (mode + "_resource_usage.json")
        self.peak_rss = 0
        self.minimum_available = psutil.virtual_memory().available
        self.samples = 0
        self.finished = threading.Event()
        self.started = time.perf_counter()
        self.thread = threading.Thread(target=self.run, daemon=True)
        self.thread.start()

    def snapshot(self):
        return {"peak_process_rss_bytes": self.peak_rss, "peak_process_rss_gib": self.peak_rss / 1024**3,
                "minimum_system_available_bytes": self.minimum_available,
                "minimum_system_available_gib": self.minimum_available / 1024**3,
                "sampling_interval_seconds": 0.1, "samples": self.samples,
                "elapsed_seconds": time.perf_counter() - self.started}

    def run(self):
        last_write = -5.0
        while not self.finished.is_set():
            self.peak_rss = max(self.peak_rss, self.process.memory_info().rss)
            self.minimum_available = min(self.minimum_available, self.psutil.virtual_memory().available)
            self.samples += 1
            elapsed = time.perf_counter() - self.started
            if elapsed - last_write >= 5.0:
                json_write(self.path, {"status": "sampling", **self.snapshot()})
                last_write = elapsed
            self.finished.wait(0.1)

    def stop(self):
        self.finished.set()
        self.thread.join(timeout=2.0)
        result = {"status": "complete", **self.snapshot()}
        json_write(self.path, result)
        return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=["inspect", "smoke", "full"], required=True)
    parser.add_argument("--minimum-free-gib", type=float, default=6.0)
    parser.add_argument("--memory-wait-seconds", type=float, default=600.0)
    args = parser.parse_args()
    entries, metadata = prepare_metadata()
    selected = selected_entries(entries, args.mode)
    (HERE / "inspect_protocol.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    with (HERE / (args.mode + "_inputs.csv")).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(selected[0]))
        writer.writeheader()
        writer.writerows(selected)
    if args.mode == "inspect":
        print(json.dumps({"status": "static_inspection_only", "metadata": metadata,
                          "selected_smoke_chains": [{key: entry[key] for key in ("input_id", "sample_id", "chain", "n_residues", "dataset")} for entry in selected]}, indent=2), flush=True)
        return
    metadata["wrapper_sha256"] = sha256(__file__)
    metadata["esm3_weights_sha256"] = sha256(ESM_WEIGHTS)
    metadata["structure_encoder_weights_sha256"] = sha256(STRUCTURE_WEIGHTS)
    proof_path = HERE / "smoke_audit.json"
    if args.mode == "full":
        proof = json.loads(proof_path.read_text(encoding="utf-8"))
        assert proof["status"] == "complete" and proof["heads_verified_against_official_forward"] == ["full", "seq", "mini"]
        assert proof["chains_requested"] == 3 and proof["official_forward_equivalence_checks"] == 9
        assert proof["sources_verified"] == ["atlas", "dynarepo", "itrace"]
        assert proof["metadata"] == metadata, "Run smoke again after any source, weight, wrapper, or cohort change"
    metadata_path = HERE / "inference_protocol.json"
    if metadata_path.exists():
        assert json.loads(metadata_path.read_text(encoding="utf-8")) == metadata
    else:
        json_write(metadata_path, metadata)
    for name in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        os.environ[name] = "2"
    os.environ["HF_HOME"] = str(HF_HOME)
    os.environ["HF_HUB_OFFLINE"] = "1"
    os.environ["TRANSFORMERS_OFFLINE"] = "1"
    os.environ["HF_HUB_DISABLE_TELEMETRY"] = "1"
    memory_gate = wait_for_ram(args.minimum_free_gib, args.memory_wait_seconds, args.mode)
    if memory_gate["status"].startswith("blocked"):
        json_write(HERE / (args.mode + "_audit.json"), {
            "status": memory_gate["status"], "metadata": metadata, "memory_gate": memory_gate,
            "chains_requested": len(selected), "chains_attempted": 0,
            "heads_verified_against_official_forward": [], "official_forward_equivalence_checks": 0,
        })
        return 3
    resources = ResourceSampler(args.mode)
    sys.path.insert(0, str(SOURCE))
    import numpy as np
    import torch
    from biotite.structure.io import pdb
    from planet_md import PlanetMD
    from planet_md.esm3 import get_structure_vae, get_tokenizers
    from planet_md.features import esm3_sequence, esm3_vqvae
    from planet_md.structure.protein_chain import ProteinChain
    import esm
    assert Path(esm.__file__).resolve().parent == ESM_PACKAGE.resolve()
    torch.set_num_threads(2)
    torch.set_num_interop_threads(1)
    torch.set_float32_matmul_precision("highest")
    startup = time.perf_counter()
    models = {variant: load_head(path, torch, PlanetMD) for variant, path in CHECKPOINTS.items()}
    tokenizers = get_tokenizers()
    esm_model = load_esm(torch, tokenizers)
    structure_model = get_structure_vae(device="cpu").eval()
    assert all(parameter.dtype == torch.float32 and parameter.device.type == "cpu" for parameter in structure_model.parameters())
    runtime_environment = {"torch": torch.__version__, "numpy": np.__version__,
                           "threads": torch.get_num_threads(), "interop_threads": torch.get_num_interop_threads(),
                           "float32_matmul_precision": torch.get_float32_matmul_precision(), "python": sys.version}
    runtime_environment_path = HERE / "runtime_environment.json"
    if args.mode == "full":
        assert json.loads(runtime_environment_path.read_text(encoding="utf-8")) == runtime_environment, "Runtime changed since smoke"
    else:
        json_write(runtime_environment_path, runtime_environment)
    print(json.dumps({"status": "models_loaded", "model_load_seconds": time.perf_counter() - startup,
                      "torch": torch.__version__, "parameters_by_head": {variant: sum(p.numel() for p in model.parameters()) for variant, model in models.items()}}), flush=True)
    cache = HERE / "sequence_cache"
    cache.mkdir(exist_ok=True)
    cache_protocol = {key: metadata[key] for key in ("esm3_weights_sha256", "upstream_source_sha256", "esm_implementation_sha256", "dtype", "cpu_threads")}
    cache_protocol_path = cache / "protocol.json"
    if cache_protocol_path.exists():
        assert json.loads(cache_protocol_path.read_text(encoding="utf-8")) == cache_protocol
    else:
        json_write(cache_protocol_path, cache_protocol)
    results = []
    for position, entry in enumerate(selected, 1):
        statuses = {}
        started = time.perf_counter()
        try:
            identity = validate_input(entry, np)
            pending = []
            for variant in CHECKPOINTS:
                status_path = HERE / variant / (entry["input_id"] + ".json")
                output = status_path.with_suffix(".npz")
                if status_path.exists() and output.exists():
                    old = json.loads(status_path.read_text(encoding="utf-8"))
                    if old.get("status") == "complete" and old.get("input_sha256") == entry["structure_sha256"]:
                        assert old["output_sha256"] == sha256(output)
                        if args.mode == "full" or old.get("official_forward_bitwise_equal"):
                            statuses[variant] = old
                            continue
                pending.append(variant)
            if pending:
                sequence = entry["sequence"]
                sequence_hash = hashlib.sha256(sequence.encode()).hexdigest()
                embedding_path = cache / (sequence_hash + ".npy")
                embedding_status_path = cache / (sequence_hash + ".json")
                cache_hit = embedding_path.exists() and embedding_status_path.exists()
                if cache_hit:
                    embedded = json.loads(embedding_status_path.read_text(encoding="utf-8"))
                    assert embedded["sequence"] == sequence and embedded["embedding_sha256"] == sha256(embedding_path)
                    array = np.load(embedding_path, allow_pickle=False)
                    assert array.dtype == np.float32 and array.shape == (1, len(sequence), 1536)
                    seq_features = torch.from_numpy(array)
                else:
                    embedding_started = time.perf_counter()
                    with torch.inference_mode():
                        seq_features = esm3_sequence(sequence, esm_model, tokenizers)[:, 1:-1, :]
                    assert seq_features.shape == (1, len(sequence), 1536) and seq_features.dtype == torch.float32
                    assert torch.isfinite(seq_features).all()
                    array = seq_features.cpu().numpy()
                    temporary = cache / (sequence_hash + ".tmp.npy")
                    np.save(temporary, array, allow_pickle=False)
                    temporary.replace(embedding_path)
                    embedded = {"sequence": sequence, "sequence_sha256": sequence_hash,
                                "embedding_sha256": sha256(embedding_path), "compute_seconds": time.perf_counter() - embedding_started,
                                "shape": list(array.shape), "dtype": "float32"}
                    json_write(embedding_status_path, embedded)
                structure_features = {}
                structure_seconds = 0.0
                needed_stages = {metadata["checkpoint_metadata"][v]["hyper_parameters"]["struct_stage"] for v in pending if v == "full"}
                if needed_stages:
                    structure_started = time.perf_counter()
                    structure = pdb.PDBFile.read(entry["structure_path"]).get_structure(model=1)
                    chain = ProteinChain.from_atomarray(structure)
                    assert chain.sequence == sequence
                    with torch.inference_mode():
                        for stage in needed_stages:
                            structure_features[stage] = esm3_vqvae(chain, structure_model, stage=stage)
                            assert structure_features[stage].shape[1] == len(sequence)
                    structure_seconds = time.perf_counter() - structure_started
                for variant in pending:
                    variant_started = time.perf_counter()
                    hp = metadata["checkpoint_metadata"][variant]["hyper_parameters"]
                    features = {"seq_feats": seq_features, "struct_feats": structure_features.get(hp["struct_stage"]) if hp["struct_features"] else None}
                    with torch.inference_mode():
                        prediction = rmsf_only(models[variant], features, torch)
                        head_compute_seconds = time.perf_counter() - variant_started
                        delta = None
                        bitwise = None
                        if args.mode == "smoke":
                            official = models[variant](features)["rmsf"]
                            delta = float(torch.max(torch.abs(official - prediction)))
                            bitwise = bool(torch.equal(official, prediction))
                            assert bitwise, (variant, "RMSF-only path differs from official forward", delta)
                    native = prediction.reshape(-1).cpu().numpy()
                    assert native.shape == (len(sequence),) and np.isfinite(native).all() and native.dtype == np.float32
                    legacy_delta = None
                    legacy = ROOT.parent / "baselines/planet_seq_cpu_smoke_20261010"
                    legacy_output = legacy / (entry["input_id"] + ".npz")
                    if args.mode == "smoke" and variant == "seq" and legacy_output.exists():
                        legacy_protocol = json.loads((legacy / "metadata.json").read_text(encoding="utf-8"))
                        legacy_status = json.loads(legacy_output.with_suffix(".json").read_text(encoding="utf-8"))
                        assert legacy_protocol["weights_sha256"] == metadata["checkpoint_metadata"]["seq"]["sha256"]
                        assert legacy_status["output_sha256"] == sha256(legacy_output)
                        with np.load(legacy_output, allow_pickle=False) as old:
                            assert old["sequence"].item() == sequence
                            legacy_delta = float(np.max(np.abs(old["rmsf_angstrom"] - 10.0 * native)))
                            assert np.allclose(old["rmsf_angstrom"], 10.0 * native, rtol=1e-5, atol=1e-4), ("legacy official CPU sequence smoke differs", legacy_delta)
                    directory = HERE / variant
                    directory.mkdir(exist_ok=True)
                    output = directory / (entry["input_id"] + ".npz")
                    temporary = directory / (entry["input_id"] + ".tmp.npz")
                    np.savez_compressed(temporary, rmsf_native_nm=native, rmsf_angstrom=10.0 * native,
                                        sequence=np.asarray(sequence), input_id=np.asarray(entry["input_id"]), **identity)
                    temporary.replace(output)
                    status = {"status": "complete", "variant": variant, "input_id": entry["input_id"], "sample_id": entry["sample_id"],
                              "chain": entry["chain"], "input_sha256": entry["structure_sha256"], "record_structure_sha256": entry["record_structure_sha256"],
                              "n_residues": len(sequence), "output_sha256": sha256(output), "sequence_embedding_cache_hit": cache_hit,
                              "sequence_embedding_compute_seconds": embedded["compute_seconds"], "structure_features_seconds": structure_seconds if hp["struct_features"] else 0.0,
                              "head_RMSF_path_seconds": head_compute_seconds, "head_validation_and_write_wall_seconds": time.perf_counter() - variant_started, "temperature_K": 300.0,
                              "structure_feature_stage": hp["struct_stage"] if hp["struct_features"] else None,
                              "all_original_heavy_atom_identities_and_coordinates_verified": True,
                              "official_forward_bitwise_equal": bitwise, "official_forward_max_abs_delta_nm": delta,
                              "legacy_official_CPU_sequence_max_abs_delta_Angstrom": legacy_delta,
                              "native_unit": "nm", "output_unit": "Angstrom", "Hchi_supported": False}
                    json_write(directory / (entry["input_id"] + ".json"), status)
                    statuses[variant] = status
        except Exception as error:
            failure = {"status": "failed", "input_id": entry["input_id"], "sample_id": entry["sample_id"],
                       "error": repr(error), "traceback": traceback.format_exc()}
            for variant in CHECKPOINTS:
                if variant not in statuses:
                    json_write(HERE / variant / (entry["input_id"] + ".json"), failure)
                    statuses[variant] = failure
        result = {"input_id": entry["input_id"], "sample_id": entry["sample_id"], "seconds": time.perf_counter() - started,
                  "heads": {variant: statuses[variant]["status"] for variant in CHECKPOINTS}}
        results.append(result)
        json_write(HERE / (args.mode + "_progress.json"), {"position": position, "requested": len(selected), "latest": result})
        print(json.dumps({"progress": f"{position}/{len(selected)}", **result}), flush=True)
        gc.collect()
    complete = all(all(status == "complete" for status in row["heads"].values()) for row in results)
    json_write(HERE / (args.mode + "_audit.json"), {
        "status": "complete" if complete else "some_failures", "metadata": metadata,
        "chains_requested": len(selected), "heads_verified_against_official_forward": list(CHECKPOINTS) if args.mode == "smoke" and complete else [],
        "sources_verified": [entry["dataset"] for entry in selected] if args.mode == "smoke" and complete else [],
        "official_forward_equivalence_checks": 9 if args.mode == "smoke" and complete else 0,
        "memory_gate": memory_gate, "resources": resources.stop(),
        "results": results,
    })
    return 0 if complete else 2


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as error:
        json_write(HERE / "fatal_error.json", {"status": "failed", "time_utc": datetime.now(timezone.utc).isoformat(),
                   "error": repr(error), "traceback": traceback.format_exc()})
        raise
