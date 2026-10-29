"""Run public ProAffinity-GNN weights with one documented upstream variable fix."""
from pathlib import Path
import hashlib
import json
import os
import runpy
import shutil
import sys
import time
from types import SimpleNamespace

import torch
import torch_geometric
from transformers import AutoTokenizer, EsmModel

ROOT = Path(__file__).resolve().parent
PDB = sys.argv[1].lower()
WORK = ROOT / "raw_outputs" / "proaffinity" / PDB
WORK.mkdir(parents=True, exist_ok=True)
upstream_source = ROOT / "vendor" / "ProAffinity-GNN_inference.py"
source = ROOT / "vendor" / "ProAffinity-GNN_inference_fixed_sequence_output.py"
source_text = upstream_source.read_text(encoding="utf-8")
bug = "for fasta in seq2:\n    input2 = tokenizer(fasta, return_tensors=\"pt\")\n    output1 = model(**input2)\n    last_hidden_state2 = output2.last_hidden_state"
fix = bug.replace("output1 = model(**input2)", "output2 = model(**input2)")
assert source_text.count(bug) == 1
source.write_text(source_text.replace(bug, fix), encoding="utf-8")
checkpoint = ROOT / "vendor" / "model.pkl"
input_file = ROOT / "structures" / f"{PDB}_crystal_protein.pdbqt"
shutil.copy2(input_file, WORK / f"{PDB}.pdbqt")
weight_link = WORK / "model.pkl"
if not weight_link.exists():
    os.link(checkpoint, weight_link)
torch.set_num_threads(2)
torch.manual_seed(0)
torch_load = torch.load
ESM_DIR = ROOT / "vendor" / "esm2_t33_650M_UR50D"
ESM_REVISION = "08e4846e537177426273712802403f7ba8261b6c"
tokenizer_load = AutoTokenizer.from_pretrained
esm_load = EsmModel.from_pretrained
esm_forward = EsmModel.forward
esm_cache_audit = []


def local_tokenizer(name, *args, **kwargs):
    if name == "facebook/esm2_t33_650M_UR50D":
        name = str(ESM_DIR)
        kwargs["local_files_only"] = True
    return tokenizer_load(name, *args, **kwargs)


def local_esm(name, *args, **kwargs):
    if name == "facebook/esm2_t33_650M_UR50D":
        name = str(ESM_DIR)
        kwargs["local_files_only"] = True
    return esm_load(name, *args, **kwargs)


def cached_esm_forward(self, *args, **kwargs):
    assert not self.training, "ESM must be in evaluation mode for deterministic inference caching"
    assert not args
    digest = hashlib.sha256(ESM_REVISION.encode())
    for key, value in sorted(kwargs.items()):
        digest.update(key.encode())
        digest.update(value.detach().cpu().numpy().tobytes())
    cache = ROOT / "esm_feature_cache" / f"{digest.hexdigest()}.pt"
    cache.parent.mkdir(exist_ok=True)
    hit = cache.exists()
    if hit:
        hidden = torch_load(cache, map_location="cpu", weights_only=True)
    else:
        hidden = esm_forward(self, **kwargs).last_hidden_state.detach().cpu()
        torch.save(hidden, cache)
    esm_cache_audit.append({"cache": str(cache), "hit": hit, "shape": list(hidden.shape)})
    return SimpleNamespace(last_hidden_state=hidden)


AutoTokenizer.from_pretrained = staticmethod(local_tokenizer)
EsmModel.from_pretrained = staticmethod(local_esm)
EsmModel.forward = cached_esm_forward


def cpu_weights_load(*args, **kwargs):
    kwargs["map_location"] = "cpu"
    kwargs["weights_only"] = True
    return torch_load(*args, **kwargs)


torch.load = cpu_weights_load
original_load = torch.nn.Module.load_state_dict
load_audit = []


def audited_state_dict_load(self, state_dict, *args, **kwargs):
    # The official script uses strict=True by default. Keep that behavior.
    result = original_load(self, state_dict, *args, **kwargs)
    load_audit.append({"missing_keys": result.missing_keys, "unexpected_keys": result.unexpected_keys})
    assert not result.missing_keys and not result.unexpected_keys
    return result


torch.nn.Module.load_state_dict = audited_state_dict_load
from torch_geometric.data import Data
original_data_validate = Data.validate
graph_checks = []
original_module_call = torch.nn.Module._call_impl

def validate_graph_network_call(self, *args, **kwargs):
    if self.__class__.__name__ == "GraphNetwork":
        for graph_name, graph in zip(["inter", "TCR", "pMHC"], args):
            assert graph.validate(raise_on_error=True)
            expected_nodes = sum(len(s) for s in (result_sequences[0] + result_sequences[1]
                                                  if graph_name == "inter" else
                                                  result_sequences[0] if graph_name == "TCR" else result_sequences[1]))
            assert graph.num_nodes == expected_nodes
            graph_checks.append({"name": graph_name, "nodes": graph.num_nodes,
                                 "edges": graph.num_edges, "valid": True})
    return original_module_call(self, *args, **kwargs)

# Node counts are fixed from the audited crystal chains before executing upstream code.
expected = next(r for r in json.loads((ROOT / "structure_audit.json").read_text(encoding="utf-8"))["records"] if r["pdb"] == PDB)
result_sequences = [[expected["chains"][c]["sequence"] for c in "DE"],
                    [expected["chains"][c]["sequence"] for c in "ABC"]]
torch.nn.Module._call_impl = validate_graph_network_call
os.chdir(WORK)
sys.argv = [str(source), "-f", f"{PDB}.pdbqt", "-c", "DE,ABC"]
started = time.perf_counter()
with torch.no_grad():
    result = runpy.run_path(str(source), run_name="__main__")
expected = next(r for r in json.loads((ROOT / "structure_audit.json").read_text(encoding="utf-8"))["records"] if r["pdb"] == PDB)
observed = dict(zip("DEABC", result["seqA"] + result["seqB"]))
for chain in "ABCDE":
    assert observed[chain] == expected["chains"][chain]["sequence"], (PDB, chain)
record = {
    "pdb": PDB, "model": "ProAffinity-GNN authors' simplified inference release",
    "upstream_commit": "0dd07bb0fb0da9ef1578ee7ed9587a536fbe4c7e",
    "torch_version": torch.__version__, "torch_geometric_version": torch_geometric.__version__,
    "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
    "upstream_source_sha256": hashlib.sha256(upstream_source.read_bytes()).hexdigest(),
    "upstream_bug_fix": "In the seq2 ESM loop, output1=model(**input2) replaced by output2=model(**input2), so pMHC uses its own chain ESM outputs rather than stale last-peptide output. One-variable correction; architecture and checkpoint unchanged.",
    "checkpoint_sha256": hashlib.sha256(checkpoint.read_bytes()).hexdigest(),
    "pdbqt_sha256": hashlib.sha256(input_file.read_bytes()).hexdigest(),
    "chain_groups": ["DE", "ABC"], "parsed_chain_lengths": {c: len(s) for c, s in observed.items()},
    "all_chain_sequences_match_input": True, "state_dict_loading": load_audit,
    "runtime_adapter": "CPU and weights_only=True checkpoint loading; frozen local ESM revision; inference under no_grad; reuse identical ESM token features in evaluation mode. Correct one upstream seq2 output variable; all graph node counts and edge indices validated. No architecture or learned parameters changed.",
    "graph_checks": graph_checks,
    "esm_revision": ESM_REVISION, "esm_feature_cache": esm_cache_audit,
    "prediction_as_reported_pKa": float(result["all_predictions"][0]),
    "elapsed_seconds": time.perf_counter() - started,
}
(WORK / "audit.json").write_text(json.dumps(record, indent=2), encoding="utf-8")
print("AUDIT_JSON=" + json.dumps(record), flush=True)
