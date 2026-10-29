"""Label-free CPU adapters for the frozen, separately measured deployment pass.

API
---
Set FLEXPPI_DEPLOYMENT_PAYLOAD_ROOT before importing this module when the
primary frozen payload is not ``training_payload`` beside this file. Import
timing belongs to the worker, not this adapter. The worker must verify input,
code, schema, config, and weight hashes outside its deployment timer.

context = prepare_context(method, checkpoint_paths=[...], config_paths=[...])
context = prepare_context(method, static_paths=[{"rmsf": ..., "chi_entropy": ...}],
                          schema_path=...)
adapter = NeuralStaticAdapter(method, context)  # worker times model loading
arrays = adapter.predict_record(record)        # worker times this + NPZ save

Methods: full_single, full_five, intrachain_single, strict_full_single,
strict_intrachain_single, static0, static5, strict_static. Five-member lists
must be in fold order 0..4. Every call reparses the PDB; static calls rebuild
all 47 features. No labels, trajectory files, cached features, cached
predictions, hashes, audit files, output files, or accuracy scores are read or
written by predict_record. Model loading reads only declared checkpoints.

record is a dict with sample_id, absolute/local structure_path, target_chains,
and optional metadata.residue_count. Other record fields are ignored. A
record may additionally supply expected_identity arrays for a cheap equality
check; the worker's frozen 32-record identity preflight is authoritative.

Static regressors supply means only. They do not supply an uncertainty model;
model_variance is omitted rather than fabricated as zeros. Neural variance is
uncalibrated: exp(log_variance) for one member; mean within-member variance
plus population variance of the member means for the five-fold mixture.
"""

from pathlib import Path
import copy
import json
import os
import sys

# Fresh worker processes are the intended calling contract. No CUDA operation.
os.environ["CUDA_VISIBLE_DEVICES"] = ""
os.environ["OMP_NUM_THREADS"] = "2"
os.environ["MKL_NUM_THREADS"] = "2"
os.environ["OPENBLAS_NUM_THREADS"] = "2"
sys.dont_write_bytecode = True

import joblib
import numpy as np
import sklearn
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.utils.validation import check_is_fitted
import torch
from threadpoolctl import threadpool_info, threadpool_limits
import yaml

HERE = Path(__file__).resolve().parent
PAYLOAD_ROOT = Path(os.environ.get(
    "FLEXPPI_DEPLOYMENT_PAYLOAD_ROOT", str(HERE / "training_payload"))).resolve()
SOURCE_ROOT = PAYLOAD_ROOT / "src"
if not SOURCE_ROOT.is_dir():
    raise FileNotFoundError("Set FLEXPPI_DEPLOYMENT_PAYLOAD_ROOT to the frozen primary payload")
sys.path.insert(0, str(SOURCE_ROOT))
sys.path.insert(1, str(HERE))
from flexppi.config import ModelConfig
from flexppi.constants import METRIC_NAMES
from flexppi.data.dataset import collate_states
from flexppi.data.structures import parse_pdb, target_mask
from flexppi.model.information_ablation import build_information_ablation
import flexppi.data.structures as _structure_module
import static_feature_definitions as _feature_module

if not Path(_structure_module.__file__).resolve().is_relative_to(SOURCE_ROOT):
    raise RuntimeError("A different flexppi package was already imported; use a fresh worker")
if Path(_feature_module.__file__).resolve() != HERE / "static_feature_definitions.py":
    raise RuntimeError("The frozen local 47-feature module was shadowed")

torch.set_num_threads(2)
try:
    torch.set_num_interop_threads(1)
except RuntimeError:
    if torch.get_num_interop_threads() != 1:
        raise RuntimeError("Use a fresh CPU worker to set PyTorch interop threads to 1")

METRICS = ("rmsf", "chi_entropy")
if tuple(METRIC_NAMES[:2]) != METRICS:
    raise RuntimeError("Frozen dataset metric order is incompatible")

METHODS = {
    "full_single": dict(kind="neural", variant="full", folds=[0], strict=False),
    "full_five": dict(kind="neural", variant="full", folds=list(range(5)), strict=False),
    "intrachain_single": dict(kind="neural", variant="intrachain", folds=[0], strict=False),
    "strict_full_single": dict(kind="neural", variant="full", folds=[0], strict=True),
    "strict_intrachain_single": dict(kind="neural", variant="intrachain", folds=[0], strict=True),
    "static0": dict(kind="static", folds=[0], strict=False),
    "static5": dict(kind="static", folds=list(range(5)), strict=False),
    "strict_static": dict(kind="static", folds=[0], strict=True),
}

STATIC_CONFIGURATION = dict(
    loss="squared_error", learning_rate=0.07, max_iter=120,
    max_leaf_nodes=31, max_depth=None, min_samples_leaf=50,
    l2_regularization=1.0, max_bins=63, early_stopping=False, random_state=7)


def _require(condition, message):
    if not condition:
        raise ValueError(message)


def _path_list(values, count, name):
    _require(values is not None and not isinstance(values, (str, Path)),
             name + " must be an ordered list of paths")
    paths = [Path(p).resolve() for p in values]
    _require(len(paths) == count, name + " has wrong member count")
    _require(len(set(paths)) == count, name + " repeats a member")
    for path in paths:
        _require(path.is_file(), "Missing " + name + ": " + str(path))
    return paths


def _check_neural_config(config, specification, fold):
    _require(isinstance(config, dict), "Checkpoint/config must be a mapping")
    _require(config.get("stage") == "pretrain", "Only absolute-state pretrained members are accepted")
    _require(config.get("seed") == 7, "Cost protocol fixes seed7")
    _require(config.get("init_checkpoint") is None, "A from-scratch member is required")
    model = config["model"]
    _require(tuple(model.get("metrics", [])) == METRICS, "Model metrics must be RMSF then Hchi")
    _require(model.get("information_ablation") == specification["variant"], "Wrong neural arm")
    _require(model.get("sequence_backend") == "learned", "Unexpected sequence backend")
    _require(model.get("use_condition_features") is False, "Condition features are outside the frozen protocol")
    _require(model.get("enable_response_head") is False, "Response heads are outside the frozen protocol")
    _require(config["train"]["epochs"] == 40, "Cost member must use the fixed 40-epoch protocol")
    data = config["data"]
    _require(data["batch_size"] == 1 and data["max_residues"] == 1024,
             "Expected batch1 and the frozen 1024-residue limit")
    stem = "strict_" + specification["variant"] if specification["strict"] else specification["variant"]
    output = str(config["output_dir"]).replace("\\", "/").rstrip("/")
    _require(output.endswith("runs/%s/fold_%d" % (stem, fold)), "Wrong checkpoint fold/output identity")
    train_name = str(data["state_manifest"]).replace("\\", "/").split("/")[-1]
    expected_train = "strict_external_train.jsonl" if specification["strict"] else "fold_%d_train.jsonl" % fold
    _require(train_name == expected_train, "Wrong declared training-member manifest")


def prepare_context(method, *, checkpoint_paths=None, config_paths=None,
                    static_paths=None, schema_path=None):
    """Prepare config/schema contracts outside model-load and prediction timers.

    This function does not hash files, load labels, or load model weights.
    The worker must bind these paths and source files to its frozen hash plan.
    """
    _require(method in METHODS, "Unknown deployment method: " + str(method))
    specification = copy.deepcopy(METHODS[method])
    count = len(specification["folds"])
    context = dict(method=method, specification=specification, max_residues=1024,
                   metrics=METRICS, payload_root=str(PAYLOAD_ROOT))
    if specification["kind"] == "neural":
        _require(static_paths is None and schema_path is None, "Neural context must not include static files")
        checkpoints = _path_list(checkpoint_paths, count, "checkpoint_paths")
        configs = _path_list(config_paths, count, "config_paths")
        config_data = []
        for fold, path in zip(specification["folds"], configs):
            value = yaml.safe_load(path.read_text(encoding="utf-8"))
            _check_neural_config(value, specification, fold)
            config_data.append(value)
        # Fold members must have identical architecture, optimizer and masks.
        if count > 1:
            for config in config_data[1:]:
                for key in ["model", "train", "seed", "stage", "init_checkpoint"]:
                    _require(config[key] == config_data[0][key], "Five-member contracts differ: " + key)
        context.update(checkpoint_paths=checkpoints, config_paths=configs, configs=config_data)
    else:
        _require(checkpoint_paths is None and config_paths is None, "Static context must not include neural files")
        _require(isinstance(static_paths, (list, tuple)) and len(static_paths) == count,
                 "static_paths must have one metric-path dict per declared fold")
        _require(schema_path is not None, "Static feature schema is required")
        schema_file = Path(schema_path).resolve()
        schema = json.loads(schema_file.read_text(encoding="utf-8"))
        names = list(schema["names"])
        _require(len(names) == 47 and len(set(names)) == 47, "Static schema must define 47 distinct features")
        _require(list(schema["sets"]["sequence_complex"]) == names, "Wrong static feature set/order")
        _require(list(schema["categorical"]) == names[:9], "Wrong categorical feature order")
        members = []
        for fold, member in zip(specification["folds"], static_paths):
            _require(set(member) == set(METRICS), "Every static member needs exactly two named metric models")
            paths = {metric: Path(member[metric]).resolve() for metric in METRICS}
            for metric, path in paths.items():
                _require(path.is_file(), "Missing static model: " + str(path))
                expected = "strict_static_%s.joblib" % metric if specification["strict"] else "static_fold%d_%s.joblib" % (fold, metric)
                _require(path.name == expected, "Wrong declared static fold/metric filename")
            members.append(paths)
        _require(len({p for member in members for p in member.values()}) == 2 * count,
                 "Static model paths repeat a member")
        context.update(static_paths=members, schema_path=schema_file,
                       feature_names=names, categorical=list(schema["categorical"]))
    return context


class NeuralStaticAdapter:
    """Constructor loads declared fitted models only; predict never uses labels."""

    def __init__(self, method, context):
        _require(method in METHODS and context.get("method") == method,
                 "Use prepare_context for this exact method")
        _require(context["specification"] == METHODS[method], "Method specification changed")
        _require(Path(context["payload_root"]).resolve() == PAYLOAD_ROOT, "Frozen source root differs")
        self.method = method
        self.context = context
        self.specification = context["specification"]
        self.members = []
        # Also constrain already-loaded NumPy/sklearn native thread pools.
        self._thread_limit = threadpool_limits(limits=2)
        torch.set_num_threads(2)
        _require(torch.get_num_threads() == 2 and torch.get_num_interop_threads() == 1,
                 "CPU thread contract was not established")
        if self.specification["kind"] == "neural":
            for fold, path, expected in zip(self.specification["folds"],
                                            context["checkpoint_paths"], context["configs"]):
                checkpoint = torch.load(path, map_location="cpu", weights_only=False)
                _require(checkpoint["epoch"] == 40, "Wrong final checkpoint epoch")
                actual = checkpoint["config"]
                _check_neural_config(actual, self.specification, fold)
                for key in ["model", "data", "train", "seed", "stage", "init_checkpoint", "output_dir"]:
                    _require(actual[key] == expected[key], "Checkpoint differs from frozen config: " + key)
                for tensor in checkpoint["model"].values():
                    _require(not tensor.is_floating_point() or tensor.dtype == torch.float32,
                             "Neural weights must use the frozen float32 dtype")
                model = build_information_ablation(ModelConfig(**actual["model"]),
                                                     self.specification["variant"]).cpu().eval()
                model.load_state_dict(checkpoint["model"], strict=True)
                model.requires_grad_(False)
                _require(all(p.device.type == "cpu" for p in model.parameters()), "Unexpected model device")
                self.members.append(model)
                del checkpoint
        else:
            _require(sklearn.__version__ == "1.6.1", "Load the fitted static controls with their scikit-learn 1.6.1 runtime")
            expected_categorical = np.asarray([name in context["categorical"]
                                               for name in context["feature_names"]])
            for paths in context["static_paths"]:
                models = {}
                for metric in METRICS:
                    model = joblib.load(paths[metric])
                    _require(isinstance(model, HistGradientBoostingRegressor), "Wrong static estimator type")
                    check_is_fitted(model)
                    _require(model.n_features_in_ == 47 and model.n_iter_ == 120,
                             "Wrong fitted static feature count or iteration budget")
                    for key, value in STATIC_CONFIGURATION.items():
                        _require(model.get_params()[key] == value, "Static model contract differs: " + key)
                    _require(np.array_equal(model.is_categorical_, expected_categorical),
                             "Static fitted categorical features differ")
                    models[metric] = model
                self.members.append(models)
        _require(len(self.members) == len(self.specification["folds"]), "Incomplete fitted members")

    def metadata(self):
        """Runtime facts; no accuracy or benchmark claim and no file hashing."""
        neural = self.specification["kind"] == "neural"
        return dict(method=self.method, device="cpu", torch_threads=torch.get_num_threads(),
                    torch_interop_threads=torch.get_num_interop_threads(),
                    native_thread_pools=threadpool_info(), member_count=len(self.members),
                    metric_names=list(METRICS), float_dtype="float32 neural; float64 static/ensemble conversion",
                    torch_version=torch.__version__, numpy_version=np.__version__, sklearn_version=sklearn.__version__,
                    parameters_per_member=[sum(p.numel() for p in model.parameters())
                                           for model in self.members] if neural else None,
                    parameter_bytes_per_member=[sum(p.numel() * p.element_size() for p in model.parameters())
                                                for model in self.members] if neural else None,
                    variance_available=neural, payload_root=str(PAYLOAD_ROOT),
                    variance_scope="Uncalibrated fitted-head variance plus population between-member variance" if neural else "No fitted variance output from static regressors",
                    scope="Per-call PDB parse, input preparation, fresh forward, output conversion; worker includes validation and serialization")

    def _parse_record(self, record):
        _require(isinstance(record, dict), "record must be a dict with a localized PDB path")
        _require(bool(record.get("sample_id")), "Missing sample_id")
        path = Path(record["structure_path"])
        _require(path.is_absolute(), "Worker must localize structure_path to an absolute path")
        # This read and parse occur on EVERY call, including repeated cost rounds.
        structure = parse_pdb(str(path))
        n = structure.num_residues
        _require(0 < n <= self.context["max_residues"], "Structure violates the frozen length limit")
        expected_n = record.get("metadata", {}).get("residue_count")
        if expected_n is not None:
            _require(int(expected_n) == n, "Declared and parsed residue counts differ")
        _require(np.isfinite(structure.backbone).all(), "Nonfinite input backbone")
        keys = list(zip(structure.chain_ids, structure.author_residue_ids))
        _require(len(set(keys)) == n, "Duplicate input chain/author residue identity")
        target_chains = record.get("target_chains", [])
        _require(isinstance(target_chains, (list, tuple)), "target_chains must be a list")
        wanted = {str(chain).strip() or "_" for chain in target_chains}
        _require(wanted.issubset(set(structure.chain_ids)), "Declared target chain absent from parsed input")
        target = target_mask(structure, target_chains)
        _require(target.any(), "No target residues")
        identity = dict(
            sample_id=np.asarray(record["sample_id"]), metric_names=np.asarray(METRICS),
            amino_acid=np.asarray(list(structure.sequence)), residue_name=np.asarray(structure.residue_names),
            chain=np.asarray(structure.chain_ids), author_residue_id=np.asarray(structure.author_residue_ids),
            residue_index=structure.residue_index.copy(), chain_index=structure.chain_index.copy(),
            aatype=structure.aatype.copy(), geometry_mask=structure.frame_mask.copy(), target_mask=target)
        for name, expected in record.get("expected_identity", {}).items():
            _require(name in identity, "Unknown expected identity field: " + name)
            _require(np.array_equal(identity[name], np.asarray(expected)), "Input identity differs: " + name)
        return structure, identity

    def predict_record(self, record):
        """Fresh label-free prediction; caller owns timing, validation and NPZ IO."""
        structure, output = self._parse_record(record)
        n = structure.num_residues
        if self.specification["kind"] == "neural":
            # Frozen collate_states needs supervision-shaped fields even though
            # forward_state does not read them. Empty masks/NaN placeholders
            # satisfy that interface without loading or inventing any labels.
            columns = len(METRIC_NAMES)
            example = dict(
                sample_id=record["sample_id"], dataset="deployment", source_group="deployment", metadata={},
                aatype=torch.from_numpy(structure.aatype), residue_index=torch.from_numpy(structure.residue_index),
                chain_index=torch.from_numpy(structure.chain_index), backbone=torch.from_numpy(structure.backbone),
                backbone_mask=torch.from_numpy(structure.backbone_mask), geometry_mask=torch.from_numpy(structure.frame_mask),
                target_mask=torch.from_numpy(output["target_mask"]), condition=torch.zeros(3, dtype=torch.float32),
                label_values=torch.full((n, columns), float("nan"), dtype=torch.float32),
                label_variances=torch.zeros((n, columns), dtype=torch.float32),
                label_mask=torch.zeros((n, columns), dtype=torch.bool))
            batch = collate_states([example])
            means, variances = [], []
            with torch.inference_mode():
                for model in self.members:
                    result = model.forward_state(batch)
                    _require(torch.equal(result["mask"][0], example["geometry_mask"]),
                             "Model changed the prespecified geometry population")
                    means.append(result["mean"][0].cpu().numpy().copy())
                    variances.append(result["log_variance"][0].exp().cpu().numpy().copy())
            member_means = np.stack(means).astype(np.float64)
            member_variances = np.stack(variances).astype(np.float64)
            _require(member_means.shape == member_variances.shape == (len(self.members), n, 2),
                     "Wrong neural prediction dimensions")
            _require(np.isfinite(member_means).all() and np.isfinite(member_variances).all()
                     and (member_variances >= 0).all(), "Nonfinite neural prediction")
            output["mean"] = member_means.mean(axis=0)
            within = member_variances.mean(axis=0)
            between = member_means.var(axis=0, ddof=0)
            output["model_variance"] = within + between
            output["within_model_variance"] = within
            output["between_member_variance"] = between
        else:
            x, names, sets, categorical = _feature_module.features(structure)
            _require(x.shape == (n, 47) and names == self.context["feature_names"],
                     "Static feature identity/order changed")
            _require(sets["sequence_complex"] == names and categorical == self.context["categorical"],
                     "Static feature set/categorical identity changed")
            # NaN phi/psi at chain ends is intentional and handled by fitted HGB.
            _require(not np.isinf(x).any(), "Infinite static features")
            means = [np.stack([models[metric].predict(x) for metric in METRICS], axis=1)
                     for models in self.members]
            member_means = np.stack(means).astype(np.float64)
            _require(member_means.shape == (len(self.members), n, 2) and np.isfinite(member_means).all(),
                     "Wrong/nonfinite static prediction")
            output["mean"] = member_means.mean(axis=0)
        output["output_mask"] = output["geometry_mask"].copy()
        _require(output["mean"].shape == (n, 2) and np.isfinite(output["mean"]).all(),
                 "Invalid deployment output")
        return output
