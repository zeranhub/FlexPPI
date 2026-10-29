"""Official BackFlip2 chain calls reassembled into the frozen complete record.

This adapter never reads labels, changes upstream code, fits, or calibrates.
Input/model/source hashes and chain identities are checked outside the timed
call by the separate collector. Every timed call parses every chain again.
"""
from pathlib import Path
import sys
import time

import numpy as np
import torch


def import_dependencies(source_root):
    sys.path.insert(0, str(Path(source_root)))
    from backflip.deployment.inference_class import BackFlip
    from backflip.data.flexibility_utils import batched_rmsf_from_covar
    return BackFlip, batched_rmsf_from_covar


class BackFlipAdapter:
    def __init__(self, context):
        BackFlip, self.converter = context["dependencies"]
        self.context = context
        self.model = BackFlip(Path(context["checkpoint_path"]), device="cpu",
                              progress_bar=False, confidence_intervals=False)
        self.metadata = {
            "parameter_count": sum(p.numel() for p in self.model.model.parameters()),
            "parameter_bytes": sum(p.numel() * p.element_size() for p in self.model.model.parameters()),
            "device": "cpu", "cpu_threads": torch.get_num_threads(),
            "checkpoint_tag": "backflip-2.1",
            "source_commit": "887c6a73b216e3fe2ff9f551ae6e2c983c21b292",
            "output_metrics": ["rmsf"],
            "context": "Independent frozen isolated chains, reassembled once per record",
        }
        if self.metadata["cpu_threads"] != 2:
            raise ValueError("BackFlip2 must use exactly two CPU threads")
        if self.metadata["parameter_count"] != 965180:
            raise ValueError("Official BackFlip2 2.1 parameter count changed")
        self.last_chain_timings = []

    def predict_record(self, record):
        frozen = self.context["records"][record["sample_id"]]
        n = len(frozen["amino_acid"])
        assembled = np.full(n, np.nan, dtype=np.float64)
        visits = np.zeros(n, dtype=np.int32)
        chain_times = []
        # These are identity maps, never cached features or predictions.
        for chain in frozen["chains"]:
            ids = np.asarray(chain["record_indices"], dtype=np.int64)
            started = time.perf_counter()
            with torch.inference_mode():
                result = self.model.predict_from_pdb(chain["structure_path"], batch_size=1)
                covariance = result["per_res_covariance"]
                rmsf = np.asarray(self.converter(covariance), dtype=np.float64).reshape(-1)
                covariance = covariance.detach().cpu().numpy()
            seconds = time.perf_counter() - started
            if rmsf.shape != (len(ids),) or covariance.shape != (len(ids), 3, 3):
                raise ValueError("Official chain output lost residue coverage")
            if not np.isfinite(rmsf).all() or not np.isfinite(covariance).all():
                raise ValueError("Nonfinite official prediction")
            trace = np.trace(covariance, axis1=-2, axis2=-1)
            if (trace < 0).any() or not np.allclose(rmsf, 10 * np.sqrt(trace),
                                                 rtol=2e-6, atol=2e-6):
                raise ValueError("Covariance nm2 to RMSF Angstrom conversion changed")
            assembled[ids] = rmsf
            visits[ids] += 1
            chain_times.append({"input_id": chain["input_id"], "chain": chain["chain"],
                                "residue_count": len(ids), "parse_forward_convert_seconds": seconds})
        if not (visits == 1).all() or not np.isfinite(assembled).all():
            raise ValueError("Every complete-record residue must be covered exactly once")
        self.last_chain_timings = chain_times
        return {
            "mean": assembled[:, None],
            "amino_acid": np.asarray(frozen["amino_acid"]),
            "chain": np.asarray(frozen["chain"]),
            "author_residue_id": np.asarray(frozen["author_residue_id"]),
            "geometry_mask": np.asarray(frozen["geometry_mask"], dtype=bool),
            "target_mask": np.asarray(frozen["target_mask"], dtype=bool),
            "output_mask": np.ones(n, dtype=bool),
        }
