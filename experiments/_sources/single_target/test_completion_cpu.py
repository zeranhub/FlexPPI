"""Small CPU checks for the two consequential new conventions."""
from pathlib import Path
import tempfile

import numpy as np
import pandas as pd

from completion_common import read_cache, fold_summary
from analyze_completion_oof import component_side


def main():
    with tempfile.TemporaryDirectory(prefix="flexppi_completion_test_") as folder:
        root = Path(folder)
        path = root / "chi_entropy_only"
        path.mkdir()
        target = np.array([[10., 1., np.nan], [20., 2., np.nan], [30., 3., np.nan]])
        valid = np.array([[True, True, False], [True, False, False], [True, True, False]])
        np.savez_compressed(path / "synthetic.npz", mean=np.array([[1.1], [2.2], [3.3]]), model_variance=np.ones((3, 1)),
                            target=target, label_variance=np.zeros((3, 3)), valid=valid,
                            metric_names=np.array(["chi_entropy"]), label_names=np.array(["rmsf", "chi_entropy", "sasa"]))
        cache = read_cache(root, dict(variant="chi_entropy_only", sample_id="synthetic"), "chi_entropy")
        assert np.array_equal(cache["target"], [1., 2., 3.])
        assert np.array_equal(cache["valid"], [True, False, True])
        assert np.allclose(cache["prediction"], [1.1, 2.2, 3.3])
    data = pd.DataFrame([dict(variant="full", metric="rmsf", dataset="atlas", fold=0, mae=1., count=1000),
                         dict(variant="full", metric="rmsf", dataset="atlas", fold=0, mae=3., count=1)])
    folds, _ = fold_summary(data, ["mae"], ["variant", "metric"])
    assert np.array_equal(folds.mae, [2., 2.]), "Records must have equal weight despite different residue counts"
    assert component_side("component_same", 1729) == component_side("component_same", 1729)
    print("PASS: Hchi output0 maps to label1/mask1; equal-record macro; deterministic component split")


if __name__ == "__main__":
    main()
