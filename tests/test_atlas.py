from pathlib import Path

import numpy as np

from flexppi.data.atlas import prepare_atlas
from flexppi.data.dataset import StateDataset
from flexppi.data.labels import load_label_file
from flexppi.data.records import load_state_records
from test_structures_and_pairing import write_chain


def test_prepare_atlas_builds_replica_mean_and_variance(tmp_path):
    root = tmp_path / "atlas"
    entry = root / "demo_A"
    entry.mkdir(parents=True)
    write_chain(entry / "demo_A.pdb", [("A", "AGSV", 0.0)])
    (entry / "demo_A_RMSF.tsv").write_text(
        "seq\tRMSF_R1\tRMSF_R2\tRMSF_R3\n"
        "A\t1.0\t2.0\t3.0\n"
        "G\t2.0\t2.0\t2.0\n"
        "S\t3.0\t4.0\t5.0\n"
        "V\t1.0\t1.5\t2.0\n",
        encoding="utf-8",
    )
    summary = prepare_atlas(str(root), str(tmp_path / "out"))
    assert summary["indexed"] == 1
    record = load_state_records(summary["manifest"])[0]
    values, variances, mask = load_label_file(record.label_path)
    assert np.isclose(values[0, 0], 2.0)
    assert np.isclose(variances[0, 0], 1.0)
    assert mask[:, 0].all()
    assert not mask[:, 1:].any()
    combined = StateDataset([summary["manifest"], summary["manifest"]])
    assert len(combined) == 2
