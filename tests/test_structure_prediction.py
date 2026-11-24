from pathlib import Path

import pytest

from flexppi.inference import predict_structure


def test_released_ensemble_reproduces_reference_prediction():
    pdb = Path(__file__).resolve().parents[1] / "examples" / "2qho_ab.pdb"
    rows = predict_structure(str(pdb), target_chains=["A"], device_name="cpu")
    assert {row["chain"] for row in rows} == {"A"}
    first = rows[0]
    assert (first["author_residue_id"], first["amino_acid"]) == ("1", "M")
    # Reference means from the five original system-split training checkpoints.
    assert first["rmsf_angstrom"] == pytest.approx(0.55992086, abs=2e-5)
    assert first["chi_entropy_nats"] == pytest.approx(2.01496234, abs=2e-5)
    assert all((row["chi_entropy_nats"] is None) == (row["amino_acid"] in "AG")
               for row in rows)
