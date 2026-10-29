"""Small deterministic regression checks for the MD bookkeeping and analysis.

These tests do not simulate or fabricate validation measurements.
"""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

BASE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("kix_analysis", BASE / "analyze_job.py")
analysis = importlib.util.module_from_spec(spec)
spec.loader.exec_module(analysis)
pilot_spec = importlib.util.spec_from_file_location("kix_pilot_comparison", BASE / "compare_pilot.py")
pilot = importlib.util.module_from_spec(pilot_spec)
pilot_spec.loader.exec_module(pilot)


class KixWorkflowTests(unittest.TestCase):
    def test_frozen_chi_definition(self):
        self.assertEqual(analysis.chi_hash(), analysis.EXPECTED_CHI_HASH)

    def test_reject_segment_gap(self):
        with tempfile.TemporaryDirectory(prefix="kix-test-") as folder:
            path = Path(folder)
            (path / "segment_0000.dcd").touch()
            (path / "segment_0002.dcd").touch()
            with self.assertRaises(ValueError):
                analysis.segment_paths(path)

    def test_numeric_order_not_lexical(self):
        with tempfile.TemporaryDirectory(prefix="kix-test-") as folder:
            path = Path(folder)
            for number in range(12):
                (path / f"segment_{number}.dcd").touch()
            self.assertEqual([p.name for p in analysis.segment_paths(path)],
                             [f"segment_{number}.dcd" for number in range(12)])

    def test_input_construct_mapping(self):
        inputs = BASE.parents[1] / "artifacts/kix_md_20260914/inputs"
        for name, count in [("apo", 87), ("myb", 112), ("mll", 106), ("both", 131)]:
            residues = {(line[21], int(line[22:26]))
                        for line in (inputs / f"{name}.pdb").read_text().splitlines()
                        if line.startswith("ATOM")}
            self.assertEqual(len(residues), count)
            if name in ("mll", "both"):
                self.assertEqual(sorted(n for c, n in residues if c == "C"), list(range(839, 858)))

    def test_twenty_unique_velocity_seeds(self):
        seeds = [20260914 + repeat * 100 + index for repeat in range(1, 6) for index in range(4)]
        self.assertEqual(len(set(seeds)), 20)

    def test_pilot_policy_only_rep01(self):
        policy = json.loads((BASE.parents[1] / "artifacts/kix_md_20260914/PILOT_POLICY.json").read_text())
        self.assertEqual(policy["allowed_repeats"], [1])
        self.assertEqual(policy["states"], ["apo", "myb", "mll", "both"])
        self.assertTrue(policy["additional_repeats_require_user_approval"])

    def test_descriptive_comparison_no_pvalue(self):
        result = pilot.description([1., 2., 3.], [2., 3., 4.])
        self.assertAlmostEqual(result["mae"], 1.)
        self.assertAlmostEqual(result["spearman_r"], 1.)
        self.assertNotIn("pvalue", result)
        flat = pilot.description([0., 0., 0.], [1., 2., 3.])
        self.assertIsNone(flat["spearman_r"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
