# FastProtFlex local CPU benchmark

This directory contains new local inference with the unchanged official model in
`E:/Project/protein_flexibility_models/FastProtFlex`. Exact source SHA-256 hashes
and R package versions are recorded in `run_protocol.json` and each record's
`R_session_info.txt`. No training or target-based calibration is performed.

The input is the complete frozen structure (3,063 records), including all chains
in their original coordinates. The model builds an atom-contact graph, so it can
use inter-chain contacts. The official 7 Å graph threshold, 2.1 Å smoothing
threshold, 15 coefficients, normalizations, and 10,000-atom branch rule are
unchanged. The value assigned to the Cα atom is used as the residue prediction.

**Output is normalized relative flexibility, not an Å-valued RMSF prediction.**
Use only common ranking/standardized-shape metrics. There is no Hχ output.

## Files used by evaluation

- `records/{sample_id}/validation.json`: authoritative success marker; require
  `status == "complete"` and verify the saved `prediction_sha256`.
- `records/{sample_id}/prediction.npz`: `prediction` (float64), `target`, `valid`,
  `chain`, `author_residue_id`, `amino_acid`, and `residue_index`. Residue order is
  exactly the frozen complete-structure truth order. The full-precision values
  are captured before rounding to the official output PDB B-factor column.
- `atom_predictions_exact.rds` and `atom_predictions_full_precision.csv` retain
  all atom outputs. Every atom identity and coordinate, plus every Cα identity,
  is checked against the frozen input before an NPZ is marked complete.
- `full_results.csv` is updated after each record. `full_audit.json` is written
  only after all requested records have been attempted. Inspect failures/coverage
  explicitly; no missing prediction is replaced with another method.

## Scheduling and smoke measurements

The three source-median smoke structures passed all mapping/finite-value checks:
ATLAS `2o70_F` (174 residues, 1,385 atoms, 1.28 s, 276 MiB), DynaRepo `A007V`
(389 residues / 2 chains, 3,072 atoms, 3.00 s, 662 MiB), and ITrace `itrace_1qrn`
(826 residues / 5 chains, 6,620 atoms, 8.12 s, 2,945 MiB). Times include launching
R for each structure; import time in the Python wrapper is excluded. The early
three-point size extrapolation is about 2.3 hours for all 3,063 records, before
waiting for available RAM or contention from other authorized local workloads.
It is an estimate, not a measured full-run runtime.

The complete label-free preflight verified all 3,063 structure hashes and
961,321 Cα counts, found zero alternate-location atoms, and found a maximum of
8,058 protein heavy atoms. Every structure therefore uses the official
all-at-once branch. The atom-count version of the three-point estimate is
2.28 hours. `preflight_audit.json` and `preflight_sizes.csv` retain the details.

One R process tree runs at a time, pinned to one logical CPU, with numerical
library thread caps of one. The R runtime may create idle helper threads; the
CPU affinity prevents additional CPU parallelism. Before the next structure is
launched, the wrapper waits for a smoke-based quadratic memory estimate times
1.15 plus 512 MiB reserve. Waiting affects scheduling only and never changes the
official graph or chunking rules. `waiting_for_memory.json` records such waits.

After 22 initial records, a larger next structure was waiting for RAM. The
wrapper was safely stopped while no R child existed, then resumed in ascending
heavy-atom count order so small structures could continue first. Every one of
the 3,063 structures remains scheduled. `full_schedule.csv` records this order
and its SHA-256 is in `run_protocol.json`; completed predictions are verified
and reused without a second forward pass. The first log is retained as
`full_local_before_size_schedule.log`. During a resumed run, `full_results.csv`
reflects the new traversal so its early row count may be lower than the number
of existing complete `validation.json` records.

The first attempt exposed a wrapper-only pandas naming collision (`raw.insert`
was a method rather than the CSV column). It was changed to `raw["insert"]`;
all three smoke structures were then re-run successfully. Earlier exception
records are preserved as `exception_before_success.json` for audit.
