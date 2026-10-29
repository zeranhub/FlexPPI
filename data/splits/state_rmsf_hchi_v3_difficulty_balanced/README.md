# state_rmsf_hchi_v3_difficulty_balanced

This directory contains 3,068 ready ATLAS, DynaRepo, and ITrace systems in a
five-fold, difficulty-balanced, group-safe split.  Existing v2 biological and
homology `cluster_id` values are indivisible.  Fold optimisation uses dataset,
source group, target-chain count, label availability, protein length, and
ground-truth RMSF distribution summaries.  It never uses FlexPPI or baseline
prediction scores.

Training consumes `fold_N_train.jsonl`; evaluation consumes
`fold_N_test.jsonl`.  `fold_assignments.tsv`, `balance_audit.tsv`,
`summary.json`, and `leakage_audit.json` are audit artifacts and should be kept
with the split.  Paths inside all JSONL files are normalised to the remote
`/ranze/...` layout.

Because this split was designed after the v2 result was inspected, report it as
a new, frozen v3 protocol rather than silently replacing the original v2
evaluation.
