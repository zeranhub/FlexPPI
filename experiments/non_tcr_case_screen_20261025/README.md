# Non-TCR case screen (2026-10-25)

This folder is a stable copy of the exploratory screen produced in
`C:\Users\DELL\.codex\worktrees\2a52\FlexPPI\output\non_tcr_case_screen_20261025`.
It covers 616 non-ITrace protein-complex records from the frozen cohort.

Files:

- `all_616_record_screen.csv` — per-record overall Full/Intra metrics, interface gains, density correlations.
- `all_chain_metrics.csv` — per-chain metrics (1,535 rows).
- `all_region_metrics.csv` — interface / outside-interface / remote metrics (1,848 rows).
- `candidate_table_enriched.csv` — record table merged with region gains (used for shortlisting).
- `pair_delta_metrics.csv`, `pair_delta_audit.json` — paired isolated/complex state-difference metrics.
- `independent_numeric_review.md/.json` — independent recomputation of the 616 tables; no errors found.
- `screen_cases.py`, `screen_pairs.py`, `independent_numeric_review.py` — original scripts.

Selection was exploratory: records were ranked by balanced MAE gain and interface-specific
gain, then inspected for biological/structural clarity. The screen is not an unbiased
performance estimate, and the paired Full/Intra models are separately trained
architectures, not a causal deletion of the partner.
