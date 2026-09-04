# state_rmsf_hchi_v2_performance_first

This directory is the versioned candidate inventory and five-fold split for
the structure-only RMSF/H-chi model. `master_manifest.tsv` contains ready,
pending, and excluded records; `summary.json` contains exact counts, source
hashes, fold balance, and leakage results. Each `fold_N_train.jsonl` contains
four folds and each `fold_N_test.jsonl` contains the held-out fold. There is no
permanent random validation split and no outer test record is used for model
selection.

ITrace training NPZ files and the ThermoPCD candidate manifest are generated
below `artifacts/data` and intentionally ignored by Git. Recreate all tracked
inventory products from the repository root with:

```powershell
python scripts\integrate_itrace_hchi_labels.py
python scripts\audit_pdbbind_ppi.py
python scripts\prepare_thermopcd_candidates.py
python scripts\prepare_performance_first_v2.py `
  --manifest artifacts\data\atlas_states.jsonl `
  --manifest artifacts\data\dynarepo_states.jsonl `
  --manifest artifacts\data\itrace_performance_first_states.jsonl `
  --manifest artifacts\data\thermopcd_performance_first_candidates.jsonl `
  --pdbbind-audit data\splits\state_rmsf_hchi_v2_performance_first\pdbbind_ppi_audit.tsv `
  --mmseqs artifacts\tools\mmseqs\bin\mmseqs
```

The exact scientific protocol and readiness boundary are frozen in
`docs/performance_first_v2_training_protocol.md`.
