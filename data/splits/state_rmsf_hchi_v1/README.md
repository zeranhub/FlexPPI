# state_rmsf_hchi_v1 fixed split

This directory is the frozen split for the absolute-state retraining run.  The
model targets are RMSF and side-chain chi entropy (`chi_entropy`) only.  SASA
labels may still exist in the source NPZ files, but the training configuration
does not select that column.  Apo/holo pair manifests and delta labels were not
read when this split was created.

ATLAS and DynaRepo state records were filtered at 1,024 residues.  Every PDB
identifier is indivisible.  All protein chains were clustered jointly with
MMseqs2 at minimum sequence identity 0.30, minimum bidirectional coverage 0.70,
sensitivity 7.5, and set-cover cluster mode 0.  Records sharing any chain
cluster were joined into one biological component.  The 244 ITrace state
records participated only as external homology guards; connected training
components were excluded rather than moved to another split.

The retained records are 1,977 train, 246 validation, and 246 untouched internal
test.  Another 391 records are explicit exclusions: 36 exceed the length limit
and 355 are homologous to the ITrace external benchmark.  The generated audit
found no biological-group, exact-PDB, or sequence-cluster overlap among the
three retained splits.  Exact counts and source-manifest SHA-256 hashes are in
`summary.json`; `master_cluster_split.tsv` records the decision for every source
and external record.

Recreate the split from the repository root with:

```powershell
python scripts\prepare_state_splits.py `
  --manifest artifacts\data\atlas_states.jsonl `
  --manifest artifacts\data\dynarepo_states.jsonl `
  --external-manifest artifacts\data\itrace_context_diagnostic\itrace_states.jsonl `
  --output-dir data\splits\state_rmsf_hchi_v1 `
  --mmseqs artifacts\tools\mmseqs\bin\mmseqs `
  --max-residues 1024 `
  --minimum-identity 0.30 `
  --minimum-coverage 0.70 `
  --train-fraction 0.80 `
  --validation-fraction 0.10 `
  --seed 7
```

On Windows the downloaded Linux MMseqs2 binary is invoked through WSL.  The
binary and its working databases live below `artifacts/` and are intentionally
not committed.
