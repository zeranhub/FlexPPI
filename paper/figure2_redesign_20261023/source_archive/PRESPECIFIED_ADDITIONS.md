# Figure 2 density revision: fixed statistical additions

The v44 evaluation protocol, reference q90 thresholds with complete ties, original
out-of-fold predictions, 3,063 record weighting, original AP/AUROC values and
record-level source-stratified resampling are unchanged. No model is retrained,
no prediction is calibrated, no observation is trimmed, and no threshold is tuned.

1. Prediction/reference conditional trends use the exact 20 reference-bin edges
   from v44, with every valid original position (961,321 RMSF and 824,988 Hchi).
   Each position receives weight inverse to its record's valid target count;
   hence each record has equal total mass. All four methods use identical target
   masks and bins. Full prediction means and 95% intervals must reproduce v44.
   Intervals use the first 1,000 of the same 10,000 source-stratified record
   bootstrap draws (seed 20261018); the full 10,000 draw construction is retained
   to reproduce the exact generator stream for later source strata.
2. Within-amino-acid Hchi PR curves use each eligible record/type's own MD q90
   threshold (at least four positions and both reference classes). Scores are
   processed in whole tie groups. Precision at recall grid 0, .01, ..., 1 is
   sampled with v44's first threshold reaching that recall. Types are averaged
   equally within record, then the 3,063 records are averaged equally. AP and
   AUROC remain averages of the original complete threshold curves, rather than
   areas computed from the sampled macro curve. The true all-selected endpoint
   at recall 1 is appended separately because recall can first reach 1 before
   all positions have been selected. Eligibility and all AP/AUROC values must
   reproduce v44's 42,108 eligible record/types. Complete type threshold arrays
   are archived for every method.
3. The source AP matrix contains the same per-record AP values, averaged within
   ATLAS, DynaRepo or ITrace, for RMSF, Hchi and within-type Hchi. It contains
   3 tasks x 3 sources x 4 methods = 36 cells. Averaging source values with their
   record counts must reproduce the original overall values. No source receives
   an extra weight in the overall score.

These are complementary displays of continuous prediction agreement and fixed
high-dynamic-site retrieval. Source contrasts are descriptive; no significance
claim is added. The original v44 plots, category results and data remain archived.
