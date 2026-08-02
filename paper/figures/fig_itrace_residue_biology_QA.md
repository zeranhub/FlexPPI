# ITrace residue-biology figure QA

## Figure contract

- Core conclusion: the frozen model partially recovers the held-out MD peptide position pattern, but misses the P6 mobility peak and remains only moderately accurate within individual peptides.
- Evidence chain: panel a directly compares complete-complex predictions with ITrace MD truth; panel b is an unvalidated cohort context response plus two same-peptide cases; panel c asks what residue-type hypotheses the model generates beyond length and position.
- Export: Nature-family double-column width (183 mm), vector PDF/SVG and 300 dpi PNG.
- Main review boundary: panels b/c are fixed-coordinate model context responses, not observed apo-to-holo changes; ITrace has no MD $H_\chi$/SASA truth.

## Data validation

```json
{
  "md_9mer_positions": 9,
  "fidelity_9mer_positions": 9,
  "context_9mer_positions": 9,
  "nine_mer_complex_count": 154,
  "external_clean_complex_count": 244,
  "peptide_residue_count": 2278,
  "amino_acid_type_count": 20,
  "case_complex_count": 2,
  "case_residue_count": 18,
  "case_sequence_count": 2
}
```

No rows were downsampled. Panel a/b use all 154 external-clean 9-mer complexes; panel c uses all 2,278 peptide residues from all 244 external-clean complexes. Ala and Gly remain visible in panel c but are marked NA for $H_\chi$.

## Academic Figure Skill QA

Pass 0 — anti-patterns: AP-0 to AP-7 PASS. Verbatim typography, palette and export blocks are present; no default/rainbow palette, four-sided border, internal data-obscuring legend or raster-only export is used.

Pass 1 — code compliance: CL-1 to CL-7 PASS. Minimum declared font is 6 pt; canvas is 183 x 190 mm; PDF fonttype 42, editable SVG text and 300 dpi PNG are emitted; spines are 0.6 pt and ticks point outward.

Pass 2 — visual logic: VI-1 to VI-7 PASS. Direct model-versus-MD prediction fidelity is the first panel, colors keep the same metric meaning, zero baselines and CI definitions are explicit, and panel labels use one placement function.

Pass 3 — visual verification: VV-1 to VV-5 PASS after three render/inspection cycles. No legend, label or confidence interval obscures data; panel edges and margins align; all labels remain legible at the 183 mm canvas; black/blue/purple/orange/grey encodings remain distinct; and every panel contains non-flat finite signal. Quantitative checks pass because all nine position profiles vary and all three paired contact-effect bootstrap CIs exclude zero.
