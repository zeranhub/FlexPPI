# Case-study figures QA report

## Figure contract

- 6UON core message: central peptide contacts explain a localized static partner-context response, while the ITrace MD target identifies p6/p8 as the flexible ridge.
- 2PTT core message: stable electrostatic anchors coexist with a conformationally ordered yet still compliant CC-prime loop.
- Archetype: asymmetric mixed-modality; structure views are the hero evidence and residue plots are the quantitative support.
- Export: 183 mm double-column figures, vector PDF/SVG masters and 300 dpi PNG previews.

## Data integrity

- 6UON: all 826 residues are retained in the source-data table and whole-structure color mapping; the local profile contains all 10 peptide residues.
- 2PTT: all 206 mapped DynaRepo residues are retained; native PDB numbering is restored by documented chain-specific offsets.
- 2B4 H42 is present in the coordinates but absent from the mapped DynaRepo table; it is masked and never interpreted as zero mobility.
- The 6UON full-minus-pMHC difference holds target-chain coordinates fixed and is labeled as a static context response, not an apo/holo trajectory.

## Pass 0 - anti-pattern scan

- AP-0 to AP-7: PASS. Mandatory typography, restrained custom palette, font embedding, vector export, clean spines and explicit font family are present.
- Molecular panels have white backgrounds and semantic chain/residue colors; no rainbow/jet maps are used.

## Pass 1 - code compliance

- CL-1 to CL-7: PASS. Text is at least 5.5 pt at final size; width is 183 mm; PDF/SVG and 300 dpi PNG outputs are generated; PDF text uses TrueType embedding.
- Structure panels are rendered at 2x device scale before composition; quantitative layers remain vector.

## Pass 2 - visual logic and reproducibility

- VI-1 to VI-7: PASS. Each panel contributes unique evidence; the same absolute RMSF scale is used for MD and FlexPPI; zero is explicit on the context-delta panel; selected residues are defined by published contacts/loops rather than post hoc score extremes.
- No inferential p-value is reported from residues within one complex because residues are not independent biological replicates.

## Pass 3 - rendered verification

- VV-1: PASS. Panel titles, legends, annotations and molecular labels do not overlap or clip after the final title/legend adjustments.
- VV-2: PASS. Axes, panel gutters and structure crops remain aligned at the 183 mm final figure width.
- VV-3: PASS. Residue identities, RMSF scales and interpretive callouts remain legible in both standalone figures and the compiled manuscript.
- VV-4: PASS. Structure colors, residue highlights and quantitative traces preserve the same semantic mapping across panels.
- VV-5: PASS. PDF/SVG masters and 300 dpi PNG previews open correctly; all 30 rendered manuscript pages were visually inspected.

## Traceability

- 6UON source data: `E:\Project\FlexPPI\paper\figures\fig_6uon_case_study_source_data.csv`.
- 2PTT source data: `E:\Project\FlexPPI\paper\figures\fig_2ptt_case_study_source_data.csv`.
- Full numerical and alignment record: `E:\Project\FlexPPI\paper\figures\fig_case_studies_statistics.json`.
