# 1MXE / A001K data audit

All model values are from the CaM-family-held-out fold 1. The 25-residue CaMKI peptide is 294–318; the structure-first anchor span is W303–M316, 14 residues, with the remaining 11 as the comparison region. All valid sidechains and all repeated amino-acid groups are retained. Full versus Intra compares two models against the same bound-state MD target; Intra is not an experimental unbound state.

SASA compares identical model/cohort coordinates with and without CaM, using 1.4 Angstrom probe and 960 points. The difference measures static solvent shielding by CaM, not a free-peptide ensemble. Heavy-atom minimum distances and all contacting CaM residues are exported. CA8 absence must not be called no contact. CaM labels use mature native numbering (local +3), peptide labels native CaMKI numbering (local +149).

MD five 100 ns windows come from a single 500 ns trajectory. They describe temporal consistency, not independent experimental or simulation repeats. Hchi is fixed-bin joint torsional entropy, not rotation speed. Cross-amino-acid magnitudes mix residue type and environment, so same-type examples and all groups are provided.

Regenerate with the standalone audit_1mxe_v85.py script in this folder using the project numerical runtime.
