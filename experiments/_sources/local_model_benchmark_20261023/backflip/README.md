# Fresh local BackFlip comparisons, 2026-10-23

The frozen main cohort contains all 3,063 corrected records (ATLAS 1,929,
DynaRepo 890, ITrace 244), 4,973 chains and 961,321 valid RMSF reference
positions. Existing labels and structures are read only. No MD is rerun and no
model is retrained, calibrated or selected by its errors.

`prepare_inputs.py` writes complete-backbone canonical PDB inputs and an explicit
original-chain / original-author-number / insertion-code identity map. Every
N/CA/C and optional O coordinate is checked against the existing corrected
structure, with exact float32 equality after parsing. Chain names become A-Z;
integer numbering prevents the official loader from deleting insertion-code
positions. Incomplete N/CA/C geometry remains represented in the original
record mapping but is excluded from prediction and evaluation by the frozen
geometry mask. No chain is removed for being short.

`run_benchmark.py` uses the official frame conversion and
`BackFlip.predict_from_frames()` with released weights, CPU, batch size 1. The
whole-complex context includes all valid input chains in their original relative
geometry. The independent-chain context makes separate calls and maps their
outputs back to the original complete record. Single-chain records share the
same fresh forward pass between these equivalent contexts; the status explicitly
records this within-run reuse. Historical predictions are never substituted.

The current BackFlip-2.1 source is the immutable official commit
`887c6a73b216e3fe2ff9f551ae6e2c983c21b292`; its downloaded checkpoint SHA-256 is
`1af6e13a0e0161b2b92fa6fb9cddd2dee8f3e71c284f5059c77d7a30b9c527b9`.
Per-residue predicted covariance has nm² units; RMSF is `10*sqrt(trace(C))` in Å.
BackFlip-1.0 and 1.0-seq use the existing read-only official commit
`60b8c83511b965925bbe990bb2be578bf3529460` and local released weights.
Their required GAFL and GATr sources are fixed to
`018247ade22812407b6dfb9093944f442a606aa0` and
`87fa7d0c5cec5d874b13814816c628843e348866`, with the two patches from GAFL's
official `install_gatr.sh` applied only inside this task directory.

The primary target remains the manuscript's existing per-chain Kabsch RMSF:
each chain's rigid translation/rotation is removed independently. This is not
identical to aligning a multichain trajectory as one rigid object. BackFlip-2
was trained on globally aligned monomer data; fresh complete-complex inference
here is an extension tested against the fixed manuscript target, not evidence
that the upstream authors validated TCR-pMHC complexes. Official default residue
position encodings are concatenated continuous indices across chains; chain
identity is not added to the model.

The shared native Python 3.12 runtime supplies torch 2.12 CPU. `torch_scatter`
is only imported by unused MD-label/training helpers; because its Windows CPU
binary is absent, the adapter provides import guards that raise immediately if
those helpers are invoked. No scatter operation, network layer, weight or model
forward computation is replaced. The complete native smoke forward succeeded,
and the new independent-chain 1AO7 prediction agrees with the frozen earlier
official computation to max 4.30e-6 Å and mean 4.78e-7 Å over 708 residues.

`smoke_results/` is the original four-thread BackFlip-2.1 smoke and any separately
labelled original-model smoke runs. The formal new queue uses exactly two CPU
threads. `results/<tag>/status.jsonl` records every success/failure, exact coverage,
timing and original residue identities; each NPZ includes `prediction`, `target`,
`valid`, `sequence`, `chain` and `author_residue_id`. Timings include canonical
PDB parsing, official frame conversion, forward and unit conversion, and exclude
dependency import, checkpoint load and input preparation. Shared single-chain
contexts must not be counted as an extra inference. `summarize_results.py`
reports available coverage explicitly and never treats an unfinished run as
complete. Read `run_metadata.json` and `status.jsonl` for actual execution state.

The confirmed mutually exclusive reporting groups are single-chain proteins
(2,193), other multichain systems (626), and ITrace TCR-pMHC complexes (244).
The 870 multichain records are additionally summarized together. Group names
describe the evaluation cohort and do not add inputs to any released model.

`dependency_patch_audit.json` freezes the two official GAFL/GATr patches, both
BackFlip-1 checkpoint hashes, and the exact opt_einsum source required by GATr.
This old opt_einsum commit needs one compatibility edit with NumPy 1.26:
testing whether a shape is empty uses `len(shape) == 0`, because comparing a
nonempty NumPy shape array with `()` now raises. Its contraction equations,
optimization, contraction order, and all model calculations remain unchanged.
The original and edited parser hashes are recorded. BackFlip-1's released
`global_rmsf` output is already in Angstrom, whereas the BackFlip-2 covariance
output uses the conversion stated above. The `-seq` model adds one-hot residue
identity features; it does not load a protein language model.

`three_group_smoke_manifest.json` fixes the input checks to 16pk_A (415 residues,
one chain), A000A (310 residues, two chains), and 1AO7 (708 residues, five chains).
`validate_outputs.py` checks every completed NPZ's hash, reference mask and target,
original residue order/sequence/chain/author identifier, coverage, and equality
of the reused single-chain contexts. Its `completion_audit.json` is authoritative
only when `complete` is true. `run_remaining_variants.ps1` is a bounded local
queue for BackFlip-1 and BackFlip-1-seq after BackFlip-2 exits successfully.
