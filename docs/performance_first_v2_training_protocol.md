# FlexPPI performance-first v2 frozen protocol

## Scope and model inputs

The release task predicts only per-residue absolute RMSF and side-chain
conformational entropy (`chi_entropy`, H-chi). SASA, apo/holo response targets,
delta RMSF, and delta H-chi are outside this run. The v2 model is
structure-only: temperature, trajectory duration, and dataset identity remain
in provenance metadata but are not injected into the network. Consequently a
biological structure receives only one main-condition target.

## Label definition

RMSF is internal C-alpha fluctuation after fitting every protein chain
independently, preventing partner translation or rotation from being counted as
local flexibility. H-chi uses every defined side-chain chi angle, three
120-degree bins per angle, a joint base-3 rotamer code, and natural-log Shannon
entropy. Residues missing the required atoms remain masked. Independent
replicates are averaged; sample variance across available replicates or windows
is stored as label uncertainty. Losses use RMSF and H-chi masks separately.

ATLAS contributes its three-replica 100 ns RMSF labels and no H-chi. DynaRepo
contributes the existing chain-fitted, 1,000-frame-window RMSF/H-chi labels.
ITrace contributes the existing chain-fitted RMSF plus H-chi calculated over
all readable frames of each of three 200 ns replicas, followed by replicate
mean and variance. These source-specific observation spans are recorded and
must be acknowledged as a label-domain limitation; they are not silently
treated as identical physical time scales. ThermoPCD will use the same
chain-fitting and H-chi definitions, and must be recomputed before it becomes a
ready training source. A future harmonized-time-window ablation is required if
the paper claims time-scale-independent dynamics.

## ThermoPCD condition policy

There is one main record per independent PDB or reviewed biological construct.
Use an exact 310 K condition when available; otherwise use exact 311 K. Merge
compatible random seeds at that selected temperature into one mean and sample
variance. Do not merge directories with ambiguous combined-temperature
metadata. Reserve 312/313 K for temperature-shift/OOD evaluation. A candidate
must have a production XTC, readable protein PDB, consistent sequence across
replicates, at most 1,024 residues, and a successful structure/trajectory atom
mapping before label extraction.

## Five-fold evidence protocol

Every exact PDB, reviewed biological system, exact complete-complex sequence
signature, and duplicate label fingerprint is confined to one fold. Single
chain proteins are additionally grouped by MMseqs2 chain family at 30% identity
and 70% bidirectional coverage. Multi-chain systems are grouped by their full
set of chain-family clusters; sharing one common partner is insufficient to
merge complete complexes. ITrace primary evaluation holds out complete
TCR-clonotype/peptide/MHC combinations and therefore does not claim that common
MHC or antibody-like chain families are unseen. Strict family-unseen and
leave-one-source-out evaluations are separate claims.

For outer fold `k`, `fold_k_test.jsonl` is untouched test data and
`fold_k_train.jsonl` contains the other four folds. Epoch count and optimizer
settings are frozen in advance; the outer test fold is never used for checkpoint
selection. Hyperparameter changes after viewing outer-fold results require a
new split/config version.

## Weighting and metrics

The loss first averages valid residues within each metric and biological
system, then averages available RMSF/H-chi tasks within that system. Sampling
gives equal expected mass to each data source, then to each biological system
inside the source; duplicate records or longer proteins therefore receive no
automatic extra weight. Evaluation reports residue-pooled results, equal-system
macro results, source-stratified results, uncertainty calibration, and
high-flexibility retrieval.

H-chi additionally reports an amino-acid-corrected analysis. For every outer
fold, the mean H-chi of each amino-acid type is fitted on that fold's training
records only. The same type mean is subtracted from both prediction and target,
and residual metrics plus within-amino-acid-type metrics are reported. This
tests whether the model ranks residues using three-dimensional context rather
than obtaining a high score primarily from residue identity.

## Readiness gate

Formal training may begin only when the baseline Git tag exists, all ready NPZ
files pass structure-length and finite-mask checks, every fold passes PDB,
system, complex-signature, duplicate-label, and configured homology audits, the
five configurations load successfully, and a one-system forward/backward smoke
test completes without NaN or out-of-memory errors. ThermoPCD remains excluded
from ready manifests until uniform labels have actually been generated.
