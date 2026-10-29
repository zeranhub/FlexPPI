# PROTAC input and output fit pre-screen, 2026-10-26

This is a structural and source-code pre-screen. No PROTAC model predictions or MD were run. Contact counts are not prediction accuracy, binding strength, or causal attribution.

## Decision

Keep 1MXE as the currently validated manuscript candidate. PROTAC examples should pass input/output suitability and matched-reference checks before narrative strength is considered. The present results do not establish that PROTAC prediction fails; they establish reasons not to prioritize these unvalidated examples over an existing validated protein-interface case.

## Exact model implementation

The source audited is the frozen family-held-out source used for the current CaM analysis, rather than the older worktree source:

`E:\Project\FlexPPI\paper\revision_v14_20261010\section33_expansion_20261025\biology_story_first_20261025\current_model_cam_predictions_20261025\family_heldout_model\src`

- `flexppi/data/structures.py:83–85` filters out residues outside the amino-acid dictionary. N/CA/C/O alone are stored at residue level.
- `flexppi/model/geometry.py` builds N/CA/C local frames and CA-distance/orientation features. PROTAC atoms are not inputs.
- `flexppi/model/information_ablation.py:178–217` defines the separately trained intrachain model, encoding each chain independently. Full/Intra tests access to other protein chains in the same structural state; it is neither a ligand-deletion experiment nor an observed apo/holo transition.
- `flexppi/data/labels.py:77–85` fits each protein chain independently before RMSF calculation. Relative rigid-body motion between proteins is removed. Protein-internal RMSF is not whole-ternary orientational flexibility.
- Hchi describes protein side-chain torsional-state occupancy. It does not measure linker dynamics, ternary residence time, or degradation rate.

The current family-fold configuration disables condition features and the response head. No paired-response output was used in this audit.

## Downloaded experimental coordinates

The RCSB author-assigned biological assemblies were used, excluding crystal neighbors. `download_manifest.json` records URLs and SHA256 checksums. Protein-chain identities:

- 1MXE assembly 1: CaM A, CaMKI peptide E.
- 5T35 assembly 1: BRD4 A, VHL D; assembly 2: BRD4 E, VHL H. Elongin B/C remain in the downloaded complete assemblies. This contact comparison specifically counts target BRD4 versus VHL, not the internal VCB support interfaces.
- 6W8I assemblies 1–3: BTK/cIAP1 A/D, B/E, C/F.
- 6W7O assemblies 1–2: BTK/cIAP1 A/C, B/D.

6W8I is 3.80 A resolution; 6W7O is 2.17 A; 5T35 is 2.70 A. Missing modeled residues and crystallographic disorder require attention in any subsequent matched-region comparison. All retained residues in these selected assemblies have complete N/CA/C frames; this does not imply that all residues in the experimental constructs are resolved.

## Measured structural contacts

A residue is counted when any of its retained heavy atoms is within 4.5 A of a heavy atom in the specified partner protein or PROTAC. Alternate locations are handled as in the model parser; the first model of each biological assembly is used. Counts are static geometric descriptions, not measures of how strongly an atom constrains motion. Total target lengths differ, so raw counts should not be interpreted as comparative binding strength.

| Structure | Assembly | Target modeled residues | Target residues contacting partner protein | Target residues contacting PROTAC | PROTAC-contact residues without partner-protein contact |
|---|---:|---:|---:|---:|---:|
| 1MXE, CaMKI peptide | 1 | 25 | 23 | n/a | n/a |
| 5T35, BRD4 | 1 | 111 | 6 | 16 | 11 |
| 5T35, BRD4 | 2 | 111 | 5 | 16 | 12 |
| 6W8I, BTK | 1 | 254 | 2 | 23 | 22 |
| 6W8I, BTK | 2 | 251 | 6 | 24 | 24 |
| 6W8I, BTK | 3 | 263 | 9 | 27 | 26 |
| 6W7O, BTK | 1 | 261 | 16 | 22 | 21 |
| 6W7O, BTK | 2 | 263 | 13 | 21 | 21 |

The qualitative disparity in 6W8I remains at 4.0 and 5.0 A: target–protein contact residues are 2–10, while target–PROTAC contacts are 20–27 across the three conformations and tested cutoffs. The complete three-cutoff results and residue distances are in the CSV/JSON outputs.

For all eight assemblies, the frozen parser produces identical amino-acid, residue-index, chain-index, backbone, and backbone-mask arrays before and after removing all non-amino-acid records. Hash equality verifies input invariance to that removal. No neural-network inference was needed or performed to establish this input identity. Under identical inference settings, a deterministic model receiving these same arrays cannot distinguish ligand presence from absence.

## Interpretation and subsequent eligibility

1MXE has an extensive protein-peptide interface that is explicitly represented by the current protein inputs, and already has model/reference analysis. CaM calcium ions are also absent from this parser: its advantage is empirical validation and a represented peptide interface, not an assertion that every physical determinant is encoded.

5T35 remains a possible protein-interface candidate, particularly for BRD4–VHL contact residues away from ligand-dominated pocket contacts. However, it has no local PROTAC prediction/reference validation in this screen, and matched experimental BRD3 ternary coordinates are not established by the previously cited BRD3 comparison.

6W8I/6W7O provide real ternary structural alternatives, but ligand-contacting regions contain many contacts absent from the input, and published whole-complex flexibility is not the same quantity as per-chain-aligned RMSF. Their present role is exploratory evaluation, not a preferred manuscript replacement.

RIPK1–VHL public trajectories may permit a useful technical test, but the four discussed states are modeled metastable conformations, and the broad trajectory collection includes linker-free encounters. It requires a separate dataset/label audit before inference; neither the public data listing nor cell degradation results are model validation.

For any subsequent case, require protein geometry that represents the proposed local constraint, reference RMSF/Hchi computed with the existing definitions, successful model prediction on matched regions, and robustness across available conformations/replicates. Only then assess whether independent functional experiments strengthen its biological interpretation.

## Reproduction and sources

Run `audit_contacts.py` with the numerical Python runtime recorded in the current project. Outputs are `structural_input_audit.csv`, `structural_input_audit.json`, and `residue_contact_details.csv`; downloaded sources and generated protein-only files are retained.

- https://www.rcsb.org/structure/1MXE
- https://www.rcsb.org/structure/5T35
- https://www.rcsb.org/structure/6W8I
- https://www.rcsb.org/structure/6W7O
- https://www.nature.com/articles/s41589-020-00686-2
- https://pmc.ncbi.nlm.nih.gov/articles/PMC6423499/
- https://zenodo.org/records/18155129
