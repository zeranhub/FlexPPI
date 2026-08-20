# Local baseline deployments

This directory records reproducible local deployments added for the frozen
FlexPPI RMSF comparisons. Upstream repositories and model artifacts under
`tmp/third_party/` are intentionally ignored by Git; the pinned revisions,
checksums, adapters, and evaluation outputs remain recorded here or under
`artifacts/evaluation/additional_baselines/`.

## RMSF-net PDB-only

- Upstream: `XintSong/RMSF-net` at commit
  `eb353b0d9d909a6c46d2808467a9d0f16e2bbda7`.
- Local image: `flexppi-baseline/rmsf-net:pdb-eb353b0`.
- Image ID: `sha256:e25d665b8e6c2f81ba4886db347a939db782a1d11fa631656e54216f1d5fb4f7`.
- Adapter: `scripts/run_rmsf_net_baseline_batch.py`.
- Checkpoint: upstream `model_cpu/sim_model.pth` (the official only-PDB CPU
  branch). No cryo-EM map intensities enter the comparison.

The adapter reproduces the documented UCSF Chimera `molmap` settings in NumPy:
4 Angstrom resolution, 1.5 Angstrom spacing, atomic-number Gaussian amplitudes,
and the default 12 Angstrom edge padding. It then runs the official network and
maps atom-voxel predictions to author-numbered residues. The released program
normalizes predicted values before writing the PDB B-factor field, so RMSF-net
is treated as a relative flexibility score; only rank, standardized-error, and
high-flexibility detection metrics are assigned to it.

Build from the repository root:

```powershell
docker build -f baselines/rmsf_net/Dockerfile -t flexppi-baseline/rmsf-net:pdb-eb353b0 .
```

## AlphaFlow MD+Templates

- Upstream: `bjing2016/alphaflow` at commit
  `0408d7c89dac444a43a9089d7427ce470b0a5e67`.
- OpenFold revision: `103d037`.
- Selected checkpoint: `alphaflow_12l_md_templates_distilled_202406.pt`.
- Checkpoint SHA-256:
  `E481DA4432CC0118E3A1CD51FFBDFB0ABECD50CE6642D9CD96FA5515451867AA`.
- Official ensemble archive SHA-256:
  `73859D4277D5CA65BDF16910BCE807232ACE44722A437EA29CDE11D93CBFBAE1`.
- Ensemble evaluator: `scripts/evaluate_alphaflow_ensembles.py`.

The official AlphaFlow entry point hard-codes CUDA transfers and its OpenFold
stack requires CUDA 11. The current workstation has no NVIDIA adapter exposed
to Docker, so source, weights, and a pinned GPU image recipe are installed
locally, but fresh official inference cannot run on this host. The completed
comparison therefore uses the authors' released 250-frame trajectories for the
six proteins overlapping the frozen ATLAS-191 cohort. This is reported as a
small ensemble-generation case study, never as a full-cohort or ITrace result.

Build the GPU image on a CUDA-capable Docker host from the repository root:

```powershell
docker build -f baselines/alphaflow/Dockerfile -t flexppi-baseline/alphaflow:0408d7c .
```
