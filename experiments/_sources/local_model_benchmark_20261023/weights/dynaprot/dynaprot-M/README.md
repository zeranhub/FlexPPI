# dynaprot-M

DynaProt-M (marginal 3x3 covariances), AlphaFlow/ATLAS split

## Original files (still present in the source repo; nothing was moved)

- checkpoint: `ckpt/marginal_logfrob_lripa_afsplit_NOsinusoidal_randomcropping_1layerreadout-step=054500.ckpt`
- model config: `configs/model/dynaprot_marginal_nosinusoidal_1layer.yaml`
- data config: `configs/data/atlas_config_randomcropping.yaml`
- data config (alternate, eval-equivalent): `configs/data/atlas_config_randomcropping_256res.yaml`

## Referenced by

- experiments/eval_afsplit.ipynb cell 5 (with data cfg atlas_config_randomcropping.yaml)
- experiments/eval_pairwise.ipynb cell 4 (with data cfg atlas_config_randomcropping_256res.yaml)
- experiments/pocketfinding.ipynb cell 5

## Produces

- Table: anisotropic blob prediction (RMWD Var / Sym. KL Var), ATLAS test split
- Table: RMSF Pearson r (also reported on FlexPert split by dynaprot-M-flexpert)
- Figure: marginal_gaussians A/B (rendered Gaussians, joint distribution vs AFMD+T/NMA)
- Figure: rmsf_grid (appendix RMSF correlations)
- Figure: crypticpocket (1ADE apo marginals), cryptic_pocket_alignment (BioEmu cases)
- Ensemble generation (with dynaprot-J): Table dynamics_comparison, Figure ensembles_magma, BPTI appendix

## Loading

```python
import yaml, torch
from dynaprot.model.architecture import DynaProt

model_cfg = yaml.safe_load(open('model_config.yaml'))
model_cfg['data_config'] = yaml.safe_load(open('data_config.yaml'))
model = DynaProt.load_from_checkpoint('model.ckpt', cfg=model_cfg).eval()
```

`save_hyperparameters()` was never called during training, so the config above is
not recoverable from the checkpoint itself -- it must be passed explicitly.
