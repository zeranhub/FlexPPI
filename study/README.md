# Study data and training protocol

The study evaluates FlexPPI on held-out protein structures. Each of the 3,063 records contributes one structure and its molecular-dynamics (MD) reference measurements. The cohort draws 1,929 records from [ATLAS](https://www.dsimb.inserm.fr/ATLAS/), 890 from [DynaRepo](https://dynarepo.inria.fr/) and 244 from [ITrace](https://itrace.life/docs.html).

## Cohort and cross-validation

[Fold membership](fold_membership.csv) gives the study record identifier (`sample_id`), dataset, PDB accession, target chains and source group. An empty chain selection includes all chains. The last two columns give the held-out fold (0–4) under two separate five-fold partitions. `system_fold` groups records by biological system and supplies the predictions for the internal residue-ranking comparison. `family_fold` gives the assignment for the sequence-family evaluation, which considers every constituent chain. When training fold 0, records assigned 0 supply the test set and records assigned 1–4 supply the training set.

## Training from local study inputs

Use a directory of study structures named `<sample_id>.pdb` and a directory of aligned MD labels named `<sample_id>.npz`. Each label file contains `metric_names` (including `rmsf` and `chi_entropy`), `values`, `variances` and a Boolean `mask` marking valid observations. The three arrays have one row per residue in the order parsed from the PDB file and one column per metric. The structure and label collections are separate inputs to this repository. The [input checksums](input_checksums.csv) identify the exact prepared PDB and label files used for each record by SHA-256.

Run these commands from the repository root after installing FlexPPI and a CUDA-enabled PyTorch build. Replace the two directory arguments with your input directories.

```bash
python study/prepare_folds.py --structures data/structures --labels data/labels
flexppi train --config study/configs/full_fold0.yaml
flexppi evaluate --config study/configs/full_fold0.yaml --checkpoint runs/full/fold_0/epoch_040.pt --split test --output runs/full/fold_0/test_metrics.json
```

The preparation step writes the five train/test manifest pairs under `data/splits/system/`. The [Full](configs/full_fold0.yaml) and [Intra](configs/intrachain_fold0.yaml) configurations train separate models for fold 0, using learned amino-acid embeddings, seed 7 and 40 epochs. Full processes the complete structure. Intra processes one chain at a time. Both use a sampling distribution that gives each dataset equal weight, then gives each biological system equal weight within its dataset. To train another fold, set the fold number in the two manifest paths and in `output_dir`. For the family partition, pass `--scheme family`, select manifests under `data/splits/family/` and use a distinct output directory for the new run.

## Benchmark tables

The [absolute-error table](benchmarks/absolute_errors.csv) reports the mean absolute error (MAE) in ångströms for RMSF and in nats for Hχ. For each held-out structure, errors are averaged over its valid residues. The table then gives each of the 3,063 structures equal weight. Full achieves a mean MAE of 0.577 Å for RMSF and 0.166 nats for Hχ under the system-based partition.

[Residue ranking](benchmarks/residue_ranking.csv) reports average precision (AP) and area under the receiver operating characteristic curve (AUROC) for Full, Intra, a static-feature baseline and an amino-acid-only baseline. A positive residue has an MD value at or above the 90th percentile for that record and target. The within-type Hχ analysis ranks residues of the same amino-acid type within a record, then averages over eligible types. Each record has equal weight in the final mean. The interval columns give 95% bootstrap intervals from 10,000 resamples of records stratified by dataset.

[Protein-complex RMSF](benchmarks/complex_rmsf.csv) compares seven methods on the same 260,340 residues from 615 DynaRepo complexes in 261 families. Each row gives the record-mean AP or AUROC for one method. Its 95% interval comes from 2,000 bootstrap resamples of families. This comparison uses the DynaRepo subset of the study cohort.

[Ubiquitin residue profiles](ubiquitin/profiles.csv) give the plotted MD and model values by position. The accompanying [solvent-accessible surface areas](ubiquitin/solvent_accessibility.csv) report side-chain and whole-residue exposure in square ångströms.

## Sensitivity to the MD observation window

For 83 records (30 ATLAS, 23 DynaRepo and 30 ITrace), reference values were also calculated over the first 100 ns, sampling every 200 ps. The [label changes by record](label_sensitivity/label_changes_per_record.csv) and [source means](label_sensitivity/label_change_macro.csv) compare the original and common-window references. The [control results by record](label_sensitivity/controls_per_record_metrics.csv), [control means](label_sensitivity/controls_record_macro.csv) and [prediction errors under both windows](label_sensitivity/prediction_effect_record_mean_sample_sd.csv) show how the choice of reference affects evaluation. In these tables, `uniform_rep1` uses the first trajectory and `uniform_available_replicas` uses all available trajectories for the record.
