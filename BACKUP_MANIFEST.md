# FlexPPI GitHub backup manifest

This repository is a source-level backup intended to restore the FlexPPI code, configuration, tests, documentation, and paper on another computer.

## Included in Git

- Python package source under `src/`
- training, evaluation, and utility scripts under `scripts/`
- configuration files under `configs/`
- tests and documentation
- manuscript source, tables, figures, figure-generation scripts, and the rendered `paper/paper.pdf`
- project metadata such as `pyproject.toml`, `README.md`, and `NOTICE`

## Intentionally excluded

The following local directories are generated, machine-specific, or too large for ordinary Git history:

- `artifacts/` — datasets, checkpoints, and evaluation products
- `output/` — generated benchmark and experiment outputs
- `.venv-md/` and other `.venv-*` directories — local Python environments
- `tmp/` and `paper/tmp/` — temporary files and rendered-page checks
- `ideaspark_run/` — generated literature-search and idea-evaluation caches
- LaTeX auxiliary files such as `*.aux`, `*.log`, and `*.out`

At the time this backup was prepared, the excluded local content was approximately:

- `artifacts/`: 1.36 GB
- `output/`: 285 MB
- `.venv-md/`: 105 MB
- `ideaspark_run/`: 9 MB
- `tmp/`: 3 MB

These files remain on the original computer and are not deleted by the Git backup.

## External data required after restore

Several commands refer to datasets outside this repository. Their paths must be supplied or remapped on the destination computer:

- ATLAS molecular-dynamics data
- DynaRepo molecular-dynamics data
- TCR3D structures
- ITrace benchmark data
- trained FlexPPI checkpoints when evaluation or continued training is required

See `README.md`, `configs/evaluate_state_local.yaml`, and `configs/evaluate_pair_local.yaml` for the currently configured local paths.

## Basic restore check

After cloning, create a new Python environment and install the project rather than copying `.venv-md/`:

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -e ".[dev]"
pytest
```

The manuscript can be rebuilt from `paper/paper.tex` with XeLaTeX when a compatible TeX distribution is installed.
