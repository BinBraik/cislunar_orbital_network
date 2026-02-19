# CR3BP Reachable Sets (rs3 / rs4)

A MATLAB codebase for generating and analyzing reachable-set atlases in the Circular Restricted Three-Body Problem (CR3BP), with overlap analysis utilities for transfer discovery.

## What this repository contains

- `src/`: reusable core functions for configuration, propagation, voxel binning, caching, overlap analysis, and plotting.
- `scripts/`: top-level runners for typical workflows.
- `tests/`: MATLAB unit/smoke tests.
- `docs/`: architecture and logging notes.
- `example_outputs/`: example CSV artifacts from completed runs.

## Quick start

### 1) Requirements

- MATLAB (recommended R2021b+)
- Toolboxes as required by your local MATLAB setup (parallel features are optional and can be disabled in config)

### 2) Setup paths

From repository root:

```matlab
rs3_setup
```

### 3) Run a single-family atlas build and plots

```matlab
run('scripts/run_rs3_one_family_atlas_and_plots.m')
```

### 4) Run overlap workflow

```matlab
run('scripts/run_rs4_overlap_and_visuals.m')
```

### 5) Run tests

```matlab
run('tests/runtests_rs3.m')
```

## Configuration entry points

- Start from `rs3_cfg_defaults`.
- Validate with `rs3_cfg_validate` before long runs.
- The grid builder/validator pair is `rs3_grid_make` + `rs3_grid_validate`.
- Family atlas caching helpers are `rs3_prepare_or_load_family`, `rs3_cache_save_family`, and `rs3_cache_try_load_family`.

## Typical output locations

Runner scripts write results under the configured output root and tag (see `cfg.io.out_root`, `cfg.io.tag`).

## Repository hygiene

This repository now includes:

- A focused `.gitignore` for MATLAB/editor/cache artifacts
- This root README for onboarding and execution flow
- A `CONTRIBUTING.md` file with coding and testing guidance

## Notes

- `docs/architecture.md` summarizes the step-based pipeline.
- `docs/logging.md` documents optional logger usage if logger files are present in your working tree.
