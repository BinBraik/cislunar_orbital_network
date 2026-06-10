# Cislunar Orbital Network

MATLAB codebase accompanying the paper:

> **Orbital Networks in the Three-Body Problem**  
> Abdullah Braik and Shane D. Ross  
> Aerospace and Ocean Engineering, Virginia Tech  
> *Preprint, 2026*

This repository implements a reachable-set-based framework for characterizing pairwise accessibility across thirteen representative cislunar periodic orbit families at a common Jacobi constant in the planar Earth–Moon CR3BP. Forward and backward reachable-set atlases are built for each family, pairwise accessibility is inferred from reachable-set overlap, and the resulting cost matrix is assembled into a weighted cislunar orbital network whose hub, gateway, relay, and bottleneck roles are analyzed over varying maneuver and time budgets.

## Repository layout

```
src/              core functions (atlas construction, CR3BP dynamics, overlap, network)
  network/        network analysis functions (centrality, Floyd-Warshall, articulation)
scripts/          runner scripts for all workflows
tests/            MATLAB unit and smoke tests
docs/             architecture and data-structure reference
example_outputs/  example CSV artifacts from completed runs
web/              interactive D3 network visualization template
```

## Requirements

- **MATLAB R2022b or later** (uses `parenAssign`, `struct` patterns from modern releases)
- **Parallel Computing Toolbox** — optional; set `cfg.par.enable = false` to run serially

No additional toolboxes or external packages are required.

## Quick start

```matlab
% 1. From the repository root, add all paths
setup

% 2. Build the atlas for one family (Lyapunov L1 by default)
run('scripts/run_atlas_one_family.m')

% 3. Compute overlap and visualizations for one pair
run('scripts/run_overlap_pair_visuals.m')

% 4. Build the full 13×13 DV proxy matrix
run('scripts/run_overlap_all_pairs.m')

% 5. Run the test suite
run('tests/runtests.m')
```

## Execution order for full paper reproduction

All atlases must exist in the cache before running overlap or network scripts.
Run `run_atlas_one_family.m` once per family (13 families listed in Table 2 of the paper),
or adapt it to loop over `cfg.families.list`.

```
STEP 1 — Build family atlases                          [paper §3]
  run_atlas_one_family.m  (×13, one per family)
  └─ outputs: atlas_cache/<family>.mat (cached automatically)

STEP 2 — Baseline DV proxy matrix                     [paper §4]
  run_overlap_all_pairs.m
  └─ outputs: minDVproxy_matrix.csv, pair_winners_top1.csv

STEP 3 — Parametric sweep over (DV_cap, Tmax) budget plane  [paper §6]
  run_overlap_dv_tmax_sweep.m
  └─ outputs: atlas_sweep_results/sweep_DVmatrix_results.mat
  ┆            (takes several hours; checkpointed automatically)

STEP 4 — Network analysis and centrality maps          [paper §5–6]
  run_network_centrality_sweep.m     (requires STEP 3)
  └─ outputs: atlas_network_results/  (CSVs + PDF figures)

STEP 5 — Differential-correction trajectory examples   [paper §7]
  run_overlap_dc_sweep.m             (requires STEP 1)
  └─ run_overlap_dc_sweep_plot.m     (post-process: requires dc_sweep output)

STEP 6 — Relay / betweenness examples (optional)      [paper §5.2, §7.3]
  run_betweenness_explainer.m        (requires STEP 3)
  └─ run_betweenness_dc.m
       └─ run_betweenness_video.m

STEP 7 — Numerical robustness verification (optional)  [paper Appendix A]
  run_verification_grid_sweep.m      (×7 grid configurations)
  └─ run_robustness_analysis.m

ANYTIME — utility runners
  run_atlas_cache_manager.m          inspect / prune the atlas cache
  run_overlap_pair_visuals.m         deep-dive overlap plots for one pair
  run_overlap_voxel_traj.m           extract and correct a single voxel transfer
```

## Representative orbit families (paper Table 2, CJ = 3.1294)

| Abbreviation | Family | Period [days] |
|---|---|---|
| LL1 | L1 Lyapunov | 12.811 |
| LL2 | L2 Lyapunov | 15.117 |
| C11a | (1,1)a-cycler | 42.140 |
| C11b | (1,1)b-cycler | 55.995 |
| C21 | (2,1)-cycler | 84.533 |
| C32 | (3,2)-cycler | 78.613 |
| R21-S | 2:1 stable resonant | 26.500 |
| R21-U | 2:1 unstable resonant | 31.039 |
| R31-S | 3:1 stable resonant | 27.252 |
| R31-U | 3:1 unstable resonant | 28.066 |
| R52-S | 5:2 stable resonant | 54.802 |
| R52-U | 5:2 unstable resonant | 56.436 |
| DPO | Distant prograde orbit | 11.184 |

Family initial conditions are stored in `src/cr3bp_family_ic.m` (paper Table 2).

## Configuration

All parameters are centralized in `atlas_cfg_defaults`. Override individual fields in
runner scripts before calling `atlas_cfg_validate`. Key parameters:

| Field | Default | Meaning |
|---|---|---|
| `cfg.grid.dx` / `cfg.grid.dy` | 0.01 nd | Spatial voxel width |
| `cfg.grid.dtheta` | 2° | Heading voxel width |
| `cfg.fan.DV_cap_nd` | 0.2 nd | Max heading-offset ΔV budget |
| `cfg.propag.Tmax` | π nd | Max integration time per direction |
| `cfg.par.enable` | true | Enable PARFOR over seed×heading jobs |
| `cfg.cache.rebuild` | false | Force atlas recomputation |

## Outputs

- **Atlases** — cached in `atlas_cache/` (keyed by a deterministic parameter fingerprint)
- **Results** — written to `atlas_results/<timestamp>/` as CSV files and PNG figures
- **Sweep data** — written to `atlas_sweep_results/` (DV/TOF matrices over budget plane)
- **Network results** — written to `atlas_network_results/` (centrality CSVs and figures)

## Citing this work

See `CITATION.cff` for the BibTeX entry. GitHub renders a "Cite this repository" button
from that file automatically.

## Further reading

- `docs/architecture.md` — detailed pipeline and data-structure reference
- `docs/logging.md` — optional logger integration
- `CHANGELOG.md` — physics bug-fix history
