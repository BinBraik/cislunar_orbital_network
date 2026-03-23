# CR3BP Reachable Sets

MATLAB pipeline for computing reachable-set atlases in the planar Earth-Moon Circular Restricted Three-Body Problem (CR3BP), finding low-energy transfer candidates between periodic orbit families, and ranking families by network centrality.

---

## What this does

Given a set of periodic orbit families, the pipeline:

1. **Voxelises** the (x, y, θ) state space into a 3D grid.
2. **Seeds** trajectories from points along each periodic orbit with a fan of heading angles within a ΔV budget.
3. **Integrates** forward and backward from each seed, logging which voxels each trajectory passes through — the Forward Reachable Set (FRS) and Backward Reachable Set (BRS).
4. **Caches** each family atlas to disk, keyed on all pipeline parameters via an MD5 fingerprint.
5. **Overlaps** family A's FRS with family B's BRS to find shared voxels — candidates for a 3-impulse transfer.
6. **Corrects** the best candidate arcs with differential correction to obtain a continuous, fuel-optimal transfer.
7. **Analyses** the network of all family pairs: Floyd-Warshall shortest paths, LCC, strength, harmonic closeness, betweenness, and PageRank.

---

## Repository layout

```
cr3bp_reachable_sets/
├── rs3_setup.m           — add repo paths (call once per session)
├── scripts/              — 15 runner scripts (entry points, see table below)
├── src/                  — ~85 library functions (rs3_*, rs4_*, rs5_*, net_*)
│   └── network/          — 9 graph-analysis functions (net_*)
├── tests/                — unit and smoke tests; run via runtests_rs3.m
├── docs/                 — architecture and logging reference
├── example_outputs/      — example CSV artifacts from completed runs
├── web/                  — interactive DV map template (dv_map_template.html)
├── rs3_cache/            — generated atlas cache (gitignored)
├── rs3_results/          — timestamped run output (gitignored)
├── rs3_sweep_results/    — parametric sweep matrices (gitignored)
├── rs3_network_results/  — network analysis output (gitignored)
└── rs4_verification_runs/ — robustness verification runs (gitignored)
```

---

## Pipeline at a glance

| Stage | Script | What it does |
|---|---|---|
| 1 – Setup | `rs3_setup` | Add `src/` and `scripts/` to the MATLAB path |
| 2–5 – Single family | `run_rs3_one_family_atlas_and_plots` | Build / load one family atlas; plot FRS, BRS, seeds |
| 5 – Overlap | `run_rs4_overlap_and_visuals` | Overlap two atlases; generate detailed figures |
| 5 – Batch overlap | `run_rs4_all_pairs_summary` | All 78 pairs at one configuration; DV matrix + CSV |
| 5 – DC single | `run_rs4_voxel_trajectories` | Re-integrate best voxel; differential correction |
| 5 – DC batch | `run_rs4_dc_sweep` | DC for all 78 pairs at fixed parameters |
| 5 – Sweep | `run_rs4_dv_tmax_sweep` | DV-proxy over (DV_cap, Tmax) grid |
| 5 – DC sweep | `run_rs4_dc_sweep_tmax` | DC sweep over (DV_cap, Tmax) grid |
| 5 – Full traj | `run_rs5_single_family_sweep` | Complete 3-impulse trajectory from one origin |
| 6 – Network | `run_network_centrality_sweep` | Centrality metrics over the sweep grid |
| 6 – Web export | `export_sweep_to_json` | Self-contained interactive HTML map |
| 6 – Betweenness | `run_betweenness_explainer` | Demonstrate bridge routing savings |
| Util – Cache | `run_rs3_cache_manager` | Inspect, review, and derive subset caches |
| Util – Plot | `run_rs4_dc_sweep_plot` | Replay figures from a completed DC sweep |
| Util – Verify | `run_verification_grid_sweep` | Robustness check over 7 grid configurations |
| Util – Robust | `run_robustness_analysis` | Cross-run rank stability analysis |

---

## Quick start

### 1. Requirements

- MATLAB R2021b or later
- Parallel Computing Toolbox (optional; set `cfg.par.enable = false` to disable)

### 2. Add paths

```matlab
rs3_setup
```

### 3. Build a single-family atlas

```matlab
run('scripts/run_rs3_one_family_atlas_and_plots.m')
```

Edit `familyName` and the grid/fan/propag knobs at the top of the script.
Figures and a `.mat` file are written to `rs3_results/<timestamp>/`.

### 4. Overlap two families

```matlab
run('scripts/run_rs4_overlap_and_visuals.m')
```

Edit `famA` and `famB`. Both atlases must be cached first.

### 5. Run the all-pairs DV-proxy sweep

```matlab
run('scripts/run_rs4_dv_tmax_sweep.m')
```

Sweeps a 20×20 (DV_cap, Tmax) grid. Output: `rs3_sweep_results/sweep_DVmatrix_results.mat`.

### 6. Run network centrality analysis

```matlab
run('scripts/run_network_centrality_sweep.m')
```

Reads the sweep .mat and produces centrality contour maps and winner maps.

### 7. Run tests

```matlab
run('tests/runtests_rs3.m')
```

---

## Scripts reference

| Script | What it does | Needs | Produces |
|---|---|---|---|
| `run_rs3_one_family_atlas_and_plots` | Build / load atlas for one family; plot FRS, BRS, seeds, grid | Nothing (builds cache) | `rs3_cache/*.mat`, figures in `rs3_results/` |
| `run_rs4_overlap_and_visuals` | Overlap two atlases; full figure set | Cached atlases for famA and famB | Overlap .mat files + figures |
| `run_rs4_all_pairs_summary` | DV-proxy for all 78 pairs at one configuration | All 13 cached atlases | `minDVproxy_matrix.csv`, `pair_winners_top1.csv` |
| `run_rs4_voxel_trajectories` | Re-integrate best voxel + differential correction | Cached atlases for famA and famB | Trajectory figures + DC comparison figure |
| `run_rs4_dc_sweep` | DC for all 78 pairs at fixed parameters | All 13 cached atlases | `rs4_dc_sweep_results.mat` |
| `run_rs4_dc_sweep_plot` | Replay figures from a completed DC sweep | `rs4_dc_sweep_results.mat` | Before/after DC figures, `transfer_summary.csv` |
| `run_rs4_dv_tmax_sweep` | DV-proxy over (DV_cap, Tmax) grid | All 13 cached atlases | `sweep_DVmatrix_results.mat` |
| `run_rs4_dc_sweep_tmax` | DC sweep over (DV_cap, Tmax) grid | All 13 cached atlases | `sweep_dc_results.mat` |
| `run_rs5_single_family_sweep` | Full 3-impulse trajectory for one origin family | Cached atlases | Trajectory figures + `result.mat` per pair |
| `run_network_centrality_sweep` | Centrality metrics over the sweep grid | `sweep_DVmatrix_results.mat` | `network_results.mat`, contour maps |
| `export_sweep_to_json` | Interactive HTML DV map | `sweep_DVmatrix_results.mat` | `web/dv_map.html` |
| `run_betweenness_explainer` | Bridge-routing savings demo | Sweep + network .mat + atlases | Trajectory figures, GIFs, `*_data.mat` |
| `run_rs3_cache_manager` | Inspect and derive subset caches | Existing cache | Derived cache files, config printout |
| `run_verification_grid_sweep` | Robustness check over 7 grid configs | Atlases at each resolution | `rs4_verification_runs/*/result.mat` |
| `run_robustness_analysis` | Cross-run rank stability analysis | Verification run results | Rank heatmaps, edge-persistence figures |

---

## Function library

| Prefix | Role | Key functions |
|---|---|---|
| `rs3_core_*` | CR3BP dynamics and propagation | `rs3_core_cr3bp_U_and_derivs`, `rs3_core_reduced_cr3bp_model`, `rs3_core_family_ic`, `rs3_core_integrate_reduced_with_fallback` |
| `rs3_cfg_*` | Configuration | `rs3_cfg_defaults` (full schema), `rs3_cfg_validate` |
| `rs3_grid_*` | Voxel grid construction | `rs3_grid_make`, `rs3_grid_validate`, `rs3_bin_xyth`, `rs3_keep_mask_xy` |
| `rs3_family_*` | Periodic orbit seed preparation | `rs3_family_prepare_seeds`, `rs3_family_build_hits`, `rs3_seed_upper_segments_by_arclength` |
| `rs3_hits_*` | Voxel hit logging | `rs3_hits_log_from_traj`, `rs3_ev_stop_reduced`, `rs3_ev_stop_full4d` |
| `rs3_rows_*` | Packed row struct operations | `rs3_rows_empty`, `rs3_rows_vcat`, `rs3_rows_subset`, `rs3_rows_mirror_lower`, `rs3_rows_to_matrix` |
| `rs3_cache_*` | Cache read/write/inspect | `rs3_cache_save_family`, `rs3_cache_try_load_family`, `rs3_cache_fingerprint_family`, `rs3_cache_inspect` |
| `rs3_*` (util) | Utilities | `rs3_prepare_or_load_family`, `rs3_atlas_derive_subset`, `rs3_fan_delta_lists`, `rs3_io_save_figure`, `rs3_md5`, `rs3_progress_tick` |
| `rs4_overlap_*` | FRS/BRS overlap computation | `rs4_overlap_pair`, `rs4_overlap_extract_voxel_info`, `rs4_overlap_visualize`, `rs4_overlap_visualize_bounds` |
| `rs4_diffcorr*` | Differential correction | `rs4_diffcorr` (main solver), `rs4_diffcorr_visualize` |
| `rs4_voxel_*` | True-DV arc extraction | `rs4_voxel_traj_extract`, `rs4_voxel_traj_visualize_single` |
| `rs5_*` | Full trajectory assembly | `rs5_build_full_traj`, `rs5_visualize_full_traj` |
| `net_*` | Graph / network analysis | `net_build_graph`, `net_floyd_warshall`, `net_lcc`, `net_centrality`, `net_articulation`, `net_plot_baseline`, `net_plot_winner_map` |

---

## Supported periodic orbit families

| Family | Short tag |
|---|---|
| Lyapunov L1 | LyapL1 |
| Lyapunov L2 | LyapL2 |
| Cycler 2:1 | Cyc21 |
| Cycler 1:1a | Cyc11a |
| Cycler 1:1b | Cyc11b |
| Cycler 3:2 | Cyc32 |
| Resonant 2:1 Stable | R21S |
| Resonant 2:1 Unstable | R21U |
| Resonant 3:1 Stable | R31S |
| Resonant 3:1 Unstable | R31U |
| Resonant 5:2 Stable | R52S |
| Resonant 5:2 Unstable | R52U |
| Distant Prograde Orbit | DPO |

---

## Configuration reference

All knobs live in `rs3_cfg_defaults`. Override fields in runner scripts before calling `rs3_cfg_validate`.

### System / units (`cfg.sys`, `cfg.units`)

| Field | Default | Meaning |
|---|---|---|
| `cfg.sys.RE_nd` | 6378/384400 | Non-dimensional Earth radius |
| `cfg.sys.RM_nd` | 1737/384400 | Non-dimensional Moon radius |
| `cfg.units.LU_m` | 384400e3 | Earth-Moon distance (m) |
| `cfg.units.VU_mps` | ~1023 m/s | Non-dimensional velocity unit |
| `cfg.units.TU_days` | ~4.343 days | Non-dimensional time unit |

### Seed / PO sampling (`cfg.seed`)

| Field | Default | Meaning |
|---|---|---|
| `cfg.seed.ds_seed` | 0.02 nd | Arclength spacing between seed points along the PO |
| `cfg.seed.N_dense` | 2001 | Points for dense PO sampling |
| `cfg.seed.y_eps` | 0 | y-threshold for upper-half seed selection |
| `cfg.seed.minSegPts` | 5 | Minimum points per upper-half arc segment |

### Grid / voxelisation (`cfg.grid`)

| Field | Default | Meaning |
|---|---|---|
| `cfg.grid.dx` / `dy` | 0.01 nd | Spatial voxel width; finer → more memory and runtime |
| `cfg.grid.dtheta` | 2° | Heading-angle voxel width |
| `cfg.grid.Rdom` | 1.2 nd | Domain radius for grid construction |
| `cfg.grid.enforce_y0_edge` | true | Force y = 0 to lie on a grid edge |

### Steering fan (`cfg.fan`)

| Field | Default | Meaning |
|---|---|---|
| `cfg.fan.DV_cap_nd` | 0.2 nd | Max heading-offset ΔV budget at each seed |
| `cfg.fan.dtheta_fan` | 1° | Angular resolution of the heading fan |

### Propagation (`cfg.propag`)

| Field | Default | Meaning |
|---|---|---|
| `cfg.propag.Tmax` | π nd | Max integration time per arc (forward and backward) |
| `cfg.propag.absTol` | 1e-9 | ODE absolute tolerance |
| `cfg.propag.relTol` | 1e-9 | ODE relative tolerance |
| `cfg.propag.v2tol` | 1e-8 | If v² < v2tol at a seed, fall back to full 4D model |

### Logging / segment-walk (`cfg.log`)

| Field | Default | Meaning |
|---|---|---|
| `cfg.log.step_len_factor` | 0.5 | ODE step_len = factor × min(dx, dy) |
| `cfg.log.maxstep_factor` | 0.5 | ODE MaxStep = factor × step_len |
| `cfg.log.segwalk.enable` | true | Enable segment-walk subsampling |
| `cfg.log.segwalk.frac` | 0.25 | Substep ≤ frac × voxel width |

### Caching (`cfg.cache`)

| Field | Default | Meaning |
|---|---|---|
| `cfg.cache.enable` | true | Cache atlases to disk |
| `cfg.cache.dir` | `rs3_cache/` | Cache directory |
| `cfg.cache.rebuild` | false | Recompute even if cache exists |
| `cfg.cache.version_tag` | `rs3_v2_keep_masked` | Bump to invalidate all existing caches |

### Differential correction (`cfg.diffcorr`)

| Field | Default | Meaning |
|---|---|---|
| `cfg.diffcorr.tol_patch` | 1e-4 | Normalised convergence tolerance for fmincon |
| `cfg.diffcorr.tol_converged` | 1e-4 | Threshold for the CONVERGED status label |
| `cfg.diffcorr.display` | `'iter'` | fmincon verbosity: `'off'` / `'iter'` / `'final'` |
| `cfg.diffcorr.MaxIterations` | 300 | fmincon iteration budget |
| `cfg.diffcorr.MaxFunEvals` | 8000 | fmincon function-evaluation budget |
| `cfg.diffcorr.N_po_dt` | 0.003 nd | Target time between PO spline knots |

### Parallelism (`cfg.par`)

| Field | Default | Meaning |
|---|---|---|
| `cfg.par.enable` | true | Enable parfor over (seed × heading) jobs |
| `cfg.par.progress_every` | 50 | Print progress every N jobs in parfor |

### Output (`cfg.io`)

| Field | Default | Meaning |
|---|---|---|
| `cfg.io.out_root` | `rs3_results/` | Root output directory |
| `cfg.io.tag` | timestamp | Subdirectory name for each run |
| `cfg.io.save_figs` | true | Save PNG figures |
| `cfg.io.save_fig` | true | Also save `.fig` for interactive use |
| `cfg.io.fig_visible` | `'off'` | `'on'` to see figures during run |
| `cfg.io.fig_resolution` | 220 | PNG resolution (DPI) |

### Diagnostics (`cfg.diag`)

| Field | Default | Meaning |
|---|---|---|
| `cfg.diag.zoom.enable` | true | Apply zoom window to trajectory figures |
| `cfg.diag.zoom.xlim` | [0.70 1.25] | Zoom x-limits (nd), centred near Moon |
| `cfg.diag.zoom.ylim` | [−0.30 0.30] | Zoom y-limits (nd) |

---

## Outputs guide

| Directory | Contents |
|---|---|
| `rs3_cache/` | MD5-keyed family atlas `.mat` files (~5–500 MB each); regenerated on demand |
| `rs3_results/<timestamp>/` | Figures and `.mat` files from single-family and overlap runs |
| `rs3_sweep_results/` | `sweep_DVmatrix_results.mat` — DV/TOF matrices for every (DV_cap, Tmax) cell |
| `rs3_dc_sweep_results/` | `sweep_dc_results.mat` — post-DC DV/TOF matrices |
| `rs3_network_results/` | `network_results.mat` — centrality metrics; contour and winner map figures |
| `rs4_verification_runs/` | Per-run `result.mat` files + robustness analysis figures |
| `rs3_betweenness_explainer/` | Per-example static figures, GIFs, and `*_data.mat` files |
| `web/dv_map.html` | Self-contained interactive DV map (generated by `export_sweep_to_json`) |
| `example_outputs/` | Example CSV artifacts checked into the repository |

---

## Running the tests

```matlab
run('tests/runtests_rs3.m')
```

The test suite covers: packed row struct operations (`test_rows`), CR3BP physics (`test_physics`), grid construction (`test_grid`), voxel hit logging (`test_hits_log`), utility functions (`test_utilities`), and integration smoke tests (`test_smoke`).

---

## Further reading

- `docs/architecture.md` — detailed pipeline and data structure reference
- `docs/logging.md` — optional logger integration
- `CONTRIBUTING.md` — development workflow, naming conventions, testing requirements
