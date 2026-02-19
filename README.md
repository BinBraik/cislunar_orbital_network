# CR3BP Reachable Sets (rs3 / rs4)

A MATLAB codebase for computing and analyzing reachable-set atlases in the planar Earth-Moon Circular Restricted Three-Body Problem (CR3BP), with overlap analysis for identifying low-energy transfer candidates between periodic orbit families.

## What this does

Given two periodic orbit families (e.g., Lyapunov L1 and a Resonant orbit), the pipeline:

1. **Voxelizes** the (x, y, θ) state space into a 3D grid.
2. **Seeds** trajectories from points along each periodic orbit with a fan of heading angles within a DV budget.
3. **Integrates** forward and backward from each seed in the reduced CR3BP model, logging which voxels each trajectory passes through (the Forward/Backward Reachable Sets).
4. **Caches** each family atlas to disk keyed on all pipeline parameters.
5. **Overlaps** family A's FRS with family B's BRS to find voxels reachable from both, which are candidates for a 3-impulse transfer.

## Repository structure

- `src/` — core functions for configuration, propagation, voxel binning, caching, overlap, and plotting
- `scripts/` — top-level runner scripts for typical workflows
- `tests/` — MATLAB unit and smoke tests
- `docs/` — architecture and logging notes
- `example_outputs/` — example CSV artifacts from completed runs

## Supported periodic orbit families

| Family | Notes |
|---|---|
| Lyapunov L1 | L1 Lyapunov family |
| Lyapunov L2 | L2 Lyapunov family |
| Cycler 21 | 2:1 resonant cycler |
| Cycler 11a / 11b | 1:1 resonant cyclers |
| Cycler 32 | 3:2 resonant cycler |
| Resonant 2to1 Stable / Unstable | 2:1 resonance |
| Resonant 3to1 Stable / Unstable | 3:1 resonance |
| Resonant 5to2 Stable / Unstable | 5:2 resonance |
| Distant Prograde Orbit | DPO |

## Quick start

### 1) Requirements

- MATLAB R2021b or later
- Parallel Computing Toolbox (optional; disable with `cfg.par.enable = false`)

### 2) Add paths

From the repository root:

```matlab
rs3_setup
```

### 3) Build a single-family atlas

```matlab
run('scripts/run_rs3_one_family_atlas_and_plots.m')
```

Edit the script to change the family name and grid/propagation parameters.

### 4) Compute overlap between two families

```matlab
run('scripts/run_rs4_overlap_and_visuals.m')
```

Edit `famA` / `famB` at the top of the script to choose families.

### 5) Run all-pairs summary

```matlab
run('scripts/run_rs4_all_pairs_summary.m')
```

### 6) Run tests

```matlab
run('tests/runtests_rs3.m')
```

## Configuration

All knobs live in `rs3_cfg_defaults`. Override fields in runner scripts before calling `rs3_cfg_validate`. Key parameters:

| Field | Default | Meaning |
|---|---|---|
| `cfg.grid.dx/dy` | 0.01 nd | Spatial voxel width |
| `cfg.grid.dtheta` | 4° | Heading voxel width |
| `cfg.fan.DV_cap_nd` | 0.1 nd | Max heading offset DV budget |
| `cfg.propag.Tmax` | π/2 nd | Max integration time (each direction) |
| `cfg.par.enable` | true | Enable PARFOR over seed×heading jobs |
| `cfg.cache.rebuild` | false | Force atlas recomputation |

## Outputs

Results are written under `cfg.io.out_root` (default: `rs3_results/`) in timestamped subdirectories. Each run produces PNG figures. Cached atlases are stored in `rs3_cache/`.

## Further reading

- `docs/architecture.md` — detailed pipeline and data structure reference
- `docs/logging.md` — optional logger integration
