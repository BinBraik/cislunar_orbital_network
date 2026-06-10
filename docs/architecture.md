# Pipeline architecture

## Overview

The pipeline has four stages that map directly to sections of the paper:

| Stage | Paper section | Key functions |
|---|---|---|
| Atlas construction | §3 | `atlas_grid_make`, `atlas_family_prepare_seeds`, `atlas_family_build_hits`, `atlas_cache_save` |
| Overlap & accessibility | §4 | `overlap_pair`, `overlap_extract_voxel_info`, `overlap_visualize_*` |
| Network analysis | §5–6 | `net_build_graph`, `net_floyd_warshall`, `net_centrality`, `net_articulation` |
| Trajectory realization | §7, App. C | `overlap_voxel_traj_extract`, `traj_diffcorr`, `traj_build_full` |

---

## Stage 1 — Atlas construction (§3)

**Step 1 — Grid** (`atlas_grid_make`, `atlas_grid_validate`)

Build the 3D (x, y, θ) voxel grid. The grid is symmetric about both axes and
enforces y = 0 on a bin edge so that upper-half and lower-half rows can be
mirrored exactly. Grid parameters live in `cfg.grid`.

**Step 2 — Seeds** (`atlas_family_prepare_seeds`)

Sample a periodic orbit densely along its arc length. Extract points on the
upper half (y ≥ 0) as seeds. Apply the admissible-domain Keep mask
(`atlas_keep_mask_xy`) to exclude seeds near Earth/Moon or outside the domain.

**Step 3 — Integration** (`atlas_family_build_hits`, `atlas_hits_log_from_traj`)

For each (seed × heading) job, apply a constant-CJ heading-change maneuver
(`atlas_fan_delta_lists`) and integrate forward in the reduced (x, y, θ) model
(`cr3bp_reduced_ode`). If the reduced model becomes singular near the zero-
velocity curve, fall back to the full 4D model (`cr3bp_integrate`). Log every
voxel hit into a packed row struct using the segment-walk algorithm
(`atlas_hits_log_from_traj`). PARFOR over jobs if `cfg.par.enable = true`.

**Step 4 — Cache** (`atlas_cache_save`, `atlas_cache_try_load`)

Cache each family atlas to a `.mat` file keyed on an MD5 fingerprint of all
pipeline parameters (`atlas_cache_fingerprint`). A cache hit skips Steps 2–3
entirely. Inspect the cache with `atlas_cache_inspect`.

---

## Stage 2 — Overlap & accessibility (§4)

**Overlap computation** (`overlap_pair`)

Compute FRS(A) ∩ BRS(B) on the shared voxel grid. The backward reachable set
of family B is derived from its forward rows via the y-axis time-reversal
symmetry (§2.3 of the paper): `R(x, y, θ) = (x, −y, π−θ)`. Lower-half rows
for both families are reconstructed on-the-fly from upper-half storage via
`atlas_rows_mirror_lower`. A pre-built θ-mirror look-up table makes the
symmetry mapping O(1) per row.

**Voxel metadata** (`overlap_extract_voxel_info`)

For each overlap voxel, extract the minimum-DV seed/heading pair from each
family, compute the heading-turn cost (ΔV_turn), the patch cost (ΔV_patch),
and the proxy total (§4.3 of the paper).

**Visualization** (`overlap_visualize`, `overlap_visualize_combo`, `overlap_visualize_bounds`)

Render FRS/BRS occupancy scatter plots, overlap voxels, and DV-proxy bound
annotations. Optional plot types are controlled by `cfg.plot.overlap.*` flags
(all enabled by default).

---

## Stage 3 — Network analysis (§5–6)

**Graph construction** (`net_build_graph`)

Assemble the 13×13 DV proxy matrix into a weighted undirected graph. Apply
feasibility threshold (DV_cap, Tmax) to determine which edges exist.

**Shortest paths** (`net_floyd_warshall`)

All-pairs shortest weighted paths over the family network. Used for harmonic
closeness centrality and multileg routing.

**Centrality** (`net_centrality`, `net_articulation`)

Harmonic closeness centrality (cost-aware, §5.2) and articulation-point
detection (gateway roles, §5.2). Short family names for plotting come from
`net_family_short_names`.

---

## Stage 4 — Trajectory realization (§7, Appendix C)

**Voxel extraction** (`overlap_voxel_traj_extract`)

From a chosen overlap voxel, extract the full arc pair: the A-side arc
(forward from family A seed) and B-side arc (backward from family B seed,
via symmetry). Reconstruct the 3-impulse proxy trajectory.

**Differential correction** (`traj_diffcorr`)

Refine the proxy into a continuous trajectory using `fmincon` with position-
continuity constraints. Formulation in Appendix C of the paper. Convergence
tolerance controlled by `cfg.diffcorr.tol_converged`.

**Visualization** (`traj_diffcorr_visualize`, `traj_build_full`, `traj_visualize_full`)

Plot before/after DC comparisons and render the full multi-leg trajectory.

---

## Key data structures

### `cfg` — configuration struct
All parameters in one place. See `atlas_cfg_defaults` for the full schema
with inline comments. Key sub-structs:

| Sub-struct | Controls |
|---|---|
| `cfg.grid` | Voxel widths, domain radius |
| `cfg.fan` | DV budget, heading resolution |
| `cfg.propag` | Integration time, ODE tolerances |
| `cfg.cache` | Cache directory, version tag, rebuild flag |
| `cfg.overlap` | Primary buffer fraction, parallel extraction flag |
| `cfg.diffcorr` | DC tolerances, solver settings |
| `cfg.par` | PARFOR enable, progress interval |
| `cfg.io` | Output root, figure settings |

### `grid3` — voxel grid
Bin edges, bin centers, Keep mask (logical Ny×Nx), and the `bin_xyth`
function handle for fast (x, y, θ) → (ix, iy, it) lookup.

### `S` — family atlas struct
Output of Steps 2–3, stored in the cache:

| Field | Description |
|---|---|
| `S.name` | Family name string |
| `S.mu`, `S.CJ`, `S.Tf_PO` | Dynamics parameters |
| `S.SeedsUpper`, `S.SeedsLower` | Seed positions (Nx3) |
| `S.Step4.rows_FRS_upper` | Packed hit rows for forward reachable set |
| `S.grid3` | Grid with per-family Keep mask applied |

### `O` — overlap struct
Output of `overlap_pair`:

| Field | Description |
|---|---|
| `O.ids` | Linear voxel IDs of overlapping voxels |
| `O.x`, `O.y`, `O.th` | Voxel centre coordinates |
| `O.countA`, `O.countB` | Hit counts per family per voxel |

---

## Packed row struct schema

Hit rows are stored as a struct-of-arrays rather than a double matrix
(~4× memory reduction at scale).

| Field | Type | Description |
|---|---|---|
| `n` | uint32 | Number of valid rows |
| `iSeed` | uint16 | Seed index |
| `iHead` | uint16 | Heading index within fan |
| `leg` | uint8 | 1 = forward, 2 = backward |
| `halfFlag` | int8 | +1 = upper half, −1 = lower half |
| `t` | single | Time of voxel hit |
| `ix`, `iy`, `it` | uint16 | Voxel bin indices |

Helper functions: `atlas_rows_empty`, `atlas_rows_vcat`, `atlas_rows_count`,
`atlas_rows_subset`, `atlas_rows_mirror_lower`, `atlas_rows_to_matrix`.

---

## Symmetry convention

The CR3BP has y = 0 symmetry (§2.3): a trajectory `(x, y, θ, t)` in the upper
half corresponds to `(x, −y, π−θ, −t)` in the lower half. The pipeline
exploits this to:

- Integrate only from upper-half seeds (Stage 1)
- Store only upper-half hit rows in the atlas
- Reconstruct lower-half rows on-the-fly via `atlas_rows_mirror_lower` (Stage 2)
- Convert forward rows to backward rows via the same symmetry in `overlap_pair`
