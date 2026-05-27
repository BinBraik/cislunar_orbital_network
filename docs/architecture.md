# Pipeline architecture (rs3 / rs4)

## Step summary

The pipeline builds reachable-set atlases for periodic orbit families in the Earth-Moon CR3BP and then finds overlap between two family atlases to identify transfer candidates.

- **Step 2** — Build the 3D (x, y, θ) voxel grid (`atlas_grid_make`, `atlas_grid_validate`).
- **Step 3** — Sample a periodic orbit densely, extract upper-half seeds, build per-family Keep mask (`atlas_family_prepare_seeds`).
- **Step 4** — For each (seed × heading) job, integrate forward and backward, log voxel hits into packed row structs, cache the atlas (`atlas_family_build_hits`, `atlas_cache_save`).
- **Step 5 / Overlap** — Compute the intersection of family A's FRS with family B's BRS. Lower-half rows are reconstructed on-the-fly from upper-half via y-axis symmetry (`overlap_pair`, `atlas_rows_mirror_lower`).
- **Visualization** — Render FRS/BRS occupancy, overlap voxels, zoom insets, and DV-proxy bounds (`overlap_visualize`, `overlap_visualize_combo`, `overlap_visualize_bounds`).

## Key data structures

- `cfg` — configuration struct (see `atlas_cfg_defaults` for the full schema with inline comments)
- `grid3` — grid definition: bin edges, voxel centers, Keep mask, `bin_xyth` function handle
- `S` — family atlas struct (output of Steps 3–4, stored in cache):
  - `S.SeedsUpper`, `S.SeedsLower` — seed positions (Nx3)
  - `S.Step4.rows_FRS_upper`, `S.Step4.rows_BRS_upper` — packed hit row structs
  - `S.grid3` — grid with per-family Keep mask applied
- `O` — overlap struct (output of Step 5):
  - `O.ids` — linear voxel IDs of overlapping voxels
  - `O.x`, `O.y`, `O.th` — voxel center coordinates
  - `O.countA`, `O.countB` — hit counts per family per voxel

## Packed row struct schema

Hit rows are stored as a struct-of-arrays rather than a double matrix (~4× memory reduction).

| Field | Type | Description |
|---|---|---|
| `n` | uint32 | Number of valid rows |
| `iSeed` | uint16 | Seed index |
| `iHead` | uint16 | Heading index within fan |
| `leg` | uint8 | 1 = forward, 2 = backward |
| `halfFlag` | int8 | +1 = upper half, -1 = lower half |
| `t` | single | Time of voxel hit |
| `ix`, `iy`, `it` | uint16 | Voxel bin indices |

## Symmetry convention

The CR3BP has y=0 symmetry: a trajectory `(x, y, θ, t)` in the upper half corresponds to a mirrored trajectory `(x, -y, π−θ, −t)` in the lower half. The pipeline exploits this to:
- Integrate only from upper-half seeds (Steps 3–4)
- Store only upper-half hit rows
- Reconstruct lower-half rows on-the-fly via `atlas_rows_mirror_lower` during overlap (Step 5)

## Caching

Cache keys are MD5 hashes of a fingerprint string that encodes all pipeline parameters (grid, seed, fan, propagation, version tag). A cache hit skips Steps 3–4 entirely. See `atlas_cache_fingerprint`, `atlas_cache_try_load`, `atlas_cache_save`.
