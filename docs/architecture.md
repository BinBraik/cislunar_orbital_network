# Pipeline architecture (rs3 / rs4)

## Step summary

The pipeline builds reachable-set atlases for periodic orbit families in the Earth-Moon CR3BP and then finds overlap between two family atlases to identify transfer candidates.

- **Step 2** — Build the 3D (x, y, θ) voxel grid (`rs3_grid_make`, `rs3_grid_validate`).
- **Step 3** — Sample a periodic orbit densely, extract upper-half seeds, build per-family Keep mask (`rs3_family_prepare_seeds`).
- **Step 4** — For each (seed × heading) job, integrate forward and backward, log voxel hits into packed row structs, cache the atlas (`rs3_family_build_hits`, `rs3_cache_save_family`).
- **Step 5 / Overlap** — Compute the intersection of family A's FRS with family B's BRS. Lower-half rows are reconstructed on-the-fly from upper-half via y-axis symmetry (`rs4_overlap_pair`, `rs3_rows_mirror_lower`).
- **Visualization** — Render FRS/BRS occupancy, overlap voxels, zoom insets, and DV-proxy bounds (`rs4_overlap_visualize`, `rs4_overlap_visualize_combo`, `rs4_overlap_visualize_bounds`).

## Key data structures

- `cfg` — configuration struct (see `rs3_cfg_defaults` for the full schema with inline comments)
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
- Reconstruct lower-half rows on-the-fly via `rs3_rows_mirror_lower` during overlap (Step 5)

## Caching

Cache keys are MD5 hashes of a fingerprint string that encodes all pipeline parameters (grid, seed, fan, propagation, version tag). A cache hit skips Steps 3–4 entirely. See `rs3_cache_fingerprint_family`, `rs3_cache_try_load_family`, `rs3_cache_save_family`.

---

## Voxel metadata and differential correction (rs4)

After the overlap step, `rs4_overlap_extract_voxel_info` computes per-voxel metadata:

- `V` struct — per-voxel DV bounds and estimated TOF:
  - `V.dv_lb_mps` — lower bound: `DV_turn_A_min + DV_turn_B_min`
  - `V.dv_patch_ub_mps` — patch upper bound: `2·v_box·sin(|Δθ|/2)`
  - `V.dv_proxy_mps` — proxy: lb + patch estimate
  - `V.tof_days` — mean transfer time estimate

`rs4_voxel_traj_extract` re-integrates the two argmin-DV arcs for the best voxel and computes the true `DV_total`:

```
DV_total_true = DV_turn_A_min + DV_patch_true + DV_turn_B_min
```

where `DV_patch_true` is measured at the closest points of each re-integrated arc to the voxel centre.

`rs4_diffcorr` then solves for continuous patch continuity using `fmincon` with two constraint modes:

- **Standard**: enforce position + velocity continuity at the patch point.
- **relax_heading**: enforce only spatial continuity; allow heading mismatch to be absorbed into the DV budget (useful for hard cases).

The resulting `Tc` struct replaces the proxy with a tightened, physically continuous transfer.

---

## Full trajectory assembly (rs5)

`rs5_build_full_traj` assembles three segments into a single continuous-time trajectory struct:

1. **Departure arc** — FRS arc from origin PO to patch point (forward integration).
2. **Patch gap** — the DV_patch impulse (represented as a gap in position).
3. **Arrival arc** — BRS arc from patch point to destination PO (time-reversed BRS, run forward).

`rs5_visualize_full_traj` plots all three segments with direction-of-motion arrows on each leg and labels showing the three impulse magnitudes (DV_turn_A, DV_patch, DV_turn_B).

Output `result.mat` fields:
- `T` — raw (pre-DC) trajectory struct
- `Tc` — differential-corrected trajectory struct
- `traj_raw`, `traj_dc` — arc arrays for re-plotting without re-running the pipeline

---

## Network analysis (net_*)

Inputs: an N×N DV matrix from `run_rs4_dv_tmax_sweep` or `run_rs4_dc_sweep_tmax`.

**Graph construction** (`net_build_graph`): directed weighted graph where edge (i, j) exists iff `DV(i,j) ≤ DV_budget`. Edge weight = DV in m/s.

**All-pairs shortest paths** (`net_floyd_warshall`): Floyd-Warshall with path reconstruction. Returns distance matrix `D` and predecessor matrix `P` for tracing optimal multi-hop routes.

**Largest Connected Component** (`net_lcc`): extracts the LCC of the undirected version of the graph.

**Centrality metrics** (`net_centrality`): five metrics are computed at each (DV_cap, Tmax) snapshot:

| Metric | Definition |
|---|---|
| Strength | Sum of outgoing edge weights (total reachable DV) |
| Harmonic closeness | Mean reciprocal shortest-path distance to all reachable nodes |
| Betweenness | Fraction of shortest paths (across all pairs) that pass through this node |
| PageRank | Stationary distribution of a random walk with teleportation |
| Articulation | Flag: node whose removal disconnects the graph |

**Winner maps** (`net_plot_winner_map`): for each (DV_cap, Tmax) cell, the family with the highest score per metric is recorded. Plotted as a colour-coded contour map over the parameter space.

---

## Betweenness explainer

Purpose: demonstrate that a high-betweenness orbit acts as a routing hub — routing through it is cheaper than flying direct.

For a bridge family B and all (origin A, destination C) pairs not involving B, the savings are:

```
savings(A, C) = DV_direct(A, C) − [DV(A, B) + DV(B, C)]
```

The top-N pairs by savings are selected. For each, the script:

1. Extracts Leg-1 arcs (A → B via FRS/BRS overlap) and Leg-2 arcs (B → C) using `rs4_voxel_traj_extract`.
2. Computes a coast arc along the bridge PO connecting the two patch seeds.
3. Plots the full three-segment path: origin PO (red) → Transfer A → bridge PO coast (green) → Transfer B → destination PO (blue).
4. Optionally animates the spacecraft along the path as a GIF with arc-length-uniform speed.

The coast arc is computed by integrating the bridge PO from the Leg-1 patch seed to the Leg-2 patch seed along the PO in the correct direction.
