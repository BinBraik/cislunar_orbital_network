# Pipeline architecture (rs3)

The code is organized around the numbered pipeline steps used in the `rs3_step*` functions.

## Step summary

- **Step 2** — Build the 3D grid (`rs3_grid_make`) and validate the Keep mask.
- **Step 4/5** — For each family, integrate seed trajectories, bin hits into voxels, and cache an atlas.
- **Step 6** — Compute overlap between two family atlases (`rs3_step6_overlap_pair`).
- **Step 7** — Refine candidate seed pairs (new refinement: `rs3_step7_new_refine_pair`).
- **Step 8** — Score candidates and reconstruct best transfer (`rs3_step8_score_pair`).
- **Step 9** — Run many pairs, write matrix summaries.

## Key data structures

- `cfg` — configuration struct (grid, propagation, fan, refinement, scoring, IO)
- `grid3` — grid definition (bin edges, Keep mask, voxel indexing)
- `A`, `B` — family atlas structs (cached output of Step 4/5)
- `O` — overlap struct (Step 6)
- `R` — refinement struct (Step 7)
- `P` — scoring packet (Step 8), including `P.best`
