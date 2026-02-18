# Cleanup audit for current runner goal (two-family overlap visualization)

This audit assumes the **current primary goal** is exactly what
`scripts/run_rs4_overlap_and_visuals.m` does: build/load two atlases,
compute overlap, and render overlap visuals.

## Scope used for this audit

- Entry runner: `scripts/run_rs4_overlap_and_visuals.m`
- Transitive calls through atlas build/cache, overlap, and visualization paths.
- We intentionally did **not** include legacy refinement/scoring pipeline stages,
  because they are not invoked by this runner.

## Files that appear unnecessary for this goal

These source files are present but not on the transitive path of
`run_rs4_overlap_and_visuals`:

- `src/rs3_cache_visual_validate.m`
- `src/rs3_family_visual_validate.m`
- `src/rs3_grid_visual_validate.m`
- `src/rs3_hits_visual_validate.m`
- `src/rs3_plot_pair_from_disk.m`
- `src/rs3_rows_vcat.m`

Interpretation:
- First four are diagnostics/plots used by the **single-family atlas plot runner**.
- `rs3_plot_pair_from_disk` is an extra utility, not needed by rs4 overlap runner.
- `rs3_rows_vcat` looks like a legacy helper superseded by packed-row concat.

## Configuration fields required by the rs4 overlap runner path

The following config fields are referenced by the runner + transitive function path
(excluding defaults/validation self-references):

- `cache.dir`, `cache.enable`, `cache.rebuild`, `cache.store_dense_po`, `cache.version_tag`
- `diag.po_stride`, `diag.progress`, `diag.smoke_checks`, `diag.zoom.enable`, `diag.zoom.xlim`, `diag.zoom.ylim`
- `fan.DV_cap_nd`, `fan.dtheta_fan`
- `grid.Rdom`, `grid.dtheta`, `grid.dx`, `grid.dy`, `grid.enforce_xy_symmetry`, `grid.enforce_y0_edge`
- `io.fig_resolution`, `io.fig_subdir`, `io.fig_visible`, `io.out_root`, `io.save_fig`, `io.save_figs`, `io.tag`, `io.verbose`
- `log.maxstep_factor`, `log.segwalk.enable`, `log.segwalk.frac`, `log.step_len_factor`
- `overlap.primary_buffer_frac`
- `par.enable`, `par.progress_every`
- `propag.Tmax`, `propag.absTol`, `propag.relTol`, `propag.v2tol`
- `seed.N_dense`, `seed.Tf_scale`, `seed.ds_seed`, `seed.minSegPts`, `seed.y_eps`
- `sys.RE_nd`, `sys.RM_nd`

## High-confidence config cleanup candidates (not needed for rs4 overlap runner)

If your repo goal is now strictly up to overlap visualization, these groups are
safe candidates to remove **after adjusting validation/defaults accordingly**:

- Entire legacy/new refinement blocks:
  - `cfg.refine.*`
  - `cfg.refine7.*`
- Step-8 scoring block:
  - `cfg.score.*`
- Candidate-retention block:
  - `cfg.cand.*`
- Family-selection defaults not used by this runner:
  - `cfg.families.*` (runner sets `famA/famB` directly)
- Unit-reporting block if you don't use printed unit conversions:
  - `cfg.units.*`
- Extra system metadata not used in overlap path:
  - `cfg.sys.model`, `cfg.sys.mu`, `cfg.sys.CJ_policy`, `cfg.sys.CJ_fixed`
- Back-compat aliases that can be dropped after code cleanup:
  - `cfg.fan.dtheta`
  - `cfg.par.use`
  - `cfg.log.maxStep_factor`
- Output toggles not used by this runner path:
  - `cfg.io.save_mat`, `cfg.io.save_csv`
- Unused parallel knobs for this path:
  - `cfg.par.mode`, `cfg.par.pool_size`, `cfg.par.batch_size`
- Misc diagnostics not used by this path:
  - `cfg.diag.enable_asserts`, `cfg.diag.estimate_memory`, `cfg.diag.print_cfg`,
    `cfg.diag.plot_each_level`, `cfg.diag.maxPlotFootprint`, `cfg.diag.maxPlotOverlap`,
    `cfg.diag.plot_style`, `cfg.diag.story_plots`, `cfg.diag.show_orbits`,
    `cfg.diag.run_step3`, `cfg.diag.run_step4`, `cfg.diag.zoom.show_gridlines`,
    `cfg.diag.zoom.max_gridlines`
- Cache knob currently not consumed on this path:
  - `cfg.cache.store_entry_state`

## Recommended cleanup order

1. **Freeze target scope**: keep only rs4 overlap workflow + atlas build/cache.
2. Remove unneeded **files** listed above.
3. Trim `rs3_cfg_defaults` to only required fields for this scope.
4. Update `rs3_cfg_validate` to validate only remaining schema.
5. Run MATLAB tests + one rs4 overlap runner smoke execution.
6. Then optionally reintroduce minimal docs for any future stages.

## Caveat

This is a usage-based audit for the current runner goal.
If you plan to revive Step 7/8 refinement/scoring later, keep those config blocks
in a separate `rsX_experimental_cfg_defaults` rather than deleting permanently.
