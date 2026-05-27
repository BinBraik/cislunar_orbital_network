# Changelog

## v1.0.0 — 2026-05-27 (initial public release)

### Physics fixes (carried forward from development history)

**CR3BP potential derivative (r^5 computation)**
- File: `src/cr3bp_potential.m`
- The term multiplying the cubic distance factor in the effective-potential gradient used
  `r^7` in an intermediate expression where `r^5` is correct. This produced wrong
  gradient magnitudes and caused trajectory drift over long integration times. Fixed.

**Voxel hit logger — trajectory transpose detection (BUG-4)**
- File: `src/atlas_hits_log_from_traj.m`
- When the ODE solver returned a transposed state matrix `X` (Nx4 instead of 4xN), the
  segment-walk extracted wrong column indices, silently logging all zeros. Fixed by
  detecting transpose via `numel(t)` and reshaping before indexing.

### Repository changes for public release

- Removed large animation artifact (`cycler_bus_animation.gif`, 46 MB); regenerate
  with `run_betweenness_video.m`.
- Renamed all internal pipeline prefixes from version-number style (`rs3_`, `rs4_`,
  `rs5_`) to semantic prefixes aligned with paper sections:
  - `rs3_*` → `atlas_*` / `cr3bp_*` (§3 atlas construction, §2 dynamics)
  - `rs4_*` → `overlap_*` / `traj_diffcorr*` (§4 overlap analysis, §7 correction)
  - `rs5_*` → `traj_*` (§7 trajectory realization)
- Added `LICENSE` (MIT), `CITATION.cff`, and `CHANGELOG.md`.
- Rewrote `README.md` with execution-order diagram, orbit family table, and
  paper section cross-references.
