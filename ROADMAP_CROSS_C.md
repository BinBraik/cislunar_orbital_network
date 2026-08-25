# Cross-C Extension Roadmap

Plan for extending the same-Jacobi-constant impulsive transfer network to
**varying Jacobi constants**: each family `F_i` contributes a representative
`γ_i(C_i)` at its own energy level, producing the family of weighted graphs

> G(C_1, C_2, …, C_13)

with minimum impulsive transfer costs as edge weights. This document is the
working plan for the successor repository (cloned from
`cislunar_orbital_network` with history intact). The paper version of the
original repository stays frozen.

Pipeline chain being built toward:

> same-C impulsive transfers → **cross-C impulsive transfers (this work)** → low-thrust transfers

Scientific question: how much of the current network structure (e.g., the
prominence of C32) is an artifact of the common-Jacobi-constant restriction,
and how much persists when each family may use a more favorable
representative?

---

## 0. Repository migration (one-time, manual)

Clone with history so provenance and future diffs against the paper version
are preserved — do **not** copy files into a fresh `git init`:

```bash
# 1. Freeze the paper version in the original repo (immutable reference)
cd cislunar_orbital_network
git tag -a v1.0-paper -m "Frozen version accompanying the 2026 paper"
git push origin v1.0-paper

# 2. Create the new (empty) repository on GitHub, e.g. cislunar_crossC_network
#    (no README/license — it must be empty to accept a full-history push)

# 3. Clone and repoint
git clone https://github.com/BinBraik/cislunar_orbital_network.git cislunar_crossC_network
cd cislunar_crossC_network
git remote rename origin upstream          # keep a link to the paper repo
git remote add origin https://github.com/BinBraik/cislunar_crossC_network.git
git push -u origin main
```

Bug fixes made in either repo can then be cherry-picked across via the
`upstream` remote.

---

## 1. Core formulation

### 1.1 Separability — this is NOT a 13-dimensional problem

Each edge weight `w_ij` of `G(C_1,…,C_13)` depends **only on the pair
(C_i, C_j)**. Therefore the object to compute is 78 two-parameter maps

```
w_ij : (C_i, C_j) → min transfer ΔV        (one map per unordered family pair)
```

after which **any** graph in the 13-parameter family — including the
best-representative graph `argmin` over all C — is assembled by table lookup
at zero cost. Time-reversal (y-axis) symmetry gives
`w_ij(a, b) = w_ji(b, a)`, so the graph remains undirected and one map per
unordered pair suffices; Floyd–Warshall, centrality, and articulation
analysis carry over unchanged.

Cost scaling with `N_C` energy levels per family:
- `13 × N_C` atlas builds — linear, embarrassingly parallel, cached.
- `78 × N_C²` overlap computations — but overlap is a cheap voxel-ID
  intersection on the shared grid (minutes, not hours). The expensive stage
  stays linear in `N_C`.

### 1.2 Transfer model: 3-impulse cross-C scheme

Atlases stay **single-C** (the reduced (x, y, θ) model cannot change C
mid-coast — speed is slaved to C). The entire energy change is concentrated
in the one patch impulse:

1. **Turn impulse** on γ_i — constant-C_i heading change (existing fan,
   unchanged: `dv_turn = 2 v sin(δ/2)`).
2. **Coast** on the C_i manifold — FRS of family i built at C_i.
3. **Patch impulse** at the overlap voxel — changes heading *and* speed:

   ```
   ΔV_patch = |v_B − v_A| (as vectors)
            = sqrt( vA² + vB² − 2 vA vB cos Δθ )
   vA = sqrt(2U(x,y) − C_i),   vB = sqrt(2U(x,y) − C_j)
   ```

4. **Coast** on the C_j manifold — BRS of family j built at C_j.
5. **Turn impulse** onto γ_j — constant-C_j heading change.

Key properties:
- **Closed form, no search.** The (x, y, θ) grid plus a family's C determine
  the full 4D state at every voxel; the patch cost is one vectorized
  `cr3bp_potential` call plus a few flops — computationally identical to the
  current `dv_patch_ub` in `overlap_visualize_bounds.m`, with two speeds
  instead of one.
- **No C-targeting.** Matching the full arrival velocity *vector* lands
  exactly on the C_j manifold with the required heading; energy matching is
  implied by vector matching.
- **Degenerate case = current code.** At C_i = C_j the formula collapses to
  `2 v sin(Δθ/2)`.
- **Proxy bias washes out at DC.** `traj_diffcorr` (fmincon) has no
  constraint on C, so the corrector is free to redistribute the energy
  change across impulses; concentrating ΔC at the patch only biases the
  proxy slightly upward, same relationship as today.
- **Interpretation shift.** Cross-C, even a perfect (x, y, θ) match carries
  an irreducible cost floor `|vA − vB| = |C_i − C_j| / (vA + vB)`; overlap
  no longer certifies *cheap*, the minimum over the overlap region does.

### 1.3 Physics identities (used for screening and sanity checks)

- Impulse–energy relation: `ΔC = −(2 v·Δv + |Δv|²)` per impulse.
- Lower bound to bridge a gap `|ΔC|` at local speed v:
  `|Δv| ≥ sqrt(v² + |ΔC|) − v` (Oberth: energy changes are cheap where U,
  hence v, is large — near the Moon).
- **Sanity prediction:** cross-C minimum-cost voxels should migrate toward
  low-altitude (high-speed) regions relative to the same-C winners. If the
  first cross-C maps do not show this tendency, audit the patch-cost
  implementation before trusting network results.

---

## 2. Work phases

### Phase 1 — Family continuation and IC tables (cheap; do first, completely)

Every family already has a converged member at C = 3.1294 (`cr3bp_family_ic`)
— continuation starts from known solutions; nothing has to be re-found.

- **Pseudo-arclength continuation** along each family with C as an *output*
  (not the continuation parameter — C-targeting fails at folds where C
  reverses along the family, likely the cause of the families that could not
  be converged to exactly 3.1294; see note in `cr3bp_family_ic.m`).
- **Multiple shooting** (10–20 patch points) for the long-period cyclers
  (42–56 days) — single-shooting monodromy correction is ill-conditioned at
  these periods.
- **Exploit x-axis symmetry** where present: perpendicular-crossing
  conditions at y = 0 halve the integration span and the unknowns.
  Asymmetric families fall back to full multiple shooting.
- Small arclength steps; log C, period, stability, fold locations along the
  way.
- Harvest members near target C levels; polish each with one final
  differential correction constrained to C = C_target exactly.

**Deliverables**
- Per-family continuation dataset: table of `(X0, Tf, C)` along the family.
- Generalized `cr3bp_family_ic(name, C)`: data-backed (interpolate + polish
  DC), replacing the hardcoded switch.
- **Family existence map**: 13 curves in the (C, family) plane — the domain
  of `G(C_1,…,C_13)`, and a candidate first figure of the paper. Note the
  formulation only needs each family valid on its *own* range; no common
  level is required (that is exactly the restriction being lifted).
- Mind ZVC topology transitions (Earth–Moon: C_L1 ≈ 3.188, C_L2 ≈ 3.172,
  C_L3 ≈ 3.012): crossing C_L3 opens the exterior realm and changes the
  network qualitatively.

### Phase 2 — Zero-integration screening (uses existing C = 3.1294 atlases)

- **Oberth lower bounds** per pair: for candidate (C_i, C_j), evaluate
  `sqrt(v² + |ΔC|) − v` voxel-wise over the existing overlap region and add
  the existing turn-cost LB → lower bound on the cross-C edge. Identifies
  which pairs / C-gaps could possibly beat the same-C edge before any new
  atlas is built. Pairs provably C-insensitive need no new atlases at all.
- **Family-ladder upper bounds**: continuation along family i from C_i to
  C* (small tangential impulses) + existing same-C edge at C* → upper
  bracket.
- **Sensitivities** `∂w/∂C` from stored voxel metadata (`∂v/∂C = −1/(2v)`).
- Output: per-pair brackets that direct where Phase-3 compute is spent.

### Phase 3 — Coarse shared C-grid atlases

- **Fixed voxel grid across all C** (non-negotiable — overlap must remain a
  pure index intersection). Only the Keep mask / ZVC forbidden region varies
  per (family, C); that is already per-family via `S.grid3`.
- 3–5 shared levels per family within its existence range → ~50–65 atlas
  builds ≈ 4–5× the current STEP-1 cost, fully parallel; the existing cache
  fingerprint (cfg includes CJ) makes this incremental and restartable.
- Seeds come from Phase-1 `cr3bp_family_ic(name, C)`.

### Phase 4 — Cross-C overlap and the 78 edge maps

- `overlap_pair` intersection is unchanged (same grid for all C).
- Generalize the patch cost to the two-speed law-of-cosines formula
  (§1.2); currently `CJstar = min(SA.CJ, SB.CJ)` in
  `overlap_visualize_bounds.m` already tolerates slightly mismatched CJ —
  replace with the exact two-speed form.
- Produce `w_ij(C_i, C_j)` heatmaps for all 78 pairs. **These maps are
  themselves the scientific deliverable**: they answer directly how much of
  C32's prominence is a same-C artifact.
- **Adaptive refinement on shared levels only**: refine the C grid where a
  map's minimum is sharp or the network ranking is undecided, and refine on
  levels shared across pairs so each new atlas serves all 12 edges of its
  family. Avoid per-pair independent 1-D optimization — it fragments atlas
  reuse.

### Phase 5 — Network assembly and analysis

- Assemble `G(C_1,…,C_13)` by lookup; compute the best-representative graph
  (per-map argmin) and compare against the same-C baseline network
  (hubs, gateways, relays, bottlenecks; C32's role).
- Run the expensive (DV_cap, Tmax) budget-plane sweep **only on the final
  selected graph(s)**, not across the C grid.

### Phase 6 — Trajectory realization

- `overlap_voxel_traj_extract` + `traj_diffcorr` on selected cross-C edges;
  the DC formulation needs no change (no C constraint). Verify the corrected
  ΔV against the proxy and the Oberth floor.

---

## 3. Code-touch inventory

Existing machinery reused unchanged: grid, atlas builder, cache, packed
rows, symmetry mirroring, overlap intersection, network stage, DC stage.

| Change | Where | Size |
|---|---|---|
| Two-speed patch cost | `overlap_visualize_bounds`, `overlap_extract_voxel_info` | small |
| Per-(family, C) seeds | `cr3bp_family_ic(name, C)` + continuation module | **new module** (main work) |
| Per-C Keep mask / ZVC | `atlas_keep_mask_xy` usage per (family, C) | small |
| Cross-C runners | `run_family_continuation.m`, `run_atlas_c_grid.m`, `run_overlap_cross_c.m`, `run_crossc_edge_maps.m` | new scripts |
| Screening | `run_crossc_screening.m` (Phase 2 bounds on existing atlases) | new script |

---

## 4. Validation checklist

- [ ] Continuation reproduces the 13 stored ICs at C = 3.1294 to tolerance.
- [ ] Cross-C patch formula reduces to current `dv_patch_ub` when C_i = C_j.
- [ ] Cross-C edge maps satisfy the Phase-2 brackets (LB ≤ w ≤ UB) everywhere.
- [ ] Diagonal of each `w_ij` map at the common level matches the published
      same-C matrix (`minDVproxy_matrix.csv`).
- [ ] Oberth sanity check: cross-C minimum-cost voxels concentrate at low
      altitude (high v) as |C_i − C_j| grows.
- [ ] Symmetry check: `w_ij(a, b) = w_ji(b, a)` numerically.

---

## 5. Open questions

- Repo name and whether the web visualization should be carried over or
  redesigned around the (C_i, C_j) edge maps.
- Per-family C ranges: truncate at family terminations/bifurcations, or also
  at ZVC realm transitions for interpretability?
- Whether to expose direction-dependent costs later (drop the undirected
  assumption) in preparation for the low-thrust extension.
