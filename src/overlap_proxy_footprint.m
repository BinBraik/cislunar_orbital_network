function F = overlap_proxy_footprint(S, grid3, VU_mps, TU_days)
%OVERLAP_PROXY_FOOTPRINT  Build compact per-voxel summary for one atlas family.
%
% Computes, for FRS and BRS separately:
%   uid      — sorted unique voxel IDs  (linear index into [Ny,Nx,Nt])
%   dv_min   — min dv_turn (m/s) over all rows in that voxel
%   t_min    — min |TOF| (days) over all rows in that voxel
%
% FRS  = direct rows (FRS_upper + FRS_lower).
% BRS  = R(FRS): mirror of FRS_upper + mirror of FRS_lower.
%        Mirror: iy → Ny-iy+1,  it → it_lut(it)  (same as overlap_pair).
%
% DV reuse: since U(x,y)=U(x,-y) in CR3BP, the DV computed at a seed and at
% its y-mirror are identical.  BRS voxels therefore reuse the DV values from
% the corresponding FRS rows — no extra potential evaluation needed. TOF is
% likewise reused (time of flight to a seed does not depend on the mirror).
%
% NOTE: per-voxel TOF is aggregated with MIN (not mean) — a voxel's proxy
% TOF is the fastest row landing in it, matching how dv_min already picks
% the cheapest row. This keeps the DV-proxy and TOF-proxy selections
% consistent optimistic bounds, and lets a min-TOF pair-winner search
% (overlap_proxy_pair) reuse the same per-voxel fields as the min-DV search.

Ny = numel(grid3.y_centers);
Nx = numel(grid3.x_centers);
Nt = numel(grid3.th_centers);

% ── Theta mirror LUT (identical formula to overlap_pair) ───────────────
thm = wrap_to_pi(pi - grid3.th_centers(:));
lut = discretize(thm, grid3.th_edges);
lut(isnan(lut)) = 0;
it_lut = uint16(lut);

% ── Pre-build delta-angle lookup matrix (vectorised, avoids cell loop) ──────
dlists = S.Step4.delta_lists;
Ns     = numel(dlists);
max_h  = max(1, max(cellfun(@numel, dlists)));
delta_mat = zeros(Ns, max_h);
for s = 1:Ns
    v = double(dlists{s});
    delta_mat(s, 1:numel(v)) = v;
end

% ── Pre-compute v0 per unique seed (avoids per-row potential evaluation) ────
pot_u = cr3bp_potential(S.SeedsUpper(:,1), S.SeedsUpper(:,2), S.mu);
v0_upper = sqrt(max(2*pot_u.U(:) - S.CJ, 0));   % [Nseeds_upper, 1]

if isfield(S,'SeedsLower') && ~isempty(S.SeedsLower)
    pot_l = cr3bp_potential(S.SeedsLower(:,1), S.SeedsLower(:,2), S.mu);
    v0_lower = sqrt(max(2*pot_l.U(:) - S.CJ, 0));
else
    v0_lower = v0_upper;
end

% ── Process FRS_upper rows ─────────────────────────────────────────────────
nu = double(S.Step4.rows_FRS_upper.n);
if nu > 0
    [ids_u, dv_u, t_u, ix_u, iy_u, it_u] = local_fp_rows( ...
        S.Step4.rows_FRS_upper, nu, v0_upper, ...
        delta_mat, Ns, max_h, Ny, Nx, Nt, VU_mps, TU_days);
else
    ids_u=zeros(0,1); dv_u=zeros(0,1); t_u=zeros(0,1);
    ix_u=zeros(0,1);  iy_u=zeros(0,1); it_u=zeros(0,1);
end

% ── Process FRS_lower rows ─────────────────────────────────────────────────
nl = double(S.Step4.rows_FRS_lower.n);
if nl > 0
    [ids_l, dv_l, t_l, ix_l, iy_l, it_l] = local_fp_rows( ...
        S.Step4.rows_FRS_lower, nl, v0_lower, ...
        delta_mat, Ns, max_h, Ny, Nx, Nt, VU_mps, TU_days);
else
    ids_l=zeros(0,1); dv_l=zeros(0,1); t_l=zeros(0,1);
    ix_l=zeros(0,1);  iy_l=zeros(0,1); it_l=zeros(0,1);
end

% ── Aggregate FRS voxels ───────────────────────────────────────────────────
[F.uid_frs, F.dv_min_frs, F.t_min_frs] = local_fp_agg( ...
    [ids_u; ids_l], [dv_u; dv_l], [t_u; t_l]);

% ── BRS: mirror FRS_upper ──────────────────────────────────────────────────
if ~isempty(ix_u)
    biy_u = Ny - iy_u + 1;
    bit_u = double(it_lut(it_u));
    ok_u  = bit_u > 0;
    ids_brs_u = sub2ind([Ny,Nx,Nt], biy_u(ok_u), ix_u(ok_u), ...
        max(1, min(Nt, bit_u(ok_u))));
    dv_bu = dv_u(ok_u);  t_bu = t_u(ok_u);
else
    ids_brs_u = zeros(0,1);  dv_bu = zeros(0,1);  t_bu = zeros(0,1);
end

% ── BRS: mirror FRS_lower ──────────────────────────────────────────────────
if ~isempty(ix_l)
    biy_l = Ny - iy_l + 1;
    bit_l = double(it_lut(it_l));
    ok_l  = bit_l > 0;
    ids_brs_l = sub2ind([Ny,Nx,Nt], biy_l(ok_l), ix_l(ok_l), ...
        max(1, min(Nt, bit_l(ok_l))));
    dv_bl = dv_l(ok_l);  t_bl = t_l(ok_l);
else
    ids_brs_l = zeros(0,1);  dv_bl = zeros(0,1);  t_bl = zeros(0,1);
end

% ── Aggregate BRS voxels ───────────────────────────────────────────────────
[F.uid_brs, F.dv_min_brs, F.t_min_brs] = local_fp_agg( ...
    [ids_brs_u; ids_brs_l], [dv_bu; dv_bl], [t_bu; t_bl]);

F.CJ   = S.CJ;
F.mu   = S.mu;
F.name = S.name;
end

% ─────────────────────────────────────────────────────────────────────────────
function [ids, dv_mps, t_days, ix_out, iy_out, it_out] = local_fp_rows( ...
        rows, n, v0_per_seed, delta_mat, Ns, max_h, Ny, Nx, Nt, VU_mps, TU_days)
%LOCAL_FP_ROWS  Extract voxel IDs, DV (m/s), and |TOF| (days) for n packed rows.
% v0_per_seed  [Nseeds,1] — pre-computed sqrt(max(2U-CJ,0)) per seed position.
%   Caller computes this once with cr3bp_potential on the seeds
%   matrix (~hundreds of evals) so this function avoids per-row pot evaluation.
ix_out = double(rows.ix(1:n));
iy_out = double(rows.iy(1:n));
it_out = double(rows.it(1:n));
ids    = sub2ind([Ny, Nx, Nt], iy_out, ix_out, it_out);

iSeed = double(rows.iSeed(1:n));
iHead = double(rows.iHead(1:n));
t_nd  = double(rows.t(1:n));

% Vectorised delta-angle lookup
lin   = sub2ind([Ns, max_h], iSeed, iHead);
delta = delta_mat(lin);

% Per-row v0 via seed lookup (O(n) index, no potential evaluation)
v0     = v0_per_seed(iSeed);
dv_mps = 2 * v0(:) .* sin(abs(delta(:)) / 2) * VU_mps;
t_days = abs(t_nd(:)) * TU_days;
end

% ─────────────────────────────────────────────────────────────────────────────
function [uid, dv_min, t_min] = local_fp_agg(ids, dv, t)
%LOCAL_FP_AGG  Aggregate per-voxel min DV and min TOF from raw row arrays.
if isempty(ids)
    uid = zeros(0,1);  dv_min = zeros(0,1);  t_min = zeros(0,1);
    return;
end
[uid, ~, ic] = unique(ids(:));
dv_min = accumarray(ic, dv(:), [], @min);
t_min  = accumarray(ic, t(:),  [], @min);
end
