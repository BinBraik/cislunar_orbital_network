function T = overlap_voxel_traj_extract(SA, SB, V, B, cfg)
%RS4_VOXEL_TRAJ_EXTRACT  Re-integrate the argmin-DV trajectory pair for the
% best overlap voxel and compute the true DV_patch.
%
% Takes the best voxel (B.imin into V), finds the single A-row and B-row
% with minimum DV_turn, densely re-integrates both arcs, finds where each
% arc passes closest to the voxel center (xc,yc), measures the true heading
% difference delta_th, and computes:
%
%   DV_patch_true = 2 * v_box_center * sin(|delta_th| / 2)
%
% where v_box_center uses min(CJ_A, CJ_B) — consistent with the existing
% proxy in overlap_visualize_bounds.
%
% Inputs
%   SA, SB : family structs from atlas_prepare_or_load
%   V      : voxel metadata struct from overlap_extract_voxel_info
%   B      : bounds struct from overlap_visualize_bounds (needs B.imin)
%   cfg    : config struct
%
% Output T (struct) — see field comments below

VU_mps  = local_cfg_get(cfg, 'units.VU_mps',  1.0);
TU_days = local_cfg_get(cfg, 'units.TU_days', 1.0);
relTol  = local_cfg_get(cfg, 'propag.relTol', 1e-9);
absTol  = local_cfg_get(cfg, 'propag.absTol', 1e-9);

grid3 = SA.grid3;
Ny    = grid3.Ny;
Nx    = grid3.Nx;
Nth   = grid3.Nth;

% =========================================================================
% Step 1: Identify best voxel
% =========================================================================
k   = double(B.imin);
vid = double(V.ids(k));
xc  = V.x(k);
yc  = V.y(k);
thc = V.th(k);

fprintf('[overlap_extract] Best voxel k=%d, vid=%d\n', k, vid);
fprintf('[overlap_extract]   center (x,y,th) = (%.6f, %.6f, %.4f deg)\n', ...
    xc, yc, rad2deg(thc));
fprintf('[overlap_extract]   DVproxy = %.3f m/s  (DVturn_A=%.3f + DVpatch_ub=%.3f + DVturn_B=%.3f)\n', ...
    B.min_dvproxy, V.dv_turn_mps_min_A(k), B.dv_patch_ub(k), V.dv_turn_mps_min_B(k));

% =========================================================================
% Step 2: Argmin A-row — find row in FRS with min DV_turn that hits vid
% =========================================================================
rowsA_F = local_rows_cat(SA.Step4.rows_FRS_upper, SA.Step4.rows_FRS_lower);
nA      = double(rowsA_F.n);
idsA    = double(atlas_rows_to_voxelid(rowsA_F, grid3));
idxA    = find(idsA(1:nA) == vid);

if isempty(idxA)
    error('[overlap_extract] No A-side rows found for voxel vid=%d.', vid);
end

[iSeed_A, iHead_A, halfFlag_A, t_A, delta_A] = ...
    local_argmin_dv_row(rowsA_F, idxA, SA);

fprintf('[overlap_extract] A argmin: iSeed=%d iHead=%d hf=%d t=%.5f nd (%.2f days) delta=%.5f rad\n', ...
    iSeed_A, iHead_A, halfFlag_A, t_A, t_A * TU_days, delta_A);

% =========================================================================
% Step 3: Argmin B-row via inverse mirror of vid
% =========================================================================
% BRS = R(FRS): the target voxel in BRS space maps back to FRS voxel vid_frs
[iy_t, ix_t, it_t] = ind2sub([Ny, Nx, Nth], vid);
iy_frs  = Ny - iy_t + 1;                                % y-mirror (self-inverse)
th_t    = grid3.th_centers(it_t);
th_frs  = wrap_to_pi(pi - th_t);                      % theta-mirror (self-inverse)
it_frs  = discretize(th_frs, grid3.th_edges);
vid_frs = sub2ind([Ny, Nx, Nth], iy_frs, ix_t, it_frs);

% FRS_lower rows → become BRS_upper after R; FRS_upper rows → become BRS_lower
nBlo   = double(SB.Step4.rows_FRS_lower.n);
nBup   = double(SB.Step4.rows_FRS_upper.n);
idsB_lo = double(atlas_rows_to_voxelid(SB.Step4.rows_FRS_lower, grid3));
idsB_up = double(atlas_rows_to_voxelid(SB.Step4.rows_FRS_upper, grid3));
idxBlo  = find(idsB_lo(1:nBlo) == vid_frs);
idxBup  = find(idsB_up(1:nBup) == vid_frs);

if isempty(idxBlo) && isempty(idxBup)
    error('[overlap_extract] No B-side rows map to inverse-mirror voxel vid_frs=%d.', vid_frs);
end

[iSeed_B, iHead_B, halfFlag_B_frs, t_B, delta_B, from_lower_B] = ...
    local_argmin_dv_row_B(SB, idxBlo, idxBup);

fprintf('[overlap_extract] B argmin: iSeed=%d iHead=%d hf_frs=%d t=%.5f nd (%.2f days) delta=%.5f rad\n', ...
    iSeed_B, iHead_B, halfFlag_B_frs, t_B, t_B * TU_days, delta_B);

% =========================================================================
% Step 4: Reconstruct initial conditions
% =========================================================================
if halfFlag_A == 1
    seed_A = SA.SeedsUpper(iSeed_A, :);
else
    seed_A = SA.SeedsLower(iSeed_A, :);
end
IC_A = [seed_A(1); seed_A(2); wrap_to_pi(seed_A(3) + delta_A)];

if from_lower_B
    seed_B_frs = SB.SeedsLower(iSeed_B, :);
else
    seed_B_frs = SB.SeedsUpper(iSeed_B, :);
end
IC_B_frs = [seed_B_frs(1); seed_B_frs(2); wrap_to_pi(seed_B_frs(3) + delta_B)];

% =========================================================================
% Step 5: Dense re-integration (no event stops — integrate to exact t)
% =========================================================================
odeOpts = odeset('RelTol', relTol, 'AbsTol', absTol);

[tA_vec, XA] = ode113( ...
    @(tt,XX) cr3bp_reduced_ode(tt, XX, SA.CJ, SA.mu, false), ...
    [0, t_A], IC_A, odeOpts);

[tB_vec, XB_frs] = ode113( ...
    @(tt,XX) cr3bp_reduced_ode(tt, XX, SB.CJ, SB.mu, false), ...
    [0, t_B], IC_B_frs, odeOpts);

% Apply R-transform to B arc: (x, y, theta) -> (x, -y, pi-theta)
x_B  =  XB_frs(:, 1);
y_B  = -XB_frs(:, 2);
th_B =  wrap_to_pi(pi - XB_frs(:, 3));

% =========================================================================
% Step 6: Closest point on each arc to voxel center (xc, yc)
% =========================================================================
[dA, i_star] = min(hypot(XA(:,1) - xc, XA(:,2) - yc));
[dB, j_star] = min(hypot(x_B        - xc, y_B       - yc));

th_A_star = XA(i_star, 3);
th_B_star = th_B(j_star);

% =========================================================================
% Step 7: True DV_patch
% =========================================================================
pot_center   = cr3bp_potential(xc, yc, SA.mu);
v_box_center = sqrt(max(2 * pot_center.U - min(SA.CJ, SB.CJ), 0));
delta_th     = abs(wrap_to_pi(th_A_star - th_B_star));
DV_patch_nd  = 2 * v_box_center * sin(delta_th / 2);
DV_patch_mps = DV_patch_nd * VU_mps;

% =========================================================================
% Step 8: DV_total_true
% =========================================================================
DV_turn_A_mps     = V.dv_turn_mps_min_A(k);
DV_turn_B_mps     = V.dv_turn_mps_min_B(k);
DV_total_true_mps = DV_turn_A_mps + DV_patch_mps + DV_turn_B_mps;

fprintf('[overlap_extract] DVturn_A=%.3f  DVpatch_true=%.3f  DVturn_B=%.3f  => DVtotal=%.3f m/s\n', ...
    DV_turn_A_mps, DV_patch_mps, DV_turn_B_mps, DV_total_true_mps);
fprintf('[overlap_extract] DVproxy was %.3f m/s  => tightening: %.3f m/s\n', ...
    B.min_dvproxy, B.min_dvproxy - DV_total_true_mps);
fprintf('[overlap_extract] delta_th=%.4f deg  miss dA=%.6f nd  dB=%.6f nd\n', ...
    rad2deg(delta_th), dA, dB);

% =========================================================================
% Pack output
% =========================================================================
T = struct();
T.vid              = vid;
T.k                = k;
T.xc               = xc;
T.yc               = yc;
T.thc              = thc;

T.iSeed_A          = iSeed_A;
T.iHead_A          = iHead_A;
T.halfFlag_A       = halfFlag_A;
T.t_A              = t_A;
T.IC_A             = IC_A;
T.seed_A           = seed_A;    % [x, y, th_nominal] on PO_A

T.iSeed_B          = iSeed_B;
T.iHead_B          = iHead_B;
T.halfFlag_B_frs   = halfFlag_B_frs;
T.from_lower_B     = from_lower_B;
T.t_B              = t_B;
T.IC_B_frs         = IC_B_frs;
T.seed_B_frs       = seed_B_frs;   % [x, y, th_nominal] on PO_B (FRS side)

T.tA_vec           = tA_vec;
T.XA               = XA;           % [n x 3]: x, y, theta along A arc

T.tB_vec           = tB_vec;
T.x_B              = x_B;          % BRS arc after R-transform
T.y_B              = y_B;
T.th_B             = th_B;

T.i_star           = i_star;       % index on A arc closest to voxel center
T.j_star           = j_star;       % index on B arc closest to voxel center
T.dA_nd            = dA;           % miss distance A (nd)
T.dB_nd            = dB;           % miss distance B (nd)
T.th_A_star        = th_A_star;    % heading of A at closest point
T.th_B_star        = th_B_star;    % heading of B (BRS) at closest point
T.delta_th_rad     = delta_th;     % |theta_A - theta_B| at closest points
T.v_box_center_nd  = v_box_center;

T.DV_patch_nd      = DV_patch_nd;
T.DV_patch_mps     = DV_patch_mps;
T.DV_turn_A_mps    = DV_turn_A_mps;
T.DV_turn_B_mps    = DV_turn_B_mps;
T.DV_total_true_mps = DV_total_true_mps;
T.DV_proxy_mps     = B.min_dvproxy;

T.tof_A_days       = t_A * TU_days;
T.tof_B_days       = t_B * TU_days;
end

% =========================================================================
% Local helpers
% =========================================================================

function [iSeed, iHead, halfFlag, t, delta] = ...
        local_argmin_dv_row(rows, idxInRows, S)
% Find the row in rows(idxInRows) with the minimum DV_turn.
% Returns the metadata of that argmin row.

iSeed_v   = double(rows.iSeed(idxInRows));
iHead_v   = double(rows.iHead(idxInRows));
halfFlag_v = double(rows.halfFlag(idxInRows));
t_v        = double(rows.t(idxInRows));

n = numel(idxInRows);
delta_v = zeros(n, 1);
seed_x  = zeros(n, 1);
seed_y  = zeros(n, 1);

for i = 1:n
    s = iSeed_v(i);
    h = iHead_v(i);
    delta_v(i) = double(S.Step4.delta_lists{s}(h));

    if halfFlag_v(i) >= 0
        seed_x(i) = S.SeedsUpper(s, 1);
        seed_y(i) = S.SeedsUpper(s, 2);
    else
        seed_x(i) = S.SeedsLower(s, 1);
        seed_y(i) = S.SeedsLower(s, 2);
    end
end

pot  = cr3bp_potential(seed_x, seed_y, S.mu);
v0   = sqrt(max(2 * pot.U - S.CJ, 0));
dv   = 2 * v0 .* sin(abs(delta_v) / 2);

[~, iBest] = min(dv);

iSeed    = iSeed_v(iBest);
iHead    = iHead_v(iBest);
halfFlag = halfFlag_v(iBest);
t        = t_v(iBest);
delta    = delta_v(iBest);
end

% -------------------------------------------------------------------------

function [iSeed, iHead, halfFlag_frs, t, delta, from_lower] = ...
        local_argmin_dv_row_B(SB, idxBlo, idxBup)
% Find the argmin DV_turn row across both FRS_lower (idxBlo) and
% FRS_upper (idxBup) of SB.  from_lower=true means the winning row
% came from rows_FRS_lower (halfFlag_frs = -1), false means FRS_upper (+1).

candidates = struct('iSeed',{},'iHead',{},'halfFlag_frs',{},'t',{},'delta',{},'dv',{});

if ~isempty(idxBlo)
    rows = SB.Step4.rows_FRS_lower;
    [iS,iH,hf,tv,delt] = local_argmin_dv_row(rows, idxBlo, SB);
    pot = cr3bp_potential(SB.SeedsLower(iS,1), SB.SeedsLower(iS,2), SB.mu);
    v0  = sqrt(max(2*pot.U - SB.CJ, 0));
    dv0 = 2 * v0 * sin(abs(delt)/2);
    candidates(end+1) = struct('iSeed',iS,'iHead',iH,'halfFlag_frs',hf,'t',tv,'delta',delt,'dv',dv0);
end

if ~isempty(idxBup)
    rows = SB.Step4.rows_FRS_upper;
    [iS,iH,hf,tv,delt] = local_argmin_dv_row(rows, idxBup, SB);
    pot = cr3bp_potential(SB.SeedsUpper(iS,1), SB.SeedsUpper(iS,2), SB.mu);
    v0  = sqrt(max(2*pot.U - SB.CJ, 0));
    dv0 = 2 * v0 * sin(abs(delt)/2);
    candidates(end+1) = struct('iSeed',iS,'iHead',iH,'halfFlag_frs',hf,'t',tv,'delta',delt,'dv',dv0);
end

dvs = [candidates.dv];
[~, iBest] = min(dvs);

iSeed        = candidates(iBest).iSeed;
iHead        = candidates(iBest).iHead;
halfFlag_frs = candidates(iBest).halfFlag_frs;
t            = candidates(iBest).t;
delta        = candidates(iBest).delta;
from_lower   = (iBest == 1) && ~isempty(idxBlo);
end

% -------------------------------------------------------------------------

function r = local_rows_cat(a, b)
if isstruct(a) && isstruct(b)
    nA = double(a.n); nB = double(b.n);
    if nA == 0, r = b; return; end
    if nB == 0, r = a; return; end
    r = atlas_rows_empty(nA + nB);
    r.n        = uint32(nA + nB);
    r.iSeed(1:nA)    = a.iSeed(1:nA);    r.iSeed(nA+1:end)    = b.iSeed(1:nB);
    r.iHead(1:nA)    = a.iHead(1:nA);    r.iHead(nA+1:end)    = b.iHead(1:nB);
    r.leg(1:nA)      = a.leg(1:nA);      r.leg(nA+1:end)      = b.leg(1:nB);
    r.halfFlag(1:nA) = a.halfFlag(1:nA); r.halfFlag(nA+1:end) = b.halfFlag(1:nB);
    r.t(1:nA)        = a.t(1:nA);        r.t(nA+1:end)        = b.t(1:nB);
    r.ix(1:nA)       = a.ix(1:nA);       r.ix(nA+1:end)       = b.ix(1:nB);
    r.iy(1:nA)       = a.iy(1:nA);       r.iy(nA+1:end)       = b.iy(1:nB);
    r.it(1:nA)       = a.it(1:nA);       r.it(nA+1:end)       = b.it(1:nB);
    return;
end
if isempty(a), r = b; return; end
if isempty(b), r = a; return; end
r = [a; b];
end

% -------------------------------------------------------------------------

function v = local_cfg_get(cfg, path, defaultVal)
v = defaultVal;
try
    parts = strsplit(path, '.');
    cur = cfg;
    for i = 1:numel(parts)
        k = parts{i};
        if ~isstruct(cur) || ~isfield(cur, k), return; end
        cur = cur.(k);
    end
    if ~isempty(cur), v = cur; end
catch
    v = defaultVal;
end
end
