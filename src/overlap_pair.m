function O = overlap_pair(SA, SB, cfg)
%OVERLAP_PAIR  Overlap voxels between A.FRS_full and B.BRS_full,
% with conservative filtering:
%   - forbidden region removed via KeepA & KeepB
%   - buffers around primaries: r > (1+bufFrac)*R
%
% Output O:
%   ids, ix,iy,it, x,y,th, countA,countB, Nx,Ny,Nt
%
% Performance rewrite:
%   - Mirror + cat eliminated: IDs extracted directly from ix/iy/it fields
%     of the upper-half rows; no 8-field packed struct is allocated for
%     mirrored or concatenated row sets.
%   - Theta mirror uses a precomputed LUT (one uint16 array lookup per row)
%     instead of per-row wrapToPi + discretize calls.
%   - unique() called once per set; uid/cnt reused for both intersect and
%     count queries.
%   - local_counts_via_unique replaces ismember(N_rows, K_overlap): query
%     is O(K log N_unique) vs old O(N_rows log K).

if nargin < 3, cfg = struct(); end

% ---------------- grid dims (robust) ----------------
grid3 = SA.grid3;

Nx = numel(grid3.x_centers);
Ny = numel(grid3.y_centers);
Nt = numel(grid3.th_centers);

gridB = SB.grid3;
assert(Nx==numel(gridB.x_centers) && Ny==numel(gridB.y_centers) && Nt==numel(gridB.th_centers), ...
    'A and B must share the same grid dims (x_centers/y_centers/th_centers).');

% ---------------- theta mirror LUT ----------------
% Precompute: for each theta bin index it (1..Nt), the bin index of
% wrapToPi(pi - th_center(it)).  LUT value 0 = no valid mirror bin.
it_lut = local_build_theta_mirror_lut(grid3, Nt);

% ---------------- voxel IDs (BRS = R(FRS), no backward rows stored) ----
% A.FRS_full = FRS_upper (direct) + FRS_lower (direct, from lower-seed integration)
idsA = local_frs_vid(SA.Step4.rows_FRS_upper, SA.Step4.rows_FRS_lower, Ny, Nx, Nt);

% B.BRS_full = R(B.FRS_full) = mirror(FRS_upper) + mirror(FRS_lower)
idsB = local_brs_vid(SB.Step4.rows_FRS_upper, SB.Step4.rows_FRS_lower, Ny, Nx, Nt, it_lut);

% ---------------- unique once; reuse for intersect and counts ----------------
[uid_A, ~, ic_A] = unique(idsA);   cnt_A = accumarray(ic_A, 1);
[uid_B, ~, ic_B] = unique(idsB);   cnt_B = accumarray(ic_B, 1);

idsO = intersect(uid_A, uid_B);    % uid_A/B already sorted — fast merge
idsO = idsO(:);

O = struct('ids',idsO,'ix',[],'iy',[],'it',[],'x',[],'y',[],'th',[], ...
           'countA',[],'countB',[],'Nx',Nx,'Ny',Ny,'Nt',Nt);

if isempty(idsO)
    fprintf('[overlap] overlap voxels: 0\n');
    return;
end

% ---------------- unpack overlap IDs ----------------
[iy, ix, it] = ind2sub([Ny, Nx, Nt], idsO);
ix = ix(:); iy = iy(:); it = it(:);

% ---------------- FILTER: Keep + buffered primaries ----------------
bufFrac = 0.05;
if isfield(cfg,'overlap') && isfield(cfg.overlap,'primary_buffer_frac') && ~isempty(cfg.overlap.primary_buffer_frac)
    bufFrac = cfg.overlap.primary_buffer_frac;
end

keepA = local_get_keep_xy(SA, Ny, Nx);
keepB = local_get_keep_xy(SB, Ny, Nx);
keepXY = keepA & keepB;

if ~(isfield(cfg,'sys') && isfield(cfg.sys,'RE_nd') && isfield(cfg.sys,'RM_nd'))
    error('cfg.sys.RE_nd and cfg.sys.RM_nd are required for primary buffer filtering.');
end
RE = cfg.sys.RE_nd;
RM = cfg.sys.RM_nd;
mu = SA.mu;

x = grid3.x_centers(ix);
y = grid3.y_centers(iy);

rE = hypot(x + mu, y);
rM = hypot(x - (1-mu), y);

okKeep  = keepXY(sub2ind([Ny, Nx], iy, ix));
okEarth = rE > (1 + bufFrac)*RE;
okMoon  = rM > (1 + bufFrac)*RM;

ok = okKeep(:) & okEarth(:) & okMoon(:);

idsO = idsO(ok);
ix   = ix(ok);
iy   = iy(ok);
it   = it(ok);

O.ids = idsO;

if isempty(idsO)
    O.ix=[]; O.iy=[]; O.it=[]; O.x=[]; O.y=[]; O.th=[];
    O.countA=[]; O.countB=[];
    fprintf('[overlap] overlap voxels after Keep+buffer: 0\n');
    return;
end

O.ix = ix; O.iy = iy; O.it = it;
O.x  = grid3.x_centers(ix);
O.y  = grid3.y_centers(iy);
O.th = grid3.th_centers(it);

% ---------------- counts per kept overlap voxel ----------------
% uid_A/B and cnt_A/B already computed above — no re-sorting needed.
O.countA = local_counts_via_unique(uid_A, cnt_A, idsO);
O.countB = local_counts_via_unique(uid_B, cnt_B, idsO);

fprintf('[overlap] overlap voxels after Keep+buffer: %d\n', numel(O.ids));

end

% =======================================================================
% Helpers
% =======================================================================

function it_lut = local_build_theta_mirror_lut(grid3, Nt)
%LOCAL_BUILD_THETA_MIRROR_LUT
% Precompute the mirrored theta bin for each of the Nt theta bins.
% Mirror formula: thm = wrapToPi(pi - th_center(it)).
% LUT value 0 means "invalid mirror bin" (defensive; should not occur for
% a well-formed symmetric theta grid).
thm = wrap_to_pi(pi - grid3.th_centers(:));   % (Nt x 1)
lut = discretize(thm, grid3.th_edges);           % (Nt x 1), NaN if out of range
lut(isnan(lut)) = 0;
it_lut = uint16(lut);
end


function ids = local_frs_vid(rows_upper, rows_lower, Ny, Nx, Nt)
%LOCAL_FRS_VID  FRS full voxel IDs = upper rows (direct) + lower rows (direct).
% No mirror needed: both halves are stored from actual forward integration.

nu = double(rows_upper.n);
if nu > 0
    iu_x = max(1, min(Nx, double(rows_upper.ix(1:nu))));
    iu_y = max(1, min(Ny, double(rows_upper.iy(1:nu))));
    iu_t = max(1, min(Nt, double(rows_upper.it(1:nu))));
    ids_u = sub2ind([Ny, Nx, Nt], iu_y, iu_x, iu_t);
else
    ids_u = zeros(0, 1);
end

nl = double(rows_lower.n);
if nl > 0
    il_x = max(1, min(Nx, double(rows_lower.ix(1:nl))));
    il_y = max(1, min(Ny, double(rows_lower.iy(1:nl))));
    il_t = max(1, min(Nt, double(rows_lower.it(1:nl))));
    ids_l = sub2ind([Ny, Nx, Nt], il_y, il_x, il_t);
else
    ids_l = zeros(0, 1);
end

ids = [ids_u(:); ids_l(:)];
end


function ids = local_brs_vid(rows_frs_upper, rows_frs_lower, Ny, Nx, Nt, it_lut)
%LOCAL_BRS_VID  BRS full voxel IDs = R(FRS_full) via CR3BP time-reversal.
% R(FRS_upper) contributes to BRS_lower; R(FRS_lower) contributes to BRS_upper.

ids_a = local_mirror_rows_to_ids(rows_frs_upper, Ny, Nx, Nt, it_lut);
ids_b = local_mirror_rows_to_ids(rows_frs_lower, Ny, Nx, Nt, it_lut);
ids = [ids_a(:); ids_b(:)];
end


function ids = local_mirror_rows_to_ids(rows, Ny, Nx, Nt, it_lut)
%LOCAL_MIRROR_ROWS_TO_IDS  Apply R map to rows and return voxel IDs.
nm = double(rows.n);
if nm == 0
    ids = zeros(0, 1);
    return;
end
im_x = double(rows.ix(1:nm));
im_y = Ny - double(rows.iy(1:nm)) + 1;         % y mirror
im_t = double(it_lut(double(rows.it(1:nm))));   % theta LUT lookup
ok   = im_t > 0 & im_y >= 1 & im_y <= Ny;      % skip invalid mirrors
im_x = max(1, min(Nx, im_x(ok)));
im_y = im_y(ok);
im_t = max(1, min(Nt, im_t(ok)));
ids  = sub2ind([Ny, Nx, Nt], im_y, im_x, im_t);
end


function c = local_counts_via_unique(uid_all, cnt_all, idsO)
%LOCAL_COUNTS_VIA_UNIQUE  Count how many raw rows in the full set fall into
% each voxel of idsO, using a pre-computed unique decomposition.
%
% Inputs:
%   uid_all  — sorted unique voxel IDs from the full set
%   cnt_all  — hit counts per unique voxel (from accumarray)
%   idsO     — query IDs (the filtered overlap voxels)
%
% This is O(K log N') where K = numel(idsO) and N' = numel(uid_all),
% replacing the old ismember(N_rows, K) pattern which was O(N log K)
% with N being the full (non-unique) row count.
if isempty(idsO)
    c = zeros(0, 1);
    return;
end
[tf, loc] = ismember(idsO(:), uid_all);
c = zeros(numel(idsO), 1);
c(tf) = cnt_all(loc(tf));
end


function keep = local_get_keep_xy(S, Ny, Nx)
assert(isfield(S,'grid3') && isfield(S.grid3,'Keep') && ~isempty(S.grid3.Keep), ...
    'Keep mask missing: expected S.grid3.Keep (built in Step 3).');

keep = logical(S.grid3.Keep);
sz = size(keep);

if isequal(sz, [Ny, Nx])
    return;
elseif isequal(sz, [Nx, Ny])
    keep = keep.';
    return;
else
    error('Keep mask size [%d,%d], expected [%d,%d] (or transposed).', sz(1), sz(2), Ny, Nx);
end
end
