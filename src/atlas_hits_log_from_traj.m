function rows = atlas_hits_log_from_traj(iSeed, iHead, leg, halfFlag, t, X, grid3, cfg)
%ATLAS_HITS_LOG_FROM_TRAJ  Log (x,y,theta) voxel hits along a trajectory.
%
% Returns a PACKED ROW STRUCT (atlas_rows_empty schema) for memory efficiency.
% Legacy callers expecting [N x 8] double should use atlas_rows_to_matrix().
%
% Logging uses segment-wise oversampling ("segment-walk") to avoid missing
% voxels when ODE steps are large.
%
% Phase 2 rewrite:
%   - nsub and theta-unwrap pre-computed for all segments before the k-loop
%     (one vectorized circ_diff call instead of one per segment).
%   - Inner j-loop eliminated: all sub-step operations use column-vector
%     arithmetic — wrapToPi, discretize, keepXY lookup (sub2ind), and
%     run-length deduplication (diff) are each one vectorized call per segment.
%   - Buffer append is a bulk range assignment (not element-by-element).

assert(isvector(t) && ~isempty(t), 't must be a nonempty vector');
t = t(:);
nPts = numel(t);

% --- BUG-4 FIX: Robust transpose detection ---
% Use numel(t) as the authoritative row count.
if ~isempty(X) && isnumeric(X)
    if size(X, 1) ~= nPts && size(X, 2) == nPts
        X = X.';
    end
end

% Convert states to reduced [x y theta]
ncols = size(X, 2);
if ncols == 3
    Xr = X;
elseif ncols == 4
    xd = X(:,3); yd = X(:,4);
    Xr = [X(:,1), X(:,2), wrap_to_pi(atan2(yd, xd))];
else
    error('X must have 3 (reduced) or 4 (full) columns (got %dx%d).', size(X,1), size(X,2));
end

% --- Pre-allocate buffer ---
bufSize = max(64, (nPts - 1) * 4);
buf_iSeed    = zeros(bufSize, 1, 'uint16');
buf_iHead    = zeros(bufSize, 1, 'uint16');
buf_leg      = zeros(bufSize, 1, 'uint8');
buf_halfFlag = zeros(bufSize, 1, 'int8');
buf_t        = zeros(bufSize, 1, 'single');
buf_ix       = zeros(bufSize, 1, 'uint16');
buf_iy       = zeros(bufSize, 1, 'uint16');
buf_it       = zeros(bufSize, 1, 'uint16');

nRow = 0;

% Cast metadata once
u_iSeed    = uint16(iSeed);
u_iHead    = uint16(iHead);
u_leg      = uint8(leg);
s_halfFlag = int8(halfFlag);

frac = cfg.log.segwalk.frac;
dx = grid3.dx; dy = grid3.dy; dth = grid3.dtheta;
doSegwalk = cfg.log.segwalk.enable;

% Optional center-based admissibility mask (Ny x Nx)
useKeep = false;
keepXY = [];
if isfield(grid3,'Keep') && ~isempty(grid3.Keep)
    keepXY = logical(grid3.Keep);
    if isequal(size(keepXY), [grid3.Ny, grid3.Nx])
        useKeep = true;
    elseif isequal(size(keepXY), [grid3.Nx, grid3.Ny])
        keepXY = keepXY.';
        useKeep = true;
    end
end

if nPts < 2
    rows = atlas_rows_empty();
    rows.n = uint32(0);
    return;
end

% ==========================================================
% Phase 2: Pre-compute per-segment quantities (all k at once)
% ==========================================================
p0_mat = Xr(1:end-1, :);    % (nSeg x 3)
p1_mat = Xr(2:end,   :);    % (nSeg x 3)

dx_segs  = p1_mat(:,1) - p0_mat(:,1);              % (nSeg x 1)
dy_segs  = p1_mat(:,2) - p0_mat(:,2);
dth_segs = circ_diff(p1_mat(:,3), p0_mat(:,3)); % vectorized, (nSeg x 1)
dt_segs  = t(2:end) - t(1:end-1);

if doSegwalk
    nsub_arr = max(1, ceil(max([ ...
        abs(dx_segs)  / max(frac*dx,  eps), ...
        abs(dy_segs)  / max(frac*dy,  eps), ...
        abs(dth_segs) / max(frac*dth, eps) ], [], 2)));
else
    nsub_arr = ones(nPts-1, 1);
end

% ==========================================================
% Segment loop — inner j-loop replaced by column-vector ops
% ==========================================================
prev_ix = uint16(0); prev_iy = uint16(0); prev_it = uint16(0);
prevSet = false;

for k = 1:(nPts-1)
    nsub  = nsub_arr(k);
    f_arr = (0:nsub)' / nsub;    % (nsub+1) x 1

    % Sub-step positions — vectorized over j
    xj  = p0_mat(k,1) + f_arr * dx_segs(k);
    yj  = p0_mat(k,2) + f_arr * dy_segs(k);
    thj = wrap_to_pi(p0_mat(k,3) + f_arr * dth_segs(k));  % (nsub+1) x 1
    tj  = single(t(k) + f_arr * dt_segs(k));

    % Vectorized binning (discretize accepts column vectors)
    % thj already in [-pi,pi) — no second wrapToPi needed
    cix_d = discretize(xj,  grid3.x_edges);
    ciy_d = discretize(yj,  grid3.y_edges);
    cit_d = discretize(thj, grid3.th_edges);

    % In-bounds mask: discretize returns NaN for out-of-domain inputs
    valid = ~isnan(cix_d) & ~isnan(ciy_d) & ~isnan(cit_d);

    % Keep mask: vectorized lookup via sub2ind (avoids per-element branching)
    if useKeep && any(valid)
        vi      = find(valid);
        lin_idx = sub2ind([grid3.Ny, grid3.Nx], ciy_d(vi), cix_d(vi));
        valid(vi) = keepXY(lin_idx);
    end

    if ~any(valid)
        continue;
    end

    % Cast to uint16 after NaN/keep filtering (uint16 cannot represent NaN)
    cix_v = uint16(cix_d(valid));
    ciy_v = uint16(ciy_d(valid));
    cit_v = uint16(cit_d(valid));
    tj_v  = tj(valid);
    nv    = numel(cix_v);

    % Run-length deduplication: drop sub-steps that repeat the previous voxel
    if prevSet
        first_new = cix_v(1) ~= prev_ix || ciy_v(1) ~= prev_iy || cit_v(1) ~= prev_it;
    else
        first_new = true;
    end

    if nv > 1
        same  = (diff(uint32(cix_v)) == 0) & ...
                (diff(uint32(ciy_v)) == 0) & ...
                (diff(uint32(cit_v)) == 0);
        is_new = [first_new; ~same(:)];
    else
        is_new = first_new;
    end

    cix_new = cix_v(is_new);
    ciy_new = ciy_v(is_new);
    cit_new = cit_v(is_new);
    tj_new  = tj_v(is_new);
    n_new   = numel(cix_new);

    if n_new == 0
        continue;
    end

    % Grow buffer if needed (doubling strategy)
    needed = nRow + n_new;
    if needed > bufSize
        newSize = max(bufSize * 2, needed + 64);
        buf_iSeed(newSize)    = uint16(0);
        buf_iHead(newSize)    = uint16(0);
        buf_leg(newSize)      = uint8(0);
        buf_halfFlag(newSize) = int8(0);
        buf_t(newSize)        = single(0);
        buf_ix(newSize)       = uint16(0);
        buf_iy(newSize)       = uint16(0);
        buf_it(newSize)       = uint16(0);
        bufSize = newSize;
    end

    % Bulk append (one range assignment per field, not per element)
    range = nRow + (1:n_new);
    buf_iSeed(range)    = u_iSeed;
    buf_iHead(range)    = u_iHead;
    buf_leg(range)      = u_leg;
    buf_halfFlag(range) = s_halfFlag;
    buf_t(range)        = tj_new;
    buf_ix(range)       = cix_new;
    buf_iy(range)       = ciy_new;
    buf_it(range)       = cit_new;
    nRow = nRow + n_new;

    % Update previous voxel state
    prev_ix = cix_new(end);
    prev_iy = ciy_new(end);
    prev_it = cit_new(end);
    prevSet = true;
end

% Trim to actual size and pack into struct
rows = atlas_rows_empty();
rows.n        = uint32(nRow);
rows.iSeed    = buf_iSeed(1:nRow);
rows.iHead    = buf_iHead(1:nRow);
rows.leg      = buf_leg(1:nRow);
rows.halfFlag = buf_halfFlag(1:nRow);
rows.t        = buf_t(1:nRow);
rows.ix       = buf_ix(1:nRow);
rows.iy       = buf_iy(1:nRow);
rows.it       = buf_it(1:nRow);
end
