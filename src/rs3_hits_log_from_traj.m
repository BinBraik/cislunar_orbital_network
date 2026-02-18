function rows = rs3_hits_log_from_traj(iSeed, iHead, leg, halfFlag, t, X, grid3, cfg)
%RS3_HITS_LOG_FROM_TRAJ  Log (x,y,theta) voxel hits along a trajectory.
%
% Returns a PACKED ROW STRUCT (rs3_rows_empty schema) for memory efficiency.
% Legacy callers expecting [N x 8] double should use rs3_rows_to_matrix().
%
% Logging uses segment-wise oversampling ("segment-walk") to avoid missing
% voxels when ODE steps are large.
%
% CHANGES from original:
%   - Output is packed struct-of-arrays (~4x less memory)
%   - Pre-allocated buffer (avoids O(n^2) row-by-row growth)
%   - BUG FIX: Transpose detection uses numel(t) to disambiguate 3xN vs Nx4

assert(isvector(t) && ~isempty(t), 't must be a nonempty vector');
t = t(:);
nPts = numel(t);

% --- BUG-4 FIX: Robust transpose detection ---
% Use numel(t) as the authoritative row count. If X has nPts along dim 1,
% it's already row-major (nPts x cols). Otherwise transpose.
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
    Xr = [X(:,1), X(:,2), rs3_wrapToPi(atan2(yd, xd))];
else
    error('X must have 3 (reduced) or 4 (full) columns (got %dx%d).', size(X,1), size(X,2));
end

% --- Pre-allocate buffer (Phase 5) ---
% Estimate: worst case ~10 voxels per ODE step (generous)
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

prev_ix = uint16(0); prev_iy = uint16(0); prev_it = uint16(0);
prevSet = false;

for k = 1:(nPts - 1)
    p0 = Xr(k,:); p1 = Xr(k+1,:);
    t0 = t(k);    t1 = t(k+1);

    % Unwrap theta across seam along shortest direction
    dth_c = rs3_circ_diff(p1(3), p0(3));
    th1u = p0(3) + dth_c;

    dx_seg = p1(1) - p0(1);
    dy_seg = p1(2) - p0(2);
    dth_seg = th1u - p0(3);

    nsub = 1;
    if doSegwalk
        nsub = ceil(max([abs(dx_seg)/(max(frac*dx, eps)), ...
                         abs(dy_seg)/(max(frac*dy, eps)), ...
                         abs(dth_seg)/(max(frac*dth, eps))]));
        nsub = max(nsub, 1);
    end

    for j = 0:nsub
        f = j / nsub;
        xj  = p0(1) + f * dx_seg;
        yj  = p0(2) + f * dy_seg;
        thj = rs3_wrapToPi(p0(3) + f * dth_seg);
        tj  = t0 + f * (t1 - t0);

        [cix, ciy, cit] = grid3.bin_xyth(xj, yj, thj);
        if isnan(cix) || isnan(ciy) || isnan(cit)
            continue;
        end
        cix = uint16(cix); ciy = uint16(ciy); cit = uint16(cit);

        if ~prevSet || cix ~= prev_ix || ciy ~= prev_iy || cit ~= prev_it
            nRow = nRow + 1;
            % Grow buffer if needed (doubling strategy)
            if nRow > bufSize
                bufSize = bufSize * 2;
                buf_iSeed(bufSize)    = 0;
                buf_iHead(bufSize)    = 0;
                buf_leg(bufSize)      = 0;
                buf_halfFlag(bufSize) = 0;
                buf_t(bufSize)        = 0;
                buf_ix(bufSize)       = 0;
                buf_iy(bufSize)       = 0;
                buf_it(bufSize)       = 0;
            end
            buf_iSeed(nRow)    = u_iSeed;
            buf_iHead(nRow)    = u_iHead;
            buf_leg(nRow)      = u_leg;
            buf_halfFlag(nRow) = s_halfFlag;
            buf_t(nRow)        = single(tj);
            buf_ix(nRow)       = cix;
            buf_iy(nRow)       = ciy;
            buf_it(nRow)       = cit;

            prev_ix = cix; prev_iy = ciy; prev_it = cit;
            prevSet = true;
        end
    end
end

% Trim to actual size and pack into struct
rows = rs3_rows_empty();
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
